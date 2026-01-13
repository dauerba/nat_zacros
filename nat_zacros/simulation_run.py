"""
SimulationRun class for managing Zacros simulation runs with multiple trajectories.

This module provides a high-level interface for loading, caching, and analyzing
collections of trajectories from a single Zacros simulation run.
"""

import json
import pickle
import numpy as np
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from .lattice import lattice
from .trajectory import trajectory

class SimulationRun:
    """
    Manages a Zacros simulation run with multiple trajectories.
    
    This class provides a high-level interface for:
    - Loading multiple trajectories with equilibration cutoff
    - Automatic caching of parsed trajectories
    - Ensemble-averaged analysis (RDF, energy statistics)
    - Metadata extraction from jobs.log
    
    Attributes
    ----------
    eq_method : str
		Equilibration detection method ('fraction' by default)
    fraction : float
        Fraction of trajectory to keep from end:
            - 1.0 = keep full trajectory (default)
            - 0.7 = keep last 70% (discard first 30%)
            - 0.5 = keep last 50% (discard first 50%)
    lattice : lattice
        Shared lattice object for all trajectories
    metadata : dict
        Simulation metadata (temperature, coverage, interactions, etc.)
    results_dir : Path
        Directory for storing cache files
    run_dir : Path
        Directory containing trajectory folders (traj_1, traj_2, ...)
    trajectories : list of trajectory
        Loaded trajectory objects
    traj_dirs : list of Path
        Paths to individual trajectory directories
    
    Examples
    --------
    >>> # Typical workflow
    >>> run = SimulationRun('fn_3leed/jobs/1')
    >>> run.load_trajectories()  # Uses cache if available
    >>> r, g, g_std = run.get_ensemble_rdf(r_max=40.0, dr=0.1)
    >>> times, energies, energies_std = run.get_ensemble_energy_vs_time()
    """
    
    def __init__(self, run_dir):
        """
        Initialize a SimulationRun.
        
        Parameters
        ----------
        run_dir : str or Path
            Path to simulation run directory (e.g., 'fn_3leed/jobs/1')
            This directory should contain traj_1, traj_2, ... subdirectories
        
        Notes
        -----
        The fraction should be determined from exploratory energy-only analysis
        """

        self.eq_method = 'fraction'  # Default equilibration method
        self.run_dir = Path(run_dir)
        self.fraction = 1.0  # Default to full trajectory
        
        # Validate run directory exists
        if not self.run_dir.exists():
            raise FileNotFoundError(f"Run directory not found: {self.run_dir}")
        
        # Auto-detect trajectory directories
        self.traj_dirs = sorted([
            d for d in self.run_dir.iterdir() 
            if d.is_dir() and d.name.startswith('traj_')
        ])
        
        if len(self.traj_dirs) == 0:
            raise ValueError(
                f"No trajectory directories (traj_*) found in {self.run_dir}\n"
                f"Expected directories like traj_1, traj_2, etc."
            )
        
        # Create lattice from first trajectory
        self.lattice = lattice(self.traj_dirs[0])
        
        # Initialize trajectory list (filled by load_trajectories)
        self.trajectories = []
        
        # Set up results directory (../../results/ from run_dir)
        jobs_dir = self.run_dir.parent
        interaction_dir = jobs_dir.parent
        self.results_dir = interaction_dir / 'results'
        self.results_dir.mkdir(exist_ok=True)
        
        # Load metadata from jobs.log
        self._load_metadata()
        
    def _load_metadata(self):
        """
        Load simulation metadata from jobs.log file.
        
        Parses the jobs.log file to extract temperature, coverage, interactions,
        and lattice dimensions for this specific run.
        
        Raises
        ------
        FileNotFoundError
            If jobs.log is not found in parent directory
        ValueError
            If run number is not found in jobs.log
        """
        jobs_log = self.run_dir.parent / 'jobs.log'
        
        if not jobs_log.exists():
            raise FileNotFoundError(
                f"jobs.log not found at: {jobs_log}\n"
                f"SimulationRun requires jobs.log in the parent directory."
            )
        
        # Extract run number from directory name (e.g., '1' from 'fn_3leed/jobs/1')
        run_number = int(self.run_dir.name)
        
        # Parse jobs.log
        with open(jobs_log, 'r') as f:
            header = f.readline().split()  # Read header line
            log_entries = [json.loads(line) for line in f if line.strip()]
        
        # Find entry matching this run number
        matching_entry = None
        for entry in log_entries:
            if entry[0] == run_number:
                matching_entry = entry
                break
        
        if matching_entry is None:
            raise ValueError(
                f"Run number {run_number} not found in {jobs_log}\n"
                f"Available run numbers: {[e[0] for e in log_entries]}"
            )
        
        # Extract metadata from log entry
        # Format: [run_num, job_name, [nx, ny], [n_ads], temp, interaction_info, ...]
        self.metadata = {
            'run_number': matching_entry[0],
            'job_name': matching_entry[1],
            'lattice_dimensions': matching_entry[2],  # [nx, ny]
            'n_cells': matching_entry[2][0] * matching_entry[2][1],
            'n_adsorbates': matching_entry[3][0],
            'temperature': matching_entry[4],  # K
            'coverage': matching_entry[3][0] / (matching_entry[2][0] * matching_entry[2][1]),
            'interactions': matching_entry[5][1:]
        }

    def _load_single_trajectory(self, traj_dir, energy_only):
        """
        Helper function for parallel trajectory loading.
        Parameters
        ----------
        traj_dir : Path
            Directory containing trajectory data
        energy_only : bool
            If True, only load times and energies (much faster).
            If False, load full state configurations.
        Returns
        -------
        trajectory
            Trajectory with equilibrated states loaded
        """
        traj = trajectory(self.lattice, traj_dir)
        traj.load_trajectory(fraction=self.fraction, load_energy=True, energy_only=energy_only)
        return traj

    def load_trajectories_parallel(self, energy_only=False, n_workers=None):
        """
        Load multiple trajectories in parallel.
        Parameters
        ----------
        energy_only : bool, default False
            If True, only load times and energies (much faster).
            If False, load full state configurations.
        n_workers : int, optional
            Number of parallel workers. If None, uses all available cores.
        Returns
        -------
        list of trajectories
            Loaded trajectories with equilibrated states
        """
        try:
            from tqdm import tqdm
            use_tqdm = True
        except ImportError:
            use_tqdm = False
        if len(self.traj_dirs) == 0:
            return []
        if n_workers is None:
            n_workers = mp.cpu_count()
        print(f"Loading {len(self.traj_dirs)} trajectories in parallel using {n_workers} workers...")
        energy_onlys = [energy_only] * len(self.traj_dirs)
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            if use_tqdm:
                trajs = list(tqdm(executor.map(self._load_single_trajectory, self.traj_dirs, energy_onlys),
                                total=len(self.traj_dirs),
                                desc="Loading trajectories",
                                unit="traj"))
            else:
                trajs = list(executor.map(self._load_single_trajectory, self.traj_dirs, energy_onlys))
        print(f"Successfully loaded {len(trajs)} trajectories")

        return trajs


    def load_trajectories(self, use_cache=True, parallel=True, energy_only=False, verbose=False):
        """
        Load trajectory data with caching support.

        Parameters
        ----------
        use_cache : bool, default True
            If True, load from cache file if available, otherwise parse and cache.
            If False, always parse from history_output.txt files.
        parallel : bool, default True
            If True, use parallel loading (recommended for full-state data).
            If False, use sequential loading.
        energy_only : bool, default False
            If True, only load energy data from trajectories (faster, less memory).
            If False, load full trajectory data.
        verbose : bool, default False
            If True, print detailed loading information.

        Notes
        -----
        Cache files are stored as:
        - Trajectories: results/{run_number}_trajs_eq.pkl

        The loaded trajectories use the fraction specified during initialization.
        """
        cache_file = self.results_dir / f"{self.metadata['run_number']}_trajs_eq.pkl"

        # Try loading from cache
        if use_cache and cache_file.exists():
            if verbose: print(f"Loading trajectories from cache: {cache_file.name}")
            with open(cache_file, 'rb') as f:
                self.trajectories = pickle.load(f)
            if verbose: print(f"Loaded {len(self.trajectories)} cached trajectories")
            return

        # Load trajectories from files
        if verbose:
            print(f"Loading {len(self.traj_dirs)} trajectories...")
            print(f"  Equilibration fraction: {self.fraction} (keeping last {self.fraction*100:.0f}%)")
            print(f"  Loading mode: {'parallel' if parallel else 'sequential'}")
            print(f"  Energy only: {energy_only}")

        if parallel:
            # Use parallel loading
            self.trajectories = self.load_trajectories_parallel(
                energy_only=energy_only,
                n_workers=None
            )
        else:
            # Sequential loading
            self.trajectories = []
            for traj_dir in self.traj_dirs:
                self.trajectories.append(self._load_single_trajectory(traj_dir, energy_only=energy_only))
        if verbose:
            print(f"Loaded {len(self.trajectories)} trajectories")
            print(f"  States per trajectory: {len(self.trajectories[0].states)}")
            print(f"  Total states: {sum(len(t.states) for t in self.trajectories)}")

        # Save to cache
        if use_cache:
            if verbose: print(f"Saving to cache: {cache_file.name}")
            with open(cache_file, 'wb') as f:
                pickle.dump(self.trajectories, f)
            size_mb = cache_file.stat().st_size / 1024**2
            if verbose: print(f"Cache saved: {size_mb:.1f} MB")

    def get_ensemble_rdf(self, r_max=40.0, dr=0.1):
        """
        Compute ensemble-averaged radial distribution function.
        
        Automatically handles g_ref caching and computation.
        
        Parameters
        ----------
        r_max : float, default 40.0
            Maximum distance for RDF (Angstroms)
        dr : float, default 0.1
            Bin width for RDF (Angstroms)
        
        Returns
        -------
        r : ndarray
            Distance values (Angstroms)
        g_avg : ndarray
            Ensemble-averaged RDF
        g_std : ndarray
            Standard deviation of RDF across trajectories
        
        Notes
        -----
        g_ref is cached at the interaction level (results/gref.pkl) since it
        depends only on lattice geometry, not on temperature or coverage.
        
        Raises
        ------
        RuntimeError
            If trajectories have not been loaded yet
        """
        if len(self.trajectories) == 0:
            raise RuntimeError(
                "No trajectories loaded. Call load_trajectories() first."
            )
        
        # Check for cached g_ref
        gref_cache_file = self.results_dir / 'gref.pkl'
        
        if gref_cache_file.exists():
            print(f"Loading g_ref from cache: {gref_cache_file.name}")
            with open(gref_cache_file, 'rb') as f:
                r_ref, g_ref = pickle.load(f)
        else:
            print(f"Computing g_ref (one-time calculation)...")
            r_ref, g_ref = self.trajectories[0].get_g_ref(r_max=r_max, dr=dr)
            print(f"Saving g_ref to cache: {gref_cache_file.name}")
            with open(gref_cache_file, 'wb') as f:
                pickle.dump((r_ref, g_ref), f)
        
        # Compute RDF for each trajectory
        print(f"Computing RDF for {len(self.trajectories)} trajectories...")
        rdfs = []
        for i, traj in enumerate(self.trajectories):
            r, g = traj.get_rdf(r_max=r_max, dr=dr, g_ref=g_ref, vectorized=True)
            rdfs.append(g)
        
        # Ensemble average
        g_avg = np.mean(rdfs, axis=0)
        g_std = np.std(rdfs, axis=0)
        
        print(f"RDF computation complete")
        return r, g_avg, g_std
    
    def get_ensemble_energy_vs_time(self, n_bins=100):
        """
        Compute ensemble-averaged energy as function of time.
        
        Uses a two-stage averaging approach for robustness:
        1. Intra-trajectory: Average all energy measurements within each time bin
           for each trajectory independently
        2. Inter-trajectory: Average the binned results across all trajectories
        
        This approach:
        - Ensures equal weighting of trajectories (each contributes one value per bin)
        - Uses all available data points (no interpolation artifacts)
        - Handles uneven sampling naturally (bins with more data get better statistics)
        - Preserves measured values without artificial smoothing
        
        Parameters
        ----------
        n_bins : int, default 100
            Number of time bins for averaging
        
        Returns
        -------
        time_centers : ndarray
            Time bin centers
        energy_avg : ndarray
            Ensemble-averaged energy at each time
        energy_std : ndarray
            Standard deviation of energy across trajectories (trajectory-to-trajectory variation)
        
        Raises
        ------
        RuntimeError
            If trajectories have not been loaded yet
        
        Notes
        -----
        Alternative approaches and their trade-offs:
        
        1. **Interpolation**: Interpolate each trajectory to common time points, then average.
           - Pros: Simple, guaranteed uniform sampling
           - Cons: Creates artificial values, smooths real fluctuations, wastes data
        
        2. **Global binning**: Bin all trajectories together with count tracking.
           - Pros: Uses all data, no interpolation
           - Cons: Trajectories with more samples get higher weight (unequal weighting)
        
        3. **Two-stage (this method)**: Bin within each trajectory, then ensemble average.
           - Pros: Equal trajectory weighting + uses all data + no artificial smoothing
           - Cons: Slightly more complex implementation
        
        The two-stage approach is preferred for ensemble statistics as it combines the
        benefits of both interpolation (equal trajectory weights) and binning (uses all
        available data without artificial smoothing).
        """
        if len(self.trajectories) == 0:
            raise RuntimeError(
                "No trajectories loaded. Call load_trajectories() first."
            )
        
        # Find common time range across all trajectories
        end_time = min([traj.times[-1] for traj in self.trajectories])
        start_time = max([traj.times[0] for traj in self.trajectories])
        
        # Create time bins for discretization
        time_bins = np.linspace(start_time, end_time, n_bins + 1)
        time_centers = 0.5 * (time_bins[:-1] + time_bins[1:])
        
        # STAGE 1: Intra-trajectory averaging
        # For each trajectory, bin its energy measurements and average within bins
        energy_hists = []
        for traj in self.trajectories:
            times, energies = traj.get_energy_vs_time()
            
            # Initialize binned energy and sample counts for this trajectory
            energy_hist = np.zeros(n_bins)
            counts = np.zeros(n_bins)
            
            # Accumulate energy measurements into bins
            for t, energy in zip(times, energies):
                if start_time <= t <= end_time:
                    bin_idx = np.digitize(t, time_bins, right=False) - 1
                    if 0 <= bin_idx < n_bins:
                        energy_hist[bin_idx] += energy
                        counts[bin_idx] += 1  # Track number of samples per bin
            
            # Average within each bin (multiple measurements → single value per bin)
            # This gives us one representative energy value per time bin for THIS trajectory
            with np.errstate(divide='ignore', invalid='ignore'):
                energy_hist = np.where(counts > 0, energy_hist / counts, np.nan)
            
            energy_hists.append(energy_hist)
        
        # STAGE 2: Inter-trajectory (ensemble) averaging
        # Each trajectory now contributes exactly one value per time bin
        # Average across trajectories with equal weighting
        energy_avg = np.nanmean(energy_hists, axis=0)
        energy_std = np.nanstd(energy_hists, axis=0)  # Trajectory-to-trajectory variation
        
        return time_centers, energy_avg, energy_std
    
    def clear_cache(self, trajectories=True, gref=False):
        """
        Clear cached data files.
        
        Parameters
        ----------
        trajectories : bool, default True
            If True, clear trajectory cache for this run
        gref : bool, default False
            If True, clear g_ref cache (affects all runs in this interaction set)
        
        Notes
        -----
        Use gref=True with caution - this will force recomputation of g_ref
        for all simulation runs in this interaction set.
        """
        if trajectories:
            cache_file = self.results_dir / f"{self.metadata['run_number']}_trajs_eq.pkl"
            if cache_file.exists():
                cache_file.unlink()
                print(f"Cleared trajectory cache: {cache_file.name}")
            else:
                print(f"No trajectory cache to clear")
        
        if gref:
            gref_file = self.results_dir / 'gref.pkl'
            if gref_file.exists():
                gref_file.unlink()
                print(f"Cleared g_ref cache: {gref_file.name}")
                print(f"WARNING: This affects all runs in this interaction set!")
            else:
                print(f"No g_ref cache to clear")
    
    def __repr__(self):
        """String representation of SimulationRun."""
        if len(self.trajectories) > 0:
            n_states = len(self.trajectories[0].states)
            traj_info = f", {len(self.trajectories)} trajectories ({n_states} states each)"
        else:
            traj_info = f", {len(self.traj_dirs)} trajectories (not loaded)"
        
        return (
            f"SimulationRun(run={self.metadata['run_number']}, "
            f"T={self.metadata['temperature']}K, "
            f"θ={self.metadata['coverage']:.3f}, "
            f"fraction={self.fraction}{traj_info})"
        )
    
    def __len__(self):
        """Return number of trajectories."""
        return len(self.trajectories)
