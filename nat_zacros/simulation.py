"""
Simulation class for managing Zacros simulation with multiple trajectories.

This module provides a high-level interface for loading, caching, and analyzing
collections of trajectories from a single Zacros simulation.
"""

import json
import pickle
import numpy as np
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from .lattice import Lattice
from .trajectory import Trajectory

class Simulation:
    """
    Manages a Zacros simulation with multiple trajectories.
    
    This class provides a high-level interface for:
    - Loading multiple trajectories
    - Automatic caching of parsed trajectories
    - Ensemble-averaged analysis (RDF, energy statistics)
    - Metadata extraction from jobs.log
    
    Attributes
    ----------
    is_valid : bool
        True if something wrong with simulation data (corrupted, missing etc.)
    lattice : Lattice
        Shared lattice object for all trajectories
    metadata : dict
        Simulation metadata (temperature, coverage, interactions, etc.)
    results_dir : Path
        Directory for storing cache and results files
    simulation_dir : Path
        Directory containing trajectory folders (traj_1, traj_2, ...)
    trajectories : list of trajectory
        Loaded trajectory objects
    traj_dirs : list of Path
        Paths to individual trajectory directories
    
    Examples
    --------
    >>> # Typical workflow
    >>> run = Simulation('fn_3leed/jobs/1')
    >>> run.load()  # Uses cache if available
    >>> run.is_valid  # Check if run data is valid
    >>> r, g, g_std = run.get_ensemble_rdf(r_max=40.0, dr=0.1)
    >>> times, energies, energies_std = run.get_ensemble_energy_vs_time()
    """

    def __init__(self, simulation_dir, metadata=None, log_file='jobs.log', results_dirname='results'):
        """
        Initialize a simulation.
        
        Parameters
        ----------
        simulation_dir : str or Path
            Path to simulation directory (e.g., 'fn_3leed/jobs/1')
            This directory should contain traj_1, traj_2, ... subdirectories
        metadata : dict, optional
            Pre-loaded metadata (if available)
        results_dirname : str, optional
            Name of the results directory (default: 'results')

        """

        # Cache file extensions dictionary
        self._EXTENSIONS = {
            'pickle': 'pkl'
        }

        self.is_valid = True  # Assume the simulation is valid initially
        self.simulation_dir = Path(simulation_dir)
        
        # Validate simulation directory exists
        if not self.simulation_dir.exists():
            print(f"Simulation directory {self.simulation_dir} does not exist.")
            self.is_valid = False
            
        else:
        
            # Auto-detect trajectory directories
            self.traj_dirs = sorted([
                d for d in self.simulation_dir.iterdir() 
                if d.is_dir() and d.name.startswith('traj_')
            ])
        
            if len(self.traj_dirs) == 0:
                print(f"No trajectory directories (traj_*) found in {self.simulation_dir}")
                self.is_valid = False
            
            else:
                
                # Create lattice from first trajectory
                self.lattice = Lattice(self.traj_dirs[0])
                if not self.lattice.is_defined:
                    print(f"Cannot load lattice data for simulation {self.simulation_dir.name}")
                    self.is_valid = False
                else:
                
                    # Initialize trajectory list (filled by load)
                    self.trajectories = []
        
                    # Set up results directory (../../results/ from simulation_dir)
                    self.results_dir = self.simulation_dir.parent.parent / results_dirname
                    self.results_dir.mkdir(exist_ok=True)
        
                    # Load metadata from log file if not provided
                    if metadata is not None:
                        self.metadata = metadata
                    else:
                        self._load_metadata(log_file)
        
    def clear_cache(self, cache=None, verbose=False):
        """
        Clear cached trajectory data files for the specified cache type.
        
        Parameters
        ----------
        cache : str, list of str, or None, default None
            If str, clear cache of specified file format. If None, clear all cache types.
        verbose : bool, default False
            If True, print detailed cache clearing information.
        """

        if cache is None:
            formats_to_clear = self._EXTENSIONS.keys()
        elif type(cache) is list:
            formats_to_clear = cache
        elif type(cache) is str:
            formats_to_clear = [cache]

        for format in formats_to_clear:
            cache_file = self.results_dir / f"{self.metadata['simulation_number']}_trajs.{self._EXTENSIONS[format]}"
            if cache_file.exists():
                cache_file.unlink()
                if verbose:
                    print(f"Cleared trajectory cache: {cache_file.name}")


    def _load_metadata(self, log_file):
        """
        Load simulation metadata from log file.
        
        Parses the log file to extract temperature, coverage, interactions,
        and lattice dimensions for this specific simulation.

        Parameters
        ----------
        log_file : str
            Name of the log file (default: 'jobs.log')
        
        Raises
        ------
        FileNotFoundError
            If log file is not found in the parent directory
        ValueError
            If simulation number is not found in log file
        """
        jobs_log = self.simulation_dir.parent / log_file
        if not jobs_log.exists():
            print(f"log file {jobs_log} not found. Cannot load metadata.")
            self.is_valid = False
            return
        
        # Extract simulation number from directory name (e.g., '1' from 'fn_3leed/jobs/1')
        simulation_number = int(self.simulation_dir.name)
        
        # Parse jobs.log
        with open(jobs_log, 'r') as f:
            header = f.readline().split()  # Read header line
            log_entries = [json.loads(line) for line in f if line.strip()]
        
        # Find entry matching this simulation number
        matching_entry = None
        for entry in log_entries:
            if entry[0] == simulation_number:
                matching_entry = entry
                break
        
        if matching_entry is None:
            print(
                f"Simulation number {simulation_number} not found in {jobs_log}\n"
                f"Available simulation numbers: {[e[0] for e in log_entries]}"
            )
            self.is_valid = False
            return

        # Extract metadata from log entry
        # Format: [simulation_num, job_name, [nx, ny], [n_ads], temp, interaction_info, ...]
        self.metadata = {
            'simulation_number': matching_entry[0],
            'job_name': matching_entry[1],
            'lattice_dimensions': matching_entry[2],  # [nx, ny]
            'n_cells': matching_entry[2][0] * matching_entry[2][1],
            'n_adsorbates': matching_entry[3][0],
            'temperature': matching_entry[4],  # K
            'coverage': matching_entry[3][0] / (matching_entry[2][0] * matching_entry[2][1]),
            'interactions': matching_entry[5][1:]
        }

    def _load_single_trajectory(self, traj_dir):
        """
        Helper function for parallel trajectory loading.
        Parameters
        ----------
        traj_dir : Path
            Directory containing trajectory data
        Returns
        -------
        trajectory
            Trajectory object
        """

        traj = Trajectory(self.lattice, traj_dir)
        traj.load()
        return traj

    def _load_trajectories_parallel(self, workers=None):
        """
        Load multiple trajectories in parallel.
        Parameters
        ----------
        workers : int, optional
            Number of parallel workers. If None, uses all available cores.
        Returns
        -------
        list of trajectories
            Loaded trajectories
        """

        if len(self.traj_dirs) == 0:
            return []
        if workers is None:
            workers = mp.cpu_count()

        with ProcessPoolExecutor(max_workers=workers) as executor:
            trajs = list(executor.map(self._load_single_trajectory, self.traj_dirs))

        return trajs


    def load(self, cache=None, workers=mp.cpu_count(), verbose=False):
        """
        Load trajectory data with caching support.

        Parameters
        ----------
        cache : str or None, default None
            If str, use caching with specified file format. If None, do not use caching.
        workers : int or None, default mp.cpu_count()
            Number of parallel workers to use for loading.
            If None, load sequentially.
        verbose : bool, default False
            If True, print detailed loading information.
        """

        # Determine cache file path and extension
        if cache is not None:
            if cache not in self._EXTENSIONS:
                raise ValueError(f"Unsupported cache format: {cache}. Supported formats: {list(self._EXTENSIONS.keys())}")
            cache_file = self.results_dir / f"{self.metadata['simulation_number']}_trajs.{self._EXTENSIONS[cache]}"

            # Try loading from cache
            if cache_file.exists():
                if verbose: print(f"Loading trajectories from cache: {cache_file.name}")

                if cache=='pickle':
                    with open(cache_file, 'rb') as f:
                        self.trajectories = pickle.load(f)
                # too verbose, drop printing of loaded info
                # if verbose: print(f"Loaded {len(self.trajectories)} cached trajectories")
                return

        # Load trajectories from files
        if verbose:
            print(f"Loading {len(self.traj_dirs)} trajectories...")
            print(f"  Loading mode: {'sequential' if workers is None else 'parallel with ' + str(workers) + ' workers'}")

        if workers is not None:
            # Use parallel loading
            self.trajectories = self._load_trajectories_parallel(workers=workers)
        else:
            # Sequential loading
            self.trajectories = []
            for traj_dir in self.traj_dirs:
                self.trajectories.append(self._load_single_trajectory(traj_dir))
        if verbose:
            print(f"Loaded {len(self.trajectories)} trajectories")
            print(f"  States per trajectory: {len(self.trajectories[0].states)}")
            print(f"  Total states: {sum(len(t.states) for t in self.trajectories)}")

        # Save to cache
        if cache is not None:
            if verbose: print(f"Saving to cache: {cache_file}")

            if cache=='pickle':
                with open(cache_file, 'wb') as f:
                    pickle.dump(self.trajectories, f)

            size_mb = cache_file.stat().st_size / 1024**2
            if verbose: print(f"Cache saved: {size_mb:.1f} MB")

# Warning! Come back to this later
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
    
# Warning! Come back to this later
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
    
    def __repr__(self):
        """String representation of simulation class."""
        if len(self.trajectories) > 0:
            n_states = len(self.trajectories[0].states)
            traj_info = f", {len(self.trajectories)} trajectories ({n_states} states each)"
        else:
            traj_info = f", {len(self.traj_dirs)} trajectories (not loaded)"
        
        return (
            f"simulation(sim={self.metadata['simulation_number']}, "
            f"T={self.metadata['temperature']}K, "
            f"θ={self.metadata['coverage']:.3f} "
        )
    
    def __len__(self):
        """Return number of trajectories."""
        return len(self.trajectories)
