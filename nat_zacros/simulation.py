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
from .globals import nz_cache_files

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
    trajectories : list of Trajectory
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
        
    # def clear_cache(self, cache=None, verbose=False):
    #     """
    #     Clear cached trajectory data files for the specified cache type.
        
    #     Parameters
    #     ----------
    #     cache : str, list of str, or None, default None
    #         If str, clear cache of specified file format. If None, clear all cache types.
    #     verbose : bool, default False
    #         If True, print detailed cache clearing information.
    #     """

    #     if cache is None:
    #         formats_to_clear = self._EXTENSIONS.keys()
    #     elif type(cache) is list:
    #         formats_to_clear = cache
    #     elif type(cache) is str:
    #         formats_to_clear = [cache]

    #     for format in formats_to_clear:
    #         cache_file = self.results_dir / f"{self.metadata['simulation_number']}_trajs.{self._EXTENSIONS[format]}"
    #         if cache_file.exists():
    #             cache_file.unlink()
    #             if verbose:
    #                 print(f"Cleared trajectory cache: {cache_file.name}")


    # -----------------------------------------------------------------------------
    # ------            Consider whether this function is needed             ------
    # ------   will we create a simulation outside of SimulationSet ?        ------
    # -----------------------------------------------------------------------------
    

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

    # WARNING: need to revise or eliminate.
    # This function is currently called from SimulationSet do do loading
    # consider moving all loading to SimulationSet and eliminating this function from Simulation class
    #,will we create Simulation objects from SimulationSet?

    def load(self, target='trajs', workers=mp.cpu_count(), verbose=False):
        """
        Load trajectory data with caching support.

        Parameters
        ----------
        target : str or None, default 'trajs'
            Specifies the type of data to load. Currently only 'trajs' is supported.
        workers : int or None, default mp.cpu_count()
            Number of parallel workers to use for loading.
            If None, load sequentially.
        verbose : bool, default False
            If True, print detailed loading information.
        """

        # Determine cache file path and extension
        if target is not None:
            if target not in nz_cache_files:
                raise ValueError(f"Unsupported cache target: {target}. Supported targets: {list(nz_cache_files.keys())}")
            
            suffix, ext = nz_cache_files[target]
            cache_file = self.results_dir / f"{self.metadata['simulation_number']}{suffix}.{ext}"

            # Try loading from cache
            if cache_file.exists():
                if verbose: print(f"Loading {target} from cache: {cache_file.name}")

                if ext == 'pkl':
                    with open(cache_file, 'rb') as f:
                        self.trajectories = pickle.load(f)
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
        if target is not None:
            suffix, ext = nz_cache_files[target]
            if verbose: print(f"Saving to cache: {cache_file}")

            if ext == 'pkl':
                with open(cache_file, 'wb') as f:
                    pickle.dump(self.trajectories, f)

            size_mb = cache_file.stat().st_size / 1024**2
            if verbose: print(f"Cache saved: {size_mb:.1f} MB")


    def get_g_ref(self, r_max=None, dr=0.1):
        """
        Calculate reference RDF for full lattice (all sites, coverage=1).
        
        This computes the number of neighbors in each distance shell,
        used to normalize the RDF such that g(r)=1 for ideal gas.
        
        Parameters
        ----------
        r_max : float, optional
            Maximum distance for RDF
        dr : float, default 0.1
            Bin width in Angstroms
            
        Returns
        -------
        r_bins : ndarray
            Bin centers
        g_ref : ndarray
            Number of neighbors in each shell (integer counts)
        """
        if r_max is None:
            v1 = self.lattice.cell_vectors[0]
            v2 = self.lattice.cell_vectors[1]
            l1 = np.linalg.norm(v1)
            l2 = np.linalg.norm(v2)
            l3 = np.linalg.norm(v1 + v2)
            r_max = min(l1, l2, l3) / 2.0
        
        # Initialize histogram
        n_bins = int(np.ceil(r_max / dr))
        bin_edges = np.linspace(0.0, r_max, n_bins + 1)
        r_bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        # Get all lattice site coordinates
        all_coords = self.lattice.coordinates
        n_sites = len(all_coords)
        counts = np.zeros(n_bins, dtype=int)

        # Vectorized calculation using pairwise_distances_pbc
        dists_matrix = self.lattice.pairwise_distances_pbc(all_coords)
        mask = np.triu(np.ones(dists_matrix.shape, dtype=bool), k=1)
        dists = dists_matrix[mask]
        valid_dists = dists[(dists > 0) & (dists <= r_max)]
        counts, _ = np.histogram(valid_dists, bins=bin_edges)

        # Normalize: 2 * counts / n_sites (factor 2 for unordered pairs)
        g_ref = 2.0 * counts / n_sites
        return r_bins, g_ref


    def get_ensemble_rdf(self, r_max=40.0, dr=0.1, fraction=1.0, normalize=True):
        """
        Compute ensemble-averaged radial distribution function.
        
        Automatically handles g_ref caching and computation.
        
        Parameters
        ----------
        r_max : float, default 40.0
            Maximum distance for RDF (Angstroms)
        dr : float, default 0.1
            Bin width for RDF (Angstroms)
        fraction : float, default 1.0
            Fraction of trajectory data to use for RDF calculation (e.g., 0.5 for last half)
        normalize : bool, default True
            If True, normalize RDF using reference g_ref        
        Returns
        -------
        r : ndarray
            Bin centers (Angstroms)
        g_avg : ndarray
            Ensemble-averaged RDF
        g_std : ndarray
            Standard deviation of RDF across trajectories
        g_ref : ndarray, optional
            Reference RDF (only returned if normalize=True)
        
        Raises
        ------
        RuntimeError
            If trajectories have not been loaded yet
        """

        if len(self.trajectories) == 0:
            raise RuntimeError(
                "No trajectories loaded. Call load_trajectories() first."
            )
        
        if normalize:
            r_ref, g_ref = self.get_g_ref(r_max=r_max, dr=dr)
        else:
            r_ref, g_ref = None, None
        
        # Compute RDF for each trajectory
        rdfs = []
        for i, traj in enumerate(self.trajectories):
            r, g = traj.get_rdf(r_max=r_max, dr=dr, fraction=fraction, g_ref=g_ref)
            rdfs.append(g)
        
        # Ensemble average
        g_avg = np.mean(rdfs, axis=0)
        g_std = np.std(rdfs, axis=0)

        if normalize:
            return r, g_avg, g_std, g_ref
        else:
            return r, g_avg, g_std
    
    def get_ensemble_accessibility(self, fraction=1.0):
        """
        Compute ensemble-averaged site accessibility histogram.
        
        Accessibility measures how many nearest neighbor sites are vacant for occupied sites,
        which affects reactivity and diffusion rates.
        
        Parameters
        ----------
        fraction : float, default 1.0
            Fraction of trajectory data to use for accessibility calculation (e.g., 0.5 for last half)
            
        Returns
        -------
        accessibility : ndarray
            Number of vacant nearest neighbors (0 to max_coordination)
        frequency_avg : ndarray
            Ensemble-averaged frequency distribution
        frequency_std : ndarray
            Standard deviation of frequency distribution across trajectories
        
        Raises
        ------
        RuntimeError
            If trajectories have not been loaded yet
        """
        if len(self.trajectories) == 0:
            raise RuntimeError(
                "No trajectories loaded. Call load_trajectories() first."
            )
        
        # Compute accessibility for each trajectory
        accessibilities = []
        for i, traj in enumerate(self.trajectories):
            acc, freq = traj.get_accessibility_histogram()
            accessibilities.append(freq)
        
        # Get max coordination from first trajectory
        max_coord = np.max(self.trajectories[0].lattice.site_coordinations)
        accessibility = np.arange(max_coord + 1)
        
        # Pad all frequencies to same length and compute ensemble statistics
        frequencies_padded = []
        for freq in accessibilities:
            padded = np.zeros(max_coord + 1)
            padded[:len(freq)] = freq
            frequencies_padded.append(padded)
        
        frequencies_array = np.array(frequencies_padded)
        frequency_avg = np.mean(frequencies_array, axis=0)
        frequency_std = np.std(frequencies_array, axis=0)
        
        return accessibility, frequency_avg, frequency_std
    
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
