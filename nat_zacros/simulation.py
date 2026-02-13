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
        Simulation metadata (temperature, coverage, energy_terms, etc.)
    results_dir : Path
        Directory for storing cache and results files
    dir : Path
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


    def __init__(self, dir, 
                 metadata=None,
                 lattice_dims=None, n_ads=None, temperature=None, energy_terms=None):
        """
        Initialize a simulation.
        
        Parameters
        ----------
        dir : str or Path
            Path to simulation directory
            This directory should contain traj_1, traj_2, ... subdirectories
        metadata : dict, optional
            Pre-defined metadata (if available)
        name : str, optional
            Name of the simulation (default: None)
        lattice_dims : list, optional
            Lattice dimensions [nx, ny] (default: None)
        n_ads : int, optional
            Number of adsorbates (default: None)
        temperature : float, optional
            Temperature in Kelvin (default: None)
        energy_terms : list, optional
            List of energy terms (default: None)
        """

        self.dir = Path(dir)
        self.is_valid = True  # Assume the simulation is valid initially
        
        if metadata is not None:
            self.lattice_dims = metadata['lattice_dimensions']
            self.n_ads = metadata['n_adsorbates']
            self.temperature = metadata['temperature']
            self.energy_terms = metadata['energy_terms']

        args_ok = True
        if lattice_dims is not None:
            self.lattice_dims = lattice_dims
        if self.lattice_dims is None:
            print(' Lattice dimensions undefined. Use lattice_dim argument when initializing.')
            args_ok = False
        if n_ads is not None:
            self.n_ads = n_ads
        if self.lattice_dims is None:
            print(' Number of adsorbates undefined. Use n_ads argument when initializing.')
            args_ok = False
        if temperature is not None:
            self.temperature = temperature
        if self.temperature is None:
            print(' Temperature undefined. Use temperature argument when initializing.')
            args_ok = False
        if energy_terms is not None:
            self.energy_terms = energy_terms
        if self.energy_terms is None:
            print(' Energy terms undefined. Use energy_terms argument when initializing.')
            args_ok = False

        if not args_ok:
            raise Exception("Stopping execution")

        # Validate simulation directory exists
        if not self.dir.exists():
            print(f"Simulation directory {self.dir} does not exist.")
            self.is_valid = False
            
        else:
        
            # Auto-detect trajectory directories
            self.traj_dirs = sorted([
                d for d in self.dir.iterdir() 
                if d.is_dir() and d.name.startswith('traj_')
            ])
        
            if len(self.traj_dirs) == 0:
                print(f"No trajectory directories (traj_*) found in {self.dir}")
                self.is_valid = False
            
            else:
                
                # Create lattice from first trajectory
                self.lattice = Lattice(self.traj_dirs[0])
                if not self.lattice.is_defined:
                    print(f"Cannot load lattice data for simulation {self.dir.name}")
                    self.is_valid = False
                else:
                
                    # Initialize trajectory list (filled by load)
                    self.trajectories = []
        

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
        print(f"Loading simulation {self.metadata['simulation_number']} with target '{target}'")
        if target is not None:
            # if target not in nz_cache_files:
            #     raise ValueError(f"Unsupported cache target: {target}. Supported targets: {list(nz_cache_files.keys())}")
            
            # suffix, ext = nz_cache_files[target]
            cache_file = self.results_dir / f"{self.metadata['simulation_number']}{target}"

            # Try loading from cache
            if cache_file.exists():
                if verbose: print(f"Loading {target} from cache: {cache_file.name}")

                if cache_file.suffix == '.pkl':
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
            # suffix, ext = nz_cache_files[target]
            if verbose: print(f"Saving to cache: {cache_file}")

            if cache_file.suffix == '.pkl':
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
