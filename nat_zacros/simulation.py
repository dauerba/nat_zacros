"""
Simulation class for managing Zacros simulation with multiple trajectories.

This module provides a high-level interface for loading, caching, and analyzing
collections of trajectories from a single Zacros simulation.
"""

import numpy as np
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from .lattice import Lattice
from .trajectory import Trajectory

class Simulation:
    """
    Summary
    -------
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
    
    Methods
    -------
    TO DO: write this section

    Examples
    --------
    >>> # Typical workflow
    >>> run = Simulation('fn_3leed/jobs/1')
    >>> run.load()  # Uses cache if available
    >>> run.is_valid  # Check if run data is valid
    >>> r, g, g_std = run.get_ensemble_rdf(r_max=40.0, dr=0.1)
    >>> times, energies, energies_std = run.get_ensemble_energy_vs_time()

    Cache / results helpers (per-simulation)
    ----------------------------------------
    These helpers operate at either a path (static) level or an instance level.

    # Path-level helpers (no Simulation instance required)
    >>> Simulation.clear_traj_cache('fn_3leed/jobs/1', traj_dir_pfx='traj')
    >>> Simulation.clear_results_path('fn_3leed/jobs/1', target=['energy', 'gref'])

    # Instance-level helpers (convenience wrappers)
    >>> sim = Simulation('fn_3leed/jobs/1', metadata={
    ...     'lattice_dimensions':[4,4], 'n_adsorbates':2, 'temperature':300, 'energy_terms':['label','E1']
    ... })
    >>> sim.clear_traj_cache()
    >>> sim.clear_results(target='all')

    Notes
    -----
    - `SimulationSet` delegates cache/result removal to these `Simulation` helpers so
      filesystem logic is centralized and `Simulation` can be used standalone.
    """


    def __init__(self, dir, 
                 metadata=None,
                 lattice_dims=None, n_ads=None, temperature=None, energy_terms=None,
                 traj_dir_pfx='traj',
                 lattice=None):
        """
        Initialize a simulation.
        
        Parameters
        ----------
        dir : str or Path
            Path to simulation directory
            This directory should contain traj_1, traj_2, ... subdirectories
        metadata : dict, optional
            Pre-defined metadata (if available)
        lattice : Lattice, optional
            Shared lattice object (default: None).
        lattice_dims : list, optional
            Lattice dimensions [nx, ny] (default: None)
        n_ads : int, optional
            Number of adsorbates (default: None)
        temperature : float, optional
            Temperature in Kelvin (default: None)
        energy_terms : list, optional
            List of energy terms (default: None)
        traj_dir_pfx : str, optional
            Prefix for trajectory directories (default: 'traj').
             Trajectory directories should be named like 'traj_0', 'traj_1', etc

        """

        self.dir = Path(dir)
        self.is_valid = True  # Assume the simulation is valid initially
        self.traj_dir_pfx = traj_dir_pfx
        self.lattice = lattice

        if metadata is not None:
            self.metadata = metadata
        else:
            self.metadata = {}

        args_ok = True
        args = [lattice_dims, n_ads, temperature, energy_terms]
        keys = ['lattice_dimensions', 'n_adsorbates', 'temperature', 'energy_terms']

        for arg, key in zip(args, keys):
            if arg is not None:
                self.metadata[key] = arg
            if key not in self.metadata:
                print(f" {key.replace('_', ' ').capitalize()} undefined.")
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
                if d.is_dir() and d.name.startswith(self.traj_dir_pfx + '_')
            ])
        
            if len(self.traj_dirs) == 0:
                print(f"No trajectory directories ({self.traj_dir_pfx}_*) found in {self.dir}")
                self.is_valid = False
            
            else:
                # Use passed lattice if available, otherwise create from first trajectory
                if self.lattice is None:
                    self.lattice = Lattice(self.traj_dirs[0])
                
                if not self.lattice.is_defined:
                    print(f"Cannot load lattice data for simulation {self.dir.name}")
                    self.is_valid = False
                else:
                    # Initialize trajectory list
                    self.trajectories = []
        

    def _load_single_trajectory(self, traj_dir, cache=True, verbose=False):
        """
        Helper function for parallel trajectory loading.
        Parameters
        ----------
        traj_dir : Path
            Directory containing trajectory data
        cache : bool, default True
            If True, attempt to load cached trajectory data if available; cache trajectory data if not already cached.
            If False, load from raw simulation data.
        verbose : bool, optional
            If True, print verbose output.
        Returns
        -------
        trajectory
            Trajectory object
        """

        traj = Trajectory(traj_dir, self.lattice)
        traj.load(cache=cache, verbose=verbose)
        return traj

    def _load_trajectories_parallel(self, cache=True, workers=None, verbose=False):
        """
        Load multiple trajectories in parallel.
        Parameters
        ----------
        cache : bool, default True
            If True, attempt to load cached trajectory data if available; cache trajectory data if not already cached.
            If False, load from raw simulation data.
        workers : int, optional
            Number of parallel workers. If None, uses all available cores.
        verbose : bool, optional
            If True, print verbose output.
            Target data to load (default: 'trajs')
        workers : int, optional
            Number of parallel workers. If None, uses all available cores.
        verbose : bool, optional
            If True, print verbose output.
            
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
            trajs = list(executor.map(self._load_single_trajectory, self.traj_dirs, 
                                      [cache]*len(self.traj_dirs), [verbose]*len(self.traj_dirs)))

        return trajs

    def load(self, cache=True, workers=mp.cpu_count(), verbose=False):
        """
        Load trajectory data with caching support.

        Parameters
        ----------
        cache : bool, default True
            If True, load cached trajectory data if available; cache trajectory data if not already cached.
            If False, load from raw simulation data.
        workers : int or None, default mp.cpu_count()
            Number of parallel workers to use for loading.
            If None, load sequentially.
        verbose : bool, default False
            If True, print detailed loading information.
        """

        # Determine cache file path and extension
        print(f"Loading simulation {self.dir.name} with caching {'enabled' if cache else 'disabled'}...")

        # Load trajectories from files
        if verbose:
            print(f"Loading {len(self.traj_dirs)} trajectories...")
            print(f"  Loading mode: {'sequential' if workers is None else 'parallel with ' + str(workers) + ' workers'}")

        if workers is not None:
            # Use parallel loading
            self.trajectories = self._load_trajectories_parallel(cache=cache, workers=workers, verbose=verbose)
        else:
            # Sequential loading
            self.trajectories = []
            for traj_dir in self.traj_dirs:
                self.trajectories.append(self._load_single_trajectory(traj_dir, cache=cache, verbose=verbose))
        if verbose:
            print(f"Loaded {len(self.trajectories)} trajectories")
            print(f"  Total states: {sum(len(t.states) for t in self.trajectories)}")


    # ------------------------------------------------------------------
    # Cache / results helpers at Simulation level
    # ------------------------------------------------------------------
    @staticmethod
    def clear_traj_cache_path(sim_folder, traj_dir_pfx='traj', verbose=False):
        """Remove internal trajectory cache files from a folder.

        This deletes per-trajectory ``traj.pkl`` files inside each ``traj_*`` folder.

        Parameters
        ----------
        sim_folder : str or Path
            Path to the simulation folder.
        traj_dir_pfx : str, optional
            Prefix for trajectory directories (default ``'traj'``).
        verbose : bool, optional
            If True, print each deleted cache file.

        Returns
        -------
        list of Path
            List of deleted cache files.
        """
        sim_folder = Path(sim_folder)
        if not sim_folder.exists():
            return []

        deleted = []
        for d in sim_folder.iterdir():
            if d.is_dir() and d.name.startswith(traj_dir_pfx + '_'):
                cf = d / 'traj.pkl'
                if cf.exists():
                    cf.unlink()
                    deleted.append(cf)
                    if verbose:
                        print(f"Deleted trajectory cache: {cf}")
        return deleted

    def clear_traj_cache(self_or_folder, traj_dir_pfx='traj', verbose=False):
        """Clear trajectory cache files.

        This method is polymorphic so that it can be called in two ways:
        1. **Path-level** (class call) – first argument is simulation folder path.
           Example: Simulation.clear_traj_cache('/path/to/run', traj_dir_pfx='traj')
        2. **Instance-level** – invoked on a Simulation object.
           Example: sim.clear_traj_cache()
        """
        if isinstance(self_or_folder, Simulation):
            sim_folder = self_or_folder.dir
            pfx = self_or_folder.traj_dir_pfx
            v = verbose
        else:
            sim_folder = self_or_folder
            pfx = traj_dir_pfx
            v = verbose

        deleted = Simulation.clear_traj_cache_path(sim_folder, traj_dir_pfx=pfx, verbose=v)
        return len(deleted)

    @staticmethod
    def clear_results_path(sim_folder, target='all', results_files=None, verbose=False):
        """Delete user-visible result files in a simulation folder.

        Parameters
        ----------
        sim_folder : str or Path
            Path to the simulation folder.
        target : str or list
            Which result types to remove. Supported keys: 'energy', 'rdf', 'accessibility', 'gref', 'all'.
        results_files : dict or None
            Mapping from result key to filename (e.g. {'energy':'energy.dat'}). When None, defaults are used.
        verbose : bool
            Print deleted filenames when True.
        """
        default_map = {
            'energy': 'energy.dat',
            'rdf': 'rdf.dat',
            'accessibility': 'accessibility.dat',
            'gref': 'g_ref.dat',
        }
        file_map = dict(default_map if results_files is None else {**results_files, **({'gref': 'g_ref.dat'} if 'gref' not in results_files else {})})

        valid_keys = set(file_map.keys())

        # Normalize target
        if isinstance(target, list):
            formats_to_clear = target
        elif isinstance(target, str):
            formats_to_clear = list(valid_keys) if target == 'all' else [target]
        else:
            raise TypeError("target must be str or list of str")

        formats_to_clear = [f for f in formats_to_clear if f in valid_keys]

        sim_folder = Path(sim_folder)
        count = 0
        for fmt in formats_to_clear:
            cf = sim_folder / file_map[fmt]
            if cf.exists():
                cf.unlink()
                count += 1
                if verbose:
                    print(f"Deleted results file: {cf}")
        return count

    def clear_results(self, target='all', results_files=None, verbose=False, search_dir=None):
        """Instance wrapper for `clear_results_path` that operates on this simulation."""
        d = search_dir if search_dir is not None else self.dir
        return self.clear_results_path(d, target=target, results_files=results_files, verbose=verbose)

    def get_g_ref(self, r_max=None, dr=0.1):
        """
        Calculate reference RDF for full lattice (all sites, coverage=1).

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
            Number of neighbors in each shell (normalized counts)
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


    def get_ensemble_rdf(self, pars_dict=None, file=None):
        """
        Compute ensemble-averaged radial distribution function.
        
        Automatically handles g_ref caching and computation.
        
        Parameters
        ----------
        pars_dict : dict, optional
            Dictionary of parameters:
            - r_max : float, default 40.0
                Maximum distance for RDF (Angstroms)
            - dr : float, default 0.1
                Bin width for RDF (Angstroms)
            - fraction : float, default 1.0
                Fraction of trajectory data to use for RDF calculation (e.g., 0.5 for last half)
            - normalize : bool, default True
                If True, normalize RDF using reference g_ref      
        file : str or Path, optional
            Path to save RDF data (r, g_avg, g_std, g_ref if normalize=True).
            If None, data will not be saved to file (default).  

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

        # Set defaults and update from pars_dict
        pars = {'r_max': 40.0, 'dr': 0.1, 'fraction': 1.0, 'normalize': True}
        if pars_dict is not None:
            pars.update(pars_dict)
        
        r_max =     pars['r_max']
        dr =        pars['dr']
        fraction =  pars['fraction']
        normalize = pars['normalize']

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

        # Save RDF data and g_ref (if present) to per-simulation folder
        if file is not None:
            Path(file).parent.mkdir(parents=True, exist_ok=True)
            if normalize:
                np.savetxt(file, np.column_stack((r, g_avg, g_std, g_ref)),
                    header=f'Parameters: r_max={r_max} dr={dr} fraction={fraction} normalize={normalize}\n' 
                            'r_Angstrom g_r g_std g_ref_r')
            else:
                np.savetxt(file, np.column_stack((r, g_avg, g_std)),
                    header=f'Parameters: r_max={r_max} dr={dr} fraction={fraction} normalize={normalize}\n' 
                            'r_Angstrom g_r g_std')

        if normalize:
            return r, g_avg, g_std, g_ref
        else:
            return r, g_avg, g_std
    

    def get_ensemble_accessibility(self, pars_dict=None, file=None):
        """
        Compute ensemble-averaged site accessibility histogram.
        
        Accessibility measures how many nearest neighbor sites are vacant for occupied sites,
        which affects reactivity and diffusion rates.
        
        Parameters
        ----------
        pars_dict : dict, optional
            Dictionary of parameters:
            - fraction : float, default 1.0
                Fraction of trajectory data to use for accessibility calculation (e.g., 0.5 for last half)
        file : str or Path, optional
            Path to save the accessibility data (accessibility, frequency_avg, frequency_std).
            If None, data will not be saved to file (default).
            
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
        # Set defaults and update from pars_dict
        pars = {'fraction': 1.0}
        if pars_dict is not None:
            pars.update(pars_dict)
        
        fraction = pars['fraction']

        if len(self.trajectories) == 0:
            raise RuntimeError(
                "No trajectories loaded. Call load() first."
            )

        # Compute accessibility for each trajectory
        f13 = []
        f2  = []
        for i, traj in enumerate(self.trajectories):
            freq_13, freq_2 = traj.get_accessibility_histogram(fraction=fraction)
            f13.append(freq_13)
            f2.append(freq_2)
        
        frequencies_13 = np.array(f13)
        frequency_13_avg = np.mean(frequencies_13, axis=0)
        frequency_13_std = np.std(frequencies_13, axis=0)

        frequencies_2 = np.array(f2)
        frequency_2_avg = np.mean(frequencies_2, axis=0)
        frequency_2_std = np.std(frequencies_2, axis=0)

        # Save accessibility data to file
        if file is not None:
            Path(file).parent.mkdir(parents=True, exist_ok=True)
            np.savetxt(file, 
                np.column_stack((frequency_13_avg, frequency_13_std,frequency_2_avg, frequency_2_std)), 
                header=f'Parameters: fraction={fraction}\n'
                       'frequency_13_avg frequency_13_std frequency_2_avg frequency_2_std ',
                fmt='%f %f %f %f')
        
        return frequency_13_avg, frequency_13_std,frequency_2_avg, frequency_2_std
    
    def get_ensemble_energy_vs_time(self, pars_dict=None, file=None):
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
        pars_dict : dict, optional
            Dictionary of parameters:
            - n_bins : int, default 100
                Number of time bins for averaging
        file : str or Path, optional    
            Path to save the energy vs time data (time, energy_avg, energy_std).
            If None, data will not be saved to file (default).
        
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

        # Set defaults and update from pars_dict
        pars = {'n_bins': 100}
        if pars_dict is not None:
            pars.update(pars_dict)
        
        n_bins = pars['n_bins']
        if n_bins <= 0:
            raise ValueError("The value of n_bins must be a positive integer.")

        if len(self.trajectories) == 0:
            raise RuntimeError(
                "No trajectories loaded. Call load() first."
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

        # Save energy vs time data to file
        if file is not None:
            Path(file).parent.mkdir(parents=True, exist_ok=True)
            np.savetxt(file, 
                np.column_stack((time_centers, energy_avg, energy_std)), 
                header=f'Parameters: n_bins={n_bins}\nTime_s Energy_eV Energy_std_eV')

        return time_centers, energy_avg, energy_std
    
    def fraction_eq_is_invalid(self, fraction, property=None, file=None, verbose=False):
        """Check if an equilibrium fraction is invalid.
        """

        if fraction is None or fraction <= 0.0 or fraction > 1.0:
            print(
                f"Skipping {property} for simulation #{self.dir.name}:\n"
                f"  invalid equilibration fraction ({fraction}). "
                 "Run find_equilibrium_fraction_fit() or set simset.fractions_eq[<n>].")
            # Clean cache for invalid fraction
            if file is not None and file.exists():
                file.unlink()
                if verbose:
                    print(f"  {property}-cache cleared for simulation #{self.dir.name} ({file.name})")
            return True
        
        else:
            return False


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
