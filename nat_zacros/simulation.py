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


    def __init__(self, dir, metadata, traj_dir_pfx='traj', results_dir=None, lattice=None):
        """
        Initialize a simulation.
        
        Parameters
        ----------
        dir : str or Path
            Path to simulation directory
            This directory should contain traj_1, traj_2, ... subdirectories
        metadata : dict
            metadata defining conditions of simulations;
            dict must include at least 
                'lattice_dimensions', 'n_adsorbates', 'temperature', 'energy_terms' keys
        traj_dir_pfx : str, optional
            Prefix for trajectory directories (default: 'traj').
             Trajectory directories should be named like 'traj_0', 'traj_1', etc
        lattice : Lattice, optional
            Shared lattice object (default: None).
        """

        self.dir = Path(dir)
        self.results_dir = results_dir if results_dir is not None else self.dir
        self.is_valid = True  # Assume the simulation is valid initially
        self.is_loaded = False
        self.traj_dir_pfx = traj_dir_pfx
        self.lattice = lattice
        self.trajectories = dict()

        self.metadata = metadata

        self._necessary_keys = [
            'lattice_dimensions', 
            'surf_species_names',
            'n_adsorbates', 
            'temperature', 
            'energy_terms'
            ]

        args_ok = True
        for key in self._necessary_keys:
            if key not in self.metadata:
                print(f" {key} undefined.")
                args_ok = False

        if not args_ok:
            self.is_valid = False
            print('Not enough metadata for a valid simulation.')

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
                    for tdir in self.traj_dirs:
                        self.trajectories[tdir.name] = Trajectory(tdir, self.lattice, self.metadata)
        

#
# Commented are the functions for parallelizing the loading
# We consider to throw them away
#
    # def _load_single_trajectory(self, traj_dir, cache=True, verbose=False):
    #     """
    #     Helper function for parallel trajectory loading.
    #     Parameters
    #     ----------
    #     traj_dir : Path
    #         Directory containing trajectory data
    #     cache : bool, default True
    #         If True, attempt to load cached trajectory data if available; cache trajectory data if not already cached.
    #         If False, load from raw simulation data.
    #     verbose : bool, optional
    #         If True, print verbose output.
    #     Returns
    #     -------
    #     trajectory
    #         Trajectory object
    #     """

    #     traj = Trajectory(traj_dir, self.lattice, self.metadata)
    #     traj.load(cache=cache, verbose=verbose)
    #     return traj

    # def _load_trajectories_parallel(self, cache=True, workers=None, verbose=False):
    #     """
    #     Load multiple trajectories in parallel.
    #     Parameters
    #     ----------
    #     cache : bool, default True
    #         If True, attempt to load cached trajectory data if available; cache trajectory data if not already cached.
    #         If False, load from raw simulation data.
    #     workers : int, optional
    #         Number of parallel workers. If None, uses all available cores.
    #     verbose : bool, optional
    #         If True, print verbose output.
    #         Target data to load (default: 'trajs')
    #     workers : int, optional
    #         Number of parallel workers. If None, uses all available cores.
    #     verbose : bool, optional
    #         If True, print verbose output.
            
    #     Returns
    #     -------
    #     list of trajectories
    #         Loaded trajectories
    #     """

    #     if len(self.traj_dirs) == 0:
    #         return []
    #     if workers is None:
    #         workers = mp.cpu_count()

    #     with ProcessPoolExecutor(max_workers=workers) as executor:
    #         trajs = list(executor.map(self._load_single_trajectory, self.traj_dirs, 
    #                                   [cache]*len(self.traj_dirs), [verbose]*len(self.traj_dirs)))

    #     return trajs

    # def load_mp(self, cache=True, workers=mp.cpu_count(), verbose=False):
    #     """
    #     Load trajectory data with caching support.

    #     Parameters
    #     ----------
    #     cache : bool, default True
    #         If True, load cached trajectory data if available; cache trajectory data if not already cached.
    #         If False, load from raw simulation data.
    #     workers : int or None, default mp.cpu_count()
    #         Number of parallel workers to use for loading.
    #         If None, load sequentially.
    #     verbose : bool, default False
    #         If True, print detailed loading information.
    #     """

    #     # Determine cache file path and extension
    #     print(f"Loading simulation {self.dir.name} with caching {'enabled' if cache else 'disabled'}...")

    #     # Load trajectories from files
    #     if verbose:
    #         print(f"Loading {len(self.traj_dirs)} trajectories...")
    #         print(f"  Loading mode: {'sequential' if workers is None else 'parallel with ' + str(workers) + ' workers'}")

    #     if workers is not None:
    #         # Use parallel loading
    #         self.trajectories = self._load_trajectories_parallel(cache=cache, workers=workers, verbose=verbose)
    #     else:
    #         # Sequential loading
    #         self.trajectories = []
    #         for traj_dir in self.traj_dirs:
    #             self.trajectories.append(self._load_single_trajectory(traj_dir, cache=cache, verbose=verbose))
    #     if verbose:
    #         print(f"Loaded {len(self.trajectories)} trajectories")
    #         print(f"  Total states: {sum(len(t.states) for t in self.trajectories)}")

    def _parse_results_header(self, file_path, verbose=False):
        """
        Parse parameters from the first line of a results file.

        Format expected: # Parameters: key1=val1 key2=val2 ...

        Parameters
        ----------
        file_path : str or Path
            Path to the results file.
        verbose : bool, optional
            If True, prints warnings on failure.

        Returns
        -------
        dict
            Dictionary of parsed parameter keys and values.
        """
        if file_path is None or not Path(file_path).exists():
            return {}
        try:
            with open(file_path, 'r') as f:
                line = f.readline().strip()
                if not line.startswith('# Parameters:'):
                    return {}
                
                # Parseheader: remove label, convert '=' to spaces, then tokenize by whitespace
                content = line.replace('# Parameters:', '').replace('=', ' ')
                parts = content.split()
                
                params = {}
                for i in range(0, len(parts), 2):
                    if i + 1 >= len(parts):
                        break
                    key = parts[i].strip()
                    val_str = parts[i+1].strip()
                    
                    # Type conversion
                    if val_str.lower() == 'true': val = True
                    elif val_str.lower() == 'false': val = False
                    else:
                        try:
                            val = float(val_str) if '.' in val_str or 'e' in val_str.lower() else int(val_str)
                        except ValueError:
                            val = val_str
                    params[key] = val
                return params
        except Exception as e:
            if verbose:
                print(f"Warning: Failed to parse cache header in {file_path}: {e}")
            return {}


    def load(self, cache=True, verbose=False):
        """
        Load trajectory data with caching support.

        Parameters
        ----------
        cache : bool, default True
            If True, load cached trajectory data if available; cache trajectory data if not already cached.
            If False, load from raw simulation data.
        verbose : bool, default False
            If True, print detailed loading information.
        """

        # Determine cache file path and extension
        print(f"Loading simulation {self.dir.name} with caching {'enabled' if cache else 'disabled'}...")

        # Load trajectories from files
        if verbose:
            print(f"Loading {len(self.traj_dirs)} trajectories...")

        for traj in self.trajectories.values():
            traj.load(cache=cache, verbose=verbose)
        self.is_loaded = True

        if verbose:
            print(f"Loaded {len(self.trajectories)} trajectories")
            print(f"  Total states: {sum(len(t.states) for t in self.trajectories.values())}")

    def unload(self, verbose=False):
        """
        Unload trajectory data.

        Parameters
        ----------
        verbose : bool, default False
            If True, print detailed loading information.
        """

        # Unload trajectories from files
        if verbose:
            print(f"Unloading simulation {self.dir.name} from memory...")

        for traj in self.trajectories.values():
            traj.unload(verbose=verbose)
        self.is_loaded = False

        if verbose:
            print(f"Unloaded {len(self.trajectories)} trajectories")
            print(f"  Total states: {sum(len(t.states) for t in self.trajectories.values())}")



    def clear_traj_cache(self, verbose=False):
        """Clear trajectory cache files.

        """
        # Clearing cache from files
        if verbose:
            print(f"Clearing trajectory cache for simulation {self.dir.name}...")

        for traj in self.trajectories.values():
            traj.clear_cache(verbose=verbose)

        return len(self.trajectories)  # Return number of trajectories cleared
    
    def clear_results(self, fnames, verbose=False):
        """Clear cached results file.

        Parameters
        ----------
        fnames : str or Path or list of them
            Path to the results file to clear.
        verbose : bool, default False
            If True, print detailed information about cleared results.
        """

        # Normalize fnames
        if isinstance(fnames, list):
            files_to_clear = fnames
        else:
            files_to_clear = [fnames]

        counter = 0
        for file in fnames:
            file = Path(file)
            if file.exists():
                file.unlink()
                counter += 1
                if verbose:
                    print(f"Simulation {self.dir.name}: {file.name} deleted.")
            else:
                if verbose:
                    print(f"Simulation {self.dir.name}: {file.name} does not exist.")
    
        return counter
        


    def get_g_ref(self, distances, r_max=None, dr=0.1):
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

        counts = np.zeros(n_bins, dtype=int)

        valid_dists = distances[(distances > 0) & (distances <= r_max)]
        counts, _ = np.histogram(valid_dists, bins=bin_edges)
        g_ref = counts / len(self.lattice)

        return r_bins, g_ref


    def get_ensemble_rdf(self, pars_dict=None, file=None, verbose=False):
        """
        Compute ensemble-averaged radial distribution function.
        
        Automatically handles g_ref caching and computation.
        
        Parameters
        ----------
        pars_dict : dict, optional
            Dictionary of parameters:
            - species_1, species_2 : str
                Species for which rdf is to be calculated. 
                Default: species with index 0
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
        pars = {'species_1': self.metadata['surf_species_names'][0], 
                'species_2': self.metadata['surf_species_names'][0], 
                'r_max': 40.0, 'dr': 0.1, 'fraction': 1.0, 'normalize': True}
        if pars_dict is not None:
            pars.update(pars_dict)
        
        species_1 = pars['species_1']
        species_2 = pars['species_2']
        r_max =     pars['r_max']
        dr =        pars['dr']
        fraction =  pars['fraction']
        normalize = pars['normalize']

        if len(self.trajectories) == 0:
            raise RuntimeError(
                "No trajectories loaded. Call load_trajectories() first."
            )

        # Loading rdf from the results file
        if file is not None:
            file_sp = file.parent / f'{file.stem}_{str(species_1).replace('*','star')}-{str(species_1).replace('*','star')}{file.suffix}'
            data = self.load_results(file_sp, pars, verbose=verbose)
            if data is not None:
                return data

        # Get lattice sites distances look-up table
        distances = self.lattice.pairwise_distances_pbc(condensed=False)

        if normalize:
            r_ref, g_ref = self.get_g_ref(distances, r_max=r_max, dr=dr)
        else:
            r_ref, g_ref = None, None
        
        # Compute RDF for each trajectory
        rdfs = []
        for traj in self.trajectories.values():
            r, g = traj.get_rdf(species_1, species_2, distances, r_max=r_max, dr=dr, fraction=fraction, g_ref=g_ref)
            rdfs.append(g)
        
        # Ensemble average
        g_avg = np.mean(rdfs, axis=0)

        # Save RDF data and g_ref (if present) to per-simulation folder
        if file is not None:
            Path(file).parent.mkdir(parents=True, exist_ok=True)
            if normalize:
                np.savetxt(file_sp, np.column_stack((r, g_avg, g_ref)),
                    header= f'Parameters: '
                            f'species_1={species_1} species_2={species_2} '
                            f'r_max={r_max} dr={dr} fraction={fraction} normalize={normalize}\n' 
                            'r_Angstrom g_r g_ref_r')
            else:
                np.savetxt(file_sp, np.column_stack((r, g_avg)),
                    header= f'Parameters: '
                            f'species_1={species_1} species_2={species_2} '
                            f'r_max={r_max} dr={dr} fraction={fraction} normalize={normalize}\n' 
                            'r_Angstrom g_r')

        if normalize:
            return r, g_avg, g_ref
        else:
            return r, g_avg
    

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
        for traj in self.trajectories.values():
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
    

    def get_ensemble_energy_vs_time(self, pars_dict=None, file=None, verbose=False):
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
        end_time = min([traj.times[-1] for traj in self.trajectories.values()])
        start_time = max([traj.times[0] for traj in self.trajectories.values()])
        
        # Create time bins for discretization
        time_bins = np.linspace(start_time, end_time, n_bins + 1)
        time_centers = 0.5 * (time_bins[:-1] + time_bins[1:])
        
        # STAGE 1: Intra-trajectory averaging
        # For each trajectory, bin its energy measurements and average within bins
        energy_hists = []
        for traj in self.trajectories.values():
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
    

    def get_ensemble_leed_intensity_vs_time(self, pars_dict=None, file=None):
        """
        Compute ensemble-averaged leed intensity as function of time.
        
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
        intensity_avg : ndarray
            Ensemble-averaged leed intensity at each time bin
        intensity_std : ndarray
            Standard deviation of leed intensity across trajectories
        
        Raises
        ------
        RuntimeError
            If trajectories have not been loaded yet
        
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
        end_time = min([traj.times[-1] for traj in self.trajectories.values()])
        start_time = max([traj.times[0] for traj in self.trajectories.values()])
        
        # Create time bins for discretization
        time_bins = np.linspace(start_time, end_time, n_bins + 1)
        time_centers = 0.5 * (time_bins[:-1] + time_bins[1:])
        
        # STAGE 1: Intra-trajectory averaging
        # For each trajectory, bin its measurements and average within bins
        hists = []
        for traj in self.trajectories.values():
            times, intensities = traj.get_leed_intensity_vs_time()
            
            # Initialize binned energy and sample counts for this trajectory
            hist = np.zeros(n_bins)
            counts = np.zeros(n_bins)
            
            # Accumulate energy measurements into bins
            for t, intensity in zip(times, intensities):
                if start_time <= t <= end_time:
                    bin_idx = np.digitize(t, time_bins, right=False) - 1
                    if 0 <= bin_idx < n_bins:
                        hist[bin_idx] += intensity
                        counts[bin_idx] += 1  # Track number of samples per bin
            
            # Average within each bin (multiple measurements → single value per bin)
            # This gives us one representative value per time bin for THIS trajectory
            with np.errstate(divide='ignore', invalid='ignore'):
                hist = np.where(counts > 0, hist / counts, np.nan)
            
            hists.append(hist)
        
        # STAGE 2: Inter-trajectory (ensemble) averaging
        # Each trajectory now contributes exactly one value per time bin
        # Average across trajectories with equal weighting
        avgs = np.nanmean(hists, axis=0)
        stds = np.nanstd(hists, axis=0)  # Trajectory-to-trajectory variation

        # Save energy vs time data to file
        if file is not None:
            Path(file).parent.mkdir(parents=True, exist_ok=True)
            np.savetxt(file, 
                np.column_stack((time_centers, avgs, stds)), 
                header=f'Parameters: n_bins={n_bins}\nTime_s leed_int lead_std')

        return time_centers, avgs, stds
    

    def get_ensemble_clusters(self, pars_dict=None, file=None):
        """
        Compute ensemble-averaged cluster properties.
        
        Parameters
        ----------
        pars_dict : dict, optional
            Dictionary of parameters:
            - cutoff : float, default the distance to the 3rd-neighbor shell
                Distance cutoff for connectivity
            - eps : float, default 1e-4
                Tolerance parameter for the approximate search in the tree (see docs for scipy.spatial.ckdtree package)
            - fraction : float, default 1.0
                Fraction of trajectory data to use for RDF calculation (e.g., 0.5 for last half)
            - method : str, default 'ckdtree'
                Method to get clusters 
        file : str or Path, optional
            Path to save RDF data (r, g_avg, g_std, g_ref if normalize=True).
            If None, data will not be saved to file (default).  

        Returns
        -------
        freqs_avg : ndarray
            Ensemble-averaged mean cluster size frequencies
        freqs_std : ndarray
            Sample Standard deviation of cluster size frequencies across trajectories
        
        Raises
        ------
        RuntimeError
            If trajectories have not been loaded yet
        """

        # Set defaults and update from pars_dict
        # If cutoff not provided, set to 3rd nn distance for the lattice
        pars = {'cutoff': self.lattice.get_nn_distance(order=3),
                'eps'   : 1e-4, 
                'fraction': 1.0,
                'method':   'ckdtree'}
        if pars_dict is not None:
            pars.update(pars_dict)
        
        cutoff   =  pars['cutoff']
        eps      =  pars['eps']
        fraction =  pars['fraction']
        method   =  pars['method']

        n_trajs = len(self)        
        if n_trajs == 0:
            raise RuntimeError(
                "No trajectories loaded. Call load_trajectories() first."
            )
        
        
        # Compute cluster size frequencies for each trajectory
        freqs_avg = np.zeros(len(self.lattice))
        freqs_std = np.zeros(len(self.lattice))
        for traj in self.trajectories.values():
            freqs = traj.get_cluster_size_freqs(cutoff=cutoff, eps=eps, fraction=fraction, method=method)
            freqs_avg += freqs
            freqs_std += freqs*freqs

        # Gete average
        freqs_avg /= n_trajs
        
        # Get sample standard deviation (note Bessel's correction in denominator)
        freqs_std = np.sqrt((freqs_std - n_trajs*freqs_avg**2) / (n_trajs - 1))

        # Trim zeros from back side
        freqs_avg = np.trim_zeros(freqs_avg,'b')
        freqs_std = freqs_std[:len(freqs_avg)]

        # Save cluster size frequencies data
        if file is not None:
            Path(file).parent.mkdir(parents=True, exist_ok=True)
            np.savetxt(file, np.column_stack((freqs_avg, freqs_std)),
                header=f'Parameters: cutoff={cutoff} eps={eps} fraction={fraction} method={method}\n' 
                        'freqs_avg freqs_std')

        return freqs_avg, freqs_std
    

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


    def load_results(self, file, target_params, verbose=False):

        data = None

        # Check results file validity
        res_params = self._parse_results_header(file)
        res_valid = bool(res_params) and all(
            np.isclose(res_params.get(k, -1e9), v) if isinstance(v, float) else res_params.get(k) == v
            for k, v in target_params.items()
        )

        if res_valid:
            try:
                # Skiprows=1 if it doesn't have an extra header line beyond Parameters
                # Actually RDF and Energy have multiple header lines. np.loadtxt handles # comments automatically.
                if verbose:
                    print(f" Loading from {file.name}...",end='')
                data = np.loadtxt(file, unpack=True)
                if verbose:
                    print()
            except Exception:
                if verbose:
                    print(f"failed.")
                

        if verbose and data is None:
            print(f' Calculating and saving to {file.name}...')

        return data

    def __repr__(self):
        """String representation of simulation class."""

        n_trajs = len(self)
        if n_trajs > 0:
            n_states = np.sum([len(tr) for tr in self.trajectories.values()])
        else:
            n_states = 0

        traj_info = f", {n_trajs} trajectories with {n_states} states in total)"
        
        return (
            f"simulation(sim={self.dir.name}, "
            f"T={self.metadata['temperature']}K, "
            f"θ={self.metadata['coverage']:.3f})"
            + traj_info
        )
    
    def __len__(self):
        """Return number of trajectories."""
        return len(self.trajectories)
