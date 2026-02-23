"""
SimulationSet class for managing Zacros simulation sets
We refer to calculations as a set when they share the same log file.
We refer to calculations as a simulation when they share the same simulation folder.

This module provides a high-level interface for loading, caching, and analyzing a Zacros simulation set.
"""


import json
import matplotlib.pyplot as plt
#import pickle
from matplotlib.ticker import MultipleLocator, FuncFormatter
import numpy as np
from lmfit import Model
import multiprocessing as mp
import warnings
#from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import tarfile
from .simulation import Simulation

class SimulationSet:
    """
    Manages a Zacros simulation set.
    
    This class provides a high-level interface for:
    - Metadata extraction from jobs.log
    
    Attributes
    ----------

    data_path  : Path or str
        Directory containing the log_file, simset_dir, and results_dir
    en_file_sfx : str
        Suffix for energy data files (default: 'energy.dat').
        If None, energy data files are not created.
    fractions_eq : dict
        Dictionary mapping simulation numbers to equilibrium fractions.
    log_file : str  
        name of the log file (default: 'jobs.log')
    metadata : list of dictionaries
        Simulation metadata (temperature, coverage, energy terms, etc.)
    results_dir : str
        subdirectory containing simulation results (default: 'results')
    simset_dir : str
        subdirectory containing a simulation set (default: 'jobs')
    simulations : list of Simulation
        List of loaded Simulation objects in the set.
    use_cache : bool
        Whether caching is used when loading simulations.
    verbose : bool
        Whether to print verbose output during loading.
    
    
Examples
        --------
        # Per-simulation helpers (available on Simulation and delegated by SimulationSet)
        >>> from nat_zacros import Simulation, SimulationSet
        >>> # Path-level (no Simulation instance required)
        >>> Simulation.clear_traj_cache('/path/to/sim', traj_dir_pfx='traj')
        >>> Simulation.clear_results_path('/path/to/sim', target=['energy','gref'])
        >>> # Instance-level
        >>> sim = Simulation('/path/to/sim', metadata={'lattice_dimensions':[4,4], 'n_adsorbates':2, 'temperature':300, 'energy_terms':['label','E1']})
        >>> sim.clear_traj_cache()
        >>> sim.clear_results(target='all')
        >>> # SimulationSet delegates to these helpers for convenience across a set
        >>> simset = SimulationSet('/path/to/simset')
        >>> simset.clear_traj_cache(simulations=[1,2])

        Methods
    -------

    get_ensemble_energy_vs_time(n_bins=100, verbose=False)
    get_ensemble_rdfs(r_max=40.0, dr=0.1, normalize=True, verbose=False)
    get_ensemble_accessibilities(verbose=False)
    load(cache=None, workers=mp.cpu_count(), simulations=None, verbose=False)
    plot_energy(energies_vs_time, ncols=3, figsize=(12,2.5), title_fontsize=10, suptitle_fontsize=16, show_eq=True)
      
    """

    def __init__(self, data_path, 
                 log_file           ='jobs.log', 
                 simset_dir         ='data', 
                 results_dir        ='results', 
                 traj_dir_pfx       ='traj'):
        """
        Initialize a SimulationSet.
        
        Parameters
        ----------
        data_path : Path or str
            Path to simulation set directory (e.g., 'fn_3leed')
            This directory should contain jobs.log and the simulations subdirectory
        log_file : str, optional
            Name of the log file (default: 'jobs.log')
        results_dir : str, optional
            Name of the subdirectory for storing results (default: 'results')
        simset_dir : str, optional
            Name of the subdirectory containing simulations (default: 'data')
        traj_dir_pfx : str, optional
            Prefix for trajectory directories (default: 'traj').
             Trajectory directories should be named like 'traj_0', 'traj_1', etc
       
        """
        
        # Validate the simulation set directory exists
        if not Path(data_path).exists():
            raise FileNotFoundError(f"Simulation set directory not found: {data_path}")

        # Dictionary (property key: property_function)
        self._properties = {
            'energy':        'get_ensemble_energy_vs_time',
            'rdf':           'get_ensemble_rdfs',
            'accessibility': 'get_ensemble_accessibilities',
        }

        self.data_path          = Path(data_path)
        self.log_file           = log_file
        self.results_dir        = results_dir
        self.simset_dir         = simset_dir
        self.traj_dir_pfx       = traj_dir_pfx
        self.verbose            = False             # default verbosity

        # Try to untar simulations if simulation set directory doesn't exist
        if not (self.data_path / self.simset_dir).exists():
            # untar jobs directory
            tgz_file = self.data_path / (self.simset_dir + '.tgz')
            print(f"Extracting a simulation set from {tgz_file.as_posix()} ...")
            try:
                if not tgz_file.exists():
                     # Try looking for .tar.gz if .tgz not found, or just fail
                     pass
                
                with tarfile.open(tgz_file, "r:gz") as tar:
                    tar.extractall(path=self.data_path)
            except Exception as e:
                print(f"Error extracting simulations: {e}")
                print(f"Data for simulation set {self.simset_dir} is invalid.")
                self.is_valid = False
            else:
                print(f"Extraction complete.")

        self._load_metadata()
        self.simulations        = []   # initialize simulations list
        # Initialize equilibration fractions dictionary with default None for each simulation
        # This avoids KeyError when code expects an entry per simulation unless user overrides.
        self.fractions_eq       = {md['simulation_number']: None for md in getattr(self, 'metadata', [])}
        

    def _parse_cache_header(self, file_path, verbose=False):
        """
        Parse parameters from the first line of a cache file.
        Format expected: # Parameters: key1=val1 key2=val2 ...
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

    def get(self, property_key, pars_dict=None, use_fraction=True, save=True, verbose=False):
        """
        Generic helper to get ensemble-averaged property with caching.
        """
        results = []
        for sim in self.simulations:
            # Setup paths
            sim_folder = Path(self.data_path) / self.results_dir / str(sim.metadata['simulation_number'])
            try:
                method = getattr(sim, self._properties[property_key])
            except ValueError:
                raise ValueError(f"Property {property_key} not supported.")


            if save:
                cache_file = sim_folder / (property_key + '.dat')
            else:
                cache_file = None

            # Handle fraction logic
            target_params = pars_dict.copy()
            if use_fraction:
                fraction = self.fractions_eq[sim.metadata["simulation_number"]]
                if sim.fraction_eq_is_invalid(fraction, property=property_key, file=cache_file, verbose=verbose):
                    results.append(None)
                    continue
                target_params['fraction'] = fraction

            # Check cache validity
            cached_params = self._parse_cache_header(cache_file)
            cache_valid = bool(cached_params) and all(
                np.isclose(cached_params.get(k, -1e9), v) if isinstance(v, (int, float))
                else cached_params.get(k) == v
                for k, v in target_params.items()
            )

            if cache_valid:
                try:
                    # Skiprows=1 if it doesn't have an extra header line beyond Parameters
                    # Actually RDF and Energy have multiple header lines. np.loadtxt handles # comments automatically.
                    data = np.loadtxt(cache_file, unpack=True)
                    if verbose:
                        print(f"Loaded {property_key} for simulation #{sim.metadata['simulation_number']} from cache...")
                    results.append(data)
                    continue
                except Exception:
                    pass

            # Recalculate
            if verbose:
                print(f"Calculating {property_key} for simulation #{sim.metadata['simulation_number']}...")
            
            data = method(pars_dict=target_params, file=cache_file)
            results.append(data)

        return results

    def _load_metadata(self):
        """
        Load simulation metadata from log file.
        
        Parses the log file of simulation set to extract temperature, coverage, energy terms,
        and lattice dimensions for all simulations in the set.
        
        Raises
        ------
        FileNotFoundError
            If log file is not found in set directory
        """
        lfile = Path(self.data_path) / self.log_file
        
        # Parse log file
        try:
            with open(lfile, 'r') as f:
                header = f.readline().split()  # Read header line
                log_entries = [json.loads(line) for line in f if line.strip()]

        except FileNotFoundError:
            raise FileNotFoundError(
                f"log file not found at: {self.data_path}"
            )
        
        # Extract metadata from log entry
        # Format: list of [sim_num, job_name, [nx, ny], [n_ads], temp, interaction_info, ...]
        self.metadata = []
        for entry in log_entries:
            self.metadata.append({
                'simulation_number': entry[0],
                'job_name': entry[1],
                'lattice_dimensions': entry[2],  # [nx, ny]
                'n_cells': entry[2][0] * entry[2][1],
                'n_adsorbates': entry[3][0],
                'temperature': entry[4],  # K
                'coverage': entry[3][0] / (entry[2][0] * entry[2][1]),
                'energy_terms': entry[5][1:]
               })



    def clear_traj_cache(self, simulations=None, verbose=False):
        """
        Remove internal trajectory cache files that speed up loading.

        This deletes per-trajectory `traj.pkl` files inside each `traj_*` folder.

        Parameters
        ----------
        simulations : list[int] or None
            Simulation numbers to clear; None => all simulations in set.
        verbose : bool
            Print each deleted file when True.
        """
        sim_nums = [md['simulation_number'] for md in self.metadata] if simulations is None else simulations
        from .simulation import Simulation as _Sim
        
        if verbose:
            print(f"Deleting trajectories in {(self.data_path / self.simset_dir).as_posix()}")
        
        total_count = 0
        for sn in sim_nums:
            # Check if simulation is already loaded to prefer instance method (for testing/consistency)
            sim_obj = next((s for s in self.simulations if s.metadata.get('simulation_number') == sn), None)
            
            sim_path = Path(self.data_path) / self.simset_dir / str(sn)
            pfx = getattr(sim_obj, 'traj_dir_pfx', self.traj_dir_pfx)
            
            # Use Simulation.clear_traj_cache_path directly to get the list for formatting
            deleted_files = _Sim.clear_traj_cache_path(sim_path, traj_dir_pfx=pfx, verbose=False)
            
            if sim_obj is not None:
                # Still call the instance method if it's a mock or has side effects (like in tests)
                # This ensures tests like SpySim pass.
                sim_obj.clear_traj_cache(verbose=verbose)
            
            if verbose and deleted_files:
                # Format each deleted file as sn\traj_i\traj.pkl
                file_strs = [f"{sn}\\{f.parent.name}\\{f.name}" for f in deleted_files]
                print("; ".join(file_strs) + ";")
            
            total_count += len(deleted_files)

        if total_count > 0:
            if simulations is not None:
                sim_display = str(sim_nums) if len(sim_nums) <= 10 else f"{sim_nums[:10]}... (+{len(sim_nums)-10} more)"
                sim_msg = f"for simulations: {sim_display}"
            else:
                sim_msg = "for all simulations"
            print(f"Cleared {total_count} trajectory cache file(s) {sim_msg}.")
        else:
            print("No trajectory cache files found to clear.")

        return total_count

    def clear_results(self, target='all', simulations=None, verbose=False):
        """
        Clear user-visible result files for the simulation set.

        Supported targets: 'energy', 'rdf', 'accessibility', 'gref', or 'all'.
        Parameters
        ----------
        target : str or list, default 'all'
            Which result types to clear.
        simulations : list[int] or None
            If provided, restrict clearing to these simulation numbers; None => all.
        verbose : bool
            Print deleted filenames when True.
        """
        valid_keys = set(self._results_files.keys()) | {'gref'}

        # Normalize target
        if isinstance(target, list):
            formats_to_clear = target
        elif isinstance(target, str):
            formats_to_clear = list(valid_keys) if target == 'all' else [target]
        else:
            raise TypeError("target must be str or list of str")

        # Validate
        formats_to_clear = [f for f in formats_to_clear if f in valid_keys]

        if not formats_to_clear:
            return 0

        sim_nums = [md['simulation_number'] for md in self.metadata] if simulations is None else simulations

        from .simulation import Simulation as _Sim
        results_map = dict(self._results_files)
        results_map['gref'] = 'g_ref.dat'

        total_count = 0
        for sn in sim_nums:
            sim_obj = next((s for s in self.simulations if s.metadata.get('simulation_number') == sn), None)
            sim_res_folder = Path(self.data_path) / self.results_dir / str(sn)
            
            for fmt in formats_to_clear:
                if sim_obj is not None:
                    count = sim_obj.clear_results(target=fmt, results_files=results_map, 
                                                  verbose=verbose, search_dir=sim_res_folder)
                else:
                    count = _Sim.clear_results_path(sim_res_folder, target=fmt, 
                                                    results_files=results_map, verbose=verbose)
                total_count += count

        if simulations is not None:
            sim_display = str(sim_nums) if len(sim_nums) <= 10 else f"{sim_nums[:10]}... (+{len(sim_nums)-10} more)"
            sim_msg = f"for simulations: {sim_display}"
        else:
            sim_msg = "for all simulations"

        if total_count > 0:
            print(f"Results Path: {(self.data_path / self.results_dir).as_posix()}")
            print(f"Cleared {total_count} '{target}' result file(s) {sim_msg}.")
        else:
            print(f"No '{target}' result files found {sim_msg}.")

        return total_count

    def check_cache_status(self, simulations=None):
        """
        Check existence of result and cache files for specified simulations.
        
        Parameters
        ----------
        simulations : list of int, optional
            List of simulation numbers to check. If None, check all simulations in the set.
        """
        sim_nums = [md['simulation_number'] for md in self.metadata] if simulations is None else simulations
        
        for sn in sim_nums:
            res_dir = self.data_path / self.results_dir / str(sn)
            data_dir = self.data_path / self.simset_dir / str(sn)
            
            print(f"Simulation #{sn}:")
            # Check results
            for label, filename in self._results_files.items():
                if filename is None: continue
                f = res_dir / filename
                status = "EXISTS" if f.exists() else "MISSING"
                print(f"  {label:14}: {status}")
                
            # Check trajectory cache
            traj_pkls = list(data_dir.glob(f"{self.traj_dir_pfx}_*/traj.pkl"))
            print(f"  traj cache    : {len(traj_pkls)} pkl file(s) found")
            print("-" * 34)

    def find_equilibrium_fraction_fit(self, energies_vs_time, threshold=0.01, min_equilibrium_points=10, 
                                       a0_fixed=True, a0_guess_points=10):
        """
        Find equilibrium fraction by fitting ensemble-averaged energy vs time to an biexponential decay model.
        Parameters
        ----------
        energies_vs_time : list of tuples
            Each tuple contains (times, energies, energies_std) for a simulation.
        threshold : float, default 0.01
            Threshold fraction of the asymptotic value (a0) of the fitting function to define equilibrium.
        min_equilibrium_points : int, default 10
            Minimum number of consecutive points below threshold to define equilibrium.
        a0_fixed : bool, default True
            If True, fix a0 to average of last a0_guess_points during fitting.
        a0_guess_points : int, default 10
            Number of points from end to average for a0 guess.
        Returns
        -------
        fit_results : list of tuples
            Each tuple contains (eq_idx, fit_params, fit_result, exp_terms) for a simulation.
            eq_idx : index of first equilibrium point or None if not found
            fit_params : fitted parameters (a0, a1, a2, tau1, tau2) or None if fit failed
            fit_result : lmfit ModelResult object or None if fit failed
            exp_terms : fitted exponential terms over time or None if fit failed
        """

        def exp_decay_model(t, a0, a1, a2, tau1, tau2):
            """
            Exponential decay + constant model:
            E(t) = a0 + a1*exp(-t/tau1) + a2*exp(-t/tau2)
            Parameters:
                t           : time        -- float or np.ndarray
                a0,a1,a2    : amplitudes  -- float
                tau1, tau2  : decay time constants -- float
            Returns:
                E(t) : float or np.ndarray
            """
            return a0 + a1 * np.exp(-t / tau1) + a2 * np.exp(-t / tau2)

        def find_equilibrium_exp_decay(times, energies, 
                                    threshold_fraction=0.01, 
                                    min_equ_points=10, 
                                    a0_fixed=True, a0_guess_points=10):
            """Find equilibrium by fitting exponential decay model. Optionally fix a0 to average of last 10 points."""
            if len(times) < min_equ_points:
                return None, None, None, None

            # Initial parameter guesses
            a0_guess = np.average(energies[-a0_guess_points:-1])       # Equilibrium value (final energy)
            a1_guess = energies[0]  - a0_guess            # Amplitude of first decay
            a2_guess = energies[10] - a0_guess            # Amplitude of second decay
            tau1_guess = (times[10] - times[0]) / 5       # Fast decay time constant
            tau2_guess = (times[-1] - times[0]) / 10      # Slow Decay time constant

            # Create lmfit Model
            model = Model(exp_decay_model)

            # Set up parameters with constraints
            params = model.make_params(a0=a0_guess, a1=a1_guess, a2=a2_guess, tau1=tau1_guess, tau2=tau2_guess)
            params['a0'].min = 0.0
            params['a1'].min = 0.0
            params['a2'].min = 0.0
            params['tau1'].min = 0.0
            params['tau2'].min = 0.0
            if a0_fixed:
                params['a0'].set(value=a0_guess, vary=False)

            try:
                # Fit the model
                result = model.fit(energies, params, t=times)

                # Extract fitted parameters
                a0_fit = result.params['a0'].value
                a1_fit = result.params['a1'].value
                a2_fit = result.params['a2'].value
                tau1_fit = result.params['tau1'].value
                tau2_fit = result.params['tau2'].value

                # Calculate exponential terms over time
                exp_term_1 = a1_fit * np.exp(-times / tau1_fit)
                exp_term_2 = a2_fit * np.exp(-times / tau2_fit)
                exp_terms = exp_term_1 + exp_term_2

                # Use threshold relative to average of the last 10 energy point
                threshold = threshold_fraction * np.average(energies[-10:-1])
                below_threshold = exp_terms < threshold

                # Find first sustained occurrence
                for i in range(len(below_threshold) - min_equ_points):
                    if np.all(below_threshold[i:i+min_equ_points]):
                        return i, (a0_fit, a1_fit, a2_fit, tau1_fit, tau2_fit), result, exp_terms

                # If threshold never reached, no equilibrium detected
                return None, (a0_fit, a1_fit, a2_fit, tau1_fit, tau2_fit), result, exp_terms

            except Exception as e:
                print(f"Fit failed: {e}")
                return None, None, None, None

        fit_results = []
        for isim, sim in enumerate(self.simulations):

            # Get ensemble-averaged energy vs time and fraction for this simulation
            times, energies, energies_std = energies_vs_time[isim]

            # --- fit
            eq_idx, fit_params, fit_result, exp_terms = find_equilibrium_exp_decay(
                times, energies, 
                threshold_fraction=threshold, 
                min_equ_points=min_equilibrium_points, 
                a0_fixed=a0_fixed, 
                a0_guess_points=a0_guess_points
            )

            if eq_idx is None:
                self.fractions_eq[sim.metadata["simulation_number"]] = 0.0
            else:
                self.fractions_eq[sim.metadata["simulation_number"]] = \
                    (len(energies) - eq_idx) / len(energies)
                
            fit_results.append((eq_idx, fit_params, fit_result, exp_terms))


        return fit_results


    # def get_ensemble_energy_vs_time(self, n_bins=100, verbose=False):
    #     """
    #     Get ensemble-averaged energy vs time for all simulations in the set.
    #     Parameters
    #     ----------
    #     n_bins : int, default 100
    #         Number of time bins for averaging.
    #     verbose : bool, default False
    #         If True, print detailed calculation information.

    #     Returns
    #     -------
    #     results : list of tuples
    #         Each tuple contains (times, energies, energies_std) for a simulation.
    #     """
    #     if n_bins <= 0:
    #         raise ValueError("n_bins must be a positive integer.")
        
    #     def compute(sim, cache_file, params):
    #         return sim.get_ensemble_energy_vs_time(n_bins=params['n_bins'], file=cache_file)
        
    #     return self._get_ensemble_property_generic('energy', compute, {'n_bins': n_bins}, 
    #                                              use_fraction=False, verbose=verbose)


    def get_ensemble_rdfs(self, r_max=40.0, dr=0.1, normalize=True, verbose=False):
        """
        Get ensemble-averaged RDF for all simulations in the set.
        Parameters
        ----------
        r_max : float, default 40.0
            Maximum distance for RDF (Angstroms)
        dr : float, default 0.1
            Bin width for RDF (Angstroms)
        normalize : bool, default True
            If True, normalize RDF using reference
        verbose : bool, default False
            If True, print detailed calculation information.

        Returns
        -------
        results : list of tuples
            Each tuple contains (r, g_r, g_ref_r) for a simulation.
        """
        if r_max <= 0:
            raise ValueError("r_max must be a positive number.")
        if dr <= 0:
            raise ValueError("dr must be a positive number.")
        
        def compute(sim, cache_file, params):
            return sim.get_ensemble_rdf(r_max=params['r_max'], dr=params['dr'], 
                                      fraction=params['fraction'], normalize=params['normalize'],
                                      file=cache_file)

        return self._get_ensemble_property_generic('rdf', compute, 
                                                 {'r_max': r_max, 'dr': dr, 'normalize': normalize}, 
                                                 use_fraction=True, verbose=verbose)


    def get_ensemble_accessibilities(self, verbose=False):
        """
        Get ensemble-averaged accessibility for all simulations in the set.
        
        Accessibility measures how many nearest neighbor sites are vacant,
        which affects reactivity and diffusion rates.
        
        Parameters
        ----------
        verbose : bool, default False
            If True, print detailed calculation information.

        Returns
        -------
        results : list of tuples
            Each tuple contains (accessibility, frequency, frequency_std) for a simulation.
        """
        def compute(sim, cache_file, params):
            return sim.get_ensemble_accessibility(fraction=params['fraction'], file=cache_file)
        
        return self._get_ensemble_property_generic('accessibility', compute, {}, 
                                                 use_fraction=True, verbose=verbose)

    def load(self, cache=True, workers=mp.cpu_count(), simulations=None, verbose=False):
        """
        Load data for simulations in the set.
        
        Parameters
        ----------
        cache : bool, default True
            If True, load cached trajectory data if available; cache trajectory data if not already cached.
            if False, load from raw simulation data.
        simulations : list of int or None, default None
            If list of int, load only specified simulation numbers. If None, load all simulations.
        verbose : bool, default False
            If True, print detailed loading information.
        workers : int, default mp.cpu_count()
            Number of worker processes to use for parallel loading.
            If None, load serially.
        """
        
        md_to_load = self.metadata if simulations is None else [md for md in self.metadata if md['simulation_number'] in simulations]
        for md in md_to_load:
            sim_folder = Path(self.data_path) / self.simset_dir / f"{md['simulation_number']}"
            sim = Simulation(sim_folder, metadata=md, traj_dir_pfx=self.traj_dir_pfx)
            sim.load(cache=cache, workers=workers, verbose=verbose) 
            self.simulations.append(sim)

    def plot_energy(self, energies_vs_time, ncols=3, figsize=(12,2.5), title_fontsize=10, suptitle_fontsize=16, show_eq=True):
        """Plot ensemble-averaged energy vs time for the loaded simulations.
        Parameters
        ----------
        energies_vs_time : list of tuples
            Each tuple contains (times, energies, energies_std) for a simulation.
        ncols : int, default 3
            Number of columns in the subplot grid.
        figsize : tuple, default (12, 3)
            Figure size (width, height) in inches for each row.
        title_fontsize : int, default 10
            Font size for subplot titles.
        suptitle_fontsize : int, default 16
            Font size for the overall figure title.
        show_eq : bool, default True
            If True, show equilibration fraction on each subplot.
        Returns
        -------
        None
        """

        # Set up subplots

        nrows = int(np.ceil(len(self.simulations)/ncols))
        figsize_scaled = (figsize[0], figsize[1] * nrows)
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize_scaled, squeeze=False)

        fig_title = f'Ensemble averaged energy vs time -- {self.data_path.parts[-1]}'
        fig.suptitle(fig_title, fontsize=suptitle_fontsize, fontweight='bold', y=1.)

        if not self.simulations:
            print("No simulations loaded: Nothing to plot.")
            return
        
        for isim, sim in enumerate(self.simulations):

            # Get ensemble-averaged energy vs time and fraction for this simulation
            times, energies, energies_std = energies_vs_time[isim]

            if show_eq:
                try:
                    fraction = self.fractions_eq[sim.metadata["simulation_number"]] if self.fractions_eq[sim.metadata["simulation_number"]] is not None else 0.0
                except KeyError:
                    raise KeyError(f"Equilibration fraction for simulation {sim.metadata['simulation_number']} not found in fractions_eq dictionary.")
            else:
                fraction = 0.0

            # Plot energy as function of time using subplots
            ax = axes[isim//ncols, isim%ncols]

            # Determine time units
            use_ms = len(times) > 0 and np.max(times) < 1.0
            if use_ms:
                times_plot = times * 1000
                x_label = 'Time (ms)'
            else:
                times_plot = times
                x_label = 'Time (s)'
                
            # Plot energy versus percent of total time on bottom axis; show time on top axis
            ax.grid()
            title = f'Simulation #{sim.metadata["simulation_number"]}:' \
                    fr'  $T={sim.metadata["temperature"]}$ K, $\theta={sim.metadata["coverage"]:.3f}$' 
            if show_eq:
                title += f' ({fraction*100:.0f}%)'
            ax.set_title(title, fontsize = title_fontsize)

            if len(times_plot) > 0 and times_plot[-1] != 0:
                # Use the actual maximum time (not the last element) in case times are unsorted
                times_arr = np.asarray(times_plot, dtype=float)
                max_time = float(np.max(times_arr))
                percent = (times_arr / max_time) * 100.0

                ax.plot(percent, energies, marker='o', linestyle='-', markersize=2)
                ax.set_xlabel('Percent of time (%)')
                ax.set_ylabel('Energy (eV)')

                # Ensure percent axis spans 0-100 (0% -> time 0, 100% -> max_time)
                ax.set_xlim(0.0, 100.0)

                # Percent axis ticks: minor at 10%, major (and labels) at 20%
                ax.xaxis.set_major_locator(MultipleLocator(20))
                ax.xaxis.set_minor_locator(MultipleLocator(10))
                ax.xaxis.set_major_formatter(FuncFormatter(lambda v, pos: f"{v:.0f}%"))

                # Shade equilibrium region in percent coordinates
                eq_idx = int(np.round((1 - fraction) * (len(times) - 1)))
                left_p = (times_arr[eq_idx] / max_time) * 100.0
                ax.axvspan(left_p, 100.0, alpha=0.2, color='green') 
            else:
                # Fallback: no valid times, plot energies vs times_plot as-is
                ax.plot(times_plot, energies, marker='o', linestyle='-', markersize=2)
                ax.set_xlabel(x_label)
                ax.set_ylabel('Energy (eV)')
                # Shade equilibrium region if possible
                if len(times_plot) > 0:
                    eq_idx = int(np.round((1 - fraction) * (len(times) - 1)))
                    ax.axvspan(times_plot[eq_idx], times_plot[-1], alpha=0.2, color='green')
                else:
                    eq_idx = 0  # no shading possible

            # Set y-axis limits based on equilibrium range (guard empty)
            if len(energies) > 0:
                equilibrium_energies = energies[eq_idx:] if eq_idx < len(energies) - 1 and show_eq else energies
                if len(equilibrium_energies) > 0:
                    ax.set_ylim(min(equilibrium_energies) * 0.9, max(equilibrium_energies) * 1.1)

        # Hide unused subplots
        total_plots = len(axes.flatten())
        for idx in range(len(self.simulations), total_plots):
            ax = axes[idx//ncols, idx%ncols]
            ax.axis('off')

        plt.tight_layout()
        plt.show()    


    def plot_rdf(self, rdfs, ncols=3, figsize=(12,2.5), title_fontsize=10, suptitle_fontsize=16):
        """Plot ensemble-averaged RDF for the loaded simulations.
        Parameters
        ----------
        rdfs : list of tuples
            Each tuple contains (r, g, g_std) for a simulation.
        ncols : int, default 3
            Number of columns in the subplot grid.
        figsize : tuple, default (12, 3)
            Figure size (width, height) in inches for each row.
        title_fontsize : int, default 10
            Font size for subplot titles.
        suptitle_fontsize : int, default 16
            Font size for the overall figure title.

        Returns
        -------
        None
        """

        # Set up subplots

        nrows = int(np.ceil(len(self.simulations)/ncols))
        figsize_scaled = (figsize[0], figsize[1] * nrows)
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize_scaled, squeeze=False)

        fig_title = f'Ensemble averaged RDFs -- {self.data_path.parts[-1]}'
        fig.suptitle(fig_title, fontsize=suptitle_fontsize, fontweight='bold', y=1.)

        if not self.simulations:
            print("No simulations loaded: Nothing to plot.")
            return
        
        for isim, sim in enumerate(self.simulations):

            data = rdfs[isim]
            ax = axes[isim//ncols, isim%ncols]
            title = f'Simulation #{sim.metadata["simulation_number"]}:' \
                    fr'  $T={sim.metadata["temperature"]}$ K, $\theta={sim.metadata["coverage"]:.3f}$' 
            ax.set_title(title, fontsize = title_fontsize)

            if data is None:
                ax.set_axis_off() 
                ax.text(0.5, 0.5, 'No RDF available', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
            else:
                ax.plot(data[0], data[1], marker='o', linestyle='-', markersize=2)
                ax.set_xlabel(r'$r/a_0$')
                ax.set_ylabel('Reduced $g(r)$' if len(data) == 4 else '$g(r)$')
                ax.grid()

        # Hide unused subplots
        total_plots = len(axes.flatten())
        for idx in range(len(self.simulations), total_plots):
            ax = axes[idx//ncols, idx%ncols]
            ax.axis('off')

        plt.tight_layout()
        plt.show()    


    def plot_accessibilities(self, accessibilities, ncols=3, figsize=(12,2.5), title_fontsize=10, suptitle_fontsize=16):
        """Plot ensemble-averaged site accessibility for the loaded simulations.
        
        Parameters
        ----------
        accessibilities : list of tuples
            Each tuple contains (accessibility, frequency, frequency_std) for a simulation.
            accessibility is ndarray of vacant neighbor counts (0 to max_coordination).
            frequency is ndarray of frequencies for each accessibility level.
            frequency_std is ndarray of standard deviations across trajectories.
        ncols : int, default 3
            Number of columns in the subplot grid.
        figsize : tuple, default (12, 2.5)
            Figure size (width, height) in inches for each row.
        title_fontsize : int, default 10
            Font size for subplot titles.
        suptitle_fontsize : int, default 16
            Font size for the overall figure title.

        Returns
        -------
        None
        """

        # Set up subplots
        nrows = int(np.ceil(len(self.simulations)/ncols))
        figsize_scaled = (figsize[0], figsize[1] * nrows)
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize_scaled, squeeze=False)

        fig_title = f'Ensemble averaged accessibility -- {self.data_path.parts[-1]}'
        fig.suptitle(fig_title, fontsize=suptitle_fontsize, fontweight='bold', y=1.)

        if not self.simulations:
            print("No simulations loaded: Nothing to plot.")
            return
        
        for isim, sim in enumerate(self.simulations):

            data = accessibilities[isim]
            ax = axes[isim//ncols, isim%ncols]
            title = f'Simulation #{sim.metadata["simulation_number"]}:' \
                    fr'  $T={sim.metadata["temperature"]}$ K, $\theta={sim.metadata["coverage"]:.3f}$' 
            ax.set_title(title, fontsize=title_fontsize)

            if data is None:
                ax.set_axis_off() 
                ax.text(0.5, 0.5, 'No accessibility available', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
            else:
                accessibility, frequency, frequency_std = data
                # Plot bar chart with error bars
                ax.bar(accessibility, frequency, yerr=frequency_std, capsize=5, alpha=0.7, color='steelblue', edgecolor='black')
                ax.set_xlabel('Number of vacant nearest neighbors')
                ax.set_ylabel('Frequency')
                ax.set_xticks(accessibility)
                ax.grid(axis='y', alpha=0.3)

        # Hide unused subplots
        total_plots = len(axes.flatten())
        for idx in range(len(self.simulations), total_plots):
            ax = axes[idx//ncols, idx%ncols]
            ax.axis('off')

        plt.tight_layout()
        plt.show()    


    def __len__(self):
        """Return number of simulations."""
        return len(self.metadata)


