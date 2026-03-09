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
from .lattice import Lattice

class SimulationSet:
    """
    Summary
    -------
    Manages a Zacros simulation set.
    
    This class provides a high-level interface for:
    - Metadata extraction from jobs.log
    - Loading and caching of multiple simulations in the set
    - Ensemble-averaged property calculations (energy vs time, RDFs, accessibilities)
    - Equilibration fraction estimation via fitting energy vs time to an exponential decay model
    - Cache management for trajectories and results across the set
    
    Attributes
    ----------
    data_path  : Path or str
        Directory containing the log_file, simset_dir, and results_dir
    en_file_sfx : str
        Suffix for energy data files (default: 'energy.dat').
        If None, energy data files are not created.
    fractions_eq : dict
        Dictionary mapping simulation keys to equilibrium fractions.
    log_file : str  
        name of the log file (default: 'jobs.log')
    metadata : list of dictionaries
        Simulation metadata (temperature, coverage, energy terms, etc.)
    results_dir : str
        subdirectory containing simulation results (default: 'results')
    simset_dir : str
        subdirectory containing a simulation set (default: 'jobs')
    simulations : dict
        Dictionary of loaded Simulation objects {simulation_number: Simulation}.
    use_cache : bool
        Whether caching is used when loading simulations.
    verbose : bool
        Whether to print verbose output during loading.
    
    Methods
    -------
    check_cache_status(simulations=None)
    clear_results(target='all', simulations=None, verbose=False)
    clear_traj_cache(simulations=None, verbose=False)
    find_equilibrium_fraction_fit(energies_vs_time, threshold=0.01, min_equilibrium_points=10, a0_fixed=True, a0_guess_points=10)
    get(property_key, pars_dict=None, simulations=None, use_fraction=False, save=True, verbose=False)
    load(cache=None, workers=mp.cpu_count(), simulations=None, verbose=False)
    plot_accessibilities(accessibilities, ncols=3, figsize=(12,2.5), title_fontsize=10, suptitle_fontsize=16)
    plot_energy(energies_vs_time, ncols=3, figsize=(12,2.5), title_fontsize=10, suptitle_fontsize=16, show_eq=True)
    plot_rdf(rdfs, ncols=3, figsize=(12,2.5), title_fontsize=10, suptitle_fontsize=16)

    Examples
    --------
    >>> from nat_zacros import SimulationSet
    >>> simset = SimulationSet('/path/to/simset')
    >>> simset.load()
    >>>
    >>> # Get ensemble-averaged properties via the unified get() method
    >>> energies = simset.get('energy', pars_dict={'n_bins': 100})
    >>> rdfs = simset.get('rdf', pars_dict={'dr': 0.1})
    >>> accs = simset.get('accessibility')
    >>>
    >>> # Plotting
    >>> simset.plot_energy(energies)
    >>> simset.plot_rdf(rdfs)
    >>> simset.plot_accessibilities(accs)

    Notes
    -----
    The `get()` method is the central entry point for all ensemble analysis.
    Supported property keys and their default parameters:
    - 'energy': {'n_bins': 100}
    - 'rdf': {'r_max': 40.0, 'dr': 0.1, 'fraction': 1.0, 'normalize': True}
    - 'accessibility': {'fraction': 1.0}
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
            'rdf':           'get_ensemble_rdf',
            'accessibility': 'get_ensemble_accessibility',
            'cluster':       'get_ensemble_clusters',
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

        # Load metadata from the log file to a dictionary of simulation objects
        self.simulations = self._load_metadata()

        # Initialize equilibration fractions dictionary with default None for each simulation
        # This avoids KeyError when code expects an entry per simulation unless user overrides.
        self.fractions_eq = {key: None for key in self.simulations.keys()}

    def _arg_to_list(self, arg):
        """
        Helper method converts an argument to a list of simulation keys.
        """

        if arg is None:
            arg = list(self.simulations.keys())
        elif not isinstance(arg, list):
            arg = [arg]

        return arg


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

    def get(self, property_key, pars_dict=None, simulations=None, use_fraction=False, save=True,
            autoload=False, cache=True, verbose=False):
        """
        Compute an ensemble-averaged property with caching support.

        This method retrieves the specified property for specified simulations in the set,
        either by loading from a cache file or by performing the calculation if
        the cache is missing or invalid.

        Parameters
        ----------
        property_key : str
            Key for the property to compute (e.g., 'energy', 'rdf', 'accessibility').
        pars_dict : dict, optional
            Dictionary of parameters to pass to the property computation method (default is None).
        simulations : Any, list[Any], or None
            Simulation numbers to clear; None => all simulations in set.
        use_fraction : bool, optional
            Whether to use the equilibration fraction for the computation (default is False).
        save : bool, optional
            Whether to save the computed property to a cache file (default is True).
        autoload : bool, optional
            Whether to automatically load simulations that are not currently loaded (default is False).
        cache : bool, optional
            Whether to use caching (default is True).
        verbose : bool, optional
            Whether to print detailed information during computation (default is False).

        Returns
        -------
        list
            A list containing the computed property data for each simulation.
            Entries may be None if the computation was skipped.

        Raises
        ------
        ValueError
            If the property_key is not recognized.
        """


        results = {}
        for key in self._arg_to_list(simulations):

            if key in self.simulations.keys():

                sim = self.simulations[key]
                # Check if simulation is loaded; if not, either autoload or skip based on the autoload flag
                if not sim.is_loaded:
                    if autoload:
                        print(f"Simulation #{key} is not loaded. Loading now...")
                        sim.load(cache=cache, verbose=verbose)
                    else:
                        print(f"Simulation #{key} is not loaded. Skipping.")
                        continue
                
                # Get the method to compute the property from the Simulation class instance
                try:
                    method = getattr(sim, self._properties[property_key])
                except (KeyError, AttributeError):
                    raise ValueError(f"Property {property_key} not supported.")

                if save:
                    results_file = self.data_path / self.results_dir / sim.dir.name / (property_key + '.dat')
                else:
                    results_file = None

                # Handle fraction logic
                target_params = pars_dict.copy() if pars_dict is not None else {}
                if use_fraction:
                    fraction = self.fractions_eq[key]
                    if sim.fraction_eq_is_invalid(fraction, property=property_key, file=results_file, verbose=verbose):
                        continue
                    target_params['fraction'] = fraction

                # Check results file validity
                res_params = self._parse_results_header(results_file)
                res_valid = bool(res_params) and all(
                    np.isclose(res_params.get(k, -1e9), v) if isinstance(v, float) else res_params.get(k) == v
                    for k, v in target_params.items()
                )

                if res_valid:
                    try:
                        # Skiprows=1 if it doesn't have an extra header line beyond Parameters
                        # Actually RDF and Energy have multiple header lines. np.loadtxt handles # comments automatically.
                        data = np.loadtxt(results_file, unpack=True)
                        if verbose:
                            print(f"Loaded {property_key} for simulation #{key} from results...")
                        results[key] = data
                        continue
                    except Exception:
                        pass

                # Recalculate
                if verbose:
                    print(f"Calculating {property_key} for simulation #{key}...")
                
                data = method(pars_dict=target_params, file=results_file)
                results[key] = data

            else:
                print(f'Warning: invalid simulation key {key} of {type(key)}. Ignoring.')



        return results

    def _load_metadata(self):
        """
        Load simulation metadata from log file and creates a dictionary of the corresponding
        simulation objects.
        
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
        # Format: dictionary of simulation-specific dictionaries
        sim_dict = {}
        for entry in log_entries:
            sn = str(entry[0])
            sdir = Path(self.data_path) / self.simset_dir / sn
            md = {
                'job_name': entry[1],
                'lattice_dimensions': entry[2],  # [nx, ny]
                'n_cells': entry[2][0] * entry[2][1],
                'n_adsorbates': entry[3][0],
                'temperature': entry[4],  # K
                'coverage': entry[3][0] / (entry[2][0] * entry[2][1]),
                'energy_terms': entry[5][1:]
                }

            sim_dict[sn] = Simulation(sdir, metadata=md, 
                                        traj_dir_pfx=self.traj_dir_pfx,
                                        results_dir=self.data_path /self.results_dir/ sn)
            
        return sim_dict
            
            
    def clear_traj_cache(self, simulations=None, verbose=False):
        """
        Remove internal trajectory cache files that speed up loading.
        This deletes per-trajectory `traj.pkl` files inside each trajectory folder.

        Parameters
        ----------
        simulations : Any, list[Any], or None
            Simulation numbers to clear; None => all simulations in set.
        verbose : bool
            Print each deleted file when True.
        """

        # Clear cache for selected simulations
        total_count = 0
        for key in self._arg_to_list(simulations):
            if key in self.simulations.keys():
                total_count += self.simulations[key].clear_traj_cache(verbose=verbose)
            else:
                print(f'Warning: invalid simulation key {key} of {type(key)}. Ignoring.')

        print(f"Deleted {total_count} trajectory cache file(s).")

    def clear_results(self, target='all', simulations=None, verbose=False):
        """
        Delete ensemble-averaged result files (e.g., energy.dat, rdf.dat) from the results directory.
        
        This manages files on disk without needing to load simulations into RAM first.

        Parameters
        ----------
        target : str or list, default 'all'
            Which result types to remove. Supported: 'energy', 'rdf', 'accessibility', 'gref', 'all'.
        simulations : int, list[int], 'all', or None
            Specification of simulations to clear; None or 'all' => check all in set.
        verbose : bool
            Whether to print detailed information about each deleted file.
            
        Returns
        -------
        int
            Total number of files deleted.
        """

        # Normalize target
        if target == 'all':
            formats_to_clear = list(self._properties.keys())
        elif isinstance(target, list):
            formats_to_clear = target
        else:
            formats_to_clear = [target]


        counter = 0
        for fmt in formats_to_clear:
            if fmt not in self._properties.keys():
                print(f"Warning: target '{fmt}' is not a recognized property. Ignoring.")
            else:
                for key in self._arg_to_list(simulations):
                    if key in self.simulations.keys():
                        sim = self.simulations[key]
                        fname = self.data_path / self.results_dir / sim.dir.name / (fmt + '.dat')
                        sim.clear_results(fname, verbose=verbose)
                        counter += 1
                    else:
                        print(f'Warning: invalid simulation key {key} of {type(key)}. Ignoring.')


        if counter == 0:
            print("No valid targets specified for clearing. No files deleted.")
        else:
            print(f"Deleted results for {counter} simulation(s).")


    def check_cache_status(self, simulations=None):
        """
        Check existence of result and cache files for specified simulations.
        
        Parameters
        ----------
        simulations : int, list[int], 'all', or None
            Specification of simulations to check; None or 'all' => check all in set.
        """
        sim_nums = self._get_simulation_numbers(simulations)
        
        print(f"{'Sim #':<6} | {'RAM':<6} | {'Traj Cache':<14} | {'Results'}")
        print("-" * 60)
        
        for sn in sim_nums:
            res_dir = self.data_path / self.results_dir / str(sn)
            data_dir = self.data_path / self.simset_dir / str(sn)
            
            # 1. RAM Status
            in_ram = "LOADED" if sn in self.simulations else "-"
            
            # 2. Trajectory Cache Status
            traj_pkls = list(data_dir.glob(f"{self.traj_dir_pfx}_*/traj.pkl"))
            traj_status = f"{len(traj_pkls)} pkls" if traj_pkls else "MISSING"
            
            # 3. Results Status
            found_results = []
            for label, filename in self._results_files.items():
                if filename and (res_dir / filename).exists():
                    found_results.append(label)
            
            results_str = ", ".join(found_results) if found_results else "NONE"
            
            print(f"{sn:<6} | {in_ram:<6} | {traj_status:<14} | {results_str}")

    def find_equilibrium_fraction_fit(self, energies_vs_time, threshold=0.01, min_equilibrium_points=10, 
                                       a0_fixed=True, a0_guess_points=10):
        """
        Find equilibrium fraction by fitting ensemble-averaged energy vs time to an biexponential decay model.
        Parameters
        ----------
        energies_vs_time : dict
            A dictionary where keys are simulation identifiers and values are tuples containing (times, energies, energies_std).
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

        fit_results = {}

        for key in energies_vs_time.keys():

            # Get ensemble-averaged energy vs time and fraction for this simulation
            times, energies, energies_std = energies_vs_time[key]

            # --- fit
            eq_idx, fit_params, fit_result, exp_terms = find_equilibrium_exp_decay(
                times, energies, 
                threshold_fraction=threshold, 
                min_equ_points=min_equilibrium_points, 
                a0_fixed=a0_fixed, 
                a0_guess_points=a0_guess_points
            )

            if eq_idx is None:
                self.fractions_eq[key] = 0.0
            else:
                self.fractions_eq[key] = \
                    (len(energies) - eq_idx) / len(energies)
                
            fit_results[key] = (eq_idx, fit_params, fit_result, exp_terms)


        return fit_results


    def load(self, cache=True, workers=mp.cpu_count(), simulations=None, verbose=False):
        """
        Load data for simulations in the set.
        
        Parameters
        ----------
        cache : bool, default True
            If True, load cached trajectory data if available; cache trajectory data if not already cached.
            if False, load from raw simulation data.
        simulations : Any, list[Any] or None
            Specification of simulations to load; None => load all in set.
        verbose : bool, default False
            If True, print detailed loading information.
        workers : int, default mp.cpu_count()
            Number of worker processes to use for parallel loading.
            If None, load serially.
        """

        # Load simulations selected
        for key in self._arg_to_list(simulations):
            if key in self.simulations.keys():
                self.simulations[key].load(cache=cache, verbose=verbose)
            else:
                print(f'Warning: invalid simulation key {key} of {type(key)}. Ignoring.')


    def unload(self, simulations=None, verbose=False):
        """
        Unload specified simulations from memory.
        
        Parameters
        ----------
        simulations : Any, list[Any] or None
            Specification of simulations to unload; None => unload all in set.
        verbose : bool, default False
            If True, print detailed unloading information.
        """
        # Unload simulations selected
        for key in self._arg_to_list(simulations):
            if key in self.simulations.keys():
                self.simulations[key].unload(verbose=verbose)
            else:
                print(f'Warning: invalid simulation key {key} of {type(key)}. Ignoring.')

    def plot_energies(self, energies_vs_time, ncols=3, figsize=(12,2.5), title_fontsize=10, suptitle_fontsize=16, show_eq=True):
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
        fig.suptitle(fig_title, fontsize=suptitle_fontsize, fontweight='bold', y=1.0) # slightly lower y for better spacing

        if not self.simulations:
            print("No simulations loaded: Nothing to plot.")
            return

        # Iterate over loaded simulations in numerical order
        for isim, key in enumerate(energies_vs_time.keys()):
            sim = self.simulations[key]

            # Get ensemble-averaged energy vs time and fraction for this simulation
            times, energies, energies_std = energies_vs_time[key]

            if show_eq:
                try:
                    fraction = self.fractions_eq[key] if self.fractions_eq[key] is not None else 0.0
                except KeyError:
                    raise KeyError(f"Equilibration fraction for simulation {key} not found in fractions_eq dictionary.")
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
            title = f'Simulation #{key}:' \
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
        for idx in range(len(energies_vs_time.keys()), total_plots):
            ax = axes[idx//ncols, idx%ncols]
            ax.axis('off')

        plt.tight_layout()
        plt.show()    


    def plot_rdfs(self, rdfs, ncols=3, figsize=(12,2.5), title_fontsize=10, suptitle_fontsize=16):
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
        
        # Iterate over loaded simulations in numerical order
        for isim, key in enumerate(rdfs.keys()):
            sim = self.simulations[key]

            data = rdfs[key]
            ax = axes[isim//ncols, isim%ncols]
            title = f'Simulation #{key}:' \
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
        for idx in range(len(rdfs), total_plots):
            ax = axes[idx//ncols, idx%ncols]
            ax.axis('off')

        plt.tight_layout()
        plt.show()    


    def plot_accessibilities(self, accessibilities, ncols=3, figsize=(12,2.5), title_fontsize=10, suptitle_fontsize=16):
        """Plot ensemble-averaged site accessibility for the loaded simulations.
        
        Parameters
        ----------
        accessibilities : dictionary of tuples
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
        
        # Iterate over loaded simulations in numerical order
        for isim, key in enumerate(accessibilities.keys()):
            sim = self.simulations[key]

            data = accessibilities[key]
            ax = axes[isim//ncols, isim%ncols]
            title = f'Simulation #{key}:' \
                    fr'  $T={sim.metadata["temperature"]}$ K, $\theta={sim.metadata["coverage"]:.3f}$' 
            ax.set_title(title, fontsize=title_fontsize)

            if data is None:
                ax.set_axis_off() 
                ax.text(0.5, 0.5, 'No accessibility available', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
            else:
                frequency_13_avg, frequency_13_std,frequency_2_avg, frequency_2_std = data
                # Plot bar chart with error bars
                accessibility = np.arange(len(frequency_13_avg))
                ax.bar(accessibility, frequency_13_avg, yerr=frequency_13_std, capsize=5, alpha=0.7, color='steelblue', edgecolor='black', label='13')
                ax.bar(accessibility, frequency_2_avg, yerr=frequency_2_std, capsize=5, alpha=0.7, color='lightcoral', edgecolor='black', label='2')
                ax.set_xlabel('Number of vacant nearest neighbors')
                ax.set_ylabel('Frequency')
                ax.set_xticks(accessibility)
                ax.grid(axis='y', alpha=0.3)
                ax.legend()

        # Hide unused subplots
        total_plots = len(axes.flatten())
        for idx in range(len(accessibilities), total_plots):
            ax = axes[idx//ncols, idx%ncols]
            ax.axis('off')

        plt.tight_layout()
        plt.show()    


    def __len__(self):
        """Return total number of simulations found in the log file."""
        return len(self.simulations)

    def __repr__(self):
        """String representation of SimulationSet class."""
        if len(self.simulations) > 0:
            n_sims = len(self.simulations)
            sim_info = f"{n_sims} simulations"
        else:
            sim_info = f"no simulation found"
        
        return (
            f"simulation_set({sim_info}) "
        )

    @property
    def loaded_ids(self):
        """Return list of simulation numbers currently loaded in memory."""
        return sorted([key for key in self.simulations.keys() if self.simulations[key].is_loaded])


