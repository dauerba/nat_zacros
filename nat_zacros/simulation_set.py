"""
SimulationSet class for managing Zacros simulation sets
We refer to calculations as a set when they share the same log file.
We refer to calculations as a simulation when they share the same simulation folder.

This module provides a high-level interface for loading, caching, and analyzing a Zacros simulation set.
"""


import json # 
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
    Manage a set of Zacros simulations that share one metadata log file.

    The class parses simulation metadata from ``jobs.log``, creates one
    :class:`Simulation` object per simulation directory, and provides a
    set-level interface for loading data, clearing caches, computing
    ensemble properties, and plotting per-simulation results.

    Attributes
    ----------
    data_path : pathlib.Path
        Root directory containing the log file, simulation directories, and
        results directory.
    log_file : str
        Name of the metadata log file.
    simset_dir : str
        Name of the subdirectory containing simulation folders.
    results_dir : str
        Name of the subdirectory used for cached property outputs.
    traj_dir_pfx : str
        Prefix used to detect trajectory directories inside each simulation.
    simulations : dict[str, Simulation]
        Simulation objects keyed by the simulation identifier read from the
        log file.
    fractions_eq : dict[str, float or None]
        Equilibrium fractions keyed by simulation identifier. Values remain
        ``None`` until estimated or assigned.
    verbose : bool
        Default verbosity flag stored on the instance.

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
    >>> simset.plot_energies(energies)
    >>> simset.plot_rdf(rdfs)
    >>> simset.plot_accessibilities(accs)

    Notes
    -----
    The :meth:`get` method is the main entry point for set-level analysis.
    Supported property keys and their default parameters are:
    - 'accessibility': {'fraction': 1.0}
    - 'cluster': {'cutoff': <3rd NN distance>, 'eps': 1e-4, 'fraction': 1.0, 'method': 'ckdtree'}
    - 'coverage': {'n_bins': 100, 'atoms_per_uc': 1}
    - 'energy': {'n_bins': 100}
    - 'leed': {'n_bins': 100, 'r': <3rd NN distance>}
    - 'rdf': {'r_max': 40.0, 'dr': 0.1, 'fraction': 1.0, 'normalize': True}
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
            'accessibility': 'get_ensemble_accessibility',
            'cluster':       'get_ensemble_clusters',
            'coverage':      'get_ensemble_coverages_vs_time',
            'energy':        'get_ensemble_energy_vs_time',
            'leed':          'get_ensemble_leed_intensity_vs_time',
            'rdf':           'get_ensemble_rdf',
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
        Convert a simulation selector to a list of internal simulation keys.

        ``None`` and ``'all'`` select all simulations. Scalar selectors are
        wrapped into a list. All explicit selectors are normalized to strings
        because simulation identifiers are stored as string keys.
        """

        if arg is None or arg == 'all':
            arg = list(self.simulations.keys())
        elif not isinstance(arg, list):
            arg = [arg]
        else:
            arg = list(arg)

        if arg is not None:
            arg = [str(key) for key in arg]

        return arg


    def get(self, property_key, pars_dict=None, simulations=None, use_fraction=False, save=True,
            autoload=False, cache=True, verbose=False):
        """
        Compute one cached ensemble property for one or more simulations.

        For each selected simulation, this method optionally loads the
        simulation, applies the requested parameter set, and dispatches to the
        corresponding :class:`Simulation` ensemble method.

        Parameters
        ----------
        property_key : str
            Key for the property to compute (e.g., 'energy', 'rdf', 'accessibility').
        pars_dict : dict, optional
            Dictionary of parameters to pass to the property computation method (default is None).
        simulations : str, int, list[str | int], or None
            Simulation identifiers to evaluate. ``None`` selects all
            simulations in the set. Integer selectors are converted to the
            string keys used internally.
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
        dict
            Mapping from simulation identifier to the computed property result.
            Simulations that are invalid, unloaded with ``autoload=False``, or
            skipped because of an invalid equilibrium fraction are omitted.

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

                # Recalculate
                if verbose:
                    print(f"Getting {property_key} for simulation #{key}...")
                
                data = method(pars_dict=target_params, file=results_file, verbose=verbose)
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
                'job_name'          : entry[1],
                'lattice_dimensions': entry[2],  # [nx, ny]
                'n_cells'           : entry[2][0] * entry[2][1],
                'surf_species_names': entry[3],
                'n_adsorbates'      : entry[4][0],
                'temperature'       : entry[5],  # K
                'coverage'          : entry[4][0] / (entry[2][0] * entry[2][1]),
                'energy_terms'      : entry[6][1:]
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
        simulations : str, int, list[str | int], or None
            Simulation identifiers to clear. ``None`` selects all simulations
            in the set. Integer selectors are converted to internal string
            keys.
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

        if total_count > 0:
            print(f"Deleted {total_count} trajectory cache file(s).")
        else:
            print(f"No trajectory cache files found to delete.")

    def clear_results(self, target='all', simulations=None, verbose=False):
        """
        Delete cached property files from the per-simulation results folders.

        This method only removes files on disk. It does not unload simulations
        or recompute any properties.

        Parameters
        ----------
        target : str or list, default 'all'
            Which property outputs to remove. ``'all'`` expands to all keys in
            ``self._properties``.
        simulations : str, int, list[str | int], 'all', or None
            Simulation identifiers to clear. ``None`` or ``'all'``
            selects all simulations in the set.
        verbose : bool
            Whether to print detailed information about each deleted file.
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
                        #fname = self.data_path / self.results_dir / sim.dir.name / (fmt + '.dat')
                        path = self.data_path / self.results_dir / sim.dir.name
                        fnames = path.glob(f'{fmt}*.dat')
                        counter += sim.clear_results(fnames, verbose=verbose)
                    else:
                        print(f'Warning: invalid simulation key {key} of {type(key)}. Ignoring.')


        if counter == 0:
            print("No valid targets specified for clearing. No files deleted.")
        else:
            print(f"Deleted results for {counter} simulation(s).")


    def check_cache_status(self, simulations=None):
        """
        Print a simple cache status table for selected simulations.

        The report includes whether the simulation object is present in the
        set, how many trajectory pickle caches are present, and which cached
        result files were detected.
        
        Parameters
        ----------
        simulations : str, int, list[str | int], 'all', or None
            Simulation identifiers to check. ``None`` or ``'all'``
            selects all simulations in the set.
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


    def load(self, cache=True, zfile='history_output', workers=mp.cpu_count(), simulations=None, verbose=False):
        """
        Load trajectory data for selected simulations in the set.
        
        Parameters
        ----------
        cache : bool, default True
            If True, load cached trajectory data if available; cache trajectory data if not already cached.
            if False, load from raw simulation data.
        zfile: str, default 'history_output'
            selects zacros output file from which to read states
        simulations : str, int, list[str | int], 'all', or None
            Simulation identifiers to load. ``None`` or ``'all'``
            selects all simulations in the set.
        verbose : bool, default False
            If True, print detailed loading information.
        workers : int, default mp.cpu_count()
            Reserved for future parallel loading support. The current
            implementation forwards sequentially to each :class:`Simulation`.
        """

        # Load simulations selected
        for key in self._arg_to_list(simulations):
            if key in self.simulations.keys():
                self.simulations[key].load(cache=cache, zfile=zfile, verbose=verbose)
            else:
                print(f'Warning: invalid simulation key {key} of {type(key)}. Ignoring.')


    def unload(self, simulations=None, verbose=False):
        """
        Unload specified simulations from memory.
        
        Parameters
        ----------
        simulations : str, int, list[str | int], 'all', or None
            Simulation identifiers to unload. ``None`` or ``'all'``
            selects all simulations in the set.
        verbose : bool, default False
            If True, print detailed unloading information.
        """
        # Unload simulations selected
        for key in self._arg_to_list(simulations):
            if key in self.simulations.keys():
                self.simulations[key].unload(verbose=verbose)
            else:
                print(f'Warning: invalid simulation key {key} of {type(key)}. Ignoring.')

    def plot_energies(self, energies_vs_time, ncols=3, figsize=(12,2.5), title_fontsize=10, suptitle_fontsize=16, 
                      show_eq=True, ylimit=None):
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
        ylimit : str or None, default None
            If 'eq', scales y-limit of plots to show approach to equilibrium from above
        Returns
        -------
        None
        """

        # Set up subplots

        nrows = int(np.ceil(len(energies_vs_time)/ncols))
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
            if np.all(np.isnan(energies)):
                sim_is_valid = False
            else:
                sim_is_valid = True

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
            title = f'Simulation #{key}:' \
                    fr'  $T={sim.metadata["temperature"]}$ K, $\theta={sim.metadata["coverage"]:.3f}$' 
            if show_eq:
                title += f' ({fraction*100:.0f}%)'
            ax.set_title(title, fontsize = title_fontsize)

            if sim_is_valid:
                ax.grid()

                if len(times_plot) > 0 and times_plot[-1] != 0:
                    # Use the actual maximum time (not the last element) in case times are unsorted
                    times_arr = np.asarray(times_plot, dtype=float)
                    max_time = float(np.max(times_arr))
                    percent = (times_arr / max_time) * 100.0

                    ax.plot(times_arr, energies, marker='o', linestyle='-', markersize=2)
                    ax.set_xlabel(x_label)
                    ax.set_ylabel('Energy (eV)')
                    ax.set_xlim(0.0, max_time)

                    ax1 = ax.twiny()  # Create a twin x-axis for time in seconds

                    # ax1.plot(percent, energies, marker='o', linestyle='-', markersize=2)
                    ax1.set_xlabel(' ')

                    # Ensure percent axis spans 0-100 (0% -> time 0, 100% -> max_time)
                    ax1.set_xlim(0.0, 100.0)

                    # Percent axis ticks: minor at 10%, major (and labels) at 20%
                    ax1.xaxis.set_major_locator(MultipleLocator(20))
                    ax1.xaxis.set_minor_locator(MultipleLocator(10))
                    ax1.xaxis.set_major_formatter(FuncFormatter(lambda v, pos: f"{v:.0f}%"))

                    # Shade equilibrium region in percent coordinates
                    eq_idx = int(np.round((1 - fraction) * (len(times) - 1)))
                    left_p = (times_arr[eq_idx] / max_time) * 100.0
                    ax1.axvspan(left_p, 100.0, alpha=0.2, color='green') 

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
                if ylimit == 'eq':
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


    def plot_leed_intensities(self, leed_intensities_vs_time, ncols=3, figsize=(12,2.5), title_fontsize=10, suptitle_fontsize=16):
        """Plot ensemble-averaged LEED intensity vs time for the loaded simulations.
        Parameters
        ----------
        LEED_intensity_vs_time : list of tuples
            Each tuple contains (times, avgs, stds) for a simulation.
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

        fig_title = f'Ensemble averaged LEED intensity vs time -- {self.data_path.parts[-1]}'
        fig.suptitle(fig_title, fontsize=suptitle_fontsize, fontweight='bold', y=1.0) # slightly lower y for better spacing

        if not self.simulations:
            print("No simulations loaded: Nothing to plot.")
            return

        # Iterate over loaded simulations in numerical order
        for isim, key in enumerate(leed_intensities_vs_time.keys()):
            sim = self.simulations[key]

            # Get ensemble-averaged energy vs time and fraction for this simulation
            times, avgs, stds = leed_intensities_vs_time[key]

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
                
            # Plot leed_intensity versus percent of total time on bottom axis; show time on top axis
            ax.grid()
            title = f'Simulation #{key}:' \
                    fr'  $T={sim.metadata["temperature"]}$ K, $\theta={sim.metadata["coverage"]:.3f}$' 
            ax.set_title(title, fontsize = title_fontsize)

            if len(times_plot) > 0 and times_plot[-1] != 0:
                # Use the actual maximum time (not the last element) in case times are unsorted
                times_arr = np.asarray(times_plot, dtype=float)
                max_time = float(np.max(times_arr))
                percent = (times_arr / max_time) * 100.0

                ax.plot(percent, avgs, marker='o', linestyle='-', markersize=2)
                ax.set_xlabel('Percent of time (%)')
                ax.set_ylabel('LEED intensity')

                # Ensure percent axis spans 0-100 (0% -> time 0, 100% -> max_time)
                ax.set_xlim(0.0, 100.0)

                # Percent axis ticks: minor at 10%, major (and labels) at 20%
                ax.xaxis.set_major_locator(MultipleLocator(20))
                ax.xaxis.set_minor_locator(MultipleLocator(10))
                ax.xaxis.set_major_formatter(FuncFormatter(lambda v, pos: f"{v:.0f}%"))

            else:
                # Fallback: no valid times, plot energies vs times_plot as-is
                ax.plot(times_plot, avgs, marker='o', linestyle='-', markersize=2)
                ax.set_xlabel(x_label)
                ax.set_ylabel('LEED intensity')

        # Hide unused subplots
        total_plots = len(axes.flatten())
        for idx in range(len(leed_intensities_vs_time.keys()), total_plots):
            ax = axes[idx//ncols, idx%ncols]
            ax.axis('off')

        plt.tight_layout()
        plt.show()    


    def plot_coverages(self, covs_vs_time, ncols=3, figsize=(12,2.5), title_fontsize=10, suptitle_fontsize=16, show_eq=True, species=None):
        """Plot ensemble-averaged coverages vs time for the loaded simulations.
        Parameters
        ----------
        covs_vs_time : dict
            Dictionary mapping simulation keys to tuples of (times, coverages, coverages_std).
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
        species : str, list of str, or None, default None
            If None, plot all species. If str, plot only that species.
            If list of str, plot only those species.
        Returns
        -------
        None
        """

        # Set up subplots

        nrows = int(np.ceil(len(covs_vs_time)/ncols))
        figsize_scaled = (figsize[0], figsize[1] * nrows)
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize_scaled, squeeze=False)

        fig_title = f'Ensemble averaged coverages vs time -- {self.data_path.parts[-1]}'
        fig.suptitle(fig_title, fontsize=suptitle_fontsize, fontweight='bold', y=1.0) # slightly lower y for better spacing

        if not self.simulations:
            print("No simulations loaded: Nothing to plot.")
            return

        # Iterate over loaded simulations in numerical order
        for isim, key in enumerate(covs_vs_time.keys()):
            sim = self.simulations[key]

            # Get ensemble-averaged energy vs time and fraction for this simulation
            times, covs, covs_std = covs_vs_time[key]

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
                
            # Plot coverage versus percent of total time on bottom axis; show time on top axis
            ax.grid()
            title = f'Simulation #{key}:' \
                    fr'  $T={sim.metadata["temperature"]}$ K, $\theta={sim.metadata["coverage"]:.3f}$' 
            if show_eq:
                title += f' ({fraction*100:.0f}%)'
            ax.set_title(title, fontsize = title_fontsize)

            # Determine which species to plot
            if species is None:
                # Plot all species
                species_indices = range(1, covs.shape[1])
            else:
                # Convert species to list if it's a string
                species_list = [species] if isinstance(species, str) else species
                # Find indices for requested species
                species_indices = []
                for sp in species_list:
                    if sp in sim.metadata["surf_species_names"]:
                        # Add 1 because covs[:,0] is total coverage
                        idx = sim.metadata["surf_species_names"].index(sp) + 1
                        species_indices.append(idx)
                    else:
                        print(f"Warning: Species '{sp}' not found in simulation {key}. Available: {sim.metadata['surf_species_names']}")

            for i in species_indices:
                ax.plot(times_plot, covs[:,i], label=f'{sim.metadata["surf_species_names"][i-1]}')
            ax.set_xlabel(x_label)
            ax.set_ylabel('Coverage (ML)')
            ax.legend()
            # Shade equilibrium region if possible
            if len(times_plot) > 0:
                eq_idx = int(np.round((1 - fraction) * (len(times) - 1)))
                ax.axvspan(times_plot[eq_idx], times_plot[-1], alpha=0.2, color='green')
            else:
                eq_idx = 0  # no shading possible

            # Set y-axis limits based on equilibrium range (guard empty)
            # if len(energies) > 0:
            #     equilibrium_energies = energies[eq_idx:] if eq_idx < len(energies) - 1 and show_eq else energies
            #     if len(equilibrium_energies) > 0:
            #         ax.set_ylim(min(equilibrium_energies) * 0.9, max(equilibrium_energies) * 1.1)

        # Hide unused subplots
        total_plots = len(axes.flatten())
        for idx in range(len(covs_vs_time.keys()), total_plots):
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


    def plot_clusters(self, clusters, scale='log', ncols=3, figsize=(12,2.5), title_fontsize=10, suptitle_fontsize=16):
        """Plot ensemble-averaged site accessibility for the loaded simulations.
        
        Parameters
        ----------
        clusters : dictionary of tuples
            Each tuple contains (frequency, frequency_std) for a simulation.
            frequency is ndarray of mean frequencies for cluster sizes.
            frequency_std is ndarray of standard deviations across trajectories.
        scale : str, default 'log'
            Y-axis scale for plotting ('log' or 'linear')
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

        fig_title = f'Ensemble averaged cluster size histogram -- {self.data_path.parts[-1]}'
        fig.suptitle(fig_title, fontsize=suptitle_fontsize, fontweight='bold', y=1.)

        if not self.simulations:
            print("No simulations loaded: Nothing to plot.")
            return
        
        # Iterate over loaded simulations in numerical order
        for isim, key in enumerate(clusters.keys()):
            sim = self.simulations[key]

            data = clusters[key]
            ax = axes[isim//ncols, isim%ncols]
            title = f'Simulation #{key}:' \
                    fr'  $T={sim.metadata["temperature"]}$ K, $\theta={sim.metadata["coverage"]:.3f}$' 
            ax.set_title(title, fontsize=title_fontsize)

            if data is None:
                ax.set_axis_off() 
                ax.text(0.5, 0.5, 'No histogram available', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
            else:
                frequency_avg, frequency_std = data
                # Plot bar chart with error bars
                cluster_sizes = np.arange(len(frequency_avg))
                ax.plot(cluster_sizes, frequency_avg)
                ax.set_xlabel('Cluster size')
                ax.set_ylabel('Frequency')
                #ax.set_xticks(cluster_sizes)
                ax.grid(axis='y', alpha=0.3)
                ax.set_yscale(scale)
                #ax.legend()

        # Hide unused subplots
        total_plots = len(axes.flatten())
        for idx in range(len(clusters), total_plots):
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
        """Return list of simulation identifiers currently loaded in memory."""
        return sorted([key for key in self.simulations.keys() if self.simulations[key].is_loaded])


