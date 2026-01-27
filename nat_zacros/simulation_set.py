"""
SimulationSet class for managing Zacros simulation sets
We refer to simulations as a set of simulation runs when they share the same log file.
We refer to simulations as a run when they share the same run folder.

This module provides a high-level interface for loading, caching, and analyzing
collections of runs from a Zacros simulation set.
"""

import json
import matplotlib.pyplot as plt
#import pickle
from matplotlib.ticker import MultipleLocator, FuncFormatter
import numpy as np
from lmfit import Model
#import multiprocessing as mp
#from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from .lattice import Lattice
from .trajectory import Trajectory
from .simulation import Simulation

class SimulationSet:
    """
    Manages a Zacros simulation set with multiple runs.
    
    This class provides a high-level interface for:
    - Metadata extraction from jobs.log
    
    Attributes
    ----------

    fractions_eq : dict
        Dictionary mapping run numbers to equilibrium fractions.
    log_file : str  
        name of the log file (default: 'jobs.log')
    metadata : list of dictionaries
        Simulation metadata (temperature, coverage, interactions, etc.)
    parallel : bool
        Whether to use parallel loading of simulations.
    results_dir : str
        subdirectory containing simulation results (default: 'results')
    runs_dir : str
        subdirectory containing simulation runs (default: 'jobs')
    set_dir  : Path or str
        Directory containing the log file, runs, and results)
    simulations : list of Simulation
        List of loaded Simulation objects for each run in the set.
    use_cache : bool
        Whether caching is used when loading simulations.
    verbose : bool
        Whether to print verbose output during loading.
    
    Examples
    --------
    >>> # Typical workflow
    >>> nzset = SimulationSet()
    """

    def __init__(self, set_dir, runs_dir='jobs', results_dir='results', log_file='jobs.log'):
        """
        Initialize a SimulationSet.
        
        Parameters
        ----------
        log_file : str, optional
            Name of the log file (default: 'jobs.log')
        results_dir : str, optional
            Name of the subdirectory for storing results (default: 'results')
        runs_dir : str, optional
            Name of the subdirectory containing simulation runs (default: 'jobs')
        set_dir : str or Path
            Path to simulation set directory (e.g., 'fn_3leed')
            This directory should contain jobs.log and the runs subdirectory
        """
        
        self.log_file           = log_file
        self.parallel           = True              # default parallel loading behavior
        self.results_dir        = results_dir
        self.runs_dir           = runs_dir
        self.set_dir            = Path(set_dir)
        self.use_cache          = False             # default caching behavior
        self.verbose            = False             # default verbosity
        self.simulations        = []
        

        self._load_metadata()
        self.simulations        = []   # initialize simulations list
        # Initialize equilibration fractions dictionary with default 1.0 for each run -- dja change 2026-01-22
        # This avoids KeyError when code expects an entry per run unless user overrides.
        self.fractions_eq       = {md['run_number']: None for md in getattr(self, 'metadata', [])}

        # Validate the set directory exists
        if not self.set_dir.exists():
            raise FileNotFoundError(f"Set directory not found: {self.set_dir}")
        


    def _load_metadata(self):
        """
        Load simulation metadata from log file.
        
        Parses the log file of simulation set to extract temperature, coverage, interactions,
        and lattice dimensions for all runs.
        
        Raises
        ------
        FileNotFoundError
            If log file is not found in set directory
        """
        lfile = Path(self.set_dir) / self.log_file
        
        # Parse log file
        try:
            with open(lfile, 'r') as f:
                header = f.readline().split()  # Read header line
                log_entries = [json.loads(line) for line in f if line.strip()]

        except FileNotFoundError:
            raise FileNotFoundError(
                f"log file not found at: {self.set_dir}"
            )
        
        # Extract metadata from log entry
        # Format: list of [run_num, job_name, [nx, ny], [n_ads], temp, interaction_info, ...]
        self.metadata = []
        for entry in log_entries:
            self.metadata.append({
                'run_number': entry[0],
                'job_name': entry[1],
                'lattice_dimensions': entry[2],  # [nx, ny]
                'n_cells': entry[2][0] * entry[2][1],
                'n_adsorbates': entry[3][0],
                'temperature': entry[4],  # K
                'coverage': entry[3][0] / (entry[2][0] * entry[2][1]),
                'interactions': entry[5][1:]
               })


    def clear_cache(self, verbose=False):
        """
        Clear cached data for all simulation runs in the set.
        
        """
        
        for md in self.metadata:
            run_folder = self.set_dir / self.runs_dir / f"{md['run_number']}"
            sim = Simulation(run_folder, metadata=md, log_file=self.log_file, results_dirname=self.results_dir)
            sim.clear_cache(verbose=verbose)

    def clear_rdf_normalization_cache(self):
        """
        Clear cached rdf normalization data.
        
        """
        gref_file = self.results_dir / 'gref.pkl'
        if gref_file.exists():
            gref_file.unlink()
            print(f"Cleared g_ref cache: {gref_file.name}")
        else:
            print(f"No g_ref cache to clear")
    

    def find_equilibrium_fraction_fit1(self, threshold=0.01, min_equilibrium_points=10):
        """
        Determine equilibrium points for all simulations in the set by fitting an exponential decay model.  
        Parameters
        ----------
        threshold : float, optional
            Relative threshold for determining equilibrium (default: 0.01)
        min_equilibrium_points : int, optional
            Minimum number of consecutive points below threshold to confirm equilibrium (default: 10)       
        Returns
        -------
        fit_results : list of tuples
            Each tuple contains (equilibrium_index, fit_parameters, fit_result, exp_term) for each simulation.
            - equilibrium_index : int or None
                Index of first equilibrium point, or None if not found.
            - fit_parameters : tuple or None
                Fitted parameters (A, tau, C) of the exponential decay model, or None if fit failed.
            - fit_result : lmfit ModelResult or None
                Full fit result object from lmfit, or None if fit failed.
            - exp_term : np.ndarray or None
                Exponential term values over time from the fitted model, or None if fit failed.
        """

        def exp_decay_model(t, A, tau, C):
            """Exponential decay model: E(t) = A*exp(-t/tau) + C"""
            return A * np.exp(-t / tau) + C

        def find_equilibrium_exp_decay(times, energies, threshold, min_eq_points):
            """Find equilibrium by fitting exponential decay model."""
            if len(times) < min_eq_points:
                return None, None, None, None
            
            # Initial parameter guesses
            C_guess = energies[-1]  # Equilibrium value ~ final energy
            A_guess = energies[0] - C_guess  # Initial amplitude
            tau_guess = times[-1] / 3  # Rough time constant
            
            # Create lmfit Model
            model = Model(exp_decay_model)
            
            # Set up parameters with constraints
            params = model.make_params(A=A_guess, tau=tau_guess, C=C_guess)
            params['A'].min = 0  # A must be positive
            params['tau'].min = 0  # tau must be positive
            
            try:
                # Fit the model
                result = model.fit(energies, params, t=times)
                
                # Extract fitted parameters
                A_fit = result.params['A'].value
                tau_fit = result.params['tau'].value
                C_fit = result.params['C'].value
                
                # Calculate exponential term over time
                exp_term = A_fit * np.exp(-times / tau_fit)
                
                below_threshold = exp_term < threshold*C_fit
                
                # Find first sustained occurrence
                for i in range(len(below_threshold) - min_eq_points):
                    if np.all(below_threshold[i:i+min_eq_points]):
                        return i, (A_fit, tau_fit, C_fit), result, exp_term
                
                # If threshold never reached, no equilibrium detected
                return None , (A_fit, tau_fit, C_fit), result, exp_term
                
            except Exception as e:
                print(f"Fit failed: {e}")
                return None, None, None, None

        fit_results = []
        for sim in self.simulations:

            # Get ensemble-averaged energy vs time and fraction for this simulation
            times, energies, energies_std = sim.get_ensemble_energy_vs_time()
            
            # Find equilibrium point
            eq_idx, fit_params, fit_result, exp_term = find_equilibrium_exp_decay(
                times, energies, threshold, min_equilibrium_points
            )

            if fit_params is None or fit_result is None:
                print(f'Run #{sim.metadata["run_number"]}: FIT FAILED\n\n')

            self.fractions_eq[sim.metadata["run_number"]] = \
                (len(energies) - eq_idx) / len(energies)
            
            fit_results.append((eq_idx, fit_params, fit_result, exp_term))

        return fit_results
            
    def find_equilibrium_fraction_fit2(self, threshold=0.01, min_equilibrium_points=10, 
                                       a0_fixed=True, a0_guess_points=10):

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

                # Use threshold relative to last energy point, not fitted amplitude
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
            times, energies, energies_std = sim.get_ensemble_energy_vs_time()

            try:
                fraction_eq = self.fractions_eq[sim.metadata["run_number"]]
            except KeyError:
                raise KeyError(f"Equilibration fraction for run {sim.metadata['run_number']} not found in fractions_eq dictionary.")

            # --- fit
            eq_idx, fit_params, fit_result, exp_terms = find_equilibrium_exp_decay(
                times, energies, 
                threshold_fraction=threshold, 
                min_equ_points=min_equilibrium_points, 
                a0_fixed=a0_fixed, 
                a0_guess_points=a0_guess_points
            )

            fit_results.append((eq_idx, fit_params, fit_result, exp_terms))
            self.fractions_eq[sim.metadata["run_number"]] = \
                (len(energies) - eq_idx) / len(energies)


        return fit_results


    def load(self, energy_only=False, parallel=True, use_cache=True, verbose=False):
        """
        Load data for all simulation runs in the set with caching support.
        
        Parameters
        ----------
        energy_only : bool, default False
            If True, only load energy data from simulation (faster, less memory).
            If False, load full simulation data.
        parallel : bool, default True
            If True, use parallel loading (recommended for full-state data).
            If False, use sequential loading.
        use_cache : bool, default True
            If True, load from cache file if available, otherwise parse and cache.
            If False, always parse from history_output.txt files.
        verbose : bool, default False
            If True, print detailed loading information.
        """

        if len(self.simulations) > 0 and energy_only:
            print("load() call is ignored: simulations set is not empty for energy only loading.")
            return

        for md in self.metadata:
            run_folder = self.set_dir / self.runs_dir / f"{md['run_number']}"
            sim = Simulation(run_folder, metadata=md, log_file=self.log_file, results_dirname=self.results_dir)
            sim._fraction_loaded = self.fractions_eq[md['run_number']] if self.fractions_eq[md['run_number']] is not None else 1.0
            sim.load(use_cache=use_cache, parallel=parallel, energy_only=energy_only, verbose=verbose)  # Load simulation data
            self.simulations.append(sim)


    def plot(self, ncols=3, figsize=(12,16), title_fontsize=10, suptitle_fontsize=16):
        """Plot ensemble-averaged energy vs time for all simulations in the set."""

        # Set up subplots
        fig, axes = plt.subplots(int(np.ceil(len(self)/ncols)), ncols, figsize=figsize)
        fig_title = f'Ensemble averaged energy vs time -- {self.set_dir.parts[-1]}'
        fig.suptitle(fig_title, fontsize=suptitle_fontsize, fontweight='bold', y=1.)

        if not self.simulations:
            print("Loading energy data automatically...")
            self.load_energy()

        for isim, sim in enumerate(self.simulations):

            # Get ensemble-averaged energy vs time and fraction for this simulation
            times, energies, energies_std = sim.get_ensemble_energy_vs_time()
            try:
                fraction = self.fractions_eq[sim.metadata["run_number"]] if self.fractions_eq[sim.metadata["run_number"]] is not None else 1.0
            except KeyError:
                raise KeyError(f"Equilibration fraction for run {sim.metadata['run_number']} not found in fractions_eq dictionary.")

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
            ax.set_title(f'Run #{sim.metadata["run_number"]}' 
                        fr'  $T={sim.metadata["temperature"]}$ K, $\theta={sim.metadata["coverage"]:.3f}$'
                        f' ({fraction*100:.0f}%)',
                        fontsize = title_fontsize)

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
                eq_idx = int(np.round((1 - fraction) * len(times)))
                left_p = (times_arr[eq_idx] / max_time) * 100.0
                ax.axvspan(left_p, 100.0, alpha=0.2, color='green')
            else:
                # Fallback: no valid times, plot energies vs times_plot as-is
                ax.plot(times_plot, energies, marker='o', linestyle='-', markersize=2)
                ax.set_xlabel(x_label)
                ax.set_ylabel('Energy (eV)')
                # Shade equilibrium region if possible
                if len(times_plot) > 0:
                    eq_idx = int(np.round((1 - fraction) * len(times)))
                    ax.axvspan(times_plot[eq_idx], times_plot[-1], alpha=0.2, color='green')

            

            # Set y-axis limits based on equilibrium range
            equilibrium_energies = energies[eq_idx:]
            ax.set_ylim(min(equilibrium_energies) * 0.9, max(equilibrium_energies) * 1.1)

        # Hide unused subplots
        total_plots = len(axes.flatten())
        for idx in range(len(self.simulations), total_plots):
            ax = axes[idx//ncols, idx%ncols]
            ax.axis('off')

        plt.tight_layout()
        plt.show()    


    def __len__(self):
        """Return number of runs."""
        return len(self.metadata)


