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
#from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import tarfile
from .lattice import Lattice
from .trajectory import Trajectory
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
        Simulation metadata (temperature, coverage, interactions, etc.)
    parallel : bool
        Whether to use parallel loading of simulations.
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
    >>> # Typical workflow
    >>> nzset = SimulationSet()
    """

    def __init__(self, data_path, log_file='jobs.log', results_dir='results', simset_dir='jobs', 
                 en_file_sfx='energy.dat', rdf_file_sfx='rdf.dat', acc_file_sfx='accessibility.dat'):
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
            Name of the subdirectory containing simulations (default: 'jobs')
        en_file_sfx : str, optional
            Suffix for energy data files (default: 'energy.dat').
            If None, energy data files are not created.
        rdf_file_sfx : str, optional
            Suffix for RDF data files (default: 'rdf.dat').
            If None, RDF data files are not created.
        acc_file_sfx : str, optional
            Suffix for accessibility data files (default: 'accessibility.dat').
            If None, accessibility data files are not created.
        """
        
        # Validate the simulation set directory exists
        if not Path(data_path).exists():
            raise FileNotFoundError(f"Simulation set directory not found: {data_path}")

        self.data_path          = Path(data_path)
        self.acc_file_sfx       = acc_file_sfx
        self.en_file_sfx        = en_file_sfx
        self.log_file           = log_file
        self.parallel           = True              # default parallel loading behavior
        self.rdf_file_sfx       = rdf_file_sfx
        self.results_dir        = results_dir
        self.simset_dir         = simset_dir
        self.verbose            = False             # default verbosity
        self.simulations        = []

        # Validate if simulation set directory exists
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
        

    def _load_metadata(self):
        """
        Load simulation metadata from log file.
        
        Parses the log file of simulation set to extract temperature, coverage, interactions,
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
                'interactions': entry[5][1:]
               })

    def clear_cache(self, cache=None, simulations=None, verbose=False):
        """
        Clear cached data for simulations in the set for the specified cache type.
        Parameters
        ----------
        cache : str or None, default None
            If str, clear cache of specified file format. If None, clear all cache types.
        simulations : list of int or None, default None
            If list of int, clear cache only for specified simulation numbers. If None, clear for all simulation set.
        verbose : bool, default False
            If True, print detailed clearing information.
        
        """
       
        md_to_clear = self.metadata if simulations is None else [md for md in self.metadata if md['simulation_number'] in simulations]
        for md in md_to_clear:
            sim_folder = Path(self.data_path) / self.simset_dir / f"{md['simulation_number']}"
            sim = Simulation(sim_folder, metadata=md, log_file=self.log_file, results_dirname=self.results_dir)
            sim.clear_cache(cache=cache, verbose=verbose)


    def clear_energy_cache(self, verbose=False):

        """
        Clear cached energy vs time data files for all simulations in the set.
        Parameters
        ----------
        verbose : bool, default False
            If True, print detailed clearing information.
        
        """
       
        for sim in self.simulations:

            if self.en_file_sfx is not None:
                en_file = Path(self.data_path) / self.results_dir / \
                            f"{sim.metadata['simulation_number']}_{self.en_file_sfx}"

                if en_file.exists():
                    en_file.unlink()
                    if verbose:
                        print(f"Energy cache cleared for simulation #{sim.metadata['simulation_number']} ({en_file.name})")


    def clear_rdf_cache(self, verbose=False):

        """
        Clear cached RDF data files for all simulations in the set.
        Parameters
        ----------
        verbose : bool, default False
            If True, print detailed clearing information.
        
        """
       
        for sim in self.simulations:

            if self.rdf_file_sfx is not None:
                rdf_file = Path(self.data_path) / self.results_dir / \
                            f"{sim.metadata['simulation_number']}_{self.rdf_file_sfx}"

                if rdf_file.exists():
                    rdf_file.unlink()
                    if verbose:
                        print(f"RDF cache cleared for simulation #{sim.metadata['simulation_number']} ({rdf_file.name})")

    def clear_accessibility_cache(self, verbose=False):
        """
        Clear cached accessibility data files for all simulations in the set.
        
        Parameters
        ----------
        verbose : bool, default False
            If True, print detailed clearing information.
        """
        for sim in self.simulations:
            if self.acc_file_sfx is not None:
                acc_file = Path(self.data_path) / self.results_dir / \
                           f"{sim.metadata['simulation_number']}_{self.acc_file_sfx}"

                if acc_file.exists():
                    acc_file.unlink()
                    if verbose:
                        print(f"Accessibility cache cleared for simulation #{sim.metadata['simulation_number']} ({acc_file.name})")

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


    def get_ensemble_energy_vs_time(self, n_bins=100, verbose=False):

        """
        Get ensemble-averaged energy vs time for all simulations in the set.
        Parameters
        ----------
        n_bins : int, default 100
            Number of time bins for averaging.
        verbose : bool, default False
            If True, print detailed calculation information.

        Returns
        -------
        results : list of tuples
            Each tuple contains (times, energies, energies_std) for a simulation.
        """

        if n_bins <= 0:
            raise ValueError("n_bins must be a positive integer.")
        
        results = []
        for sim in self.simulations:

            n_bins_file = 0
            if self.en_file_sfx is not None:
                en_file = Path(self.data_path) / self.results_dir / \
                            f"{sim.metadata['simulation_number']}_{self.en_file_sfx}"
            else:
                en_file = None

            try:
                # Load energy vs time data from file
                with open(en_file, 'r') as f:
                    n_bins_file = int(f.readline().split()[-1])
                    header = f.readline()  # skip header
                times, energies, energies_std = np.loadtxt(en_file, unpack=True)
            except: 
                pass
                
            if verbose:                
                print(f"Calculating mean energy for simulation #{sim.metadata['simulation_number']}...")
                print(f"Parameters:  n_bins={n_bins}")
                print(f"Cached vals: n_bins={n_bins_file}")

            if n_bins != n_bins_file:
                # Recalculate averages
                times, energies, energies_std = sim.get_ensemble_energy_vs_time(n_bins=n_bins)
                # Save energy vs time data to file
                try:
                    np.savetxt(en_file, 
                        np.column_stack((times, energies, energies_std)), 
                        header=f'Parameters: n_bins = {n_bins}\nTime_s Energy_eV Energy_std_eV')
                except:
                    pass
        
            results.append((times, energies, energies_std))

        return results


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
        
        results = []
        for sim in self.simulations:

            if self.rdf_file_sfx is not None:
                rdf_file   = Path(self.data_path) / self.results_dir / \
                               f"{sim.metadata['simulation_number']}_{self.rdf_file_sfx}"
            else:
                rdf_file = None

            # Get equilibrium fraction for this simulation
            fraction_eq = self.fractions_eq[sim.metadata["simulation_number"]]
            # if not or badly specified, skip and clear cache if exists
            if fraction_eq is None or fraction_eq <= 0.0 or fraction_eq > 1.0:
                # Clean cache for invalid fraction
                if rdf_file.exists():
                    rdf_file.unlink()
                    if verbose:
                        print(f"RDF cache cleared for simulation #{sim.metadata['simulation_number']} ({rdf_file.name})")
                # put placeholder None in results to maintain order
                results.append(None)  
                continue  
           
            r_max_file = 0
            dr_file = 0
            normalize_file = False
            fraction_eq_file = 0
            try:
                # Load RDF data from file
                with open(rdf_file, 'r') as f:
                    pars_line = f.readline().strip().split()
                    r_max_file = float(pars_line[3])
                    dr_file = float(pars_line[5])
                    fraction_eq_file = float(pars_line[7])
                    normalize_file = pars_line[9].lower() == 'true'
                data = np.loadtxt(rdf_file, unpack=True)
            except: 
                pass

            if not np.isclose(r_max, r_max_file) or \
                not np.isclose(dr, dr_file) or \
                not np.isclose(fraction_eq, fraction_eq_file) or \
                normalize != normalize_file:
                if verbose:                
                    print(f"Calculating RDF for simulation #{sim.metadata['simulation_number']}...")
                    print(f"  Parameters:  r_max={r_max}, dr={dr}, fraction={fraction_eq}, normalize={normalize}")
                    print(f"  Cached vals: r_max={r_max_file}, dr={dr_file}, fraction={fraction_eq_file}, normalize={normalize_file}")

                # Recalculate reference and ensemble RDFs for equilibrium fraction
                data = sim.get_ensemble_rdf( r_max=r_max, dr=dr, fraction=fraction_eq, normalize=normalize)
                # Save RDF data to file
                try:
                    np.savetxt(rdf_file, 
                        np.column_stack(data), 
                        header=f'  Parameters: r_max= {r_max} dr= {dr} fraction= {fraction_eq} normalize= {normalize}\n'
                                +'r_Angstrom g_r g_std g_ref_r' if normalize else 'r_Angstrom g_r g_std')
                except: 
                    pass

            else:
                if verbose:                
                    print(f"Using cached RDF for simulation #{sim.metadata['simulation_number']}...")
                    print(f"Parameters:  r_max={r_max}, dr={dr}, fraction={fraction_eq}, normalize={normalize}")
                
            results.append(data)

        return results


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
            accessibility is ndarray of vacant neighbor counts (0 to max_coordination).
            frequency is ndarray of frequencies for each accessibility level.
            frequency_std is ndarray of standard deviations across trajectories.
        """

        results = []
        for sim in self.simulations:

            if self.acc_file_sfx is not None:
                acc_file = Path(self.data_path) / self.results_dir / \
                           f"{sim.metadata['simulation_number']}_{self.acc_file_sfx}"
            else:
                acc_file = None

            # Try to load accessibility data from file
            acc_file_exists = False
            if acc_file is not None:
                try:
                    data = np.loadtxt(acc_file, skiprows=1, unpack=True)
                    acc_file_exists = True
                    if verbose:
                        print(f"Loaded cached accessibility for simulation #{sim.metadata['simulation_number']}")
                except:
                    pass

            if verbose and not acc_file_exists:
                print(f"Calculating accessibility for simulation #{sim.metadata['simulation_number']} ...")

            if not acc_file_exists:
                # Recalculate ensemble accessibility
                data = sim.get_ensemble_accessibility()
                
                # Save accessibility data to file
                if acc_file is not None:
                    try:
                        accessibility, frequency, frequency_std = data
                        # Save as columns: accessibility_level, frequency, frequency_std
                        np.savetxt(acc_file, 
                            np.column_stack([accessibility, frequency, frequency_std]), 
                            header='accessibility_level frequency frequency_std',
                            fmt='%d %f %f')
                        if verbose:
                            print(f"Cached accessibility to {acc_file.name}")
                    except Exception as e:
                        if verbose:
                            print(f"Warning: Could not save accessibility cache: {e}")
                
            results.append(data)

        return results


    def load(self, cache=None, workers=mp.cpu_count(), simulations=None, verbose=False):
        """
        Load data for simulations in the set.
        
        Parameters
        ----------
        cache : str or None, default None
            If str, cache to a specified file format. If None, do not use caching.
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
            sim = Simulation(sim_folder, metadata=md, log_file=self.log_file, results_dirname=self.results_dir)
            sim.load(cache=cache, workers=workers, verbose=verbose)  # Load simulation data
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


