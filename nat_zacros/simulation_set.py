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
import numpy as np
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
    trimming_method : str
        Method for trimming trajectories when loading (None, 'fit', 'avg')
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
        self.trimming_method    = None              # default trimming method
        self.use_cache          = False             # default caching behavior
        self.verbose            = False             # default verbosity
        

        self._load_metadata()
        self.simulations        = []   # initialize simulations list
        self.fractions_eq       = {}  # initialize equilibration fractions dictionary

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


    def clear_cache(self, trajectories=True, gref=False):
        """
        Clear cached data for all simulation runs in the set.
        
        Parameters
        ----------
        trajectories : bool, default True
            If True, clear cached trajectory data.
        gref : bool, default False
            If True, clear cached gref data.
        """
        for md in self.metadata:
            run_folder = self.set_dir / self.runs_dir / f"{md['run_number']}"
            sim = Simulation(run_folder, metadata=md, log_file=self.log_file, results_dirname=self.results_dir)
            sim.clear_cache(trajectories=trajectories, gref=gref)

    def load_energy(self, use_cache=True, parallel=False, verbose=False):
        """
        Load time and energy data for all simulation runs in the set with caching support.
        
        Parameters
        ----------
        use_cache : bool, default True
            If True, load from cache file if available, otherwise parse and cache.
            If False, always parse from history_output.txt files.
        parallel : bool, default True
            If True, use parallel loading (recommended for full-state data).
            If False, use sequential loading.
        verbose : bool, default False
            If True, print detailed loading information.

        """

        if len(self.simulations) > 0:
            print("load_energy() call is ignored: simulations set is not empty.")
            return

        for md in self.metadata:
            run_folder = self.set_dir / self.runs_dir / f"{md['run_number']}"
            sim = Simulation(run_folder, metadata=md, log_file=self.log_file, results_dirname=self.results_dir)
            sim.load(use_cache=use_cache, parallel=parallel, energy_only=True, verbose=verbose)  # Load energy only from simulation data
            self.simulations.append(sim)


    def load(self):
        """
        Load all simulation runs in the set.
        
        """
        self.simulations = []

        for md in self.metadata:
            run_folder = self.set_dir / self.runs_dir / f"{md['run_number']}"
            sim = Simulation(run_folder, metadata=md, log_file=self.log_file, results_dirname=self.results_dir)
            if self.trimming_method == 'fraction':
                sim._fraction_loaded = self.fractions_eq[md['run_number']] if self.fractions_eq[md['run_number']] is not None else 1.0
                sim.load(use_cache=self.use_cache, parallel=False, energy_only=True, verbose=self.verbose)  # Load energy only from simulation data



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
                fraction = self.fractions_eq[sim.metadata["run_number"]]
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
                
            ax.plot(times_plot, energies, marker='o', linestyle='-', markersize=2)
            ax.set_xlabel(x_label)
            ax.set_ylabel('Energy (eV)')
            ax.set_title(f'Run #{sim.metadata["run_number"]}' 
                        fr'  $T={sim.metadata["temperature"]}$ K, $\theta={sim.metadata["coverage"]:.3f}$'
                        f' ({fraction*100:.0f}%)',
                        fontsize = title_fontsize)
            ax.grid()

            # Shade equilibrium region
            eq_idx = int((1 - fraction)*(len(times) - 1))
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
