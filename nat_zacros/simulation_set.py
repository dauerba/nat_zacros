"""
SimulationSet class for managing Zacros simulation sets
We refer to simulations as a set of simulation runs when they share the same log file.
We refer to simulations as a run when they share the same run folder.

This module provides a high-level interface for loading, caching, and analyzing
collections of runs from a Zacros simulation set.
"""

import json
#import pickle
#import numpy as np
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
    set_dir  : Path or str
        Directory containing the log file, runs, and results)
    runs_dir : str
        subdirectory containing simulation runs (default: 'jobs')
    results_dir : str
        subdirectory containing simulation results (default: 'results')
    log_file : str  
        name of the log file (default: 'jobs.log')
    metadata : list of dictionaries
        Simulation metadata (temperature, coverage, interactions, etc.)
    
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
        set_dir : str or Path
            Path to simulation set directory (e.g., 'fn_3leed')
            This directory should contain jobs.log and the runs subdirectory
        runs_dir : str, optional
            Name of the subdirectory containing simulation runs (default: 'jobs')
        results_dir : str, optional
            Name of the subdirectory for storing results (default: 'results')
        log_file : str, optional
            Name of the log file (default: 'jobs.log')
        """
        
        self.set_dir     = Path(set_dir)
        self.run_dir     = runs_dir
        self.results_dir = results_dir
        self.log_file    = log_file
        
        # Validate the set directory exists
        if not self.set_dir.exists():
            raise FileNotFoundError(f"Set directory not found: {self.set_dir}")
        
        # Load metadata list from log file
        self._load_metadata()


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

    
    def __len__(self):
        """Return number of runs."""
        return len(self.metadata)
