"""
Module: nat_zacros
==================

This module provides classes for working with Zacros simulations:
- `lattice`: FCC(111) surface lattice
- `state`: Adsorbate configuration on the lattice
- `trajectory`: Sequence of states over time
"""

try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:
    from importlib_metadata import version, PackageNotFoundError

try:
    __version__ = version("nat-zacros")
except PackageNotFoundError:
    __version__ = "unknown"


"""
Performance Optimization Guide
-------------------------------
This module includes several performance optimizations for RDF and trajectory analysis.
Understanding when to use each approach is critical for optimal performance.

**Key Optimizations (Recommended for all use cases):**

1.  **Vectorized Distance Calculations** (50-100x speedup)
    - Automatically enabled by default in trajectory.get_rdf(vectorized=True)
    - Uses NumPy broadcasting to compute all pairwise distances at once
    - Replaces nested Python loops with compiled NumPy operations

2.  **Parallel Loading** (5-10x speedup)
    - Use load_trajectories_parallel() for loading multiple trajectories
    - Each trajectory reads its own file → good I/O parallelism
   
3.  **Binary Caching with pickle** (100x speedup for repeated analysis)
    - Save parsed trajectories to pickle files after first load
    - Subsequent loads read binary instead of parsing text (0.5s vs 60s) 
      
    - Example usage:
        cache_file = 'trajectories_eq.pkl'
        if not Path(cache_file).exists():
            trajs = [load and parse trajectories]
            with open(cache_file, 'wb') as f:
                pickle.dump(trajs, f)
        else:
            with open(cache_file, 'rb') as f:
                trajs = pickle.load(f)

    - dja comment 2025-12-19: 
        Timing needs to be rechecked. It is based on timing
        without parallelization of loading. 
        With vectorized distance calcs and parallel loading
        the time for loading + RDF calculation is now only ~2s
        ---------------------------------------------
        CONSIDER REMOVING CACHING (for simplicity).
        ---------------------------------------------

4.  **Parallel RDF Computation** - compute_rdf_parallel(), compute_rdf_parallel_states()
    
    - Parallelization of Computations was tried and not found effective for our current use case.
    - With vectorization, RDF computation is very fast (~2s for 10 trajectories)
    - Parallelization overhead (process spawn, pickle, IPC) is ~2-4 seconds
    - Therefore parallelization will only be beneficial when computation time >> 20-30 seconds
    - Use cases where there may be benefit:
        * compute_rdf_parallel(): >50 trajectories, or very long trajectories
        * compute_rdf_parallel_states(): >100 trajectories with many states each

    - For typical use (10-20 trajectories): sequential computation is FASTER

-------------------------------------------------------------------
        Performance Analysis Summary
-------------------------------------------------------------------
TODO: recheck and revise these benchmarks and computer system and problem size parameters
Performance Benchmark (typical system: 10 trajectories, ~100 states each, 14 cores):
    - Sequential loading: 64s
    - Parallel loading: 6-10s
    - RDF computation (sequential, vectorized): 2s
    - RDF computation (parallel): 3-5s (slower due to overhead!)
"""

# Import all public classes and functions
import json
from pathlib import Path
from .lattice import Lattice
from .state import State
from .trajectory import Trajectory
from .simulation import Simulation
from .simulation_set import SimulationSet

# Import parallel and RDF functions from data_processing
# from .data_processing import (
#    compute_rdf_parallel,
#    compute_rdf_parallel_states,
# )

# Define what gets imported with "from nat_zacros import *"
__all__ = [
    'Lattice',
    'State',
    'Trajectory',
    'Simulation',
    'SimulationSet'
]

import os
try:
    from setuptools_scm import get_version
    pkg_root = os.path.dirname(os.path.dirname(__file__))
    __version__ = get_version(root=pkg_root)
except Exception:
    __version__ = "unknown"
__author__ = 'akandra, dauerba'
