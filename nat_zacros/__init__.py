"""
Module: nat_zacros
==================

This module provides a heiracrhical set of classes for working with Zacros simulations:
- `lattice`: FCC(111) surface lattice
- `state`: Adsorbate configuration on the lattice
- `trajectory`: Sequence of states over time
- `simulation`: Single Zacros simulation run
- `simulation_set`: Collection of multiple simulations

Support is provided for loading and parsing Zacros output files, analyzing trajectories, and computing properties 
such as radial distribution functions (RDFs), cluster size distrbutions, and a metric of the accessibility of adsorbates 
inside clusters.
"""

def _custom_local_scheme(version):
    """Custom local scheme to format timestamp as YYYY-MM-DD-HH:MM"""
    from datetime import datetime
    if version.dirty:
        timestamp = datetime.now().strftime("%Y-%m-%d-%H:%M")
        return f"+dirty.d{timestamp}"
    elif version.distance:
        timestamp = datetime.now().strftime("%Y-%m-%d-%H:%M")
        return f"+{version.node}.d{timestamp}"
    else:
        return ""

try:
    from ._version import version as __version__
except ImportError:
    try:
        from setuptools_scm import get_version
        # Use custom local scheme for formatted timestamp
        __version__ = get_version(root='..', relative_to=__file__, local_scheme=_custom_local_scheme)
    except Exception:
        __version__ = "unknown"

"""
Performance Optimization
-------------------------------
1.  **Parallel Loading** (5-10x speedup)
    - Use load_trajectories_parallel() for loading multiple trajectories
    - Each trajectory reads its own file → good I/O parallelism
   
2.  **Binary Caching ** 
    - Save parsed trajectories to cache files after first load
    - Subsequent loads read binary instead of parsing text 
    - (0.xx s vs yys for 10 trajectories with 200 states each) 

3.  **Vectorized Distance Calculations** (50-100x speedup)
    - Uses NumPy broadcasting to compute all pairwise distances at once
    - Replaces nested Python loops with compiled NumPy operations

4.  **Parallel RDF Computation** 
    - Parallelization of Computations was tried and not found effective for our current use case.
    - With vectorization, RDF computation is very fast (~2s for 10 trajectories)
    - Parallelization overhead (process spawn, pickle, IPC) is ~2-4 seconds
    - Therefore parallelization will only be beneficial when computation time >> 10 seconds
    - Use cases where there may be benefit:
        * compute_rdf_parallel(): >50 trajectories, or very long trajectories
        * compute_rdf_parallel_states(): >100 trajectories with many states each
    - For typical use (10-20 trajectories): sequential computation is FASTER

-------------------------------------------------------------------
        Performance Analysis Summary
-------------------------------------------------------------------
TODO: recheck and revise these benchmarks and computer system and problem size parameters
Performance Benchmark (typical system: 10 trajectories, ~100 states each, 14 cores):
    - Sequential loading: ??? 64s
    - Parallel loading: 6-10s
    - RDF computation (sequential, vectorized): 2s
    - RDF computation (parallel): 3-5s (slower due to overhead!)
"""

# Import all public classes and functions
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

__author__ = 'akandra, dauerba'
