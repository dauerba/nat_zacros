# nat_zacros

> [!CAUTION]
> **EARLY DEVELOPMENT STAGE - DO NOT MODIFY DOCUMENTATION OR TESTS**
> This project is in active flux. Backward compatibility is NOT maintained.
> There are currently NO USERS. Do not spend time updating README, documentation, or unit tests.
> Focus solely on functional code development and bug fixes.

Python package for analyzing Zacros kinetic Monte Carlo simulations of surface reactions.

## Features  -- NEEDS TO BE REDONE

- Create submission scripts for Zacros Calculations and Analysis Code to Extract Results
- Implemented as a hiearchy of classes in Python
    - Lattice           geomentry of the lattice used in the simulation and pbc for distances
    - State             occupation of lattice sites by adsorbatetes
    - Trajectory        a sequence of states
        - methods to load trajectories from Zacros output files
    - Simulation        a set of trajecotries for given onditionos of temperature, coverage, ...
        - methods to calculate ensemble average propoerties like Energy, RDF, Cluster Size Distributions ...
    - Simulation_Set    a set of Simulation with high level interface
        - methods to implement high level interface to other classes, anallyze and plto results


## Installation

### From GitHub (development)
```bash
pip install git+https://github.com/dauerba/nat_zacros.git
```

### For local development
```bash
git clone https://github.com/dauerba/nat_zacros.git
cd nat_zacros
pip install -e .


## Quick Start   NEEDS REWRITE

```python
from nat_zacros import lattice, trajectory

show example of Simulation set class to produce energy, rdf, and accessibility distributions
```

## Performance Optimization

The package includes several performance optimizations:

1. **Vectorized distance calculations** (50-100x speedup)
2. **Parallel trajectory loading** (????x speedup)
3. **Binary caching with pickle** ( ????x speedup for repeated analysis)


See module docstring for detailed performance guide.

## Requirements
ADD PIP > xxx to support editable install
Consider more stringent requirements for Python
- Python >= 3.8
- NumPy >= 1.20
- SciPy >= 1.7
- Optional: tqdm (for progress bars in parallel functions)

## Project Structure

```
nat_zacros/
├── nat_zacros/
│   ├── __init__.py           # Package entry point
├── lattice.py                # Lattice geometry (currently only FCC111) is supported
├── state.py                  # Adsorbates configuration on the lattice
├── trajectory.py             # Sequence of states from a zacros kMC simulation
│   ├── simulation.py         # Simulation compriseing 1 or more trajectories with methods to 
|   |                         -  Load and cache trajectories, compute ensemble averaged properties (energy, RDF, ...)
│   ├── simulation_set.py
│
│
├── scripts/
├── tests/                    # Unit tests 
├── examples/                 # Example notebooks
│
├── pypoject.toml             # Package installation
└── README.md
```

## Contributing

This package is part of the O_Pt111 project for studying oxygen adsorption on Pt(111) surfaces.

**Current Maintainers:**
- Primary: akandra (pending repository transfer)
- Developer: dauerba (refactoring and packaging)

## Recent API changes (unreleased)

- Added `SimulationSet.clear_traj_cache()` — removes internal trajectory cache files (`traj.pkl`) used only to speed loading.
- Added `SimulationSet.clear_results()` — removes user-visible result files (`energy`, `rdf`, `accessibility`, `gref`).
- Added `Simulation.clear_traj_cache()` / `Simulation.clear_results_path()` (path-level helpers) and instance wrapper `Simulation.clear_traj_cache()` / `Simulation.clear_results()`; `SimulationSet` delegates to these helpers.
- Removed `SimulationSet.clear_cache()` (was deprecated). Use `clear_traj_cache()` or `clear_results()` instead.
- Unit tests added for the new cache/result-clearing behavior.
- Updated example notebooks to use `clear_traj_cache()` / `clear_results()` instead of the deprecated API.

Example (per-simulation helpers)

```python
from nat_zacros import Simulation, SimulationSet

# Path-level (no Simulation instance needed)
Simulation.clear_traj_cache('/path/to/run', traj_dir_pfx='traj')
Simulation.clear_results_path('/path/to/run', target=['energy','gref'])

# Instance-level
sim = Simulation('/path/to/run', metadata={'lattice_dimensions':[4,4], 'n_adsorbates':2, 'temperature':300, 'energy_terms':['label','E1']})
sim.clear_traj_cache()
sim.clear_results(target='all')

# Existing SimulationSet API still works (delegates to Simulation helpers)
simset = SimulationSet('/path/to/'); simset.clear_traj_cache(simulations=[1,2])
```

## License

MIT License - see LICENSE file for details

## Citation

If you use this package in your research, please cite:
```
[Add citation information when available]
```

## Related Projects

- [Zacros](http://zacros.org/) - Kinetic Monte Carlo software for catalysis
- [O_Pt111](https://github.com/akandra/O_Pt111) - Parent project for Pt(111) simulations
