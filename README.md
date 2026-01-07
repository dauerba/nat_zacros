# nat_zacros

Python package for analyzing Zacros kinetic Monte Carlo simulations of surface reactions.

## Features

- **Lattice Geometry**: FCC(111) surface lattice with periodic boundary conditions
- **State Management**: Parse and manipulate adsorbate configurations
- **Trajectory Analysis**: Load and analyze KMC trajectories
- **RDF Calculations**: Fast vectorized radial distribution functions with PBC
- **Parallel Processing**: Efficient loading and analysis of multiple trajectories

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
```
  If the virtual environment is set by conda then
follow recommendations from: https://github.com/conda/conda/issues/5861
```bash
git clone https://github.com/dauerba/nat_zacros.git
cd nat_zacros
VIRTUAL_ENV = $CONDA_PREFIX
pip install --src $VIRTUAL_ENV/src -e .
```

## Quick Start

```python
from nat_zacros import lattice, trajectory

# Load lattice and trajectory
lat = lattice(dirname='zacros_run_0080')
traj = trajectory(lat, dirname='zacros_run_0080')

# Load trajectory data (equilibrated portion)
traj.load_trajectory(energy_only=True)
traj.load_equilibrated_states(fraction=0.5)

# Compute radial distribution function
r, g = traj.get_rdf(r_max=40.0, dr=0.1)

# Plot
import matplotlib.pyplot as plt
plt.plot(r, g)
plt.xlabel('Distance (Å)')
plt.ylabel('g(r)')
plt.show()
```

## Performance Optimization

The package includes several performance optimizations:

1. **Vectorized distance calculations** (50-100x speedup)
2. **Parallel trajectory loading** (5-10x speedup)
3. **Binary caching with pickle** (100x speedup for repeated analysis)

See module docstring for detailed performance guide.

## Requirements

- Python >= 3.8
- NumPy >= 1.20
- SciPy >= 1.7
- Optional: tqdm (for progress bars in parallel functions)

## Project Structure

```
nat_zacros/
├── nat_zacros/
│   ├── __init__.py       # Package entry point
│   ├── lattice.py        # FCC(111) lattice geometry
│   ├── state.py          # Adsorbate configurations
│   ├── trajectory.py     # State sequences and RDF analysis
│   └── rdf.py            # Parallel computation utilities
├── tests/                # Unit tests (future)
├── examples/             # Example notebooks (future)
├── setup.py              # Package installation
└── README.md
```

## Contributing

This package is part of the O_Pt111 project for studying oxygen adsorption on Pt(111) surfaces.

**Current Maintainers:**
- Primary: akandra (pending repository transfer)
- Developer: dauerba (refactoring and packaging)

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
