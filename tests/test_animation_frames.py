import importlib.util
import sys
import types
from pathlib import Path


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_first_five_frames():
    repo = Path(__file__).resolve().parents[1]
    traj_path = repo / 'nat_zacros' / 'trajectory.py'
    state_path = repo / 'nat_zacros' / 'state.py'

    # Create a fake package so relative imports inside modules work
    pkg = types.ModuleType('nat_zacros')
    pkg.__path__ = [str(repo / 'nat_zacros')]
    sys.modules['nat_zacros'] = pkg

    state_mod = load_module('nat_zacros.state', state_path)
    sys.modules['nat_zacros.state'] = state_mod
    traj_mod = load_module('nat_zacros.trajectory', traj_path)

    Trajectory = traj_mod.Trajectory
    State = state_mod.State

    # Minimal fake lattice with attributes used by State/Trajectory
    class FakeLattice:
        def __init__(self, nsites=9):
            import numpy as np
            self._n = nsites
            self.cell_vectors = [np.array([1.0,0.0]), np.array([0.5, np.sqrt(3)/2])]
            self.size = (3,3)
            self.n_cell_sites = 1
            self.site_coordinations = np.array([3]*nsites)
            self.site_type_names = ['top']
            self.site_types = np.ones(nsites, dtype=int)
            self.coordinates = np.column_stack([np.arange(nsites,dtype=float), np.zeros(nsites)])
            self.site_nns = [list(range(nsites)) for _ in range(nsites)]
        def __len__(self):
            return self._n
        def minimum_image_distance(self, a,b):
            import numpy as np
            return np.linalg.norm(np.array(a)-np.array(b))

    # Monkeypatch State.plot to be a no-op for headless test
    State.plot = lambda self, **kwargs: None

    metadata = {'surf_species_names':['A*'], 'lattice_dimensions':None, 'n_adsorbates':0, 'temperature':300, 'energy_terms':[]}
    lat = FakeLattice(nsites=5)
    traj = Trajectory('.', lat, metadata)
    traj.states = [State(lat, metadata) for _ in range(5)]
    import numpy as np
    traj.times = np.linspace(0,1,5)
    traj.energies = np.zeros(5)

    anim = traj.animation(frames=range(5))
    seq = list(anim.new_frame_seq())
    assert seq == list(range(5))
