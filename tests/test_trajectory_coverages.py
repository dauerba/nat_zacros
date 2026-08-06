import numpy as np

from nat_zacros.trajectory import Trajectory


class FakeState:
    def __init__(self, coverages):
        self._coverages = np.asarray(coverages, dtype=float)

    def get_coverages(self, atoms_per_uc=1):
        return self._coverages


def test_get_coverages_vs_time_returns_array_of_coverages(tmp_path):
    traj = Trajectory(dir=tmp_path, lattice=object(), metadata={"lattice_dimensions": (1, 1)})
    traj.states = [FakeState([0.5, 0.3, 0.2]), FakeState([0.4, 0.4, 0.2])]
    traj.initial_state = None
    traj.state_deltas = []
    traj.times = np.array([0.0, 1.0], dtype=float)

    times, coverages = traj.get_coverages_vs_time(atoms_per_uc=1)

    assert np.array_equal(times, traj.times)
    assert coverages.shape == (2, 3)
    assert np.allclose(coverages[0], [0.5, 0.3, 0.2])
    assert np.allclose(coverages[1], [0.4, 0.4, 0.2])
