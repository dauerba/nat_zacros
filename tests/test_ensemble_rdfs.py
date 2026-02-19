import json
from pathlib import Path
import numpy as np
from nat_zacros.simulation_set import SimulationSet


def test_get_ensemble_rdfs_returns_4tuple_and_writes_gref(tmp_path):
    """Ensure normalized RDF returns (r, g_avg, g_std, g_ref) and writes g_ref.dat"""
    # minimal jobs.log with one simulation
    jobs_log = tmp_path / 'jobs.log'
    entries = [[1, "job1", [4, 4], [2], 300, ["label", "E1"]]]
    with open(jobs_log, 'w') as f:
        f.write("header\n")
        for e in entries:
            f.write(json.dumps(e) + "\n")

    # prepare per-simulation folder
    sim_dir = tmp_path / 'data' / '1'
    sim_dir.mkdir(parents=True, exist_ok=True)

    simset = SimulationSet(tmp_path)
    # set a valid equilibration fraction so calculation proceeds
    simset.fractions_eq[1] = 1.0

    # Dummy Simulation that returns a 4-tuple when normalize=True
    class DummySim:
        def __init__(self, sim_num):
            self.metadata = {'simulation_number': sim_num}

        def get_ensemble_rdf(self, r_max, dr, fraction, normalize):
            r = np.array([0.0, 1.0, 2.0])
            g_avg = np.array([1.0, 1.1, 0.9])
            g_std = np.array([0.0, 0.05, 0.05])
            g_ref = np.array([1.0, 1.0, 1.0])
            if normalize:
                return (r, g_avg, g_std, g_ref)
            return (r, g_avg, g_std)

    simset.simulations = [DummySim(1)]

    # Run the method under test
    results = simset.get_ensemble_rdfs(normalize=True, verbose=False)

    # Basic checks on returned structure
    assert isinstance(results, list)
    assert len(results) == 1
    entry = results[0]
    assert entry is not None
    assert isinstance(entry, tuple)
    assert len(entry) == 4
    for arr in entry:
        assert isinstance(arr, np.ndarray)

    # Check files on disk (per-simulation folder)
    sim_folder = tmp_path / 'data' / '1'
    gref_file = sim_folder / 'g_ref.dat'
    rdf_file = sim_folder / 'rdf.dat'

    assert gref_file.exists()
    assert rdf_file.exists()

    # rdf.dat should contain 4 columns when normalized (r, g_avg, g_std, g_ref)
    data = np.loadtxt(rdf_file, unpack=True)
    assert data.shape[0] == 4

    # Re-run to ensure cached read path also returns 4-tuple
    results2 = simset.get_ensemble_rdfs(normalize=True, verbose=False)
    entry2 = results2[0]
    assert isinstance(entry2, tuple)
    assert len(entry2) == 4
