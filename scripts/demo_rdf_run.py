import json
import tempfile
from pathlib import Path
import numpy as np
from nat_zacros.simulation_set import SimulationSet

with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    jobs_log = td / 'jobs.log'
    entries = [[1, "job1", [4, 4], [2], 300, ["label", "E1"]]]
    with open(jobs_log, 'w') as f:
        f.write("header\n")
        for e in entries:
            f.write(json.dumps(e) + "\n")

    sim_dir = td / 'data' / '1'
    sim_dir.mkdir(parents=True, exist_ok=True)

    simset = SimulationSet(td)
    simset.fractions_eq[1] = 1.0

    class DummySim:
        def __init__(self, sim_num):
            self.metadata = {'simulation_number': sim_num}
        def get_ensemble_rdf(self, r_max, dr, fraction, normalize):
            r = np.array([0.0, 1.0, 2.0])
            g_avg = np.array([1.0, 1.1, 0.9])
            g_std = np.array([0.0, 0.05, 0.05])
            if normalize:
                g_ref = np.array([1.0, 1.0, 1.0])
                return (r, g_avg, g_std, g_ref)
            return (r, g_avg, g_std)

    simset.simulations = [DummySim(1)]

    res = simset.get_ensemble_rdfs(normalize=True, verbose=True)

    print('Returned list length:', len(res))
    entry = res[0]
    print('Entry is None?:', entry is None)
    if entry is not None:
        print('Is tuple?', isinstance(entry, tuple))
        print('Tuple length:', len(entry))
        for i, arr in enumerate(entry):
            print(i, type(arr), getattr(arr, 'shape', None))

    gref_path = sim_dir / 'g_ref.dat'
    print('g_ref exists on disk?', gref_path.exists())
    if gref_path.exists():
        print('g_ref content:\n', gref_path.read_text())
