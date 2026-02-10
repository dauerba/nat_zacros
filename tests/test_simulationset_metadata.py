import json
from pathlib import Path
from nat_zacros.simulation_set import SimulationSet


def test_simulationset_parses_jobs_log(tmp_path):
    # create a minimal jobs.log
    jobs = [
        [1, 'job1', [2,2], [4], 300.0, ['meta', 'info']],
        [2, 'job2', [1,1], [1], 500.0, ['meta', 'info']]
    ]
    log_path = tmp_path / 'jobs.log'
    with open(log_path, 'w') as f:
        f.write('header\n')
        for entry in jobs:
            f.write(json.dumps(entry) + '\n')

    ss = SimulationSet(tmp_path)
    assert isinstance(ss.metadata, list)
    assert len(ss.metadata) == 2
    assert ss.metadata[0]['simulation_number'] == 1
