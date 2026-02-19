import json
from pathlib import Path
from nat_zacros.simulation_set import SimulationSet


def _make_simset_tree(tmp_path):
    """Create a minimal simulation-set layout with cache/result files."""
    root = tmp_path
    # jobs.log (header + two simulation entries)
    jobs_log = root / 'jobs.log'
    entries = [
        [1, "job1", [4, 4], [2], 300, ["label", "E1"]],
        [2, "job2", [4, 4], [2], 300, ["label", "E1"]],
    ]
    with open(jobs_log, 'w') as f:
        f.write("header\n")
        for e in entries:
            f.write(json.dumps(e) + "\n")

    # Create simulation directories and per-trajectory cache files
    for sim in (1, 2):
        sim_dir = root / 'data' / str(sim)
        traj0 = sim_dir / 'traj_0'
        traj1 = sim_dir / 'traj_1'
        traj0.mkdir(parents=True, exist_ok=True)
        traj1.mkdir(parents=True, exist_ok=True)
        # create traj.pkl files (simulate cached trajectories)
        (traj0 / 'traj.pkl').write_bytes(b'cached')
        (traj1 / 'traj.pkl').write_bytes(b'cached')

    # Sim-specific results files (new layout)
    for sim in (1, 2):
        sim_dir = root / 'data' / str(sim)
        sim_dir.mkdir(parents=True, exist_ok=True)
        # aggregated trajs cache per-sim
        (sim_dir / 'trajs_eq.pkl').write_bytes(b'agg')

    # user-visible results for simulation 1
    sim1 = root / 'data' / '1'
    (sim1 / 'energy.dat').write_text('dummy')
    (sim1 / 'rdf.dat').write_text('dummy')
    (sim1 / 'accessibility.dat').write_text('dummy')
    # g_ref cache (text format)
    (sim1 / 'g_ref.dat').write_text('0.0 1.0')

    return root


def test_clear_traj_cache_removes_traj_and_agg(tmp_path):
    root = _make_simset_tree(tmp_path)
    simset = SimulationSet(root)

    # Sanity: files exist before clearing
    assert (root / 'data' / '1' / 'traj_0' / 'traj.pkl').exists()
    assert (root / 'data' / '1' / 'trajs_eq.pkl').exists()

    # Clear only traj caches
    simset.clear_traj_cache(simulations=[1], verbose=False)

    # Per-trajectory cache for sim 1 removed, sim 2 still present
    assert not (root / 'data' / '1' / 'traj_0' / 'traj.pkl').exists()
    assert (root / 'data' / '2' / 'traj_0' / 'traj.pkl').exists()

    # Aggregated cache for sim 1 removed, sim 2 still present
    assert not (root / 'data' / '1' / 'trajs_eq.pkl').exists()
    assert (root / 'data' / '2' / 'trajs_eq.pkl').exists()


def test_clear_results_removes_user_files_and_gref(tmp_path):
    root = _make_simset_tree(tmp_path)
    simset = SimulationSet(root)

    # Sanity: files exist
    assert (root / 'data' / '1' / 'energy.dat').exists()
    assert (root / 'data' / '1' / 'rdf.dat').exists()
    assert (root / 'data' / '1' / 'accessibility.dat').exists()
    assert (root / 'data' / '1' / 'g_ref.dat').exists()

    # Clear specific results (energy + gref)
    simset.clear_results(target=['energy', 'gref'], simulations=[1], verbose=False)

    # energy for sim 1 removed, other files untouched
    assert not (root / 'data' / '1' / 'energy.dat').exists()
    assert (root / 'data' / '1' / 'rdf.dat').exists()
    assert (root / 'data' / '1' / 'accessibility.dat').exists()

    # gref removed
    assert not (root / 'data' / '1' / 'g_ref.dat').exists()

    # Now clear all remaining results
    simset.clear_results(target='all', verbose=False)
    assert not (root / 'data' / '1' / 'rdf.dat').exists()
    assert not (root / 'data' / '1' / 'accessibility.dat').exists()
    # aggregated trajs are NOT part of clear_results and should still exist for sim 2
    assert (root / 'data' / '2' / 'trajs_eq.pkl').exists()
