import json
import pytest
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
        # per-package design: do not create an aggregated `trajs_eq.pkl` cache here

    # user-visible results for simulation 1
    sim1 = root / 'data' / '1'
    (sim1 / 'energy.dat').write_text('dummy')
    (sim1 / 'rdf.dat').write_text('dummy')
    (sim1 / 'accessibility.dat').write_text('dummy')
    # g_ref cache (text format)
    (sim1 / 'g_ref.dat').write_text('0.0 1.0')

    return root


def test_clear_traj_cache_removes_traj_cache_only(tmp_path):
    root = _make_simset_tree(tmp_path)
    simset = SimulationSet(root)

    # Sanity: per-trajectory cache files exist before clearing
    assert (root / 'data' / '1' / 'traj_0' / 'traj.pkl').exists()

    # Clear only traj caches
    simset.clear_traj_cache(simulations=[1], verbose=False)

    # Per-trajectory cache for sim 1 removed, sim 2 still present
    assert not (root / 'data' / '1' / 'traj_0' / 'traj.pkl').exists()
    assert (root / 'data' / '2' / 'traj_0' / 'traj.pkl').exists()


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
    # aggregated caches are not part of clear_results (package does not create them)


def test_get_ensemble_rdfs_warns_on_invalid_fraction(tmp_path):
    root = _make_simset_tree(tmp_path)
    simset = SimulationSet(root)

    # Attach a lightweight dummy Simulation to avoid heavy I/O; fraction remains None by default
    class DummySim:
        def __init__(self, sim_num):
            self.metadata = {'simulation_number': sim_num}
    simset.simulations = [DummySim(1)]

    # Expect a UserWarning and a None entry in results
    with pytest.warns(UserWarning, match=r"Skipping RDF for simulation #1"):
        rdfs = simset.get_ensemble_rdfs(verbose=False)
    assert rdfs[0] is None

    # ------------------------------------------------------------------
    # Instance-preference tests
    # ------------------------------------------------------------------
    class SpySim:
        def __init__(self, num):
            self.metadata = {'simulation_number': num}
            self.cache_cleared = False
            self.results_cleared = []
        def clear_traj_cache(self, verbose=False):
            self.cache_cleared = True
        def clear_results(self, target='all', results_files=None, verbose=False):
            self.results_cleared.append(target)

    simset.simulations = [SpySim(1)]
    # clearing cache should hit instance method
    simset.clear_traj_cache(simulations=[1], verbose=False)
    assert simset.simulations[0].cache_cleared
    # clearing results should hit instance method for each format
    simset.clear_results(target=['energy','rdf'], simulations=[1], verbose=False)
    assert 'energy' in simset.simulations[0].results_cleared
    assert 'rdf' in simset.simulations[0].results_cleared

    # Also test numeric invalid fraction (0.0)
    simset.fractions_eq[1] = 0.0
    with pytest.warns(UserWarning, match=r"Skipping RDF for simulation #1"):
        rdfs = simset.get_ensemble_rdfs(verbose=False)
    assert rdfs[0] is None


# ------------------------------------------------------------------
# Tests for Simulation-level helpers
# ------------------------------------------------------------------

def test_simulation_clear_traj_cache_and_instance(tmp_path):
    sim_dir = tmp_path / 'simA'
    traj0 = sim_dir / 'traj_0'
    traj0.mkdir(parents=True, exist_ok=True)
    (traj0 / 'traj.pkl').write_bytes(b'cached')

    from nat_zacros.simulation import Simulation

    # Static/path-level helper (new name)
    Simulation.clear_traj_cache(sim_dir, traj_dir_pfx='traj', verbose=False)
    assert not (traj0 / 'traj.pkl').exists()


    # Recreate and test instance wrapper
    (traj0 / 'traj.pkl').write_bytes(b'cached')
    sim = Simulation(sim_dir, metadata={'lattice_dimensions': [4,4], 'n_adsorbates': 2, 'temperature': 300, 'energy_terms': ['label','E1']})
    sim.clear_traj_cache(verbose=False)
    assert not (traj0 / 'traj.pkl').exists()


def test_simulation_clear_results_path_and_instance(tmp_path):
    sim_dir = tmp_path / 'simB'
    sim_dir.mkdir(parents=True, exist_ok=True)
    # create result files
    (sim_dir / 'energy.dat').write_text('e')
    (sim_dir / 'rdf.dat').write_text('r')
    (sim_dir / 'accessibility.dat').write_text('a')
    (sim_dir / 'g_ref.dat').write_text('g')

    from nat_zacros.simulation import Simulation

    # Path-level removal
    Simulation.clear_results_path(sim_dir, target=['energy', 'gref'], results_files=None, verbose=False)
    assert not (sim_dir / 'energy.dat').exists()
    assert not (sim_dir / 'g_ref.dat').exists()
    assert (sim_dir / 'rdf.dat').exists()

    # Instance wrapper
    sim = Simulation(sim_dir, metadata={'lattice_dimensions': [4,4], 'n_adsorbates': 2, 'temperature': 300, 'energy_terms': ['label','E1']})
    sim.clear_results(target='all', results_files=None, verbose=False)
    assert not (sim_dir / 'rdf.dat').exists()
    assert not (sim_dir / 'accessibility.dat').exists()


def test_simulation_get_g_ref(tmp_path):
    from nat_zacros.simulation import Simulation
    from nat_zacros.lattice import Lattice

    sim_dir = tmp_path / 'sim_gref'
    sim_dir.mkdir(parents=True, exist_ok=True)

    # Create a Simulation object (no traj dirs needed for this test)
    sim = Simulation(sim_dir, metadata={'lattice_dimensions': [4,4], 'n_adsorbates': 2, 'temperature': 300, 'energy_terms': ['label','E1']})
    # assign a default lattice so get_g_ref can run
    sim.lattice = Lattice()

    r_bins, g_ref = sim.get_g_ref(r_max=3.0, dr=0.5)
    assert r_bins.shape == g_ref.shape
    assert (g_ref >= 0).all()
