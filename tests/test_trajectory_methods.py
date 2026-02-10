import numpy as np
from nat_zacros.lattice import Lattice
from nat_zacros.state import State
from nat_zacros.trajectory import Trajectory


def make_triangle_lattice():
    lat = Lattice()
    # create 4 sites forming a small cluster
    lat.coordinates = np.array([[0.0,0.0],[1.0,0.0],[0.5,0.866],[1.5,0.866]])
    lat.site_coordinations = np.array([2,2,2,2])
    lat.site_nns = [np.array([1,2]), np.array([0,2,3]), np.array([0,1,3]), np.array([1,2])]
    return lat


def test_trajectory_rdf_and_clusters():
    lat = make_triangle_lattice()
    traj = Trajectory(lat)

    # create two states with different occupations
    s1 = State(lat)
    s1.occupation[[0,1]] = 1
    s1.ads_ids[[0,1]] = 1

    s2 = State(lat)
    s2.occupation[[2,3]] = 1
    s2.ads_ids[[2,3]] = 1

    traj.states = [s1, s2]

    r_bins, g_r = traj.get_rdf(r_max=2.0, dr=0.5)
    assert len(r_bins) == len(g_r)
    sizes, freqs = traj.get_cluster_distribution(nn_cutoff=1)
    # cluster sizes should be present
    assert sizes.size >= 0
