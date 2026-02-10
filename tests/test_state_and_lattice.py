import numpy as np
from nat_zacros.lattice import Lattice
from nat_zacros.state import State


def test_state_basic_operations():
    lat = Lattice()  # default small lattice
    # synthesize a small 3-site lattice for testing
    lat.coordinates = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 0.866]])
    lat.site_coordinations = np.array([2, 2, 2])
    lat.site_nns = [np.array([1,2]), np.array([0,2]), np.array([0,1])]

    st = State(lat)
    # check initial coverage
    assert st.get_coverage() == 0.0
    assert st.n_ads() == 0

    # set occupations
    st.occupation[0] = 1
    st.occupation[2] = 1
    st.ads_ids[0] = 1
    st.ads_ids[2] = 1

    assert set(st.get_occupied_sites()) == {0, 2}
    assert set(st.get_empty_sites()) == {1}
    coords = st.get_occupied_coords()
    assert coords.shape[0] == 2
    assert st.n_ads() == 2
