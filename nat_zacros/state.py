"""
State class for adsorbate configurations on lattice.

This module provides the `state` class for representing snapshots of adsorbate
configurations on a surface lattice from Zacros simulations.
"""

import numpy as np
from pathlib import Path


class state:
    """
    Adsorbate configuration on a surface lattice.
    
    Represents a snapshot of which sites are occupied by which species
    at a particular moment in a Zacros simulation.
    
    Attributes
    ----------
    lattice : lattice object
        Reference to the underlying surface lattice
    folder : Path or None
        Directory containing history_output.txt
    n_gas_species : int
        Number of gas-phase species
    gas_species_names : list of str
        Names of gas-phase species
    n_surf_species : int
        Number of surface-bound species
    surf_species_names : list of str
        Names of surface species
    surf_species_dent : list
        Denticity of each surface species
    ads_ids : ndarray, shape (N,)
        Species ID at each lattice site (0 = empty)
    occupation : ndarray, shape (N,)
        Occupation status of each site (0 = empty, >0 = occupied)
    dentation : ndarray, shape (N,)
        Denticity at each site
        
    Methods
    -------
    get_state(idx=0)
        Read configuration from history_output.txt
    get_coverage()
        Calculate fraction of occupied sites
    get_occupied_sites()
        Get indices of occupied sites
    get_empty_sites()
        Get indices of empty sites
    get_occupied_coords()
        Get Cartesian coordinates of occupied sites
    n_ads()
        Get total number of adsorbates
    """
    

    def __init__(self, lattice, dirname=None):
        """
        Initialize state object.
        
        Parameters
        ----------
        lattice : lattice object
            The surface lattice on which adsorbates are placed
        dirname : str or Path, optional
            Directory containing history_output.txt. If provided,
            will automatically load the first state.
        """
        self.folder = None

        # Store reference to lattice
        self.lattice = lattice

        # default is no species
        self.n_gas_species = 0
        self.gas_species_names = []
        self.n_surf_species = 0
        self.surf_species_names = []
        self.surf_species_dent = []

        # Arrays defining the adsorbed species on the lattice
        # with indices corresponding to lattice site indices starting at 1

        nsites = len(self.lattice) 
        self.ads_ids =    np.zeros(nsites, dtype=int)
        self.occupation = np.zeros(nsites, dtype=int)
        self.dentation =  np.zeros(nsites, dtype=int)

        if dirname is not None:
            self.folder = Path(dirname) 
            self.get_state()


    def get_state(self, idx=0):
        """
        Read configuration from history_output.txt file.
        
        Parameters
        ----------
        idx : int, optional
            Index of the configuration to read (default: 0 = first state)
            
        Notes
        -----
        Reads from history_output.txt which contains sequential snapshots.
        Each snapshot lists the occupation, species ID, and denticity for
        every lattice site.
        """

        
        # Read configuration from history_output.txt file

        self.folder = Path(self.folder)

        try:
            with open(self.folder / 'history_output.txt', 'r') as f:
                content = f.readlines()    

            nsites = len(self.lattice)
            for site in range(nsites):
                parts = content[7 + idx*(nsites+1) + site].split()
                self.ads_ids[site]    = int(parts[1])
                self.occupation[site] = int(parts[2])
                self.dentation[site]  = int(parts[3])
        except:
            print(f'cannot read history_output.txt from {str(self.folder)}')

    
    def get_coverage(self):
        """
        Calculate the coverage (fraction of occupied sites).
        
        Returns
        -------
        float
            Fraction of sites that are occupied (0.0 to 1.0)
        """
        return np.count_nonzero(self.occupation) / len(self.lattice)

    
    def get_occupied_sites(self):
        """
        Get indices of all occupied sites.
        
        Returns
        -------
        ndarray
            Array of site indices where occupation > 0
        """
        return np.where(self.occupation > 0)[0]

    
    def get_empty_sites(self):
        """
        Get indices of all empty sites.
        
        Returns
        -------
        ndarray
            Array of site indices where occupation == 0
        """
        return np.where(self.occupation == 0)[0]

    
    def get_occupied_coords(self):
        """
        Get Cartesian coordinates of occupied sites.
            
        Returns
        -------
        ndarray
            (N, 2) array of coordinates for occupied sites
        """
        mask = self.occupation > 0
        return self.lattice.coordinates[mask]


    def n_ads(self):
        """
        Get total number of adsorbates on the surface.
        
        Returns
        -------
        int
            Number of occupied sites
        """
        return np.count_nonzero(self.occupation)


    def __repr__(self):
        """String representation of state"""
        coverage = self.get_coverage()
        return f"state(nsites={len(self.lattice)}, n_adsorbates={self.n_ads()}, coverage={coverage:.3f})"
