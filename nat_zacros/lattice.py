"""
Lattice class for FCC(111) surface geometry.

This module provides the `lattice` class for representing periodic FCC(111) 
surface lattices used in Zacros kinetic Monte Carlo simulations.
"""

import numpy as np
from pathlib import Path


class lattice:
    """
    FCC(111) surface lattice with periodic boundary conditions.
    
    Attributes
    ----------
    cell_vectors : ndarray, shape (2, 2)
        Supercell lattice vectors (unit_cell_vectors * size)
    coordinates : ndarray, shape (N, 2)
        Cartesian coordinates of all lattice sites
    folder : str (must be str for pickle compatibility) or None
        Directory containing lattice_input.dat and lattice_output.txt
    fractional_coordinates : ndarray
        Fractional coordinates of sites within unit cell
    is_defined : bool
        True if lattice parameters have been set
    n_cell_sites : int
        Number of sites per unit cell
    n_site_types : int
        Number of distinct site types
    neighboring_structure : list of tuples
        Neighboring connectivity pattern
    site_coordinations : ndarray
        Number of nearest neighbors for each site
    site_nns : list of ndarray
        Nearest neighbor indices for each site
    site_types : list or ndarray
        Site type for each site
    site_type_names : list of str
        Names of site types (e.g., ['fcc', 'hcp'])
    size : ndarray, shape (2,)
        Number of unit cell repetitions in each direction
    type : str
        Lattice type (e.g., 'periodic_cell')
    unit_cell_vectors : ndarray, shape (2, 2)
        Unit cell lattice vectors in Cartesian coordinates
        
    Methods
    -------
    apply_pbc(coords)
        Apply periodic boundary conditions to wrap coordinates
    get_cell_area()
        Calculate simulation cell area
    get_lattice()
        Read lattice definition from input/output files
    get_nn_distance(order=1)
        Get nearest neighbor distance for given order
    minimum_image_distance(coord1, coord2)
        Calculate minimum distance between two points with PBC
    pairwise_distances_pbc(coords)
        Calculate all pairwise distances with PBC (vectorized)
    """
    
    def __init__(self, dirname=None):
        """
        Initialize lattice object.
        
        Parameters
        ----------
        dirname : str or Path, optional
            Directory containing lattice_input.dat and lattice_output.txt.
            If None, creates default FCC(111) lattice with 1.0 Å spacing.
        """

        self.folder = None
        self.is_defined = False # True if lattice parameters have been set
        # default is an fcc(111) lattice with nearest-neighbor distance of 1.0 Angstrom
        # and a single fcc adsorption site per unit cell
        self.type = "periodic_cell"
        self.unit_cell_vectors =  np.array([ [1.0,          0.0], 
                                                [1/2, np.sqrt(3)/2] ])
        self.size = np.array([1, 1])
        self.cell_vectors = self.unit_cell_vectors * self.size[:, np.newaxis]
        self.n_site_types = 1
        self.site_type_names = ['fcc']
        self.n_cell_sites = 1
        self.site_types = ['fcc']
        self.fractional_coordinates = np.array([[1/3, 1/3]])
        self.neighboring_structure = [ (1,1, 'north'),
                                      (1,1, 'east'),
                                      (1,1, 'southeast') ]
        self.coordinates = np.array([[1/3, 1/3]])
        self.site_types = np.array([1])
        self.site_coordinations = np.array([6])
        # list of arrays of nearest neighbors for each lattice site
        self.site_nns = [ np.array([0,0,0,0,0,0]) ]

        if dirname is not None:
            self.folder = str(Path(dirname))
            self.get_lattice()

    def get_lattice(self):
        """
        Read lattice definition from input and output files.
        
        Reads:
        - lattice_input.dat: Unit cell definition and repeat pattern
        - lattice_output.txt: Full supercell coordinates and connectivity
        
        Raises
        ------
        FileNotFoundError
            If lattice files cannot be found in self.folder
        """

        if self.folder is None:
            print('nothing to get: lattice folder not defined')
            return
        
        # Set is_defined to True initially; will set to False if reading fails
        self.is_defined = True 


        folder_p = Path(self.folder)
        #
        # Read lattice input file
        #
        try:
            with open(folder_p / 'lattice_input.dat', 'r') as f:
                content = [line for line in f.readlines() if (line.strip() and not line.startswith('#'))]
            content = [line.split('#')[0] for line in content]
            for i,line in enumerate(content):
                ls = line.split()
                if ls[0] == 'lattice':
                    self.type = ls[1]
                if ls[0] == 'cell_vectors':
                    self.unit_cell_vectors = np.array([ [float(x) for x in content[i+1].split()],
                                                        [float(x) for x in content[i+2].split()] ])
                if ls[0] == 'repeat_cell':
                    self.size = np.array([ int(x) for x in ls[1:3] ], dtype=int)
                if ls[0] == 'n_site_types':
                    self.n_site_types = int(ls[1])
                if ls[0] == 'site_type_names':
                    self.site_type_names = ls[1:]
                if ls[0] == 'n_cell_sites':
                    self.n_cell_sites = int(ls[1])
                if ls[0] == 'site_types':
                    self.site_types = ls[1:]
                    # convert to names if given as indices
                    if self.site_types[0].isdigit():
                        self.site_types = [ self.site_type_names[int(x)-1] for x in self.site_types ]
                if ls[0] == 'site_coordinates':
                    self.fractional_coordinates = np.zeros((self.n_cell_sites,2), dtype=float)
                    for j in range(self.n_cell_sites):
                        self.fractional_coordinates[j,:] = np.array([float(x) for x in content[i+1+j].split()[:2]])
                if ls[0] == 'neighboring_structure':
                   self.neighboring_structure = []
                   j = 0
                   while content[i+1+j].split()[0] != 'end_neighboring_structure':
                       parts = content[i+1+j].split()
                       self.neighboring_structure.append( (int(parts[0].split('-')[0]), int(parts[0].split('-')[1]), parts[1]) )
                       j += 1
        except:
            print(f'cannot read lattice_input.dat from {self.folder}')
            self.is_defined = False

        #
        # Read lattice output file
        #
        site_coordinates = []
        site_types = []
        site_coordinations = []
        site_nns = []

        try:
            with open(folder_p / 'lattice_output.txt') as f:
                v1 = f.readline().split()[1:3]
                v2 = f.readline().split()[1:3]
                self.cell_vectors = np.array([v1, v2], dtype=float)
                for line in f:
                    ls = line.split()
                    site_coordinates.append(ls[1:3])
                    site_types.append(int(ls[3]))
                    site_coordinations.append(int(ls[4]))
                    # -1 is for python
                    site_nns.append(np.array([ int(ls[5+i])-1 for i in range(int(ls[4]))], dtype=int))

        except:
            print(f'cannot read lattice_output.txt from {self.folder}')
            self.is_defined = False

        self.coordinates = np.array(site_coordinates, dtype=float)
        self.site_types = np.array(site_types, dtype=int)
        self.site_coordinations = np.array(site_coordinations, dtype=int)
        self.site_nns = site_nns

    def __len__(self):
        """Return total number of lattice sites."""
        # If coordinates have been loaded from lattice_output.txt, use actual count
        if len(self.coordinates) > 0:
            return len(self.coordinates)
        # Otherwise use calculated size from unit cell
        return self.size[0] * self.size[1] * self.n_cell_sites

    def apply_pbc(self, coords):
        """
        Apply periodic boundary conditions to wrap coordinates into primary cell.
        
        Parameters
        ----------
        coords : array_like, shape (N, 2) or (2,)
            Cartesian coordinates to wrap
            
        Returns
        -------
        ndarray
            Wrapped coordinates in primary cell
        """

        coords = np.atleast_2d(coords)
        # Convert to fractional coordinates
        # r_cart = frac @ cell_vectors, so frac = r_cart @ inv(cell_vectors)
        cell_inv = np.linalg.inv(self.cell_vectors)
        frac_coords = coords @ cell_inv
        
        # Wrap to [0, 1)
        frac_coords = frac_coords - np.floor(frac_coords)
        
        # Convert back to Cartesian
        cart_coords = frac_coords @ self.cell_vectors
        
        return cart_coords.squeeze()

    def minimum_image_distance(self, coord1, coord2):
        """
        Calculate minimum image distance between two points with PBC.
        
        Parameters
        ----------
        coord1, coord2 : array_like, shape (2,)
            Cartesian coordinates of two points
            
        Returns
        -------
        float
            Minimum distance respecting periodic boundary conditions
        """
        # Displacement vector
        dr = np.array(coord2) - np.array(coord1)
        
        # Convert to fractional coordinates
        # cell_vectors rows are v1, v2, so: r_cart = frac @ cell_vectors
        # Therefore: frac = r_cart @ inv(cell_vectors)
        cell_inv = np.linalg.inv(self.cell_vectors)
        frac_dr = dr @ cell_inv
        
        # Apply minimum image convention: shift by -1, 0, or +1
        frac_dr = frac_dr - np.rint(frac_dr)
        
        # Convert back to Cartesian
        cart_dr = frac_dr @ self.cell_vectors
        
        return np.linalg.norm(cart_dr)

    def pairwise_distances_pbc(self, coords):
        """
        Calculate all pairwise distances with PBC (vectorized).
        
        Parameters
        ----------
        coords : ndarray, shape (N, 2)
            Cartesian coordinates of N points
            
        Returns
        -------
        distances : ndarray, shape (N, N)
            Matrix of pairwise distances with PBC
            
        Notes
        -----
        This vectorized implementation is much faster than calling 
        minimum_image_distance in nested loops (speedup: ~10-100x for N>100).
        """
        coords = np.asarray(coords)
        n = len(coords)
        
        # Compute all displacement vectors: dr[i,j] = coords[j] - coords[i]
        # Broadcasting: (N, 1, 2) - (1, N, 2) = (N, N, 2)
        dr = coords[np.newaxis, :, :] - coords[:, np.newaxis, :]
        
        # Convert to fractional coordinates
        cell_inv = np.linalg.inv(self.cell_vectors)
        # Reshape for batch matrix multiplication: (N*N, 2) @ (2, 2)
        frac_dr = dr.reshape(-1, 2) @ cell_inv
        
        # Apply minimum image convention
        frac_dr = frac_dr - np.rint(frac_dr)
        
        # Convert back to Cartesian
        cart_dr = frac_dr @ self.cell_vectors
        cart_dr = cart_dr.reshape(n, n, 2)
        
        # Compute norms
        distances = np.linalg.norm(cart_dr, axis=2)
        
        return distances

    def get_nn_distance(self, order=1):
        """
        Get nearest neighbor distance for FCC(111) lattice.
        
        Parameters
        ----------
        order : int
            Neighbor order (1=1nn, 2=2nn, etc.)
            
        Returns
        -------
        float
            Distance to nth nearest neighbor
            
        Notes
        -----
        For FCC(111) with lattice constant a:
        1nn = a, 2nn = sqrt(3)*a, 3nn = 2*a, etc.
        """
        a = np.linalg.norm(self.unit_cell_vectors[0])
        
        # Distance formulas for FCC(111)
        nn_distances = {
            1: a,
            2: np.sqrt(3) * a,
            3: 2 * a,
            4: np.sqrt(7) * a,
            5: 3 * a,
            6: np.sqrt(12) * a,
            7: np.sqrt(13) * a,
            8: 4 * a,
            9: np.sqrt(19) * a,
        }
        
        if order in nn_distances:
            return nn_distances[order]
        else:
            raise ValueError(f"Neighbor order {order} not implemented. Valid orders: 1-9")

    def get_cell_area(self):
        """
        Calculate area of the simulation cell.
        
        Returns
        -------
        float
            Area in square distance units
        """
        # 2D cross product: |v1 × v2| = v1_x * v2_y - v1_y * v2_x
        v1, v2 = self.cell_vectors
        return abs(v1[0] * v2[1] - v1[1] * v2[0])

    def __repr__(self):
        """String representation of lattice"""
        return (f"lattice(type='{self.type}', size={tuple(self.size)}, "
                f"nsites={len(self)}, area={self.get_cell_area():.2f})")
