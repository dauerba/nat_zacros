"""
Lattice class for FCC(111) surface geometry.

This module provides the Lattice class for representing periodic FCC(111) 
surface lattices used in Zacros kinetic Monte Carlo simulations.
"""

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

class Lattice:
    """
    FCC(111) surface lattice with periodic boundary conditions.
    
    Attributes
    ----------
    cell_vectors : ndarray, shape (2, 2)
        Supercell lattice vectors (unit_cell_vectors * size)
    coordinates : ndarray, shape (N, 2)
        Cartesian coordinates of all lattice sites
    folder : str (must be str for pickle cross-platform compatibility) or None
        Directory containing lattice_input.dat and lattice_output.txt. Stored as string (not Path)
    fractional_coordinates : ndarray
        Fractional coordinates of sites within unit cell
    is_defined : bool
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
    get_cell_area()
        Calculate simulation cell area
    load()
        Read lattice definition from input/output files
    get_nn_distance(order=1)
        Get nearest neighbor distance for given order
    pairwise_distances_pbc(coords)
        Calculate all pairwise distances with PBC (vectorized)
    """
    
    def __init__(self, dirname=None):
        """
        Initialize Lattice object.
        
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
            self.load()

    def load(self):
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
            # Use tolerant UTF-8 decoding because external generators may include
            # non-UTF8 bytes in comment headers that are irrelevant for parsing.
            with open(folder_p / 'lattice_input.dat', 'r', encoding='utf-8', errors='replace') as f:
                content = []

                # dja change 2026-08-20
                # correct bug that only ignored comments at the start of the line

                # read lines ignoring comments and empty lines
                for raw_line in f:
                    line = raw_line.split('#', 1)[0].strip()
                    if line:
                        content.append(line)
                # end dja change
                
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
            print(f'error reading lattice_input.dat from {self.folder}')
            self.is_defined = False

        #
        # Read lattice output file
        #
        site_coordinates = []
        site_types = []
        site_coordinations = []
        site_nns = []

        try:
            with open(folder_p / 'lattice_output.txt', 'r', encoding='utf-8', errors='replace') as f:
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

    def pairwise_distances_pbc(self, condensed=False):
        """
        Calculate pairwise distances between all lattice sites with PBC (vectorized).
        
        Parameters
        ----------
        condensed : bool, default False
            Defines the form in which distances matrix is returned
            
        Returns
        -------
        distances : ndarray, shape (N, N)
            If condensed,    ( N*(N-1)/2 ) vector of pairwise distances with PBC
            If not condensed,        (N*N) matrix of pairwise distances with PBC

            
        Notes
        -----
        This vectorized implementation is much faster than calling 
        minimum_image_distance in nested loops (speedup: ~10-100x for N>100).
        """
        coords = self.coordinates
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

        # Return results based on condensed arg
        if condensed:
            mask = np.triu(np.ones(distances.shape, dtype=bool), k=1)
            return distances[mask]
        else:
            return distances

    # def get_nn_distance(self, order=1):
    #     """
    #     Get nearest neighbor distance for FCC(111) lattice.
        
    #     Parameters
    #     ----------
    #     order : int
    #         Neighbor order (1=1nn, 2=2nn, etc.)
            
    #     Returns
    #     -------
    #     float
    #         Distance to nth nearest neighbor
            
    #     Notes
    #     -----
    #     For FCC(111) with lattice constant a:
    #     1nn = a, 2nn = sqrt(3)*a, 3nn = 2*a, etc.
    #     """
    #     a = np.linalg.norm(self.unit_cell_vectors[0])
        
    #     # Distance formulas for FCC(111)
    #     nn_distances = {
    #         1: a,
    #         2: np.sqrt(3) * a,
    #         3: 2 * a,
    #         4: np.sqrt(7) * a,
    #         5: 3 * a,
    #         6: np.sqrt(12) * a,
    #         7: np.sqrt(13) * a,
    #         8: 4 * a,
    #         9: np.sqrt(19) * a,
    #     }
        
    #     if order in nn_distances:
    #         return nn_distances[order]
    #     else:
    #         raise ValueError(f"Neighbor order {order} not implemented. Valid orders: 1-9")

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


    def plot(self, ax=None, scaling=1, markers=None, colors=None, legend=True,
                    show_axis = True, link=True, figsize = (8,6), legend_loc = 'outside center right'):
        """Plot the lattice.

        Parameters
        ----------
        ax : matplotlib.axes.Axes or None; default None
            The axes on which to plot. If None, a new figure and axes are created.
        scaling : float; default 1
            Scaling factor for the size of the markers.
        ads_scaling:

        """

        # Library of markers and colors for plotting different species and site types
        markers_default =   ["s", "o", "v", "D", "p", "^", "+", "x", "*", "P", 
                                "H", "X", "d", "h", ",", ".", "<", ">", "1", "2"]
        colors_default = ["gray", "r", "g", "b", "m", "c", "k",
                        "tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple",
                        "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan",
                        "gold", "turquoise", "lime", "indigo"]

        if markers is None:
            markers = markers_default
        elif len(markers) < len(self.site_type_names):
            markers = markers_default
            print(f"Not enough markers provided for {len(self.site_type_names)} site type(s). Using defaults.")
        if colors is None:
            colors = colors_default
        elif len(colors) < len(self.site_type_names):
            colors = colors_default
            print(f"Not enough colors provided for {len(self.site_type_names)} site type(s). Using defaults.")

        # Create a new figure and axes if none are provided
        if not ax:
            fig, ax = plt.subplots(figsize=figsize)
        
        title = "Lattice structure"

        # Set up the plot
        ax.set_title(title)
        ax.set_xlabel(r'x ($\AA$)')
        ax.set_ylabel(r'y ($\AA$)')
        ax.set_aspect(1.0)

        # Set marker size based on the number of lattice sites and scaling factor
        size = 1.5*440 / np.sqrt(len(self)) * scaling

        # Get lattice site coordinates
        x, y = self.coordinates.T

        if link:
            for i in range(len(x)):
                for xn, yn in self.coordinates[self.site_nns[i]]:
                    if np.abs(x[i] - xn) < 0.5*np.max(x) and np.abs(y[i] - yn) < 0.5*np.max(y):
                        ax.plot([x[i], xn],[y[i], yn], lw=0.5, color='lightgray', zorder=0)
        

        # Loop over site types
        for i_st in range(len(self.site_type_names)):
            # Get indices for a site type i_st+1
            ids = np.where(self.site_types == i_st+1)[0]

            ax.scatter(x[ids], y[ids], color=colors[i_st], marker= markers[i_st], s=size, 
                    label= self.site_type_names[i_st])

        if legend:
            fig.legend(loc=legend_loc, frameon=True)
        if not show_axis:
            ax.axis('off')
                    
    

    def __len__(self):
        """Return total number of Lattice sites."""
        # If coordinates have been loaded from lattice_output.txt, use actual count
        if len(self.coordinates) > 0:
            return len(self.coordinates)
        # Otherwise use calculated size from unit cell
        return self.size[0] * self.size[1] * self.n_cell_sites

    def __repr__(self):
        """String representation of Lattice class"""
        return (f"Lattice(type='{self.type}', size={tuple(self.size)}, "
                f"nsites={len(self)}, area={self.get_cell_area():.2f})")
