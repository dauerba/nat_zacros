"""
State class for adsorbate configurations on lattice.

This module provides the `state` class for representing snapshots of adsorbate
configurations on a surface lattice from Zacros simulations.
"""

import numpy as np
from pathlib import Path
from scipy.spatial import cKDTree


class State:
    """
    Adsorbate configuration on a surface lattice.
    
    Represents a snapshot of which sites are occupied and by which species
    at a particular moment in a Zacros simulation.
    
    Attributes
    ----------
    lattice : Lattice object
        Reference to the underlying surface lattice
    folder : str (must be str for pickle cross-platform compatibility) or None
        Directory containing history_output.txt. Stored as string (not Path)
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
    load(idx=0)
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
    

    def __init__(self, lattice, metadata, dirname=None):
        """
        Initialize state object.
        
        Parameters
        ----------
        lattice : Lattice object
            The surface lattice on which adsorbates are placed
        dirname : str or Path, optional
            Directory containing history_output.txt. If provided,
            will automatically load the first state.
        """
        self.dir = str(dirname) if dirname is not None else None

        # Store reference to lattice
        self.lattice = lattice

        # Get metadate and check that necessary keys are present
        self.metadata = metadata

        self._necessary_keys = [
            'surf_species_names'
            ]

        args_ok = True
        for key in self._necessary_keys:
            if key not in self.metadata:
                print(f" {key} undefined.")
                args_ok = False

        if not args_ok:
            self.is_valid = False
            print('Not enough metadata for a valid simulation.')

        # default is no species
        self.gas_species_names = []
        self.n_gas_species = len(self.gas_species_names)
        self.n_surf_species = len(self.metadata['surf_species_names'])

        self.surf_species_names = dict(zip(self.metadata['surf_species_names'],range(self.n_surf_species)))
        self.surf_species_dent = []

        # Arrays defining the adsorbed species on the lattice
        # with indices corresponding to lattice site indices starting at 1

        nsites = len(self.lattice) 
        self.ads_ids =    np.zeros(nsites, dtype=int)
        self.occupation = np.zeros(nsites, dtype=int)
        self.dentation =  np.zeros(nsites, dtype=int)


    def load(self, idx=0, verbose=False):
        """
        Read configuration from history_output.txt file.
        
        Parameters
        ----------
        idx : int, optional
            Index of the configuration to read (default: 0 = first state)
            If negative, counts from the end (e.g., -1 = last state)
        verbose : bool, optional
            If True, print verbose output.
            
        Notes
        -----
        Reads from history_output.txt which contains sequential snapshots.
        Each snapshot lists the occupation, species ID, and denticity for
        every lattice site.
        """

        if verbose:
            print(f"Loading state {idx} from trajectory {self.dir}")

        # Read configuration from history_output.txt file
        try:
            with open(Path(self.dir) / 'history_output.txt', 'r') as f:
                content = f.readlines()    

            nsites = len(self.lattice)
            for site in range(nsites):
                parts = content[7 + idx*(nsites+1) + site].split()
                self.ads_ids[site]    = int(parts[1])
                self.occupation[site] = int(parts[2])
                self.dentation[site]  = int(parts[3])
        except Exception as e:
            print(f'Error loading a state from trajectory {str(self.dir)}: {e}')
    

    def get_accessibility(self):
        """
        Calculate accessibility (number of vacant nearest neighbors) for each occupied site.
        
        Returns
        -------
        accessibility : ndarray
            Number of vacant nearest neighbors (0 to max_coordination)
            
        Notes
        -----
        Accessibility measures how many nearest neighbor sites are vacant,
        which affects reactivity and diffusion rates.
        """

        acc13_list = []
        acc2_list  = []
        
        for site_idx in self.get_occupied_sites():
            # Get nearest neighbors for this site
            shell_1 = self.lattice.site_nns[site_idx]

            # 2nd shell (angle between lines connecting site_idx with its nn site 
            #            and the nn site with a 2nd shell site is 120 degrees)
            shell_2 = [self.lattice.site_nns[s][i-1] for i, s in enumerate(shell_1)]

            # 3rd shell (the line connecting site_idx with its nn site and its 3rd shell site is straight)
            shell_3 = [self.lattice.site_nns[s][i] for i, s in enumerate(shell_1)]

            # Count vacant neighbors
            vacant_13 = np.sum((self.occupation[shell_1] + self.occupation[shell_3]) == 0)
            vacant_2  = np.sum(self.occupation[shell_2] == 0)

            acc13_list.append(vacant_13)
            acc2_list.append(vacant_2)

                
        return acc13_list, acc2_list

    def get_clusters(self, cutoff, eps=1e-4):
        """
        Cluster 2D points with periodic boundary conditions
        using kD-Tree and Union method

        Parameters
        ----------
        cutoff : float
            Distance cutoff for connectivity (e.g. 3rd NN distance).
        eps   : float, default 1e-4
            Tolerance parameter for the approximate search in the tree (see docs for scipy.spatial.ckdtree package)

        Returns
        -------
        labels : (n_ads,) int
            Cluster label for each original point (0..nclusters-1).
        clusters : list of ndarray
            Indices of points in each cluster.
        sizes : list of int 
            Sizes of clusters
        """

        class UnionFind:
            def __init__(self, n):
                self.parent = np.arange(n)
                self.rank = np.zeros(n, dtype=int)

            def find(self, a):
                p = self.parent
                while p[a] != a:
                    p[a] = p[p[a]]
                    a = p[a]
                return a

            def union(self, a, b):
                ra = self.find(a); rb = self.find(b)
                if ra == rb:
                    return
                if self.rank[ra] < self.rank[rb]:
                    self.parent[ra] = rb
                else:
                    self.parent[rb] = ra
                    if self.rank[ra] == self.rank[rb]:
                        self.rank[ra] += 1
                
        
        pts = self.get_occupied_coords()

        if pts.size == 0:
            return np.array([], dtype=int), [], []

        # Get number of occupied sites
        n_ads = len(pts)

        # Build augmented points = original points shifted by translations i*v1 + j*v2 with i,j in {-1,0,1}
        shifts = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]
        aug_pts = np.zeros((n_ads * len(shifts), 2), dtype=float)
        orig_idx = np.zeros(n_ads * len(shifts), dtype=int)

        k = 0
        for si, (i, j) in enumerate(shifts):
            shift_vec = i*self.lattice.cell_vectors[0] + j*self.lattice.cell_vectors[1]
            aug_pts[k:k+n_ads] = pts + shift_vec
            orig_idx[k:k+n_ads] = np.arange(n_ads)
            k += n_ads

        # KD-tree on augmented points
        tree = cKDTree(aug_pts)
        # Get pairs inside of cutoff*(1+eps) with eps ensuring that single-precision lattice vectors would work
        pairs = tree.query_pairs(cutoff,eps=eps, output_type='ndarray')  # array of shape (M,2)

        uf = UnionFind(n_ads)
        for a, b in pairs:
            ia = orig_idx[a]
            ib = orig_idx[b]
            if ia != ib:
                uf.union(ia, ib)

        # Extract roots and relabel to contiguous labels
        roots = np.array([uf.find(i) for i in range(n_ads)])
        unique_roots, inv = np.unique(roots, return_inverse=True)
        labels = inv
        clusters = [np.nonzero(labels == k)[0] for k in range(len(unique_roots))]
        sizes = [len(c) for c in clusters]

        return labels, clusters, sizes

    def get_clusters_bfs(self, cutoff, eps=1e-4):
        """
        Cluster 2D points with periodic boundary conditions
        using BFS (Breadth First Search) method

        Parameters
        ----------
        cutoff : float
            Distance cutoff for connectivity (e.g. 3rd NN distance).
        eps   : float, default 1e-4
            Tolerance parameter for the cutoff

        Returns
        -------
        labels : (n_ads,) int
            Cluster label for each original point (0..nclusters-1).
        clusters : list of ndarray
            Indices of points in each cluster.
        sizes : list of int 
            Sizes of clusters
        """

        cutoff = cutoff*(1 + eps)

        pts = self.get_occupied_coords()

        if pts.size == 0:
            return np.array([], dtype=int), [], []

        # Get number of occupied sites
        n_ads = len(pts)

        # Build adjacency matrix
        sizes = []
        clusters = []
        visited = np.zeros(n_ads, dtype=bool)
        
        for i in range(n_ads):
            if visited[i]:
                continue
                
            # Start new cluster with BFS
            queue   = np.zeros(n_ads, dtype=bool)
            cluster = np.zeros(n_ads, dtype=bool)
            queue[i]   = True
            visited[i] = True
            
            while np.any(queue):
                current = np.where(queue)[0][-1]
                queue[current] = False
                cluster[current] = True
                # Check neighbors
                for j in range(n_ads):
                    if not visited[j]:
                        dist = self.lattice.minimum_image_distance(pts[current], pts[j])
                        #dist = np.linalg.norm(pts[current] - pts[j])
                        if dist <= cutoff:
                            visited[j] = True
                            queue[j]   = True

            cl = np.where(cluster)[0]
            clusters.append(cl)
            sizes.append(len(cl))

        return np.arange(len(clusters)), clusters, sizes


    def get_coverage(self):
        """
        Calculate the coverage (fraction of occupied sites).
        
        Returns
        -------
        float
            Fraction of sites that are occupied (0.0 to 1.0)
        """
        return np.count_nonzero(self.occupation) / len(self.lattice)

    def get_coverage(self, species = None):
        """
        Calculate the coverage (fraction of occupied sites).

        Parameters
        ----------
        species : integer or None; default None
            Index of a species; if None, include all species
        
        Returns
        -------
        float
            Fraction of sites that are occupied (0.0 to 1.0)
        """
        if species is None:
            return np.count_nonzero(self.occupation) / len(self.lattice)
        else:
            return np.count_nonzero(self.occupation == species)

    
    def get_leed_intensity(self):
        """
        Calculate LEED intensity.
        
        Returns
        -------
        intensity : float
            LEED intensity
            
        """

        n_pairs = 0
        for site_idx in self.get_occupied_sites():

            # Get nearest neighbors for this site
            shell_1 = self.lattice.site_nns[site_idx]

            # 3rd shell (the line connecting site_idx with its nn site and its 3rd shell site is straight)
            shell_3 = [self.lattice.site_nns[s][i] for i, s in enumerate(shell_1)]

            # Count vacant neighbors
            n_pairs  += np.sum(self.occupation[shell_3] > 0)
        
        # Avoid double-counting
        n_pairs //= 2
                
        return n_pairs**2


    def get_rdf(self, species_1, species_2, distances, r_max=None, dr=0.1, g_ref=None):
        """
        Calculate rdf for two species.
        
        Parameters
        ----------
        species_1, species_2 : str
            Species for which rdf is to be calculated
        distances : (N,N) array of floats
            Matrix of lattice site distances
        r_max : float, optional
            Maximum distance for RDF (default: half of the cell diagonal)
        dr : float, default 0.1
            Bin width in Angstroms
        g_ref : ndarray, optional
            Reference RDF for normalization (from full lattice at coverage=1).
            If provided, normalizes by number of neighbors in each shell.
            
        Returns
        -------
        r_bins : ndarray
            Bin centers  
        g_r : ndarray
            RDF values normalized such that g(r)=1 for ideal gas
            
        """

        if r_max is None:
            # Default: half the minimum cell dimension
            v1 = self.lattice.cell_vectors[0]
            v2 = self.lattice.cell_vectors[1]
            l1 = np.linalg.norm(v1)
            l2 = np.linalg.norm(v2)
            l3 = np.linalg.norm(v1 + v2)
            r_max = min(l1, l2, l3) / 2.0

        # Initialize histogram
        n_bins = int(np.ceil(r_max / dr))
        bin_edges = np.linspace(0.0, r_max, n_bins + 1)
        r_bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        g_r = np.zeros(n_bins)

        if species_1 == species_2:
            occupied_sites_1 = self.get_occupied_sites(species=species_1) 
            occupied_sites_2 = occupied_sites_1 
            # Get distances between occupied sites:
            # by selecting distances submatrix using np.ix_ mesh constructor
            dist_sel    = distances[np.ix_(occupied_sites_1, occupied_sites_2)]
            # and taking the elements above diagonal
            n_occupied_sites_1 = len(occupied_sites_1)
            n_occupied_sites_2 = n_occupied_sites_1
            dist_vector = dist_sel[np.triu_indices(n_occupied_sites_1,1)]
            # Factor of 2 to account for unordered pairs
            degeneracy = 2

        else:
            occupied_sites_1 = self.get_occupied_sites(species=species_1) 
            occupied_sites_2 = self.get_occupied_sites(species=species_2) 
            # Get distances between occupied sites:
            # by selecting distances submatrix using np.ix_ mesh constructor
            dist_vector = distances[np.ix_(occupied_sites_1, occupied_sites_2)].flatten()

            n_occupied_sites_1 = len(occupied_sites_1)
            n_occupied_sites_2 = len(occupied_sites_2)
            degeneracy = 1

        # Histogram
        counts, _ = np.histogram(dist_vector[(dist_vector > 0) & (dist_vector < r_max)], bins=bin_edges)
        
        # Normalize by g_ref if provided (number of neighbors in each shell)
        if g_ref is not None:
            counts_n = np.zeros_like(counts, dtype=float)
            np.divide(counts, g_ref, out=counts_n, where=g_ref!=0)
            g_r = counts_n * degeneracy / (n_occupied_sites_1 * n_occupied_sites_2) * len(self.lattice)
        else:
            g_r = counts * degeneracy / (n_occupied_sites_1 * n_occupied_sites_2) * len(self.lattice)

        return r_bins, g_r


    def get_occupied_sites(self, species = None):
        """
        Get indices of sites occupied by a species.

        Parameters
        ----------
        species : str or None; default None
            Surface species name; if None, include all species
        
        Returns
        -------
        ndarray
            Array of occuied site indices
        """
        if species is None:
            return np.where(self.occupation > 0)[0]
        else:
            return np.where(self.occupation == self.surf_species_names[species] + 1)[0]

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
        """String representation of State class"""
        coverage = self.get_coverage()
        return f"State(nsites={len(self.lattice)}, surf_species={self.surf_species_names}, "\
               f"n_adsorbates={self.n_ads()}, coverage={coverage:.3f})"
