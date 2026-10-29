"""
State class for adsorbate configurations on lattice.

This module provides the `state` class for representing snapshots of adsorbate
configurations on a surface lattice from Zacros simulations.
"""

import numpy as np
from pathlib import Path
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt


class State:
    """
    Adsorbate configuration on a surface lattice.
    
    Represents a snapshot of which sites are occupied and by which species
    at a particular moment in a Zacros simulation.
    
    Attributes
    ----------
    lattice : Lattice object
        Reference to the underlying surface lattice
    dir : str or None
        Directory containing the source trajectory files. Stored as a string for
        simple serialization.
    n_surf_species : int
        Number of surface-bound species
    surf_species_names : list of str
        Names of surface species
    surf_species_dent : list
        Denticity of each surface species
    occupation : ndarray, shape (N,)
        Surface-species code at each site (0 = empty, 1..n = entries in
        ``surf_species_names`` plus one)
    dentation : ndarray, shape (N,)
        Denticity at each site
        
    Methods
    -------
    load(idx=0)
        Read configuration from history_output.txt
    get_accessibility()
        Get accessibilities
    get_clusters()
        Get clusters
    get_coverages()
        Calculate fraction of occupied sites
    get_species_sites()
        Get indices of occupied sites
    get_empty_sites()
        Get indices of empty sites
    get_species_coords()
        Get Cartesian coordinates of occupied sites
    n_ads()
        Get total number of adsorbates
    """
    

    def __init__(self, lattice, input_pars, traj_path=None):
        """
        Initialize state object.
        
        Parameters
        ----------
        lattice : Lattice object
            The surface lattice on which adsorbates are placed
        dirname : str or Path, optional
            Directory containing Zacros output. 
        """

        self.traj_path = Path(traj_path) if traj_path is not None else None

        # Store reference to lattice
        self.lattice = lattice

        # Set attributes
        self.input_pars = input_pars
        self.n_surf_species = int(input_pars['n_surf_species'][0])
        self.surf_species_names = dict(zip(input_pars['surf_specs_names'], range(1,self.n_surf_species+1)))
        self.surf_species_names['*'] = 0
        self.surf_species_dent = dict(zip(input_pars['surf_specs_names'], [int(d) for d in input_pars['surf_specs_dent']]))
        self.surf_species_dent['*'] = 1
        self.dent_lookup = np.array([1] + input_pars['surf_specs_dent'], dtype=int)

        # Arrays of gas species and adsorbed species on the lattice
        self.occupation = np.zeros(self.lattice.n_sites, dtype=int)
        self.dentation =  np.zeros(self.lattice.n_sites, dtype=int)


    def load(self, idx=0, verbose=False):
        """
        Read one configuration block from ``history_output.txt``.
        
        Parameters
        ----------
        idx : int, optional
            Index of the configuration to read (default: 0 = first state)
            If negative, counts from the end (e.g., -1 = last state)
        verbose : bool, optional
            If True, print verbose output.
            
        Notes
        -----
        Each configuration block supplies the adsorbate ID, surface-species
        code, and denticity for every lattice site.
        """

        if verbose:
            print(f"Loading state {idx} from trajectory {self.dir}")

        # Read configuration from history_output.txt file
        try:
            with open(Path(self.dir) / 'history_output.txt', 'r') as f:
                content = f.readlines()    

            nsites = self.lattice.n_sites
            for site in range(nsites):
                parts = content[7 + idx*(nsites+1) + site].split()
                self.occupation[site] = int(parts[2])
                self.dentation[site]  = int(parts[3])
        except Exception as e:
            print(f'Error loading a state from trajectory {str(self.dir)}: {e}')
    

    def get_accessibility(self):
        """
        Calculate shell-resolved accessibilities for each occupied site.
        
        Returns
        -------
        acc13_list : list[int]
            Number of vacant first-plus-third-shell neighbors for each occupied
            site.
        acc2_list : list[int]
            Number of vacant second-shell neighbors for each occupied site.
            
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
                        dist = np.linalg.norm(pts[current] - pts[j])
                        if dist <= cutoff:
                            visited[j] = True
                            queue[j]   = True

            cl = np.where(cluster)[0]
            clusters.append(cl)
            sizes.append(len(cl))

        return np.arange(len(clusters)), clusters, sizes


    def get_nads(self):
        """
        Calculate the coverages for surface species (the ratio of adsorbate molecules 
        to the number of structural surface molecules).

        Returns
        -------
        coverages: vector of floats (number of surface species + 1)
            partial coverages of surface species, 
            index 0 refers to the coverage of empty sites
        """
        return np.bincount(self.occupation) // self.dent_lookup

    
    def get_leed_intensity(self, distances, r):
        """
        Calculate LEED intensity.
        
        Returns
        -------
        intensity : float
            LEED intensity
            
        """
        # Get indices of occupied sites
        occupied_sites = self.get_occupied_sites()
        # Select the corresponding part of distances matrix
        d_occ    = distances[np.ix_(occupied_sites, occupied_sites)]

        n_pairs = len(d_occ[np.isclose(d_occ, r)])/2

        n_ads_in_pairs = np.sum(np.any(np.isclose(d_occ, r), axis=1))

        #return n_pairs**2
        return n_ads_in_pairs**2


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
            occupied_sites_1 = self.get_species_sites(species=species_1) 
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
            occupied_sites_1 = self.get_species_sites(species=species_1) 
            occupied_sites_2 = self.get_species_sites(species=species_2) 
            # Get distances between occupied sites:
            # by selecting distances submatrix using np.ix_ mesh constructor
            dist_vector = distances[np.ix_(occupied_sites_1, occupied_sites_2)].flatten()

            n_occupied_sites_1 = len(occupied_sites_1)
            n_occupied_sites_2 = len(occupied_sites_2)
            degeneracy = 1

        # Histogram
        counts, _ = np.histogram(dist_vector[(dist_vector > 0) & (dist_vector < r_max)], bins=bin_edges)
        
        # Normalize by g_ref if provided (number of neighbors in each shell)
        # Warning: not checked for n_site_types > 1
        if g_ref is not None:
            counts_n = np.zeros_like(counts, dtype=float)
            np.divide(counts, g_ref, out=counts_n, where=g_ref!=0)
            g_r = counts_n * degeneracy / (n_occupied_sites_1 * n_occupied_sites_2) * self.lattice.n_sites
        else:
            g_r = counts * degeneracy / (n_occupied_sites_1 * n_occupied_sites_2) * self.lattice.n_sites

        return r_bins, g_r


    def get_species_sites(self, species, site_type=None):
        """
        Get indices of sites of site_type occupied by a species.

        Parameters
        ----------
        species : str;
            Surface species name; if None, include all species
        site_type : str or None; default None
            Site type name; if None, include all site types

        Returns
        -------
        ndarray
            Array of occupied site indices
        """
        mask = self.occupation == self.surf_species_names[species]

        if site_type is not None:
            site_types = self.lattice.site_types
            mask = mask & (site_types == self.lattice.site_type_names.index(site_type))

        return np.where(mask)[0]

    def get_occupied_sites(self, site_type=None):
        """
        Get indices of occupied sites of site_type 

        Parameters
        ----------
        site_type : str or None; default None
            Site type name; if None, include all site types

        Returns
        -------
        ndarray
            Array of occupied site indices
        """
        mask = self.occupation > 0

        if site_type is not None:
            site_types = self.lattice.site_types
            mask = mask & (site_types == self.lattice.site_type_names.index(site_type) )

        return np.where(mask)[0]

    def get_species_coords(self, species, site_type=None):
        """
        Get indices of sites of site_type occupied by a species.

        Parameters
        ----------
        species : str;
            Surface species name; if None, include all species
        site_type : str or None; default None
            Site type name; if None, include all site types

        Returns
        -------
        ndarray
            Array of occupied site indices
        """
        indices = self.get_species_sites(species, site_type=site_type)
        return self.lattice.coordinates[indices]

    def get_occupied_coords(self, site_type=None):
        """
        Get indices of occupied sites of site_type 

        Parameters
        ----------
        site_type : str or None; default None
            Site type name; if None, include all site types

        Returns
        -------
        ndarray
            Array of occupied site indices
        """
        mask = self.occupation > 0

        indices = self.get_occupied_sites(site_type=site_type)
        return self.lattice.coordinates[indices]


    def plot(self, ax=None, scaling=1, ads_scaling=1.2, markers=None, colors=None, legend=True,
                    show_axis = True, coverage_norm=1):
        """Plot the state of the lattice.

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
        colors_default = ["lightgray", "r", "g", "b", "m", "c", "k",
                        "tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple",
                        "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan",
                        "gold", "turquoise", "lime", "indigo"]

        if markers is None:
            markers = markers_default
        elif len(markers) < len(self.lattice.site_type_names):
            markers = markers_default
            print(f"Not enough markers provided for {len(self.lattice.site_type_names)} site type(s). Using defaults.")
        if colors is None:
            colors = colors_default
        elif len(colors) < self.n_surf_species + 1:
            colors = colors_default
            print(f"Not enough colors provided for {self.n_surf_species + 1} species. Using defaults.")

        # Create a new figure and axes if none are provided
        if not ax:
            fig, ax = plt.subplots(figsize=(8,6))
        
        # Get temperature and coverages for the title
        Ts = float(self.input_pars['temperature'][0])
        sp_names = [sp[:-1] if sp.endswith('*') else sp for sp in list(self.surf_species_names.keys()) if sp != '*']
        covs = self.get_nads() 
        title = f'T = {Ts} K, ' + ', '.join([fr'$\theta_{{{sp}}} = {{{covs[self.surf_species_names[sp+'*']]/ (coverage_norm*self.lattice.n_cells):.3f}}}$' for sp in sp_names])

        # Set up the plot
        ax.set_title(title)
        ax.set_xlabel(r'x ($\AA$)')
        ax.set_ylabel(r'y ($\AA$)')
        ax.set_aspect(1.0)

        # Set marker size based on the number of lattice sites and scaling factor
        size = 1.5*440 / np.sqrt(self.lattice.n_sites) * scaling

        # Put empty sites on the plot
        for ist, st in enumerate(self.lattice.site_type_names):
            x, y = self.get_species_coords('*', site_type=st).T
            if len(x) > 0:  # plot only if there are points to plot
                ax.scatter(x, y, color=colors[0], marker= markers[ist], s=size, zorder=2,
                        label= st)

        # Loop over site types and surface species to plot occupied and empty sites
        for isp in range(self.n_surf_species):
            for ist, st in enumerate(self.lattice.site_type_names):

                sp_name = list(self.surf_species_names.keys())[isp]
                x, y = self.get_species_coords(species=sp_name, site_type=st).T
                marker_size = size * ads_scaling
                if len(x) > 0:  # plot only if there are points to plot
                    ax.scatter(x, y, color=colors[isp+1], marker= markers[ist], s=marker_size, zorder=isp+2,
                            label= sp_name + '-' + st)
        if legend:
            ax.legend(loc='lower right', frameon=False)
        if not show_axis:
            ax.axis('off')

        plt.show()
                    

    def __repr__(self):
        """String representation of State class"""
        occ = [ (k, len(self.get_species_sites(k))) for k in self.surf_species_names.keys() if k != '*']
        return f"State(nsites={self.lattice.n_sites}, adsorbates={occ})"
