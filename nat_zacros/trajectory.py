"""
Trajectory class for sequences of lattice states over time.

This module provides the `trajectory` class for managing and analyzing
sequences of adsorbate configurations from Zacros KMC simulations.
"""

import pickle
import numpy as np
from pathlib import Path
from .state import State

class Trajectory:
    """
    Container for a sequence of lattice states over time.
    
    Attributes
    ----------
    lattice : Lattice object
        The underlying surface lattice
    states : list of State objects
        Sequence of configurations
    times : ndarray
        Time points for each state
    energies : ndarray
        Total energy for each state
# Warning 307!
    folder : str (must be str for pickle compatibility)
        Directory containing trajectory data
        
    Methods
    -------
    load(...)
        Load states from history_output.txt
    add_state(state, time, energy)
        Add a state to the trajectory
    get_energy_vs_time()
        Get energy as function of time
    estimate_equilibration(fraction=0.5)
        Estimate equilibration index
    get_accessibility_histogram()
        Calculate histogram of site accessibility
    get_cluster_distribution(nn_cutoff)
        Calculate cluster size distribution
    get_coverage_vs_time()
        Get coverage as function of time
    get_energy_vs_time()
        Get energy as function of time
    get_g_ref(r_max, dr)
        Calculate reference RDF for normalization
    get_rdf(r_max, dr, g_ref, vectorized)
        Calculate radial distribution function
    """
    
    def __init__(self, dir, lattice, cache_file='traj.pkl'):
        """
        Initialize trajectory with lattice and optional data folder.
        
        Parameters
        ----------
        lattice : Lattice object
            The surface lattice for this trajectory
        dir : str or Path, optional
            Directory containing history_output.txt
        """

        self.dir = Path(dir)
        self.lattice = lattice
        self.states = []
        self.times = []
        self.energies = []
        self.cache_file = self.dir / cache_file
        
    def get_energy_vs_time(self):
        """
        Get energy as a function of time.

        Returns
        -------
        times : ndarray
            Time points
        energies : ndarray
            Energies at each time point
        """
        return self.times, self.energies
      
                
    def load(self, cache=True, verbose=False):
        """
        Load states from history_output.txt.

        Parameters
        ----------
        cache : bool, default True
            If True, attempt to load cached trajectory data if available; cache trajectory data if not already cached.
            If False, load from raw simulation data.
        verbose : bool, default False
            If True, print detailed loading information.
        """

        # Try loading from cache
        if cache:
            if self.cache_file.exists():
                with open(self.cache_file, 'rb') as f:
                    self.states, self.times, self.energies = pickle.load(f)
                if verbose:
                    size_mb = self.cache_file.stat().st_size / 1024**2
                    print(f"Loaded trajectory {self.dir.name} with {len(self.states)} states from cache ({size_mb:.1f} MB).")
                return

        try:
            # Single-pass streaming parse (fast)
            nsites = len(self.lattice)

            # Read in a trajectory from history_output.txt
            with open(self.dir / 'history_output.txt', 'r') as f: 
                content = f.readlines()

            # Count total number of states in trajectory
            n_states = sum(line.lstrip().startswith('configuration') for line in content)

            self.times = np.empty(n_states, dtype=float)
            self.energies = np.empty(n_states, dtype=float)
            self.states = [None] * n_states

            # Find the first configuration line
            pos = 0
            while pos < len(content) and not content[pos].lstrip().startswith('configuration'):
                pos += 1
            
            # Find the second configuration line to determine block size
            pos2 = pos + 1
            while pos2 < len(content) and not content[pos2].lstrip().startswith('configuration'):
                pos2 += 1
            
            block_size = pos2 - pos  # Lines per configuration block
            # Check consistency
            expected_block_size = 1 + nsites  # 1 header + nsites data
            if block_size != expected_block_size:
                raise ValueError(f'block size: {block_size}, expected {expected_block_size}.')


            for k in range(n_states):

                pos = 6 + k * block_size
                line = content[pos]

                try:
                    parts = line.split()
                    time = float(parts[3])
                    energy = float(parts[5])
                except (ValueError, IndexError):
                    raise ValueError(f'{str(self.dir.name)}: Failed to parse line: {line.strip()}')

                st = State(self.lattice, dirname=self.dir)

                for site in range(nsites):
                    site_line = content[pos + 1 + site]
                    p = site_line.split()
                    st.ads_ids[site] = int(p[1])
                    st.occupation[site] = int(p[2])
                    st.dentation[site] = int(p[3])

                self.states[k] = st
                self.times[k] = time
                self.energies[k] = energy

            # Save to cache
            if cache:
                with open(self.cache_file, 'wb') as f:
                    pickle.dump([self.states, self.times, self.energies], f)
                    size_mb = self.cache_file.stat().st_size / 1024**2

            if verbose:
                print(f"Loaded {f' and cached ({size_mb:.1f} MB)' if cache else ''} trajectory {self.dir.name} with {len(self.states)} states.")

        except Exception as e:
            print(f'Error loading trajectory from {str(self.dir)}: {e}')


    def unload(self, verbose=False):
        """
        Unload trajectory data.

        Parameters
        ----------
        verbose : bool, default False
            If True, print detailed loading information.
        """

        n_states = len(self.states)
        self.states = []
        self.times = np.array([])
        self.energies = np.array([])
        if verbose:
            print(f" {self.dir.name} with {n_states} states unloaded.")


    def clear_cache(self, verbose=False):
        """
        Clear cached trajectory data.

        Parameters
        ----------
        verbose : bool, default False
            If True, print detailed clearing information.
        """
        if self.cache_file.exists():
            self.cache_file.unlink()
            if verbose:
                print(f"Cleared cache for trajectory {self.dir.name}")



        
    # ==========================================================================
    # RDF (Radial Distribution Function) Analysis Methods
    # ==========================================================================

    def get_rdf(self, r_max=None, dr=0.1, fraction=1.0, g_ref=None):
        """
        Calculate radial distribution function averaged over trajectory.
        
        Parameters
        ----------
        r_max : float, optional
            Maximum distance for RDF (default: half of cell diagonal)
        dr : float, default 0.1
            Bin width in Angstroms
        fraction : float, default 1.0
            Fraction of trajectory to use for RDF calculation (e.g., 0.5 for last half)
        g_ref : ndarray, optional
            Reference RDF for normalization (from full lattice at coverage=1).
            If provided, normalizes by number of neighbors in each shell.
            
        Returns
        -------
        r_bins : ndarray
            Bin centers  
        g_r : ndarray
            RDF values normalized such that g(r)=1 for ideal gas
            
        Notes
        -----
        RDF is averaged over fraction of states.
        Normalization follows zacros_functions.py: divides counts by g_ref 
        (number of neighbors in each shell) and by coverage.
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

        # Determine which states to include based on fraction
        eq_idx = int((1.0 - fraction) * len(self.states))
        states_to_use = self.states[eq_idx:]
        
        # Get average coverage
        avg_coverage = np.mean([s.get_coverage() for s in states_to_use])
        
        # Accumulate over all states in states_to_use
        for st in states_to_use:
            occupied_coords = st.get_occupied_coords()
            n_occupied = len(occupied_coords)
            
            if n_occupied < 2:
                continue
            
            counts = np.zeros(n_bins, dtype=int)
            
            # Vectorized distance calculation (much faster for large n_occupied)
            distances = self.lattice.pairwise_distances_pbc(occupied_coords)
            # Get upper triangle (no diagonal, no double counting)
            mask = np.triu(np.ones(distances.shape, dtype=bool), k=1)
            valid_dists = distances[mask]
            valid_dists = valid_dists[(valid_dists > 0) & (valid_dists <= r_max)]
            
            # Histogram
            counts, _ = np.histogram(valid_dists, bins=bin_edges)
            
            # Normalize by g_ref if provided (number of neighbors in each shell)
            if g_ref is not None:
                counts_n = np.zeros_like(counts, dtype=float)
                np.divide(counts, g_ref, out=counts_n, where=g_ref!=0)
                g_r += counts_n / n_occupied / avg_coverage
            else:
                g_r += counts / n_occupied / avg_coverage
        
        # Normalize by number of states
        # Factor of 2 to account for unordered pairs
        if len(states_to_use) > 0:
            g_r = 2 * g_r / len(states_to_use)
                    
        return r_bins, g_r
        
    def get_cluster_size_freqs(self, cutoff=2.0, eps=1e-4, fraction=1.0, method='ckdtree'):
        """
        Calculate cluster size frequencies averaged over trajectory.
        
        Parameters
        ----------
        cutoff : float
            Cutoff distance for clustering.
            If float: explicit distance in Angstroms
        eps   : float, default 1e-4
            Tolerance parameter for the approximate search in the tree (see docs for scipy.spatial.ckdtree package)
            
        Returns
        -------
        cluster_sizes : ndarray
            Unique cluster sizes
        frequencies : ndarray
            Fraction of time each cluster size appears
            
        Notes
        -----
        Uses connected components algorithm with PBC-aware distances.
        """
            
        
        if method == 'ckdtree':
            func = 'get_clusters'
        elif method == 'bfs':
            func = 'get_clusters_bfs'
        else:
            print(f"Warning: unknown method '{method}' for cluster analysis. Defaulting to 'ckdtree'.")
            func = 'get_clusters'
       
        # Initialize cluster size counts (index is cluster size)
        cluster_size_freqs = np.zeros(len(self.lattice))
        
        # Determine number of snapshots to use based on fraction
        n_states = len(self.states)
        use_states = self.states[int(n_states * (1.0 - fraction)):]

        n_clusters = 0
        for st in use_states:
            labels, clusters, sizes = getattr(st, func)(cutoff, eps=eps)

            # Accumulate cluster sizes
            for size in sizes: 
                cluster_size_freqs[size-1] += 1
            n_clusters += len(sizes)

        cluster_size_freqs /= n_clusters  # Normalize by number of states  

        return cluster_size_freqs
        
    def get_accessibility_histogram(self, fraction=1.0):
        """
        Calculate histogram of site accessibility (number of vacant nearest neighbors).

        Parameters
        ----------
        fraction : float, default 1.0
            Fraction of snapshots to use (e.g., 0.5 for last half).
        
        Returns
        -------
        accessibility : ndarray
            Number of vacant nearest neighbors (0 to max_coordination)
        frequencies : ndarray
            Fraction of occupied sites with each accessibility
            
        Notes
        -----
        Accessibility measures how many nearest neighbor sites are vacant,
        which affects reactivity and diffusion rates.
        """
        accessibilities_13 = []
        accessibilities_2  = []
        
        # Determine number of snapshots to use based on fraction
        n_states = len(self.states)
        use_states = self.states[int(n_states * (1.0 - fraction)):]
        # Maximum site coordination number
        max_coord = np.max(self.lattice.site_coordinations)


        for st in use_states:
            acc13, acc2 = st.get_accessibility()
            accessibilities_13 += acc13
            accessibilities_2  += acc2

        counts = np.bincount(accessibilities_13, minlength=max_coord)
        frequencies_13 = counts / counts.sum()

        counts = np.bincount(accessibilities_2, minlength=max_coord)
        frequencies_2 = counts / counts.sum()

        return frequencies_13, frequencies_2
        

    def get_coverage_vs_time(self):
        """
        Get coverage as a function of time.
        
        Returns
        -------
        times : ndarray
            Time points
        coverages : ndarray
            Coverage at each time point
        """
        coverages = np.array([s.get_coverage() for s in self.states])
        return self.times, coverages
        
    def __len__(self):
        """
        Number of configurations in trajectory.
        
        Returns number of states if loaded, otherwise number of time points.
        This allows len() to work correctly for energy_only trajectories.
        """
        if len(self.states) > 0:
            return len(self.states)
        else:
            return len(self.states)
        
    def __getitem__(self, idx):
        """
        Access states by index.
        
        Parameters
        ----------
        idx : int or slice
            Index or slice for states
            
        Returns
        -------
        state or list of states
        
        Examples
        --------
        >>> traj[0]          # First state
        >>> traj[-1]         # Last state  
        >>> traj[10:20:2]    # Every other state from 10 to 20
        """
        return self.states[idx]
        
    def __repr__(self):
        """String representation of Trajectory"""
        if len(self) > 0:
            t_range = f"t=[{self.times[0]:.2f}, {self.times[-1]:.2f}]"
        else:
            t_range = "empty"
        return f"Trajectory(nstates={len(self)}, {t_range}, lattice={len(self.lattice)} sites)"
