"""
Trajectory class for sequences of lattice states over time.

This module provides the `trajectory` class for managing and analyzing
sequences of adsorbate configurations from Zacros KMC simulations.
"""

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
    load_equilibrated_states(fraction=0.5)
        Reload only equilibrated portion with full state data
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
    
    def __init__(self, lattice, dirname=None):
        """
        Initialize trajectory with lattice and optional data folder.
        
        Parameters
        ----------
        lattice : Lattice object
            The surface lattice for this trajectory
        dirname : str or Path, optional
            Directory containing history_output.txt
        """
        self.lattice = lattice
        self.states = []
        self.times = []
        self.energies = []
        #

        # self.folder = Path(dirname) if dirname else None
        # The line above caused a compatibility problem with pickle
        # Path objects are not always picklable in all environments
        # Convert to string for safety
        #
        self.folder = str(Path(dirname)) if dirname else None
        
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
      
    # def load(self, dirname=None, start=0, end=None, step=1, load_energy=True, energy_only=False, fraction=1.0):
    #     """
    #     Load states from history_output.txt.
    #     """
    #     folder_p = Path(dirname) if dirname else Path(self.folder) if self.folder else None
    #     try:
    #         if energy_only:
    #             # Fast path: scan file for configuration headers only
    #             if fraction < 1.0:
    #                 with open(folder_p / 'history_output.txt', 'r') as f:
    #                     n_total = sum(1 for line in f if line.lstrip().startswith('configuration'))
    #                 start = max(start, int((1.0 - fraction) * n_total))

    #             with open(folder_p / 'history_output.txt', 'r') as f:
    #                 idx = 0
    #                 for line in f:
    #                     if line.lstrip().startswith('configuration'):
    #                         if idx < start:
    #                             idx += 1
    #                             continue
    #                         if end is not None and idx >= end:
    #                             break
    #                         if (idx - start) % step != 0:
    #                             idx += 1
    #                             continue

    #                         parts = line.split()
    #                         time = float(parts[3])
    #                         energy = float(parts[5]) if load_energy and len(parts) > 5 else 0.0

    #                         self.times.append(time)
    #                         self.energies.append(energy)
    #                         idx += 1
    #         else:
    #             # Single-pass streaming parse (fast)
    #             nsites = len(self.lattice)

    #             if fraction < 1.0 or end is None:
    #                 with open(folder_p / 'history_output.txt', 'r') as f:
    #                     n_states = sum(1 for line in f if line.lstrip().startswith('configuration'))
    #             else:
    #                 n_states = end

    #             if fraction < 1.0:
    #                 start = max(start, int((1.0 - fraction) * n_states))
    #             if end is None:
    #                 end = n_states

    #             with open(folder_p / 'history_output.txt', 'r') as f:
    #                 idx = -1
    #                 for line in f:
    #                     if not line.lstrip().startswith('configuration'):
    #                         continue

    #                     idx += 1
    #                     parts = line.split()
    #                     time = float(parts[3])
    #                     energy = float(parts[5]) if load_energy and len(parts) > 5 else 0.0

    #                     if idx < start or idx >= end or ((idx - start) % step != 0):
    #                         for _ in range(nsites):
    #                             next(f, None)
    #                         continue

    #                     st = State(self.lattice)
    #                     st.folder = str(folder_p)

    #                     for site in range(nsites):
    #                         site_line = next(f, None)
    #                         if site_line is None:
    #                             raise RuntimeError("Unexpected end of file while reading state block.")
    #                         p = site_line.split()
    #                         st.ads_ids[site] = int(p[1])
    #                         st.occupation[site] = int(p[2])
    #                         st.dentation[site] = int(p[3])

    #                     self.states.append(st)
    #                     self.times.append(time)
    #                     self.energies.append(energy)

    #     except Exception as e:
    #         print(f'Error loading trajectory from {str(folder_p)}: {e}')

    #     self.times = np.array(self.times)
    #     self.energies = np.array(self.energies)

    # ...existing code...
                  
    def load(self, dirname=None, start=0, end=None, step=1, load_energy=True, energy_only=False, fraction=1.0):
        """
        Load states from history_output.txt.
        """
        folder_p = Path(dirname) if dirname else Path(self.folder) if self.folder else None
        try:
            if energy_only:
                # Count total configurations only if needed
                if fraction < 1.0 or end is None:
                    with open(folder_p / 'history_output.txt', 'r') as f:
                        n_total = sum(1 for line in f if line.lstrip().startswith('configuration'))
                else:
                    n_total = end

                if fraction < 1.0:
                    start = max(start, int((1.0 - fraction) * n_total))
                if end is None:
                    end = n_total

                n_keep = max(0, (end - start + (step - 1)) // step)
                self.times = np.empty(n_keep, dtype=float)
                self.energies = np.empty(n_keep, dtype=float)

                with open(folder_p / 'history_output.txt', 'r') as f:
                    idx = 0
                    k = 0
                    for line in f:
                        if line.lstrip().startswith('configuration'):
                            if idx < start:
                                idx += 1
                                continue
                            if end is not None and idx >= end:
                                break
                            if (idx - start) % step != 0:
                                idx += 1
                                continue

                            parts = line.split()
                            self.times[k] = float(parts[3])
                            self.energies[k] = float(parts[5]) if load_energy and len(parts) > 5 else 0.0
                            k += 1
                            idx += 1

                # Trim in case of mismatch
                self.times = self.times[:k]
                self.energies = self.energies[:k]

            else:
                # Single-pass streaming parse (fast)
                nsites = len(self.lattice)

                if fraction < 1.0 or end is None:
                    with open(folder_p / 'history_output.txt', 'r') as f:
                        n_states = sum(1 for line in f if line.lstrip().startswith('configuration'))
                else:
                    n_states = end

                if fraction < 1.0:
                    start = max(start, int((1.0 - fraction) * n_states))
                if end is None:
                    end = n_states

                n_keep = max(0, (end - start + (step - 1)) // step)
                self.times = np.empty(n_keep, dtype=float)
                self.energies = np.empty(n_keep, dtype=float)
                self.states = [None] * n_keep

                with open(folder_p / 'history_output.txt', 'r') as f:
                    idx = -1
                    k = 0
                    for line in f:
                        if not line.lstrip().startswith('configuration'):
                            continue

                        idx += 1
                        parts = line.split()
                        time = float(parts[3])
                        energy = float(parts[5]) if load_energy and len(parts) > 5 else 0.0

                        if idx < start or idx >= end or ((idx - start) % step != 0):
                            for _ in range(nsites):
                                next(f, None)
                            continue

                        st = State(self.lattice)
                        st.folder = str(folder_p)

                        for site in range(nsites):
                            site_line = next(f, None)
                            if site_line is None:
                                raise RuntimeError("Unexpected end of file while reading state block.")
                            p = site_line.split()
                            st.ads_ids[site] = int(p[1])
                            st.occupation[site] = int(p[2])
                            st.dentation[site] = int(p[3])

                        self.states[k] = st
                        self.times[k] = time
                        self.energies[k] = energy
                        k += 1

                # Trim in case of mismatch
                self.states = self.states[:k]
                self.times = self.times[:k]
                self.energies = self.energies[:k]

        except Exception as e:
            print(f'Error loading trajectory from {str(folder_p)}: {e}')
        
        
    def load_equilibrated_states(self, fraction=0.5, method='fraction', dirname=None):
        """
        Reload trajectory with full state data only for equilibrated portion.
        
       
        
        Parameters
        ----------
        fraction : float, default 0.5
            Fraction of trajectory to skip for equilibration
        method : str, default 'fraction'
            Method for equilibration detection
        dirname : str or Path, optional
            Override folder location
            
        Returns
        -------
        None
            Modifies self.states in place, clearing old states and loading
            only equilibrated configurations.
            
        Examples
        --------
        >>> # Phase 1: Fast energy-only loading
        >>> traj = Trajectory(lat, dirname)
        >>> traj.load(energy_only=True)
        >>> 
        >>> # Phase 2: Reload equilibrated states for analysis
        >>> traj.load_equilibrated_states(fraction=0.5)
        >>> r, g = traj.get_rdf()  # Now works with full state data
        
        Notes
        -----
        Requires that times/energies are already loaded (from energy_only mode).
        Will clear existing states and reload from file.
        """
        if len(self.times) == 0:
            raise RuntimeError("No trajectory data loaded. Run load() first.")
            
        # Determine equilibration index
        eq_idx = self.estimate_equilibration(fraction=fraction, method=method)
        
        # Clear existing states
        self.states = []
        
        # Reload only equilibrated portion with full state data
        folder = Path(dirname) if dirname else Path(self.folder) if self.folder else None
        
        if folder is None:
            raise RuntimeError('Error: folder not specified')
        
        try:
            # Reload with full state parsing, starting from equilibration point
            # Keep existing times/energies, just populate states
            for idx in range(eq_idx, len(self.times)):
                st = State(self.lattice)
                st.folder = folder
                st.load(idx=idx)
                self.states.append(st)
                
        except Exception as e:
            print(f'Error loading equilibrated states from {str(folder)}: {e}')
    
    # ==========================================================================
    # RDF (Radial Distribution Function) Analysis Methods
    # ==========================================================================
    #
    # PERFORMANCE NOTES:
    # ------------------
    # The RDF methods have been heavily optimized through vectorization.
    # Key performance characteristics:
    #
    # 1. get_rdf() with vectorized=True (default):
    #    - Uses lattice.pairwise_distances_pbc() for vectorized distance calculations
    #    - 50-100x faster than nested Python loops
    #    - Typical: 2s for 10 trajectories × 100 states
    #
    # 2. Sequential vs Parallel:
    #    - Sequential loop over trajectories: ~2s for typical datasets
    #    - Parallel RDF functions: ~3-5s (slower due to 2-4s overhead!)
    #    - Recommendation: Use simple sequential loop for RDF computation
    #
    # 3. When parallel RDF helps:
    #    - Only beneficial when computation time >> 30 seconds
    #    - Typically: >50 trajectories or very large systems
    #    - Always benchmark before using parallel functions
    #
    # See module docstring for complete performance analysis.
    # ==========================================================================
    
    def get_g_ref(self, r_max=None, dr=0.1, vectorized=True):
        """
        Calculate reference RDF for full lattice (all sites, coverage=1).
        
        This computes the number of neighbors in each distance shell,
        used to normalize the RDF such that g(r)=1 for ideal gas.
        
        Parameters
        ----------
        r_max : float, optional
            Maximum distance for RDF
        dr : float, default 0.1
            Bin width in Angstroms
        vectorized : bool, default True
            Use vectorized distance calculation for all pairs
            
        Returns
        -------
        r_bins : ndarray
            Bin centers
        g_ref : ndarray
            Number of neighbors in each shell (integer counts)
        """
        if r_max is None:
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
        # Get all lattice site coordinates
        all_coords = self.lattice.coordinates
        n_sites = len(all_coords)
        counts = np.zeros(n_bins, dtype=int)
        if vectorized:
            # Vectorized calculation using pairwise_distances_pbc
            dists_matrix = self.lattice.pairwise_distances_pbc(all_coords)
            mask = np.triu(np.ones(dists_matrix.shape, dtype=bool), k=1)
            dists = dists_matrix[mask]
            valid_dists = dists[(dists > 0) & (dists <= r_max)]
            counts, _ = np.histogram(valid_dists, bins=bin_edges)
        else:
            # Original nested loop
            for i in range(n_sites - 1):
                for j in range(i + 1, n_sites):
                    dist = self.lattice.minimum_image_distance(
                        all_coords[i], all_coords[j]
                    )
                    if 0 < dist <= r_max:
                        bin_idx = int(dist / dr)
                        if bin_idx < n_bins:
                            counts[bin_idx] += 1
        # Normalize: 2 * counts / n_sites (factor 2 for unordered pairs)
        g_ref = 2.0 * counts / n_sites
        return r_bins, g_ref
        
    def get_rdf(self, r_max=None, dr=0.1, g_ref=None, vectorized=True):
        """
        Calculate radial distribution function averaged over trajectory.
        
        Parameters
        ----------
        r_max : float, optional
            Maximum distance for RDF (default: half of cell diagonal)
        dr : float, default 0.1
            Bin width in Angstroms
        g_ref : ndarray, optional
            Reference RDF for normalization (from full lattice at coverage=1).
            If provided, normalizes by number of neighbors in each shell.
        vectorized : bool, default True
            Use vectorized distance calculations for better performance
            
        Returns
        -------
        r_bins : ndarray
            Bin centers  
        g_r : ndarray
            RDF values normalized such that g(r)=1 for ideal gas
            
        Notes
        -----
        RDF is calculated for occupied sites only and averaged over all states.
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
        
        # Get average coverage
        avg_coverage = np.mean([s.get_coverage() for s in self.states])
        
        # Accumulate over all states
        for st in self.states:
            occupied_coords = st.get_occupied_coords()
            n_occupied = len(occupied_coords)
            
            if n_occupied < 2:
                continue
            
            counts = np.zeros(n_bins, dtype=int)
            
            if vectorized:
                # Vectorized distance calculation (much faster for large n_occupied)
                distances = self.lattice.pairwise_distances_pbc(occupied_coords)
                # Get upper triangle (no diagonal, no double counting)
                mask = np.triu(np.ones(distances.shape, dtype=bool), k=1)
                valid_dists = distances[mask]
                valid_dists = valid_dists[(valid_dists > 0) & (valid_dists <= r_max)]
                
                # Histogram
                counts, _ = np.histogram(valid_dists, bins=bin_edges)
            else:
                # nested loop implementation
                for i in range(n_occupied - 1):
                    for j in range(i + 1, n_occupied):
                        dist = self.lattice.minimum_image_distance(
                            occupied_coords[i], occupied_coords[j]
                        )
                        if 0 < dist <= r_max:
                            bin_idx = int(dist / dr)
                            if bin_idx < n_bins:
                                counts[bin_idx] += 1
            
            # Normalize by g_ref if provided (number of neighbors in each shell)
            if g_ref is not None:
                counts_n = np.zeros_like(counts, dtype=float)
                np.divide(counts, g_ref, out=counts_n, where=g_ref!=0)
                g_r += counts_n / n_occupied / avg_coverage
            else:
                g_r += counts / n_occupied / avg_coverage
        
        # Normalize by number of states
        # Factor of 2 to account for unordered pairs
        if len(self.states) > 0:
            g_r = 2 * g_r / len(self.states)
                    
        return r_bins, g_r
        
    def get_cluster_distribution(self, nn_cutoff=1):
        """
        Calculate cluster size distribution averaged over trajectory.
        
        Parameters
        ----------
        nn_cutoff : int or float
            Nearest neighbor distance cutoff for clustering.
            If int: nth nearest neighbor distance
            If float: explicit distance in Angstroms
            
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
        if isinstance(nn_cutoff, int):
            cutoff_dist = self.lattice.get_nn_distance(nn_cutoff) * 1.1  # 10% tolerance
        else:
            cutoff_dist = nn_cutoff
            
        all_clusters = []
        
        for st in self.states:
            occupied_sites = st.get_occupied_sites()
            n_occupied = len(occupied_sites)
            
            if n_occupied == 0:
                continue
                
            # Build adjacency matrix
            occupied_coords = st.get_occupied_coords()
            clusters = []
            visited = np.zeros(n_occupied, dtype=bool)
            
            for i in range(n_occupied):
                if visited[i]:
                    continue
                    
                # Start new cluster with BFS
                cluster = []
                queue = [i]
                visited[i] = True
                
                while queue:
                    current = queue.pop(0)
                    cluster.append(current)
                    
                    # Check neighbors
                    for j in range(n_occupied):
                        if not visited[j]:
                            dist = self.lattice.minimum_image_distance(
                                occupied_coords[current], occupied_coords[j]
                            )
                            if dist < cutoff_dist:
                                visited[j] = True
                                queue.append(j)
                                
                clusters.append(len(cluster))
                
            all_clusters.extend(clusters)
            
        # Calculate distribution
        if len(all_clusters) > 0:
            unique_sizes, counts = np.unique(all_clusters, return_counts=True)
            frequencies = counts / counts.sum()
            return unique_sizes, frequencies
        else:
            return np.array([]), np.array([])
        
    def get_accessibility_histogram(self):
        """
        Calculate histogram of site accessibility (number of vacant nearest neighbors).
        
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
        all_accessibility = []
        
        for st in self.states:
            occupied_sites = st.get_occupied_sites()
            
            for site_idx in occupied_sites:
                # Get nearest neighbors for this site
                nn_sites = self.lattice.site_nns[site_idx]
                
                # Count vacant neighbors
                vacant_nn = np.sum(st.occupation[nn_sites] == 0)
                all_accessibility.append(vacant_nn)
                
        # Calculate histogram
        if len(all_accessibility) > 0:
            max_coord = np.max(self.lattice.site_coordinations)
            accessibility = np.arange(max_coord + 1)
            counts = np.zeros(max_coord + 1)
            
            for val in all_accessibility:
                counts[val] += 1
                
            frequencies = counts / counts.sum()
            return accessibility, frequencies
        else:
            return np.array([]), np.array([])
        
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
            return len(self.times)
        
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
