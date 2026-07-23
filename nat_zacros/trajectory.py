"""
Trajectory class for sequences of lattice states over time.

This module provides the `trajectory` class for managing and analyzing
sequences of adsorbate configurations from Zacros KMC simulations.
"""

import pickle
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import copy

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
    
    def __init__(self, dir, lattice, metadata, cache_file='traj.pkl'):
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

        self.metadata = metadata

        self._necessary_keys = [
            'lattice_dimensions', 
            'surf_species_names',
            'n_adsorbates', 
            'temperature', 
            'energy_terms'
            ]

        args_ok = True
        for key in self._necessary_keys:
            if key not in self.metadata:
                print(f" {key} undefined.")
                args_ok = False

        if not args_ok:
            self.is_valid = False
            print('Not enough metadata for a valid simulation.')



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
      
                
    def load(self, cache=True, zfile='history_output', verbose=False):
        """
        Load states from the output defined by zfile.

        Parameters
        ----------
        cache : bool, default True
            If True, attempt to load cached trajectory data if available; cache trajectory data if not already cached.
            If False, load from raw simulation data.
        zfile: str, default 'history_output'
            selects zacros output file from which to read states
        verbose : bool, default False
            If True, print detailed loading information.
        """

        if zfile=='history_output':

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

                # Get the number of gas species
                n_gas_species = len(content[0].split()) - 1

                # Count total number of states in trajectory
                n_states = sum(line.lstrip().startswith('configuration') for line in content)

                if n_states == 0: # history_output.txt does not have configurations
                    n_states = 1
                    st = State(self.lattice, self.metadata, dirname=self.dir)

                    self.times = np.zeros(n_states, dtype=float)
                    self.energies = np.zeros(n_states, dtype=float)
                    self.states = [st] * n_states

                else:

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
                    if n_gas_species == 0:
                        expected_block_size = 1 + nsites  # 1 header + nsites data
                    else:
                        expected_block_size = 2 + nsites  # 1 header + nsites data + gas species line
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

                        st = State(self.lattice, self.metadata, dirname=self.dir)

                        for site in range(nsites):
                            site_line = content[pos + 1 + site]
                            p = site_line.split()
                            st.ads_ids[site] = int(p[1])
                            st.occupation[site] = int(p[2])
                            st.dentation[site] = int(p[3])
                        
                        if n_gas_species > 1:
                            st.gas_species_change = [int(s) for s in content[pos + 1 + nsites].split()]

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
        
        elif zfile=='general_output':
            
            try:
                # Get number of lattice sites
                nsites = len(self.lattice)
                # Read in a trajectory from general_output.txt
                with open(self.dir / f'{zfile}.txt', 'r') as f: 
                    content = f.readlines()

                # Start from the end to get number of states
                for line in reversed(content):

                    if 'Events occurred' in line:
                        # Get the number of states including initial state
                        n_states = int(line.split()[-1]) + 1
                        break

                # Allocate arrays
                self.times = np.empty(n_states, dtype=float)
                self.energies = np.empty(n_states, dtype=float)
                self.states = [None] * n_states

                # Get initial state
                st = State(self.lattice, self.metadata, dirname=self.dir)
                ads_counter = 0
                for line in content:
                    
                    if 'Surface species dentation' in line:
                        dent = [int(str) for str in line.split()[3:]]

                    if 'adparticle' in line:
                        l_split = line.split()
                        sp_idx = st.surf_species_names[ l_split[0] ]
                        ads_counter += 1

                        sites = [int(s) for s in l_split[-dent[sp_idx]:]]
                        for site in sites:
                            # Account for Zacros numbering of sites and species
                            st.occupation[site-1] = sp_idx + 1
                            st.dentation[site-1] = dent[sp_idx]
                            st.ads_ids[site-1] = ads_counter
                
                    if 'Total adlayer energy:' in line:
                        en = float(line.split()[-1])
                        # We need the initial state only
                        break

                self.states[0] = st
                self.times[0] = 0.0
                self.energies[0] = en

                # Finde the 1st line with kmc step
                for iline, line in enumerate(content):
                    if 'KMC step' in line: break                        

                # Loop to get states forming the trajectory
                for k in range(1, n_states):

                    # Get the state info:
                    # kmc step
                    step = int(line.split()[-1])
                    # reaction
                    reaction_name = content[iline+1].split()[-1]
                    # time
                    time = float(content[iline+2].split()[-1])
                    # sites involved
                    sites_inv = [ int(s) for s in content[iline+3].split()[2:] ]
                    # change of the energy
                    delta_en = float(content[iline+7].split()[-1])

                    # Shift the line index to the next KMC step record
                    iline += 9

                    # Get the previous state
                    # st = deepcopy(self.states[k-1])
                    st = State(self.lattice, self.metadata, dirname=self.dir)
                    st.occupation = copy.copy(self.states[k-1].occupation)
                    st.dentation  = copy.copy(self.states[k-1].dentation)
                    st.ads_ids    = copy.copy(self.states[k-1].ads_ids)

                    # modify it
                    if 'hopping' in reaction_name:
                        dent = len(sites_inv) // 2
                        for i in range(dent):
                            idx_r = sites_inv[i] - 1
                            idx_p = sites_inv[-i-1] - 1
                            if 'rev' in reaction_name:
                                idx_p = sites_inv[i] - 1
                                idx_r = sites_inv[-i-1] - 1

                            # save reactant data
                            occ = st.occupation[idx_r]
                            ads = st.ads_ids[idx_r]

                            # Remove reactant from the lattice
                            st.occupation[idx_r] = 0
                            st.dentation[ idx_r] = 0
                            st.ads_ids[   idx_r] = 0

                            # Put product on the lattice
                            st.occupation[idx_p] = occ
                            st.dentation[ idx_p] = dent
                            st.ads_ids[   idx_p] = ads

                    elif 'desorption' in reaction_name:
                        dent = len(sites_inv)
                        for i in range(dent):
                            idx_r = sites_inv[i] - 1

                            # Remove reactant from the lattice
                            st.occupation[idx_r] = 0
                            st.dentation[ idx_r] = 0
                            st.ads_ids[   idx_r] = 0
                    
                    else:
                        print('trajectory load function: general_output reading warning: reaction unknown.')

                    self.states[k] = st
                    self.times[k] = time
                    self.energies[k] = delta_en + self.energies[k-1]


            except Exception as e:
                print(f'Error loading trajectory from {str(self.dir)}: {e}')

        else:
            print(f'Error loading trajectory: output file {zfile} unknown.')


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

    def get_rdf(self, species_1, species_2, distances, r_max=None, dr=0.1, fraction=1.0, g_ref=None):
        """
        Calculate radial distribution function averaged over trajectory.
        
        Parameters
        ----------
        species_1, species_2 : str
            Species for which rdf is to be calculated
        distances : (N,N) array of floats
            Matrix of lattice site distances
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
            
        # Determine which states to include based on fraction
        eq_idx = int((1.0 - fraction) * len(self.states))
        states_to_use = self.states[eq_idx:]
        
        # Accumulate over all states in states_to_use
        r_bins, g_r_avg = states_to_use[0].get_rdf(species_1, species_2, distances, r_max=r_max, dr=dr, g_ref=g_ref)
        for st in states_to_use[1:]:
            r_bins, g_r = st.get_rdf(species_1, species_2, distances, r_max=r_max, dr=dr, g_ref=g_ref)
            g_r_avg += g_r

        # Normalize by number of states
        if len(states_to_use) > 0:
            g_r_avg = g_r_avg / len(states_to_use)
                    
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
        

    def get_leed_intensity_vs_time(self, distances, r):
        """
        Get LEED intensity for states in a trajectory

        Returns
        -------
        intensities : ndarray
            LEED intensities vs time
        """

        intensities = np.zeros(len(self))

        for i, st in enumerate(self.states):
            intensities[i] = st.get_leed_intensity(distances, r)

        return self.times, intensities
        

    def get_coverages_vs_time(self, atoms_per_uc=1):
        """
        Get coverages of surface species as a function of time.
        
        Returns
        -------
        times : ndarray
            Time points
        coverages : 2D-array (nstates, n_surf_species+1)
            Partial coverages at each time point
        """
        
        coverages = np.array([s.get_coverages(atoms_per_uc=atoms_per_uc) for s in self.states])

        return self.times, coverages

    def get_ads_path(self, idx, trim_zeros=False):
        """
        Get the path of an adsorbate idx
        """

        ads_path = np.empty(len(self), dtype=int)

        for i,st in enumerate(self.states):
            ads_path[i] = st.get_ads_site(idx)

        if trim_zeros:
            return np.trim_zeros(ads_path, 'b')
        else:
            return ads_path


    def animation(self, frames, interval=220, precision=2,
                  # options for plot function of the state class
                  scaling=1, ads_scaling=1.2, markers=None, colors=None, legend=True, show_axis=True):

        fig, ax = plt.subplots(figsize=(8, 6))

        def to_superscript_scientific(number, precision=2):
            # 1. Format the number into standard 'e' notation
            e_string = f"{number:.{precision}e}"
            
            # 2. Split the significand (mantissa) and the exponent
            significand, exponent = e_string.split("e")
            
            # 3. Clean up the exponent (convert to int to strip leading zeros and '+' signs)
            exponent_int = int(exponent)

            str = rf'${significand}\times 10^{{{exponent_int}}}$'
        
            return str

        def update(frame):
            ax.cla()
            self.states[frame].plot(ax=ax, scaling=scaling, ads_scaling=ads_scaling, 
                                    markers=markers, colors=colors, legend=legend, show_axis=show_axis)

            str = rf'$t_{{{frame}}} = $ {to_superscript_scientific(self.times[frame], precision=precision)} s'
            ax.text(0.03, 0.9, str, transform=ax.transAxes, ha='left', fontsize=10)

            return ax


        return FuncAnimation(fig, update, frames=frames, interval=interval, blit=False)


    def __len__(self):
        """
        Number of configurations in trajectory.
        
        Returns number of states if loaded, otherwise number of time points.
        This allows len() to work correctly for energy_only trajectories.
        """
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
            t_range = f"t=[{self.times[0]:.2e}, {self.times[-1]:.2e}]"
        else:
            t_range = "empty"
        return f"Trajectory(nstates={len(self)}, {t_range}, lattice={len(self.lattice)} sites)"
