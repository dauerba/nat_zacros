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
from .lattice import Lattice
from .state import State

class Trajectory:
    """
    Store one Zacros trajectory as a time-ordered sequence of lattice states.

    Depending on how it is loaded, a trajectory can hold either full lattice
    snapshots from ``history_output.txt`` or event-style state deltas derived
    from ``general_output.txt``.

    Attributes
    ----------
    dir : pathlib.Path
        Directory containing the Zacros output files for this trajectory.
    lattice : Lattice
        Lattice shared with the parent simulation.
    metadata : dict
        Metadata required to interpret species identifiers and coverage.
    states : list[State]
        Loaded state snapshots. Populated when reading ``history_output``.
    initial_state : State or None
        Initial state reconstructed from ``general_output.txt`` when available.
    state_deltas : list[tuple[str, list[int]]]
        Event records stored as ``(reaction_name, sites_involved)`` tuples when
        reading ``general_output``.
    times : ndarray
        Simulation times corresponding to the loaded snapshots or event records.
    energies : ndarray
        Total energies aligned with ``times``.
    cache_file : pathlib.Path
        Pickle file used to cache parsed ``history_output`` snapshots.
    """
    
    def __init__(self, traj_path, lattice=None, input_pars=None, reactions=None, 
                 cache_file='traj.pkl', verbose=False):
        """
        Initialize a trajectory container.
        
        Parameters
        ----------
        dir : str or Path, optional
            Directory containing the trajectory output files.
        lattice : Lattice object
            The surface lattice used by this trajectory.
        metadata : dict
            Metadata required to interpret surface species.
        cache_file : str, default 'traj.pkl'
            File name of the pickle cache stored inside ``dir``.
        """

        self.traj_path = Path(traj_path)

        # Load lattice from lattice input and output files if necessary
        if lattice is None:
            self.lattice = Lattice(self.traj_path)
        else:
            self.lattice = lattice

        # Parse input files if necessary
        if input_pars is None:
            self.input_pars = self._get_input_pars(verbose=verbose)
        else:
            self.input_pars = input_pars

        self.n_gas_species = int(self.input_pars['n_gas_species'][0])
        self.gas_species_names = self.input_pars['gas_specs_names']
        self.n_surf_species = int(self.input_pars['n_surf_species'][0])
        self.surf_species_names = dict(zip(self.input_pars['surf_specs_names'], range(1,self.n_surf_species+1)))
        self.surf_species_names['*'] = 0
        self.surf_species_dent = dict(zip(self.input_pars['surf_specs_names'], [int(d) for d in self.input_pars['surf_specs_dent']]))
        self.surf_species_dent['*'] = 1

        # Get the initial seed from state_input
        self.nads0 = self._get_initial_nads() if len(self._get_initial_nads()) > 0 \
                                              else dict(zip(self.input_pars['surf_specs_names'], [0] * self.n_surf_species))

        if reactions is None:
            # Initialize reactions atribute
            self.reactions_surf = dict()
            self.reactions_gas =  dict()
            # Parse reaction mechanism input
            reactions_surf, reactions_gas = self._parse_reaction_mechanism(verbose=verbose)

            # Converts species names to the corresponding indices
            # and dentations to int in the reaction mechanism dict
            for r in reactions_surf:
                self.reactions_surf[r] = np.array([[self.surf_species_names[sp_name], dent] for sp_name, dent in reactions_surf[r]],
                                                    dtype=int)

            # reactions_gas attribute has the following structure: {reaction name : {gas_species_name : its change} }
            # Set values for for all gas species changes to 0
            for r in reactions_surf:
                self.reactions_gas[r] = dict(zip(self.gas_species_names, [0] * self.n_gas_species))
            # Get the actual values from the parsed information
            for r in reactions_gas:
                for sp in reactions_gas[r]:
                    self.reactions_gas[r][sp] = reactions_gas[r][sp]
        else:
            self.reactions_surf, self.reactions_gas = reactions

        if verbose:
            print(f'\n  reactions_surf attribute:')
            for k,v in self.reactions_surf.items():
                print(f'   {k}:\n {v}')
            print(f'\n  reactions_gas attribute:')
            for k,v in self.reactions_gas.items():
                print(f'   {k}:\n {v}')

        self.initial_state = None
        self.delta_states = []
        self.gas_species_changes = dict()
        self.times = []
        self.energies = []
        self.cache_file = self.traj_path / cache_file


    def create_states(self, idx):
        """
        Create a new state by applying reaction deltas to the current state.
        
        Parameters
        ----------
        deltas : a tuple or a list of tuples
            Each tuple contains (reaction_name, sites_inv) representing
            the reaction and the involved sites that change the state.
            
        Returns
        -------
        State or None
            New state after applying the supplied deltas, or ``None`` if one of
            the reaction names is not recognized so no changes are defined.

        
        """

        if isinstance(idx, int):
            start = 0
            stop = idx + 1
            step  = stop
        elif isinstance(idx, slice):
            start = idx.start or 0          
            stop  = idx.stop  or len(self)
            step  = idx.step  or 1
            if start == 0:
                states = [ self.initial_state ]
            else:
                states = []
        else:
            raise TypeError("Invalid index type. Must be int or slice.")

        # Initialize a new state object with the same lattice and metadata
        new_state = State(self.lattice, self.input_pars, self.traj_path)

        # Initialize arrays for the new state
        new_state.occupation = np.copy(self.initial_state.occupation)
        new_state.dentation  = np.copy(self.initial_state.dentation)

        # Loop over deltas up to the stop state
        for i in range(stop-1):
            sites_inv = self.delta_states[i][0,:]
            new_state.occupation[sites_inv] = self.delta_states[i][1]
            new_state.dentation[sites_inv]  = self.delta_states[i][2]
            if i + 1 >= start and (i + 1) % step == 0:
                # Initialize a new state object with the same lattice and metadata
                state_2_save = State(self.lattice, self.input_pars, self.traj_path)
                # Copy arrays from the new state
                state_2_save.occupation = np.copy(new_state.occupation)
                state_2_save.dentation  = np.copy(new_state.dentation)
                states.append(state_2_save)

        if isinstance(idx, int):
            return new_state
        else:
            return states


    def _get_initial_nads(self, verbose=False):

        # Load state_input.dat 
        # accounting for comments marked with '#'
        nads = dict()
        try:
            with open(self.traj_path / 'state_input.dat') as f:
                for line in f:
                    lsplit = line.partition('#')[0].strip().split()
                    if lsplit:
                        if lsplit[0] == 'seed_on_sites':
                            sp_name = lsplit[1]
                            nads[sp_name] = nads[sp_name] + 1 if sp_name in nads.keys() else 1
                        if lsplit[0] == 'seed_multiple':
                            sp_name = lsplit[1]
                            n = int(lsplit[2])
                            nads[sp_name] = nads[sp_name] + n if sp_name in nads.keys() else n
        except:
            if verbose:
                print(f"  Cannot get seeding info from state_input.dat")

        if verbose:
            print(f'  Number of adsorbates seeded:\n')
            print(nads)

        return nads
        

    def _get_input_pars(self, verbose=False):

        # Load simulation_input.dat as a dict
        # accounting for comments marked with '#'
        # assuming the first entry to be a key

        input_pars = dict()
        try:
            with open(self.traj_path / 'simulation_input.dat') as f:
                for line in f:
                    lsplit = line.partition('#')[0].strip().split()
                    if lsplit:
                        input_pars[lsplit[0]] = lsplit[1:]

            # Add gas species info if not present
            if input_pars.get('gas_specs_names') is None:
                input_pars['gas_specs_names'] = []
                input_pars['n_gas_species'] = [0]
        except:
            self.is_valid = False
            print(f"Cannot get simulation input parameters from simulation_input.dat")

        if verbose:
            print(f'Input parameters read for trajectory: {self.traj_path}:')
            for k, v in input_pars.items():
                print(f'  {k} : {v}')

        return input_pars
        
    def _parse_reaction_mechanism(self, verbose=False):
        # Load mechanism_input.dat

        # Initialize reaction dictionaries with reaction names as keys
        # values:
        #   surf: a list with length of number of involved sites of lists consisting of a species name and dentate as a str
        #   gas:  a dict with a gas species names as a key and its change as a value (int)
        reactions_surf = dict()
        reactions_gas  = dict()
        try:
            with open(self.traj_path /'mechanism_input.dat') as f:
                content = f.readlines()
                reaction_name_rev = None # Used to account for reversibility
                
                for iline, line in enumerate(content):
                    lsplit = line.partition('#')[0].strip().split()
                    if lsplit:
                        # Irreversible reaction
                        if lsplit[0] == 'step':
                            reaction_name = lsplit[1]

                        # Reaction is Reversible if it is marked as 'reversible step'
                        if lsplit[0] == 'reversible_step':
                            reaction_name     = lsplit[1] + '_fwd'
                            reaction_name_rev = lsplit[1] + '_rev'

                        # Reaction has gas products     
                        if lsplit[0] == 'gas_reacs_prods':
                            reactions_gas[reaction_name] = dict(zip(lsplit[1::2], [int(i) for i in lsplit[2::2]]))
                            
                            # if the gas products reaction is reversible, reaction_name_rev will exist 
                            # and changes in gas species will be inverted
                            if  reaction_name_rev is not None:
                                reactions_gas[reaction_name_rev]  = dict(zip(lsplit[1::2], [-int(i) for i in lsplit[2::2]]))

                        # Get the number of sites involved in the reaction
                        if lsplit[0] == 'sites':
                            n_sites = int(lsplit[1])

                        # for reversible case, take the initial state to be the final state
                        if lsplit[0] == 'initial' and reaction_name_rev is not None:
                                reactions_surf[reaction_name_rev] = [l.split()[1:] for l in content[iline+1:iline+1+n_sites]]
                        
                        # Get the final state
                        if lsplit[0] == 'final':
                            reactions_surf[reaction_name] = [l.split()[1:] for l in content[iline+1:iline+1+n_sites]]
                            reaction_name_rev = None
                            
        except:
            self.is_valid = False
            print(f"Cannot get reaction mechanism from mechanism_input.dat")

        if verbose:
            print(f'\nReaction mechanism parameters for trajectory {self.traj_path}:')
            print(f'  surface:')
            for k, v in reactions_surf.items():
                print(f'    {k} : {v}')
            print(f'  gas:')
            for k, v in reactions_gas.items():
                print(f'    {k} : {v}')
            print()

        return reactions_surf, reactions_gas


    def load(self, cache=True, verbose=False):
        """
        Load trajectory data from Zacros output files.

        The loader supports two modes controlled by ``zfile``:
        - ``'history_output'`` (default): reads full per-site configurations from ``history_output.txt``
          and populates ``self.states``, ``self.times``, and ``self.energies``. When ``cache=True`` the
          parsed trajectory is saved to and loaded from a pickle cache file (``traj.pkl`` by default).
        - ``'general_output'``: reads event/delta-style records from ``general_output.txt`` and
          populates ``self.times``, ``self.energies`` and ``self.state_deltas``. The initial state
          is read (if available) from ``general_output.txt`` earlier in this method.

        Parameters
        ----------
        cache : bool, default True
            If True, attempt to load a cached trajectory (pickle) when available and write the cache
            after successfully parsing ``history_output.txt``. Ignored for ``general_output`` mode.
        zfile : str, default 'history_output'
            Which Zacros output to read: either ``'history_output'`` or ``'general_output'``.
        verbose : bool, default False
            If True, print detailed loading and caching information.

        Notes
        -----
        - The method sets ``self.initial_state`` when an initial configuration can be parsed from
          ``general_output.txt``.
        - For ``history_output`` mode the per-snapshot loading uses :meth:`State.load` semantics.
                - ``general_output`` mode populates ``times``, ``energies``, and
                    ``state_deltas`` but does not populate ``states``.
        """


        if verbose:
            print(f'  Loading trajectory {str(self.traj_path)}...')

        # Try loading from cache
        if cache:
            if self.cache_file.exists():
                with open(self.cache_file, 'rb') as f:
                    self.initial_state, self.delta_states, self.gas_species_changes, self.times, self.energies = pickle.load(f)
                if verbose:
                    size_mb = self.cache_file.stat().st_size / 1024**2
                    print(f"    Loaded trajectory {self.traj_path.name}: {self} from cache ({size_mb:.1f} MB).")
                return

        # Auto-selection of the trajectory source

        # Flag indicating succesful loading
        is_loaded = False

        if self.input_pars['event_report'][0] == 'on': # Implies that the source could be 'general_output.txt'

            if verbose:
                print(f'    Trying to load from general_output.txt...')

            # Try to load states from from 'general_output.txt'
            try:
                with open(self.traj_path / 'general_output.txt', 'r') as f: 
                    content = f.readlines()

                # Initialize the initial state object
                st = State(self.lattice, self.input_pars, traj_path=self.traj_path)

                # Get the initial state
                ads_counter = 0
                for line in content:

                    if 'adparticle to be seeded on site' in line:
                        l_split = line.split()
                        sp_idx = st.surf_species_names[ l_split[0] ]
                        dent = st.surf_species_dent[ l_split[0] ]
                        sites = [int(s) for s in l_split[-dent:]]
                        for site in sites:
                            # Account for Zacros numbering of sites and species
                            st.occupation[site-1] = sp_idx
                            st.dentation[site-1]  = dent

                        ads_counter += 1

                    if 'Total adlayer energy:' in line:
                        en = float(line.split()[-1])
                        # We need the initial state only
                        break

                if ads_counter > 0:
    
                    # Initial state read in successfully
                    self.initial_state = st

                    # Proceed with reading kmc events

                    # Get number of lattice sites
                    nsites = self.lattice.n_sites

                    # Start from the end to get number of states
                    for line in reversed(content):

                        if 'Events occurred' in line:
                            # Get the number of states
                            n_states = int(line.split()[-1]) + 1
                            break

                    # Allocate arrays
                    self.times = np.empty(n_states, dtype=float)
                    self.energies = np.empty(n_states, dtype=float)
                    for gas_species in self.gas_species_names:
                        self.gas_species_changes[gas_species] = np.zeros(n_states - 1, dtype=int)

                    # Save initial values for time and energy
                    self.times[0]    = 0.0
                    self.energies[0] = en

                    # Find the 1st line with kmc step
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
                        sites_inv = [ int(s) - 1 for s in content[iline+3].split()[2:] ]
                        # change of the energy
                        delta_en = float(content[iline+7].split()[-1])

                        # Shift the line index to the next KMC step record
                        iline += 9
                        self.times[k] = time
                        self.energies[k] = delta_en + self.energies[k-1]
                        # Get the changes and append them to deltas
                        occs  = np.array(self.reactions_surf[reaction_name][:,0], dtype=int)
                        dents = np.array(self.reactions_surf[reaction_name][:,1], dtype=int)

                        # Append changes to deltas
                        self.delta_states.append( np.array([sites_inv, occs, dents]) )
                        for gas_species in self.gas_species_names:
                            self.gas_species_changes[gas_species][k-1] = self.reactions_gas[reaction_name][gas_species]

                    is_loaded = True

                else:
                    if verbose:
                        print(f'      Failed to load the initial state.')

            except FileExistsError:
                if verbose:
                    print(f"      Loading failed.")
                    print(f'         Error:{FileExistsError}')

        else: # event_report off: loading from history_output.txt

            try:
                # Single-pass streaming parse
                nsites = self.lattice.n_sites
                # Read in a trajectory from history_output.txt
                with open(self.traj_path / 'history_output.txt', 'r') as f: 
                    content = f.readlines()

                # Count total number of states in trajectory
                n_states = sum(line.lstrip().startswith('configuration') for line in content)

                for gas_species in self.gas_species_names:
                    self.gas_species_changes[gas_species] = np.zeros(n_states - 1, dtype=int)

                if n_states > 0:

                    # Allocate memory for class' attributes
                    self.times = np.empty(n_states, dtype=float)
                    self.energies = np.empty(n_states, dtype=float)

                    # Get the size of the configuration block

                    # Find the first configuration line
                    pos1 = 0
                    while pos1 < len(content) and not content[pos1].lstrip().startswith('configuration'):
                        pos1 += 1
                    # Find the second configuration line to determine block size
                    pos2 = pos1 + 1
                    while pos2 < len(content) and not content[pos2].lstrip().startswith('configuration'):
                        pos2 += 1
                    block_size = pos2 - pos1  # Lines per configuration block
                    # Check consistency
                    if self.n_gas_species == 0:
                        expected_block_size = 1 + nsites  # 1 header + nsites data
                    else:
                        expected_block_size = 2 + nsites  # 1 header + nsites data + gas species line
                    if block_size != expected_block_size:
                        raise ValueError(f'block size: {block_size}, expected {expected_block_size}.')

                    # Get initial state
                    line = content[pos1]
                    try:
                        parts = line.split()
                        time = float(parts[3])
                        energy = float(parts[5])
                    except (ValueError, IndexError):
                        raise ValueError(f'{str(self.traj_path.name)}: Failed to parse line: {line.strip()}')

                    st = State(self.lattice, self.input_pars, traj_path=self.traj_path)

                    # Get the 1st configuration as initial state
                    for site in range(nsites):
                        site_line = content[pos1 + 1 + site]
                        p = site_line.split()
                        st.occupation[site] = int(p[2])
                        st.dentation[site] = int(p[3])
                    
                    self.initial_state = st

                    self.times[0]    = time
                    self.energies[0] = energy

                    # Create arrays for state propagation
                    occupation_prev = np.copy(self.initial_state.occupation)
                    occupation      = np.empty_like(occupation_prev)
                    dentation       = np.empty_like(occupation)

                    # Get the rest of trajectory
                    for k in range(1, n_states):

                        pos = pos1 + k * block_size
                        line = content[pos]

                        try:
                            parts = line.split()
                            time = float(parts[3])
                            energy = float(parts[5])
                        except (ValueError, IndexError):
                            raise ValueError(f'{str(self.traj_path.name)}: Failed to parse line: {line.strip()}')

                        for site in range(nsites):
                            site_line = content[pos + 1 + site]
                            p = site_line.split()
                            occupation[site] = int(p[2])
                            dentation[site] = int(p[3])

                        # Get sites involved into a change and append changes to deltas
                        sites_inv = np.flatnonzero(occupation_prev == occupation)
                        self.delta_states.append( np.array([sites_inv, occupation[sites_inv], dentation[sites_inv]]) )
                        for i, gas_species in enumerate(self.gas_species_names):
                            self.gas_species_changes[gas_species][k-1] = int(content[pos + 1 + nsites].split()[i])

                        if self.n_gas_species > 0:
                            self.gas_species_changes[k-1,:] = [int(s) for s in content[pos + 1 + nsites].split()]
                        # Propagate occupation array
                        occupation_prev[sites_inv] = occupation[sites_inv]

                        self.times[k] = time
                        self.energies[k] = energy

                is_loaded = True

            except FileExistsError:
                if verbose:
                    print(f"      Loading failed.")
                    print(f'         Error:{FileExistsError}')

        # Save to cache
        if cache:
            with open(self.cache_file, 'wb') as f:
                pickle.dump([self.initial_state, self.delta_states, self.gas_species_changes, self.times, self.energies], f)
                size_mb = self.cache_file.stat().st_size / 1024**2

        if verbose:
            print(f"    Loaded {f'and cached ({size_mb:.1f} MB)' if cache else ''} {self.traj_path.name}: {self}.")

        


    def unload(self, verbose=False):
        """
        Unload trajectory data.

        Parameters
        ----------
        verbose : bool, default False
            If True, print detailed loading information.
        """

        self.initial_state = None
        self.delta_states = []
        self.gas_species_changes = dict()
        self.times = []
        self.energies = []

        if verbose:
            print(f" {self.traj_path.name}: {self} unloaded.")


    def clear_cache(self, verbose=False):
        """
        Clear cached trajectory data.

        Parameters
        ----------
        verbose : bool, default False
            If True, print detailed clearing information.
            
        Returns
        -------
        int
            1 if a cache file was deleted, 0 otherwise
        """
        if self.cache_file.exists():
            self.cache_file.unlink()
            if verbose:
                print(f"Cleared cache for trajectory {self.traj_path.name}")
            return 1
        else:
            if verbose:
                print(f"  Cache for trajectory {self.traj_path.name} not found")

        return 0


    def get_energy_vs_time(self, verbose=False):
        """
        Get energy as a function of time.

        Returns
        -------
        times : ndarray
            Time points
        energies : ndarray
            Energies at each time point
        """

        if verbose and len(self.times) == 0:
            print(f'Trajectory is empty: {self}')
            
        return self.times, self.energies
      
        
    # ==========================================================================
    # RDF (Radial Distribution Function) Analysis Methods
    # ==========================================================================

    def get_rdf(self, species_1, species_2, distances, r_max=None, dr=0.1, fraction=1.0, g_ref=None):
        """
        Calculate the RDF averaged over a selected fraction of snapshots.
        
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
        ``fraction`` selects the final part of the trajectory. For example,
        ``fraction=0.5`` averages over the last half of the loaded states.
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
        eq_idx = int((1.0 - fraction) * len(self))
        states_to_use = self[eq_idx:]
        
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
        cluster_size_freqs = np.zeros(self.lattice.n_sites)
        
        # Determine number of snapshots to use based on fraction
        n_states = len(self)
        use_states = self[int(n_states * (1.0 - fraction)):]

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
        Calculate accessibility histograms over a selected fraction of states.

        Parameters
        ----------
        fraction : float, default 1.0
            Fraction of snapshots to use (e.g., 0.5 for last half).
        
        Returns
        -------
        frequencies_13 : ndarray
            Histogram for vacant first-plus-third-shell neighbors.
        frequencies_2 : ndarray
            Histogram for vacant second-shell neighbors.
            
        Notes
        -----
        Accessibility measures how many nearest neighbor sites are vacant,
        which affects reactivity and diffusion rates.
        """
        accessibilities_13 = []
        accessibilities_2  = []
        
        # Determine number of snapshots to use based on fraction
        n_states = len(self)
        use_states = self[int(n_states * (1.0 - fraction)):]
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

        for i in range(len(self)):
            intensities[i] = self[i].get_leed_intensity(distances, r)

        return self.times, intensities
        

    def get_coverages_vs_time(self, coverage_norm=1, verbose=False):
        """
        Get coverages of surface species as a function of time.
        
        Returns
        -------
        times : ndarray
            Time points
        coverages : 2D-array (nstates, n_surf_species+1)
            Partial coverages at each time point
        """

        def cov_iter():
            for st in self[:]:
                yield st.get_nads()

        return self.times, np.vstack(list(cov_iter())) / (coverage_norm*self.lattice.n_cells)

    def get_ads_path(self, idx, trim_zeros=False, verbose=False):
        """
        Get the path of an adsorbate idx
        """

        n_states = len(self)

        def ads_site_iter():
            for st in self[:]:
                yield st.get_ads_site(idx)

        ads_path = np.fromiter(ads_site_iter(), dtype=int, count=n_states)

        if trim_zeros:
            return np.trim_zeros(ads_path, 'b')
        else:
            return ads_path

    def get_2d_lattice_indices(self, ads_id, expand = False):
        """
        Convert 1D lattice site indices to 2D matrix coordinates for a specific adsorbate.
        
        Parameters
        ----------
        ads_id : int
            The unique adsorbate ID to track along the trajectory
        expand: boolean, default: False
            if True, expands the lattice beyond periodic boundaries
        
        Returns
        -------
        cols : array of int
            Column indices (x-coordinates) of the adsorbate at each state
        rows : array of int
            Row indices (y-coordinates) of the adsorbate at each state
            
        """

        # Get sites occupied by adsorbate with ads_id
        sites = self.get_ads_path(ads_id, trim_zeros=True)

        # Convert site indices to the 2D form    
        cols = self.lattice.cols[sites]
        rows = self.lattice.rows[sites]

        if expand:

            # Get differences in cols and rows
            diffs_col = np.diff(cols)
            diffs_row = np.diff(rows)
            
            wrap_jumps_col = np.where(diffs_col == self.lattice.n_cols - 1, -self.lattice.n_cols, \
                             np.where(diffs_col == 1 - self.lattice.n_cols, self.lattice.n_cols, 0))
            wrap_jumps_row = np.where(diffs_row == self.lattice.n_rows - 1, -self.lattice.n_rows, \
                             np.where(diffs_row == 1 - self.lattice.n_rows, self.lattice.n_rows, 0))

            cumulative_shift_col = np.concatenate([[0], np.cumsum(wrap_jumps_col)])
            cumulative_shift_row = np.concatenate([[0], np.cumsum(wrap_jumps_row)])

            return (np.array(cols + cumulative_shift_col, dtype=int), 
                    np.array(rows + cumulative_shift_row, dtype=int) )
        
        else:
            return cols, rows


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
            self[frame].plot(ax=ax, scaling=scaling, ads_scaling=ads_scaling,
                    markers=markers, colors=colors, legend=legend, show_axis=show_axis)

            label = rf'$t_{{{frame}}} = $ {to_superscript_scientific(self.times[frame], precision=precision)} s'
            ax.text(0.03, 0.9, label, transform=ax.transAxes, ha='left', fontsize=10)

            return ax


        return FuncAnimation(fig, update, frames=frames, interval=interval, blit=False)


    def __len__(self):
        """
        Number of configurations in trajectory.

        Returns the number of loaded states
        """
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

        return self.create_states(idx)


    def __repr__(self):
        """String representation of Trajectory"""
        if len(self) > 0:
            t_range = f"t=[{self.times[0]:.2e}, {self.times[-1]:.2e}]"
        else:
            t_range = "empty"
        return f"Trajectory(nstates={len(self)}, {t_range}, lattice={self.lattice.n_sites} sites)"
