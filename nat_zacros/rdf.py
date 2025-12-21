"""
Parallel RDF computation and trajectory loading functions.

This module provides parallel processing utilities for RDF calculations and
trajectory loading in Zacros simulations.

⚠️  PERFORMANCE WARNING - READ BEFORE USING
--------------------------------------------
Parallel RDF functions are rarely beneficial for typical datasets.

With vectorized distance calculations, RDF computation is very fast (~2s).
Parallelization overhead (process spawn, pickling, IPC) is ~2-4 seconds.

RESULT: Sequential computation is FASTER for most use cases.

Use parallel RDF functions ONLY if:
  - You have >50 trajectories, OR
  - Sequential computation takes >30 seconds, AND
  - You have benchmarked and confirmed speedup

For typical use (10-20 trajectories):
  Sequential: ~2s
  Parallel: ~3-5s (SLOWER!)

HOWEVER: Parallel trajectory loading IS recommended and effective for all cases.
"""

import numpy as np
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
from pathlib import Path


# ==============================================================================
# Parallel RDF computation functions
# ==============================================================================

def _compute_single_rdf(args):
    """
    Helper function for parallel RDF computation (trajectory-level).
    
    Parameters
    ----------
    args : tuple
        (trajectory, r_max, dr, g_ref) tuple for computing RDF
        
    Returns
    -------
    g : ndarray
        RDF values for this trajectory
        
    Notes
    -----
    This is a module-level function to ensure it's pickle-able for multiprocessing.
    """
    traj, r_max, dr, g_ref = args
    r, g = traj.get_rdf(r_max=r_max, dr=dr, g_ref=g_ref)
    return g


def _compute_state_rdf(args):
    """
    Helper function for parallel RDF computation (state-level).
    
    Parameters
    ----------
    args : tuple
        (state, lattice, r_max, dr, g_ref, n_bins) tuple
        
    Returns
    -------
    tuple
        (counts, n_occupied, coverage) for this state
        
    Notes
    -----
    This processes a single state and returns histogram counts for aggregation.
    Used for fine-grained parallelism when you have many states.
    """
    state, lattice, r_max, dr, g_ref, n_bins = args
    
    occupied_coords = state.get_occupied_coords()
    n_occupied = len(occupied_coords)
    coverage = state.get_coverage()
    
    if n_occupied < 2:
        return np.zeros(n_bins), n_occupied, coverage
    
    # Vectorized distance calculation
    distances = lattice.pairwise_distances_pbc(occupied_coords)
    mask = np.triu(np.ones(distances.shape, dtype=bool), k=1)
    valid_dists = distances[mask]
    valid_dists = valid_dists[(valid_dists > 0) & (valid_dists <= r_max)]
    
    # Histogram
    bin_edges = np.linspace(0.0, r_max, n_bins + 1)
    counts, _ = np.histogram(valid_dists, bins=bin_edges)
    
    return counts, n_occupied, coverage


def compute_rdf_parallel(trajectories, r_max=None, dr=0.1, g_ref=None, n_workers=None):
    """
    Compute RDF averaged over multiple trajectories using parallel processing.
    
    ⚠️  WARNING: With vectorized distance calculations, RDF computation is very fast.
    Parallelization overhead (2-4s) often exceeds computation time for typical datasets.
    Only use this function when sequential computation time >> 20-30 seconds.
    
    **When to use:**
    - >50 trajectories
    - Very long trajectories (>1000 states each)
    - Large lattices (>10,000 sites)
    
    **For typical use (10-20 trajectories):** Sequential computation is FASTER.
    Use trajectory.get_rdf() in a simple loop instead.
    
    Parameters
    ----------
    trajectories : list of trajectory objects
        Trajectories to average over. Each trajectory should have states loaded.
    r_max : float, optional
        Maximum distance for RDF calculation (Angstroms). 
        If None, uses minimum cell dimension / 2.
    dr : float, optional
        Bin width for RDF histogram (default: 0.1 Angstrom)
    g_ref : ndarray, optional
        Reference RDF for normalization (from full lattice).
        If None, no normalization by coordination numbers.
    n_workers : int, optional
        Number of parallel workers. If None, uses all available cores.
        
    Returns
    -------
    r : ndarray
        Distance bin centers (Angstroms)
    g_avg : ndarray
        Average RDF over all trajectories
    g_std : ndarray
        Standard deviation of RDF across trajectories
        
    Examples
    --------
    >>> # For typical datasets, use sequential computation (FASTER):
    >>> rdfs = []
    >>> for traj in trajs:
    >>>     r, g = traj.get_rdf(r_max=40.0, dr=0.1, g_ref=g_ref)
    >>>     rdfs.append(g)
    >>> g_avg = np.mean(rdfs, axis=0)
    >>> 
    >>> # Only for very large datasets (>50 trajectories):
    >>> r, g_avg, g_std = compute_rdf_parallel(trajs, r_max=40.0, dr=0.1, 
    ...                                         g_ref=g_ref, n_workers=4)
    
    Notes
    -----
    - Uses ProcessPoolExecutor for true parallel computation
    - Parallelization overhead: ~2-4 seconds (process spawn, pickling, IPC)
    - With vectorization, RDF for 10 trajectories takes ~2s sequentially
    - Therefore parallel version is SLOWER for typical use cases
    - Benchmark your specific case before using this function
    """
    # Try to import tqdm for progress bar
    try:
        from tqdm import tqdm
        use_tqdm = True
    except ImportError:
        use_tqdm = False
    
    if len(trajectories) == 0:
        raise ValueError("No trajectories provided")
    
    # Determine number of workers
    if n_workers is None:
        n_workers = mp.cpu_count()
    
    # Prepare arguments for each trajectory
    args_list = [(traj, r_max, dr, g_ref) for traj in trajectories]
    
    # Parallel computation
    print(f"Computing RDF for {len(trajectories)} trajectories using {n_workers} workers...")
    
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        if use_tqdm:
            # With progress bar
            rdfs = list(tqdm(executor.map(_compute_single_rdf, args_list),
                            total=len(trajectories),
                            desc="Computing RDF",
                            unit="traj"))
        else:
            # Without progress bar (fallback)
            rdfs = list(executor.map(_compute_single_rdf, args_list))
    
    # Get distance axis from first trajectory
    r, _ = trajectories[0].get_rdf(r_max=r_max, dr=dr, g_ref=g_ref)
    
    # Compute statistics across trajectories
    rdfs = np.array(rdfs)
    g_avg = np.mean(rdfs, axis=0)
    g_std = np.std(rdfs, axis=0)
    
    print(f"\nSuccessfully computed RDF averaged over {len(trajectories)} trajectories")
    
    return r, g_avg, g_std


def compute_rdf_parallel_states(trajectories, r_max=None, dr=0.1, g_ref=None, n_workers=None):
    """
    Compute RDF with state-level parallelism (better for few trajectories with many states).
    
    ⚠️  WARNING: With vectorized distance calculations, RDF computation is very fast.
    Parallelization overhead (2-4s per trajectory) often exceeds computation time.
    Only use this function when sequential computation time >> 30 seconds.
    
    **When to use:**
    - >100 trajectories with many states each
    - Extremely large systems where per-state computation is slow
    - Already benchmarked and confirmed benefit
    
    **For typical use (10-20 trajectories, ~100 states each):** Sequential is FASTER.
    The overhead from spawning processes for each state dominates.
    
    Parameters
    ----------
    trajectories : list of trajectory objects
        Trajectories to average over. Each trajectory should have states loaded.
    r_max : float, optional
        Maximum distance for RDF calculation (Angstroms). 
    dr : float, optional
        Bin width for RDF histogram (default: 0.1 Angstrom)
    g_ref : ndarray, optional
        Reference RDF for normalization.
    n_workers : int, optional
        Number of parallel workers. If None, uses all available cores.
        
    Returns
    -------
    r : ndarray
        Distance bin centers (Angstroms)
    g_avg : ndarray
        Average RDF over all trajectories
    g_std : ndarray
        Standard deviation of RDF across trajectories
        
    Notes
    -----
    - Parallelizes over individual states instead of trajectories
    - Higher parallelization overhead than compute_rdf_parallel()
    - With vectorization, sequential computation is very fast (~2s for 10x100 states)
    - Parallel overhead exceeds benefit for typical datasets
    - Always benchmark before using in production code
    """
    try:
        from tqdm import tqdm
        use_tqdm = True
    except ImportError:
        use_tqdm = False
    
    if len(trajectories) == 0:
        raise ValueError("No trajectories provided")
    
    # Determine r_max if not provided
    if r_max is None:
        traj0 = trajectories[0]
        v1 = traj0.lattice.cell_vectors[0]
        v2 = traj0.lattice.cell_vectors[1]
        l1 = np.linalg.norm(v1)
        l2 = np.linalg.norm(v2)
        l3 = np.linalg.norm(v1 + v2)
        r_max = min(l1, l2, l3) / 2.0
    
    # Initialize histogram
    n_bins = int(np.ceil(r_max / dr))
    bin_edges = np.linspace(0.0, r_max, n_bins + 1)
    r_bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    
    # Determine number of workers
    if n_workers is None:
        n_workers = mp.cpu_count()
    
    # Process each trajectory
    rdfs_per_traj = []
    
    for traj_idx, traj in enumerate(trajectories):
        # Prepare arguments for each state in this trajectory
        args_list = [(state, traj.lattice, r_max, dr, g_ref, n_bins) 
                    for state in traj.states]
        
        total_states = len(args_list)
        if total_states == 0:
            continue
        
        print(f"Processing trajectory {traj_idx+1}/{len(trajectories)} ({total_states} states) with {n_workers} workers...")
        
        # Parallel computation across states
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            if use_tqdm:
                results = list(tqdm(executor.map(_compute_state_rdf, args_list, chunksize=max(1, total_states // (n_workers * 4))),
                                total=total_states,
                                desc=f"Traj {traj_idx+1}",
                                unit="state"))
            else:
                results = list(executor.map(_compute_state_rdf, args_list, chunksize=max(1, total_states // (n_workers * 4))))
        
        # Aggregate results for this trajectory
        g_r = np.zeros(n_bins)
        avg_coverage = np.mean([r[2] for r in results])
        
        for counts, n_occupied, coverage in results:
            if n_occupied < 2:
                continue
            
            # Normalize by g_ref if provided
            if g_ref is not None:
                counts_n = np.zeros_like(counts, dtype=float)
                np.divide(counts, g_ref, out=counts_n, where=g_ref!=0)
                g_r += counts_n / n_occupied / avg_coverage
            else:
                g_r += counts / n_occupied / avg_coverage
        
        # Normalize by number of states (factor of 2 for unordered pairs)
        if len(results) > 0:
            g_r = 2 * g_r / len(results)
        
        rdfs_per_traj.append(g_r)
    
    # Compute statistics across trajectories
    rdfs = np.array(rdfs_per_traj)
    g_avg = np.mean(rdfs, axis=0)
    g_std = np.std(rdfs, axis=0)
    
    print(f"\nSuccessfully computed RDF averaged over {len(trajectories)} trajectories")
    print(f"  Total states processed: {sum(len(t.states) for t in trajectories)}")
    
    return r_bins, g_avg, g_std


# ==============================================================================
# Parallel trajectory loading functions
# ==============================================================================

def _load_single_trajectory_equilibrated(args):
    """
    Helper function for parallel trajectory loading.
    
    Parameters
    ----------
    args : tuple
        (lattice, traj_dir, fraction, method) for loading
        
    Returns
    -------
    trajectory
        Trajectory with equilibrated states loaded
        
    Notes
    -----
    Module-level function for multiprocessing pickling.
    """
    from .trajectory import trajectory
    
    lattice, traj_dir, fraction, method = args
    
    # Create trajectory and do two-phase loading
    traj = trajectory(lattice, traj_dir)
    traj.load_trajectory(load_energy=True, energy_only=True)
    traj.load_equilibrated_states(fraction=fraction, method=method)
    return traj
    

def load_trajectories_parallel(lattice, traj_dirs, fraction=0.5, method='fraction', n_workers=None):
    """
    Load multiple trajectories in parallel with equilibrated states.
    
    ✅ RECOMMENDED: This optimization is effective for all use cases.
    I/O operations (file reading/parsing) are slow enough that parallel loading
    provides consistent speedup (5-10x) without excessive overhead.
    
    Unlike parallel RDF computation, loading is I/O-bound and benefits from
    parallelization even for small numbers of trajectories.
    
    Parameters
    ----------
    lattice : lattice object
        Common lattice for all trajectories
    traj_dirs : list of Path
        List of trajectory directory paths
    fraction : float, default 0.5
        Fraction of trajectory to skip for equilibration
    method : str, default 'fraction'
        Equilibration detection method
    n_workers : int, optional
        Number of parallel workers. If None, uses all available cores.
        
    Returns
    -------
    list of trajectory
        Loaded trajectories with equilibrated states
        
    Examples
    --------
    >>> # Sequential loading (slow - 60s for 10 trajectories)
    >>> trajs = []
    >>> for traj_dir in traj_dirs:
    >>>     traj = trajectory(lat, traj_dir)
    >>>     traj.load_trajectory(energy_only=True)
    >>>     traj.load_equilibrated_states(fraction=0.5)
    >>>     trajs.append(traj)
    >>> 
    >>> # Parallel loading (fast - 6-10s for 10 trajectories)
    >>> trajs = load_trajectories_parallel(lat, traj_dirs, fraction=0.5)
    
    Notes
    -----
    - Typical speedup: 5-10x with 10+ cores
    - Effective for all dataset sizes (I/O-bound operations benefit from parallelism)
    - Each trajectory reads its own file → true parallel I/O
    - Overhead is minimal compared to file parsing time
    - Combine with pickle caching for even better repeated-analysis performance
    """
    try:
        from tqdm import tqdm
        use_tqdm = True
    except ImportError:
        use_tqdm = False
    
    if len(traj_dirs) == 0:
        return []
    
    # Determine number of workers
    if n_workers is None:
        n_workers = mp.cpu_count()
    
    # Prepare arguments
    args_list = [(lattice, traj_dir, fraction, method) for traj_dir in traj_dirs]
    
    print(f"Loading {len(traj_dirs)} trajectories in parallel using {n_workers} workers...")
    
    # Parallel loading
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        if use_tqdm:
            trajs = list(tqdm(executor.map(_load_single_trajectory_equilibrated, args_list),
                            total=len(traj_dirs),
                            desc="Loading trajectories",
                            unit="traj"))
        else:
            trajs = list(executor.map(_load_single_trajectory_equilibrated, args_list))
    
    print(f"Successfully loaded {len(trajs)} trajectories")
    if len(trajs) > 0:
        if fraction < 1.0:
            print(f"  Example: {len(trajs[0])} states per trajectory (equilibrated)")
        else:
            print(f"  Example: {len(trajs[0])} states per trajectory")
    
    return trajs
