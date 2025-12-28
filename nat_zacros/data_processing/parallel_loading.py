"""
Parallel trajectory loading functions for Zacros simulations.

This module provides parallel processing utilities for loading multiple trajectories efficiently.
"""

import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor

def _load_single_trajectory_equilibrated(args):
    """
    Helper function for parallel trajectory loading.
    Parameters
    ----------
    args : tuple
        (lattice, traj_dir, fraction, method, energy_only) for loading
    Returns
    -------
    trajectory
        Trajectory with equilibrated states loaded
    Notes
    -----
    Module-level function for multiprocessing pickling.
    """
    from ..trajectory import trajectory
    lattice, traj_dir, fraction, method, energy_only = args
    traj = trajectory(lattice, traj_dir)
    traj.load_trajectory(fraction=fraction, load_energy=True, energy_only=energy_only)
    return traj

def load_trajectories_parallel(lattice, traj_dirs, fraction=0.5, method='fraction', energy_only=False, n_workers=None):
    """
    Load multiple trajectories in parallel with equilibrated states.
    Parameters
    ----------
    lattice : lattice object
        Common lattice for all trajectories
    traj_dirs : list of Path
        List of trajectory directory paths
    fraction : float, default 0.5
        Fraction of trajectory to keep from the end (for equilibration)
    method : str, default 'fraction'
        Equilibration detection method
    energy_only : bool, default False
        If True, only load times and energies (much faster).
        If False, load full state configurations.
    n_workers : int, optional
        Number of parallel workers. If None, uses all available cores.
    Returns
    -------
    list of trajectory
        Loaded trajectories with equilibrated states
    """
    try:
        from tqdm import tqdm
        use_tqdm = True
    except ImportError:
        use_tqdm = False
    if len(traj_dirs) == 0:
        return []
    if n_workers is None:
        n_workers = mp.cpu_count()
    args_list = [(lattice, traj_dir, fraction, method, energy_only) for traj_dir in traj_dirs]
    print(f"Loading {len(traj_dirs)} trajectories in parallel using {n_workers} workers...")
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