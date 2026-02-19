# Parallel Loading Performance Analysis

## Problem Summary

Performance measurements on Windows (dja X1 Carbon laptop, 14 cores) showed that **parallel trajectory loading is significantly slower than sequential loading** when using cached data:

| Load Source | Workers | Simulations | Time (s) | Speedup vs Sequential |
|:------------|:--------|:-----------|--------:|:----------------------|
| .pkl files  | None    | 19         | 3.1     | 1.0x                  |
| .pkl files  | 1       | 19         | 72.8    | 0.04x (23x slower)    |
| .pkl files  | 10      | 19         | 156.6   | 0.02x (50x slower)    |
| .pkl files  | 14 (default) | 19    | 160.9   | 0.02x (52x slower)    |

Loading from raw `history_output.txt` sequentially took 60.8s with `workers=None`.

## Root Cause: Unnecessary Pickling

The performance regression is caused by **pickling and unpickling overhead in the parallel code path**:

### Current Implementation (Problematic)

In [simulation.py](nat_zacros/simulation.py):

```python
# Line 156-185: _load_trajectories_parallel()
with ProcessPoolExecutor(max_workers=workers) as executor:
    trajs = list(executor.map(self._load_single_trajectory, self.traj_dirs, 
                              [cache]*len(self.traj_dirs), [verbose]*len(self.traj_dirs)))
```

Since `_load_single_trajectory` is a **bound method** (has `self`), `ProcessPoolExecutor` must:

1. **Pickle the entire Simulation object** (including `self.lattice` with all coordinate data) to send to each worker
2. Unpickle it in the worker process
3. Execute the method (which may just unpickle a cached .pkl file—not CPU intensive)
4. Pickle the returned Trajectory object (with large state arrays) back to main process
5. Unpickle results in main process

For cached loads, steps 1-2 and 4-5 are **completely wasted**—the function only deserializes pre-computed data from disk.

### Why Sequential Path is Faster

In [simulation.py](nat_zacros/simulation.py), lines 218-221:

```python
else:
    # Sequential loading
    self.trajectories = []
    for traj_dir in self.traj_dirs:
        self.trajectories.append(self._load_single_trajectory(traj_dir, cache=cache, verbose=verbose))
```

Here, `_load_single_trajectory` is called directly in the same process:
- ✅ `self` is already in memory (just a reference, no serialization)
- ✅ Only unpickles the actual cached data from .pkl files (unavoidable cost)
- ✅ No inter-process communication overhead

## Code Path Implications

### When is `cache=True`?
- The Trajectory data is unpickled from a pickle file on disk
- The Lattice object (`self.lattice`) is already in memory and doesn't need to be sent anywhere
- Pickling the entire Simulation object to worker processes is **wasteful**

### When is `cache=False`?
- Trajectories are parsed from `history_output.txt` (CPU-intensive)
- The Lattice object must be available in the worker to parse coordinates
- Even here, pickling the entire Simulation object is inefficient (only the Lattice is needed)

## Recommendations

### 1. **Disable Parallel Loading for Cached Data (Immediate Fix)**

Modify [simulation.py](nat_zacros/simulation.py) `load()` method:

```python
def load(self, cache=True, workers=mp.cpu_count(), verbose=False):
    """Load trajectory data with caching support..."""
    
    # Use sequential loading when cache is enabled (pickling overhead > compute savings)
    if cache:
        workers = None
    
    # Rest of method...
```

**Expected result:** 52x speedup for cached loads (160.9s → 3.1s)

### 2. **Optimize Parallel Loading (Long-term Fix)**

For non-cached loads, create a module-level worker function to avoid pickling the entire Simulation object:

```python
# Module-level function (not a method)
def _load_trajectory_worker(args):
    """Load a single trajectory in worker process.
    
    Only receives the data needed (traj_dir, lattice data), not the entire Simulation object.
    """
    traj_dir, lattice_coords, lattice_cell_vectors, cache, verbose = args
    
    # Recreate minimal Lattice object in worker
    from .lattice import Lattice
    lattice = Lattice.__new__(Lattice)
    lattice.coordinates = lattice_coords
    lattice.cell_vectors = lattice_cell_vectors
    
    traj = Trajectory(traj_dir, lattice)
    traj.load(cache=cache, verbose=verbose)
    return traj
```

Then in `_load_trajectories_parallel()`:

```python
def _load_trajectories_parallel(self, cache=True, workers=None, verbose=False):
    if len(self.traj_dirs) == 0:
        return []
    
    with ProcessPoolExecutor(max_workers=workers) as executor:
        args_list = [
            (traj_dir, self.lattice.coordinates, self.lattice.cell_vectors, 
             cache, verbose)
            for traj_dir in self.traj_dirs
        ]
        trajs = list(executor.map(_load_trajectory_worker, args_list))
    
    return trajs
```

This reduces pickling overhead while keeping parallel benefits for expensive raw-data parsing.

### 3. **Consider `threading` Instead of `ProcessPoolExecutor`**

For I/O-bound cached loads, Python's `ThreadPoolExecutor` (doesn't require pickling) might be faster:

```python
from concurrent.futures import ThreadPoolExecutor

with ThreadPoolExecutor(max_workers=workers) as executor:
    trajs = list(executor.map(self._load_single_trajectory, self.traj_dirs, ...))
```

**Caveat:** Only benefits if GIL is not contended. Profile before adopting.

## Testing & Validation

To validate improvements, benchmark with this test code:

```python
import time

# Clear cache first
simset.clear_traj_cache()

# Test sequential hidden
start = time.time()
simset.load(cache=True, workers=None)
seq_time = time.time() - start
print(f"Sequential (cache=True): {seq_time:.2f}s")

# Test parallel (before fix)
simset.simulations = []
start = time.time()
simset.load(cache=True, workers=14)
par_time = time.time() - start
print(f"Parallel with 14 workers (cache=True): {par_time:.2f}s")
print(f"Overhead ratio: {par_time / seq_time:.1f}x")
```

## References

- **Python ProcessPoolExecutor pickling:** arguments and return values are serialized using `pickle`
- **Windows process creation:** Significantly slower than Linux `fork()`, exacerbating overhead
- **GIL:** Not a factor for ProcessPoolExecutor (separate Python interpreters per process)
