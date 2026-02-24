# nat_zacros Architecture Overview

This document outlines the architectural patterns and internal data management logic of the `nat_zacros` package, designed for efficient analysis of large ensembles of Zacros kinetic Monte Carlo (KMC) simulations.

## 1. Core Hierarchy & Data Manager Pattern

The library is structured around three main components that separate geometric logic from ensemble management:

1. **Lattice**: The geometric foundation. Handles site coordinates, periodic boundary conditions (PBC), and neighbor lists.
2. **SimulationSet (The Manager)**: The central orchestrator. It manages the lifecycle of multiple simulations (e.g., across a range of temperatures/coverages).
3. **Simulation (The Worker)**: Handles data for a specific directory. It performs the heavy lifting for calculating properties (RDF, Energy, etc.) across multiple trajectory files.

### 1.1 Memory Optimization: Shared Geometry
A key architectural feature is **Shared Lattice Initialization**. When a `SimulationSet` is created, it initializes a single `Lattice` object. This object is then shared across all `Simulation` instances in the set. This prevents the memory leak and overhead of creating thousands of identical geometric models for a large ensemble.

### 1.2 ID-Keyed State Tracking
The `SimulationSet` no longer stores simulations as a list. Instead, `simset.simulations` is a dictionary keyed by the simulation ID (e.g., `1, 2, 40`). This allows for $O(1)$ lookups and prevents duplicate entries when loading data incrementally.

## 2. The "RAM World" vs. "Disk World"

To handle datasets that may exceed available RAM, `nat_zacros` operates on a strict separation between loaded objects and filesystem caches:

*   **Disk World (Persistent)**: 
    *   **Trajectory Caches**: `.pkl` files containing parsed KMC states (speeds up loading by 50-100x).
    *   **Result Files**: `.dat` files containing calculated properties (RDFs, Energies) with standardized `# Parameters:` headers for cache validation.
*   **RAM World (Transient)**: 
    *   **Simulation Objects**: Stored in `SimulationSet.simulations` (a dictionary keyed by simulation ID).
    *   **Trajectory Data**: The raw states stored inside `Simulation` objects.

### 2.1 Lazy Loading (Autoload)
The unified `get()` API implements defensive loading. If you request a property for a simulation that isn't currently in RAM, the system will:
1.  Identify the missing simulation ID.
2.  Trigger `load(simulations=[id], cache=True)` to fetch it from disk.
3.  Perform the calculation.
4.  Optionally `unload()` to free memory.

### 2.2 Explicit RAM Management
Users can explicitly prune the "RAM World" using:
- `simset.unload(simulations=[...])`: Frees RAM while keeping the "Disk World" results intact.
- `simset.check_cache_status()`: Provides a tabular dashboard showing which simulations are in RAM vs. which have persistent results on disk.

## 3. Unified Property API (`get`)

Instead of multiple disparate methods, all analysis is routed through a single dispatch system:
`results = simset.get(property_key, pars_dict={...})`

### 3.1 Parameter-Aware Caching
The system uses a "Manager-Worker" validation pattern:
1.  **Manager** checks the disk for a result file (e.g., `rdf.dat`).
2.  It reads the metadata header: `# Parameters: n_bins=100 ...`.
3.  If the requested parameters match the cached parameters, it loads the data immediately.
4.  If they mismatch (or the file is missing), it instructs the **Worker** (Simulation) to re-calculate and overwrite.

## 4. Current Limits and Future Directions

### 4.1 In-Memory Requirement (Current)
Currently, a trajectory must be fully loaded into RAM to be analyzed. While the `SimulationSet` can manage which *simulations* are loaded, individual *trajectories* cannot yet be streamed.

### 4.2 Block-wise Streaming (Future)
For datasets where a single trajectory exceeds 10GB, we are planning a **Block-wise Streaming Architecture**. This will transition `Trajectory.states` from a list to an iterator/generator, allowing analysis to happen in chunks without ever holding the full sequence in memory.

---
*For questions regarding specific implementation details, refer to `nat_zacros/simulation_set.py`.*
