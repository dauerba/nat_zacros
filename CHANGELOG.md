# Changelog

All notable changes to this project will be documented in this file.

## Unreleased

### Added
- `SimulationSet.clear_traj_cache()` — removes internal per-trajectory cache files (`traj.pkl`) and aggregated `<sim>_trajs_eq.pkl`.
- `SimulationSet.clear_results()` — removes user-visible results (`energy`, `rdf`, `accessibility`, `gref`).
- Unit tests for `clear_traj_cache()` and `clear_results()` (pytest).
- Updated example notebooks and documentation to use the new API.

### Changed
- `SimulationSet.clear_cache()` is now deprecated and delegates to the new functions; emits DeprecationWarning.

## [0.0.0] - unreleased
- Initial development (see repository history)
