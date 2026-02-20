DEVELOPMENT POLICY — nat_zacros

Status
------
This project is in active, early development. Public API and on-disk formats may change without notice.

Policy
------
- Backward compatibility is *not* guaranteed. Prioritize correctness, clarity, and rapid iteration over preserving legacy APIs.
- Remove deprecated APIs immediately where doing so simplifies the codebase and tests.
- When making breaking changes, update these three places at minimum:
  - `CHANGELOG.md` (document the change and migration notes)
  - Unit tests (add/modify tests that express the new contract)
  - Examples / notebooks (update usage examples and sample workflows)

Cache & Persistence
-------------------
- Canonical cache: per-trajectory `traj.pkl` files located inside each `traj_*` folder.
- The package does **not** create or rely on any aggregated `trajs_eq.pkl` cache; such files are considered external/user-managed artifacts and may be removed.
- Prefer storing small, stable serializable summaries (times/energies/arrays) instead of pickling full objects for long-term caches.

Developer workflow (recommended)
--------------------------------
1. Write or update unit tests that capture the intended API/behavior.
2. Update `CHANGELOG.md` and example notebooks to reflect the change.
3. Run the full test suite (`pytest -q`) and update docs/examples.
4. Open a PR with a clear description and migration notes.

Rationale
---------
- Early-stage projects benefit from fast iteration and a small, well-tested public surface.
- Avoid brittle serialized formats and duplicated cache layers.

If you want a stricter compatibility policy later, add a `STABILITY.md` and adopt semantic versioning and deprecation timelines.