# Migration to `pyproject.toml` and Versioning Overhaul

This document details the changes made to the `nat-zacros` package structure to modernize the build system and improve version control integration.

## 1. Summary of Changes

### A. Build System Migration
- **Removed `setup.py`**: The legacy setup script has been replaced.
- **Added `pyproject.toml`**: This file now governs the build configuration, dependencies, and metadata (following PEP 621).
- **Editable Installs**: We now use `pip install -e .` which is fully supported by `pyproject.toml`.

### B. Versioning (`setuptools_scm`)
- **Automated Versioning**: The package version is no longer hardcoded in `__init__.py`. It is derived automatically from Git tags.
- **Scheme**: We use the standard `guess-next-dev` scheme.
    - **Clean state**: `v0.0.4` -> `0.0.4`
    - **Development state**: `v0.0.4` + 1 commit -> `0.0.5.dev1` (indicates development towards 0.0.5).
    - **Dirty state**: Adds a suffix like `+d20260121` so you know if your local code has uncommitted changes.
- **Dynamic Runtime Version**: The `__init__.py` has been updated to fetch the version dynamically from `setuptools_scm` when running in editable mode. This means local version numbers update immediately upon commit without needing to reinstall.

### C. Code Fixes
- **`simulation.py`**: Replaced system `tar` command with Python's `tarfile` module. This fixes issues on Windows where `tar.exe` might not be in the PATH.
- **`simulation_set.py`**:
    - Fixed `AttributeError` by initializing `self.simulations = []`.
    - `plot()` now automatically calls `load_energy()` if data hasn't been loaded yet.

---

## 2. Instructions for Colleague Migration

To update your local environment to work with this new structure, follow these steps:

### Step 1: Update Repository
Pull the latest changes (or checkout the branch):
```bash
git fetch origin
git checkout migration-toml  # (or main, if already merged)
```

### Step 2: Re-install Package
You do **not** need to delete your Conda environment, but you must reinstall the package to register the new metadata system.

1.  Activate your Conda environment:
    ```bash
    conda activate natzacros
    ```

2.  Uninstall the old version (recommended to clear old egg-info):
    ```bash
    pip uninstall -y nat-zacros
    ```

### Step 3: Install
    ```bash
    pip install -e .
    ```
    *Note: Ensure you are in the root directory that contains `pyproject.toml`.*

### Step 3: Verify
Run this command in Python to ensure the version is being picked up correctly:
```bash
python -c "import nat_zacros; print(nat_zacros.__version__)"
```

---

## 3. Instructions to Merge & Release v0.0.5

When you are ready to merge the `migration-toml` branch into `main` and create the official `v0.0.5` release, follow these steps.

### Step 1: Merge into Main
```bash
# Switch to the main branch
git checkout main

# Merge the changes
git merge migration-toml
```

### Step 2: Create the Release Tag
Now that `main` has the new code, tag it. `setuptools_scm` will look at this tag to generate the version number `0.0.5` cleanly.

```bash
git tag v0.0.5
```

### Step 3: Push to Remote (GitHub)
Push both the code changes and the new tag.

```bash
git push origin main
git push origin v0.0.5
```

### Step 4: Final Verification
After pushing, anyone who installs the package from this commit (or pulls main) will see the version `0.0.5`.

```bash
python -c "import nat_zacros; print(nat_zacros.__version__)"
# Output should be: 0.0.5
```
