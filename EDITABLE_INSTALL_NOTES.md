# Editable Install Issues with Modern Python

## Question to Claude Sonnet 4.5

I was trying to use pip editable install to make an environment that included nat_zacros package. This was done in a  terminal by navigating to the root of the GIT repository and using "pip install -e ." . I had problems with the version information, derived from the GIT tag and current comitt

Are ther known problems with using editable install with modern python code?

## Claude Sonnet 4.5 Reply: 

### Known Problems with `pip install -e .`

#### 1. setuptools_scm Version Detection Failures
The package installed in "editable" mode may not properly detect git information:
- Version often shows as `0.0.0` or throws errors
- This happens because setuptools_scm looks for `.git` directory relative to the installed package location
- In editable mode, the package location can confuse the git detection

#### 2. PEP 660 (New Editable Install Standard)
Modern setuptools (≥64.0.0) uses a new editable install mechanism:
- Can break dynamic version detection and file discovery
- The old behavior used symlinks; the new one uses import hooks
- This can interfere with setuptools_scm's ability to find git metadata

### Solutions

#### Option 1: Use Legacy Editable Install
Forces the old symlink-based editable install behavior:
```bash
pip install -e . --config-settings editable_mode=compat
```

#### Option 2: Set SETUPTOOLS_SCM_PRETEND_VERSION
Override version detection with a fixed version:
```bash
# Windows (cmd)
set SETUPTOOLS_SCM_PRETEND_VERSION=0.1.0
pip install -e .

# Windows (PowerShell)
$env:SETUPTOOLS_SCM_PRETEND_VERSION="0.1.0"
pip install -e .

# Linux/Mac
export SETUPTOOLS_SCM_PRETEND_VERSION=0.1.0
pip install -e .
```

#### Option 3: Write Version to File
Update `pyproject.toml` to write version to a file during build:
```toml
[tool.setuptools_scm]
write_to = "nat_zacros/_version.py"
fallback_version = "0.0.0"
```

Then modify `nat_zacros/__init__.py` to read from this file:
```python
try:
    from ._version import version as __version__
except ImportError:
    from setuptools_scm import get_version
    __version__ = get_version(root='..', relative_to=__file__)
```

#### Option 4: Reinstall After Each Commit
If version needs to update with each commit:
```bash
pip install -e . --no-build-isolation
```

#### Option 5: Use pyproject.toml Only (Modern Approach)
Move all configuration from setup.py to pyproject.toml (recommended for new projects).

---

### Analysis of Your Current Configuration

#### Current Files

**pyproject.toml:**
```toml
[build-system]
requires = ["setuptools>=42", "setuptools_scm>=6.2"]
build-backend = "setuptools.build_meta"
```

**setup.py:**
```python
setup(
    name='nat-zacros',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    # ... other settings ...
)
```

#### Issues Found

1. **Missing setuptools_scm configuration**
   - No `[tool.setuptools_scm]` section in `pyproject.toml`
   - No `write_to` directive to persist version information
   - No fallback version specified

2. **Hybrid setup (setup.py + pyproject.toml)**
   - Configuration split between two files
   - `setup_requires` in setup.py is deprecated (should use build-system.requires in pyproject.toml)
   - This is valid but can cause confusion with modern pip versions

3. **No version fallback mechanism**
   - If git tags are not available or setuptools_scm fails, there's no fallback
   - Can cause immediate failure in editable installs

#### Recommended Fix for Your Setup

Add this to your `pyproject.toml`:

```toml
[build-system]
requires = ["setuptools>=64", "setuptools_scm>=8"]
build-backend = "setuptools.build_meta"

[tool.setuptools_scm]
write_to = "nat_zacros/_version.py"
fallback_version = "0.0.0+unknown"
version_scheme = "post-release"
local_scheme = "no-local-version"
```

Then update `nat_zacros/__init__.py` to include:

```python
try:
    from ._version import __version__
except ImportError:
    try:
        from setuptools_scm import get_version
        __version__ = get_version(root='..', relative_to=__file__)
    except (ImportError, LookupError):
        __version__ = "unknown"
```

#### Why This Helps

1. **`write_to`** - Creates a `_version.py` file that persists version info
2. **`fallback_version`** - Provides a default when git isn't available
3. **Try/except in __init__.py** - Multiple fallbacks for version detection
4. **Modern setuptools** - Better handling of editable installs
5. **`local_scheme = "no-local-version"`** - Cleaner version strings without local git hash

#### Quick Test After Installing

After running `pip install -e .`, test in Python:

```python
import nat_zacros
print(nat_zacros.__version__)  # Should show version, not error

# Also check if _version.py was created:
import nat_zacros
print(nat_zacros.__file__)  # Find install location
# Then check if _version.py exists in that directory
```

#### Additional Notes

- If you're using git tags for versioning, ensure you have at least one tag: `git tag v0.1.0`
- For editable installs with version updates, you may need to reinstall after tagging
- Consider migrating fully to pyproject.toml for cleaner configuration (Python 3.8+ standard)

---

### References

- [setuptools_scm documentation](https://github.com/pypa/setuptools_scm)
- [PEP 660 - Editable installs](https://peps.python.org/pep-0660/)
- [setuptools version documentation](https://setuptools.pypa.io/en/latest/userguide/distribution.html)
