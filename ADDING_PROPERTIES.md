# Adding New Properties to nat-zacros

This guide describes the standardized architecture for adding new ensemble-averaged properties (like RDF, Energy, or Accessibility) to the `nat-zacros` package.

## Overview: The "Manager-Worker" Pattern

The package uses a delegated architecture:
- **`Simulation` (Worker)**: Handles heavy computation for a single simulation folder. It knows how to calculate the property and save/load it to/from a specific file.
- **`SimulationSet` (Manager)**: Orchestrates calculations across multiple simulations, manages the `results/` directory structure, handles parameter-aware caching, and identifies which simulations need recalculation.

---

## Step 1: Register the Result File

In `nat_zacros/simulation_set.py`, update the `__init__` method of `SimulationSet` to include a default filename for the new property in the `_results_files` dictionary.

```python
self._results_files = {
    'energy': en_file,
    'rdf': rdf_file,
    'accessibility': acc_file,
    'new_property': 'new_property.dat', # Add your property here
}
```

---

## Step 2: Implement Calculation in `Simulation`

In `nat_zacros/simulation.py`, implement the worker method. This method must:
1. Take the necessary analysis parameters.
2. Accept an optional `file` (Path) argument.
3. Handle data saving if `file` is provided, including the **Standardized Caching Header**.

### The Caching Header
The first line of the saved file **MUST** start with `# Parameters:` followed by `key=value` pairs. This is used by `SimulationSet` for cache validation.

```python
def get_ensemble_new_property(self, param1=1.0, fraction=1.0, file=None):
    # ... calculation logic ...
    data = ... 
    
    if file is not None:
        header = f"Parameters: param1={param1} fraction={fraction}\nColumn1 Column2"
        np.savetxt(file, data, header=header)
    
    return data
```

---

## Step 3: Implement Wrapper in `SimulationSet`

In `nat_zacros/simulation_set.py`, implement the manager method using the `_get_ensemble_property_generic` helper. This helper automatically handles:
- Locating the correct file in the `results/` folder.
- Validating the cache parameters against the header.
- Handling the equilibration fraction (`fraction_eq`).

```python
def get_ensemble_new_properties(self, param1=1.0, verbose=False):
    def compute(sim, cache_file, params):
        return sim.get_ensemble_new_property(
            param1=params['param1'], 
            fraction=params['fraction'], 
            file=cache_file
        )

    return self._get_ensemble_property_generic(
        'new_property',      # Key from _results_files
        compute,             # Calculation callback
        {'param1': param1},  # Params to check in cache header
        use_fraction=True,   # Whether to include ensemble fraction_eq
        verbose=verbose
    )
```

---

## Step 4: Add Status Checking

Update `SimulationSet.check_cache_status` to include your new property in the status report.

```python
# In SimulationSet.check_cache_status
results_map = {
    'rdf': 'rdf.dat', 
    'energy': 'energy.dat', 
    'new_prop': 'new_property.dat' # Add here
}
```

---

## Step 5: (Optional) Implement Plotting

Add a visualization method to `SimulationSet` (e.g., `plot_new_property`) following the existing subplot grid patterns used in `plot_rdf` or `plot_energy`.
