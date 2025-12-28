# nat-zacros Example Notebooks

This folder contains example Jupyter notebooks demonstrating how to use the nat-zacros package for Zacros KMC analysis.

## Notebooks

- **demo_simulation_run_class.ipynb**: Demonstrates the recommended workflow using the `SimulationRun` class for loading, caching, and analyzing Zacros trajectories.

## Usage Instructions

1. **Path Setup**: The notebooks use a `user_paths` dictionary to set the path to the nat-zacros package. You must add your username and path if not already present:

    ```python
    user_paths = {
        'your_username': Path.home() / 'path' / 'to' / 'nat_zacros',
        # Add more users as needed
    }
    ```
    Replace `'your_username'` and the path as appropriate for your system.

2. **Data Directory**: The example expects Zacros output data in a `zacros_calculations` folder. Adjust the path if your data is elsewhere.

3. **Environment**: Make sure you have activated the correct Python environment and installed all dependencies (see the main README for details).

4. **Running the Notebook**: Open the notebook in Jupyter or VS Code and run the cells in order. If you encounter a `FileNotFoundError`, check that your paths are set correctly and that the required data files exist.

---

For more information, see the main nat-zacros documentation and README.
