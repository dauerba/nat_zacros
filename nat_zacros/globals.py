"""
Centralized constants and configuration for nat_zacros.
"""

# Dictionary of cache file (suffix, extension) tuples
# Used by Simulation and SimulationSet classes
nz_cache_files = {
    'trajs': ('_trajs', 'pkl'),
    'energy': ('_energy', 'dat'),
    'rdf': ('_rdf', 'dat'),
    'accessibility': ('_accessibility', 'dat'),
}
