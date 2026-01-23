from setuptools_scm import get_version

print("--- Version Debug Info ---")
try:
    version = get_version(root='.', relative_to=__file__)
    print(f"Detected Version: {version}")
except Exception as e:
    print(f"Error getting version: {e}")
