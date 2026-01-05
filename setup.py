from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='nat-zacros',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    author='dauerba, akandra',
    author_email='your-email@example.com',  # Update with actual email
    description='Analysis package for Zacros kinetic Monte Carlo simulations',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/dauerba/nat_zacros',
    packages=find_packages(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    python_requires='>=3.8',
    install_requires=[
        'numpy>=1.20',
        'scipy>=1.7',
    ],
    extras_require={
        'parallel': ['tqdm>=4.60'],
        'dev': ['pytest>=6.0', 'pytest-cov>=2.12'],
    },
)
