from setuptools import setup, find_packages

setup(
    name="metautils",
    version="0.0.1",
    description='Utilities for metagenomics data analysis.',
    author='Kenzo-Hugo Hillion',
    author_email='kehillio@pasteur.fr',
    install_requires=[
        'numpy',
        'pandas',
        'plotly',
        'jupyter'
    ],
    packages=find_packages()
)
