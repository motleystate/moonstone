from setuptools import setup, find_packages

setup(
    name="moonstone",
    version="0.0.1",
    description='Utilities for metagenomics data analysis using machine learning.',
    author='Kenzo-Hugo Hillion, Sean Kennedy',
    author_email='kehillio@pasteur.fr',
    install_requires=[
        'numpy',
        'pandas',
        'plotly',
        'jupyter'
    ],
    packages=find_packages()
)
