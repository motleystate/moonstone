from setuptools import setup, find_packages

setup(
    name="moonstone",
    version="0.0.1",
    description='Utilities for metagenomics data analysis using machine learning.',
    author='Kenzo-Hugo Hillion, Sean Kennedy',
    author_email='kehillio@pasteur.fr',
    install_requires=[
        'pandas==0.25.3',
        'scikit-learn==0.21.3',
        'matplotlib',
        'plotly',
        'statsmodels'
    ],
    packages=find_packages(),
    entry_points={'console_scripts': ['moonstone=moonstone.main:run']},
)
