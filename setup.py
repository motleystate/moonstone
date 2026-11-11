from setuptools import setup, find_packages

exec(open('moonstone/version.py').read())

setup(
    name="moonstone",
    version=__version__,  # noqa
    description='Utilities for metagenomics data analysis using machine learning.',
    author='Kenzo-Hugo Hillion, Agnès Baud, Mariela Furstenheim, Sean Kennedy',
    author_email='kehillio@pasteur.fr',
    install_requires=[
        'pandas>=3.0.0',
        'matplotlib>=3.3.0',
        'plotly>=5.17.0',
        'statsmodels>=0.14.0',
        'python-slugify>=4.0.1',
        'pyaml>=20.4.0',
        'numpy>=1.26.4',
        'scikit-bio>=0.6.0',
        'scikit-learn>=1.5.0',
        'hdmedians>=0.14.2',
        'cython>=0.29.21',
        'scipy>=1.9.0',
        'filetype>=1.2.0',
        'joblib>=1.3.0',
        'threadpoolctl>=3.1.0'
    ],
    python_requires='>=3.13,<3.14',
    packages=find_packages(),
    entry_points={'console_scripts': ['moonstone=moonstone.main:run']},
)
