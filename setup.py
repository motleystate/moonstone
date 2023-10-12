from setuptools import setup, find_packages

exec(open('moonstone/version.py').read())

setup(
    name="moonstone",
    version=__version__,  # noqa
    description='Utilities for metagenomics data analysis using machine learning.',
    author='Kenzo-Hugo Hillion, Agnès Baud, Mariela Furstenheim, Sean Kennedy',
    author_email='kehillio@pasteur.fr',
    install_requires=[
        'pandas==2.0.2',
        'matplotlib==3.3.0',
        'plotly==5.17.0',
        'statsmodels==0.13.0',
        'python-slugify==4.0.1',
        'pyaml==20.4.0',
        'numpy==1.24.3',
        'scikit-bio==0.5.9',
        'scikit-learn==1.3.1',
        'hdmedians==0.14.2',
        'cython==0.29.21',
        'scipy==1.9.0'
    ],
    packages=find_packages(),
    entry_points={'console_scripts': ['moonstone=moonstone.main:run']},
)
