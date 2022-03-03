from setuptools import setup, find_packages

exec(open('moonstone/version.py').read())

setup(
    name="moonstone",
    version=__version__,  # noqa
    description='Utilities for metagenomics data analysis using machine learning.',
    author='Kenzo-Hugo Hillion, Agnès Baud, Mariela Furstenheim, Sean Kennedy',
    author_email='kehillio@pasteur.fr',
    install_requires=[
        'pandas==1.0.1',
        'matplotlib==3.3.0',
        'plotly==5.6.0',
        'statsmodels==0.11.1',
        'python-slugify==4.0.1',
        'pyaml==20.4.0',
        'numpy==1.18.1',
        'scikit-bio==0.5.6',
        'scikit-learn==0.21.3',
        'hdmedians==0.13',
        'cython==0.29.21',
        'scipy==1.5.2'
    ],
    packages=find_packages(),
    entry_points={'console_scripts': ['moonstone=moonstone.main:run']},
)
