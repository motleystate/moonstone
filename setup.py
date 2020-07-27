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
        'scikit-learn==0.21.3',
        'matplotlib',
        'plotly',
        'statsmodels',
        'python-slugify',
        'pyaml',
        'numpy==1.18.1'
    ],
    packages=find_packages(),
    entry_points={'console_scripts': ['moonstone=moonstone.main:run']},
)
