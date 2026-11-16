"""Example input files, installed with the package. All values are simulated (see make_examples.py)."""
from importlib import resources

NAMES = ('species_counts.tsv', 'resistance_genes.csv', 'metadata.csv', 'taxa_counts.csv', 'kraken2_counts.tsv',
         'moonbeam_bracken.tsv')


def path(name: str) -> str:
    """Absolute path of an example file."""
    if name not in NAMES:
        raise FileNotFoundError(f"No example file named {name!r}. Available: {', '.join(NAMES)}.")
    return str(resources.files(__name__) / name)
