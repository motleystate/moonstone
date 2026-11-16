"""
The tools shown in the interface, and the modules of the Python API it documents.

This is the only place where a tool is declared: pages, forms, examples, checks and tests are derived from it.
Options, defaults and help texts are not written here; they are read from each tool's command-line parser.
"""
import argparse
import importlib
import inspect
from dataclasses import dataclass, field
from typing import Callable, Optional

from moonstone import examples, main
from moonstone.analysis import compositional_pipeline
from moonstone.analysis.permanova import read_taxonomy_counts
from moonstone.filtering import read_threshold_analysis
from moonstone.normalization import table as normalization_table
from moonstone.parsers import handleCounts
from moonstone.reads import config as reads_config, run as reads_run, samples as reads_samples, table as reads_table

GUIDE = 'https://moonstone.readthedocs.io/en/latest/user_guide/'

METADATA_RULES = (
    'Comma-separated, one row per sample.',
    'The first column holds the sample IDs, written exactly as in the header of the count table '
    '("01" and "1" are different samples).',
    'Every other column is a variable (disease, age, batch ...).',
)
TAXONOMY_TABLE_RULES = (
    'Tab- or comma-separated, one row per taxon and one column per sample.',
    'Taxonomy in columns named kingdom, phylum, class, order, family, genus, species (at least one of them).',
    'Every other column is a sample and must be numeric: normalized read counts, not relative abundances.',
    'No blank count cell: write 0 when the taxon was not found. An NCBI_taxonomy_ID column is ignored.',
)
FEATURE_TABLE_RULES = (
    'Tab- or comma-separated, one row per feature and one column per sample.',
    'Either a taxonomy table (columns kingdom ... species; counts are summed to the chosen level), or a plain '
    'table whose first column names the features (for example resistance genes).',
    'Sample columns must be numeric, with no blank cell.',
)


@dataclass(frozen=True)
class Tool:
    """One runnable tool: how it is run, what its form shows and what its input looks like."""
    key: str
    title: str
    group: str                    # 'Prepare' or 'Analyse'
    icon: str
    question: str                 # the question the tool answers
    module: str                   # run as: python -m <module>
    build_parser: Callable[[], argparse.ArgumentParser]
    parse: Callable[[list], argparse.Namespace]    # full check; raises SystemExit on a bad command
    fixed: tuple                  # argv tokens always sent
    options: tuple                # parser dests shown on the form, in order
    advanced: tuple               # parser dests shown under "More options"
    required: tuple               # dests whose flag is always sent
    needs_metadata: bool
    outdir_flag: Optional[str]    # None: the output directory is the third positional argument
    count_sep: Optional[str]      # separator the tool reads the count table with; None: tab or comma
    read_samples: Callable[[str], list]    # sample IDs as the tool reads them
    input_rules: tuple
    example: dict                 # {'countfile': name, 'metadata': name or None, 'values': {dest: value}}
    docs: tuple                   # dotted names whose docstrings are shown
    guide_url: str
    outputs: dict                 # file name pattern (relative to the run directory) -> what the file holds
    input_kind: str = 'table'     # 'table', or 'reads_dir': a folder of read files
    input_label: str = 'Count table'
    choice_labels: dict = field(default_factory=dict)    # {option: {choice: text shown in the menu}}
    requires: tuple = ()          # ((Python module, how to install it), ...): needed to run the tool


def _feature_samples(path):
    return list(compositional_pipeline.read_feature_table(path)[0].columns)


READS_RULES = (
    'A folder of sequencing reads (FASTQ, gzipped or not); its sub-folders are searched.',
    'Paired-end files end with _1/_2, _R1/_R2 or _R1_001/_R2_001 before the extension; single-end files have no '
    'such ending. Paired and single-end reads cannot be mixed in one run.',
    'The sample name is the file name without that ending; it may contain letters, digits, ".", "_" and "-".',
    'Also needed: a Kraken2 database with its Bracken k-mer file, and, to remove host reads, a Kraken2 database '
    'of the host genomes.',
)

TOOLS = {tool.key: tool for tool in [
    Tool(
        key='reads', title='Reads to counts', group='From reads', icon=':material/genetics:',
        question='Which taxa are in my sequencing reads, and how many reads of each?',
        module='moonstone.reads.run', build_parser=reads_run.build_parser,
        parse=lambda argv: reads_run.build_parser().parse_args(argv),
        fixed=(),
        options=('db', 'host_db', 'remove_host', 'match', 'single_end', 'keep_trimmed_reads', 'keep_kraken_output',
                 'slurm'),
        advanced=('no_trimming', 'no_bracken', 'read_len', 'bracken_threshold', 'min_reads', 'min_classified_pct',
                  'cores', 'dry_run'),
        # always sent: the database, so the command shown is complete; and --remove-host, which given with no
        # name switches host removal off
        required=('db', 'remove_host'),
        needs_metadata=False, outdir_flag=None, count_sep=None,
        read_samples=lambda path: [sample.name for sample in reads_samples.discover(path)[0]],
        input_rules=READS_RULES, example={},      # no example: reads and databases are too large to ship
        docs=('moonstone.reads.run.check_inputs', 'moonstone.reads.samples.discover', 'moonstone.reads.config.build',
              'moonstone.reads.launch.submit', 'moonstone.reads.launch.status'),
        guide_url=GUIDE + 'reads.html',
        outputs={
            'output/pipeline_report.md': 'The run report: versions, settings, reads per sample, % host, % classified, '
                                         'warnings.',
            'output/mapping/all_samples.bracken.tsv': 'Bracken-corrected raw counts of all samples: the result of '
                                                      'the run.',
            'output/mapping/all_samples.kraken.tsv': 'Kraken2 counts of all samples, before the Bracken correction.',
        },      # samples.tsv (the read files of each sample) and config.yaml (the settings) are listed below them
        input_kind='reads_dir', input_label='Reads folder',
        choice_labels={'remove_host': {name: f'{name.capitalize()} ({taxid})'
                                       for name, taxid in reads_config.HOSTS.items()}},
        requires=(('snakemake', 'pip install "moonstone[reads]", and Kraken2, Bracken, Trimmomatic and KrakenTools '
                                '(see environment-hpc.yml)'),),
    ),
    Tool(
        key='reads_table', title='Counts to Moonstone table', group='From reads', icon=':material/table_convert:',
        question="How do I get the counts of a reads run in the layout Moonstone's analyses read?",
        module='moonstone.reads.table', build_parser=reads_table.build_parser,
        parse=lambda argv: reads_table.build_parser().parse_args(argv),
        fixed=(), options=(), advanced=(), required=(),
        needs_metadata=False, outdir_flag=None, count_sep='\t',
        read_samples=lambda path: [column for column in reads_table.read_moonbeam_table(path).columns
                                   if column != reads_table.ID_COLUMN],
        input_rules=(
            'The table of a reads run: output/mapping/all_samples.bracken.tsv (or all_samples.kraken.tsv).',
            'Tab-separated, one row per taxon; a first comment line, then a header starting with #OTU ID; the last '
            'column, taxonomy, holds the lineage.',
            'The counts are written unchanged; only the layout and the sample names change.',
        ),
        example={'countfile': 'moonbeam_bracken.tsv', 'metadata': None, 'values': {}},
        docs=('moonstone.reads.table.to_moonstone_table', 'moonstone.reads.table.read_moonbeam_table'),
        guide_url=GUIDE + 'reads.html',
        outputs={'counts_raw.csv': 'The raw counts in Moonstone layout: taxonomy columns, NCBI_taxonomy_ID, one '
                                   'column per sample.',
                 'counts_raw.csv.methods.txt': 'What was done, to quote in your methods.'},
        input_label='Table of a reads run',
    ),
    Tool(
        key='read_threshold', title='Read threshold', group='Prepare', icon=':material/filter_alt:',
        question='Which taxa are prevalent and abundant enough to keep?',
        module='moonstone.filtering.read_threshold_analysis',
        build_parser=read_threshold_analysis.build_parser,
        parse=lambda argv: read_threshold_analysis.build_parser().parse_args(argv),
        fixed=('--save-plots', '--save-filtered'), options=(), advanced=(), required=(),
        needs_metadata=False, outdir_flag='--output-dir', count_sep=',',
        read_samples=lambda path: list(read_threshold_analysis.load_species_table(path).columns),
        input_rules=(
            'Comma-separated, one row per taxon and one column per sample.',
            'Columns species and NCBI_taxonomy_ID are required; kingdom ... genus are optional and ignored.',
            'Every other column is a sample: read counts (raw or normalized), not relative abundances.',
            'No blank count cell: write 0 when the taxon was not found.',
        ),
        example={'countfile': 'taxa_counts.csv', 'metadata': None, 'values': {}},
        docs=('moonstone.filtering.read_threshold_analysis',
              'moonstone.filtering.read_threshold_analysis.analyze_normalized_reads',
              'moonstone.filtering.read_threshold_analysis.rate_threshold_index'),
        guide_url=GUIDE + 'available_classes.html#choosing-a-threshold-from-the-data',
        outputs={
            'filtered_taxa.csv': 'The taxa kept, with their prevalence, mean reads and total reads.',
            'combined_score_knee.png': 'Sorted scores and the knee that sets the threshold.',
            'combined_score_derivative.png': 'The rate-of-change criterion that locates the knee.',
            'combined_score_retained.png': 'Prevalence against abundance, with the taxa kept and removed.',
        },
    ),
    Tool(
        key='normalize', title='Normalize', group='Prepare', icon=':material/balance:',
        question='How do I make the counts of samples sequenced to different depths comparable?',
        module='moonstone.normalization.table', build_parser=normalization_table.build_parser,
        parse=lambda argv: normalization_table.build_parser().parse_args(argv),
        fixed=(), options=('method',), advanced=(), required=('method',),
        needs_metadata=False, outdir_flag=None, count_sep=None,
        read_samples=lambda path: [str(column) for column in normalization_table.read_counts(path).columns],
        input_rules=(
            'Any Moonstone count table, tab- or comma-separated: one row per feature, one column per sample.',
            'The rows are identified by the taxonomy columns (kingdom ... species) and NCBI_taxonomy_ID when '
            'present, otherwise by the first column.',
            'Sample columns must be numeric with no blank cell, and no sample may be empty.',
        ),
        example={'countfile': 'taxa_counts.csv', 'metadata': None, 'values': {'method': 'geometric_mean'}},
        docs=('moonstone.normalization.table.normalize', 'moonstone.normalization.table.normalize_table',
              'moonstone.normalization.counts.geometric_mean.GeometricMeanNormalization'),
        guide_url=GUIDE + 'reads.html',
        outputs={'counts_*.csv': 'The normalized table, with the same identifying columns as the input.',
                 'counts_*.csv.methods.txt': 'What was done, with the parameters actually used, and what the table '
                                             'suits.'},
    ),
    Tool(
        key='classify', title='Classification', group='Analyse', icon=':material/category:',
        question='Can the microbial profile of a sample predict a metadata variable?',
        module='moonstone.main', build_parser=main.build_parser, parse=main.parse_arguments,
        fixed=('-os',),
        options=('classify', 'model', 'kernel', 'roc', 'group_column'),
        advanced=('explore_kernels', 'feature_selection', 'folds', 'seed', 'filtering'),
        required=('classify',),
        needs_metadata=True, outdir_flag=None, count_sep=',',
        read_samples=lambda path: [str(s) for s in handleCounts.Inputs(path).opencounts().index],
        input_rules=(
            'Comma-separated, one row per feature (taxon or gene) and one column per sample.',
            'The first column names the features; the header row holds the sample IDs.',
            'Values are normalized counts. Rows without a feature name are dropped.',
        ),
        example={'countfile': 'resistance_genes.csv', 'metadata': 'metadata.csv',
                 'values': {'classify': 'disease', 'model': ['svm', 'random_forest'], 'roc': True,
                            'group_column': 'patient'}},
        docs=('moonstone.analysis.classification.Classification',),
        guide_url=GUIDE + 'classification.html',
        outputs={
            '*_cv_scores.csv': 'Balanced accuracy, accuracy and area under the ROC curve for each test fold.',
            '*_permutation_importance.csv': 'How much the score drops when each feature is shuffled.',
            '*_permutation_importance.pdf': 'The most important features, as a plot.',
            '*_roc.pdf': 'Receiver Operating Characteristic curves on the test folds.',
        },
    ),
    Tool(
        key='permanova', title='PERMANOVA', group='Analyse', icon=':material/scatter_plot:',
        question='Do groups of samples differ in overall composition?',
        module='moonstone.main', build_parser=main.build_parser, parse=main.parse_arguments,
        fixed=('-os',),
        options=('permanova', 'covariates', 'strata', 'taxa_level', 'select', 'group_column'),
        advanced=('categorical', 'no_read_threshold', 'permutations', 'threads', 'seed'),
        required=('permanova',),
        needs_metadata=True, outdir_flag=None, count_sep=None,
        read_samples=lambda path: list(read_taxonomy_counts(path).columns),
        input_rules=TAXONOMY_TABLE_RULES,
        example={'countfile': 'species_counts.tsv', 'metadata': 'metadata.csv',
                 'values': {'permanova': ['disease'], 'covariates': ['age', 'antibiotics'], 'strata': 'batch',
                            'select': ['Type=inclusion'], 'group_column': 'patient'}},
        docs=('moonstone.analysis.permanova.run_permanova', 'moonstone.analysis.permanova.read_taxonomy_counts'),
        guide_url=GUIDE + 'permanova.html',
        outputs={'permanova.csv': 'One row per taxonomic level and model term: R2, F, p-values and the '
                                  'dispersion check (PERMDISP).'},
    ),
    Tool(
        key='correlation', title='Correlation', group='Analyse', icon=':material/hub:',
        question='Which features rise and fall together across samples?',
        module='moonstone.main', build_parser=main.build_parser, parse=main.parse_arguments,
        fixed=('-os', '--correlation'),
        options=('correlation_method', 'second_countfile', 'select', 'group_column'),
        advanced=('taxa_level', 'min_prevalence', 'top', 'correction', 'alpha', 'permutations', 'seed'),
        required=(),
        needs_metadata=True, outdir_flag=None, count_sep=None,
        read_samples=_feature_samples, input_rules=FEATURE_TABLE_RULES,
        example={'countfile': 'species_counts.tsv', 'metadata': 'metadata.csv',
                 'values': {'second_countfile': examples.path('resistance_genes.csv'),
                            'select': ['Type=inclusion'], 'group_column': 'patient'}},
        docs=('moonstone.analysis.correlation.feature_correlation',
              'moonstone.analysis.compositional_pipeline.read_feature_table'),
        guide_url=GUIDE + 'correlation.html',
        outputs={
            'correlation_edges.csv': 'One row per pair of features: r, permutation p-value and q-value.',
            'correlation_matrix.csv': 'The same r values as a feature-by-feature matrix.',
            'correlation_gephi_edges.csv': 'The significant pairs, ready to import in Gephi.',
            'methods.txt': 'A methods paragraph describing exactly what was run.',
        },
    ),
    Tool(
        key='associate', title='Association', group='Analyse', icon=':material/stacked_line_chart:',
        question='Which features are associated with a metadata variable, after adjusting for covariates?',
        module='moonstone.main', build_parser=main.build_parser, parse=main.parse_arguments,
        fixed=('-os',),
        options=('associate', 'covariates', 'reference', 'select', 'group_column'),
        advanced=('categorical', 'standardize', 'ancombc', 'taxa_level', 'min_prevalence', 'top', 'correction',
                  'alpha'),
        required=('associate',),
        needs_metadata=True, outdir_flag=None, count_sep=None,
        read_samples=_feature_samples, input_rules=FEATURE_TABLE_RULES,
        example={'countfile': 'species_counts.tsv', 'metadata': 'metadata.csv',
                 'values': {'associate': ['disease'], 'covariates': ['age', 'antibiotics'],
                            'reference': ['disease=healthy'], 'ancombc': True, 'select': ['Type=inclusion'],
                            'group_column': 'patient'}},
        docs=('moonstone.analysis.association.feature_association', 'moonstone.analysis.association.ancombc_check'),
        guide_url=GUIDE + 'correlation.html',
        outputs={
            'association.csv': 'One row per feature and model term: coefficient, confidence interval, p- and '
                               'q-value (and the ANCOM-BC columns with --ancombc).',
            'ancombc.csv': 'The ANCOM-BC sensitivity check on the same model (with --ancombc).',
            'methods.txt': 'A methods paragraph describing exactly what was run.',
        },
    ),
]}

API_MODULES = {
    'Parsers': ('moonstone.parsers',),
    'Validation': ('moonstone.validation',),
    'Filtering': ('moonstone.filtering', 'moonstone.filtering.basics_filtering',
                  'moonstone.filtering.read_threshold_analysis'),
    'Normalization': ('moonstone.normalization.counts.geometric_mean', 'moonstone.normalization.counts.total_counts',
                      'moonstone.normalization.counts.random_selection'),
    'Transformers': ('moonstone.transformers.compositional',),
    'Analysis': ('moonstone.analysis.diversity', 'moonstone.analysis.classification', 'moonstone.analysis.permanova',
                 'moonstone.analysis.correlation', 'moonstone.analysis.association',
                 'moonstone.analysis.differential_analysis'),
    'Plots': ('moonstone.plot',),
}


def api_objects() -> dict:
    """
    The public classes and functions of the Python API, by group: {group: [(dotted import path, object)]}.

    A module lists what it defines; a package also lists what its ``__init__`` re-exports. Each object is
    listed once, under the first (shortest) import path.
    """
    seen, groups = set(), {}
    for group, module_names in API_MODULES.items():
        groups[group] = []
        for module_name in module_names:
            module = importlib.import_module(module_name)
            for name, obj in vars(module).items():
                if name.startswith('_') or not (inspect.isclass(obj) or inspect.isfunction(obj)):
                    continue
                defined_in = getattr(obj, '__module__', '')
                exported = hasattr(module, '__path__') and defined_in.startswith('moonstone.')
                if (defined_in == module_name or exported) and id(obj) not in seen:
                    seen.add(id(obj))
                    groups[group].append((f'{module_name}.{name}', obj))
    return groups


def resolve(dotted: str):
    """The module, class or function named by a dotted path."""
    try:
        return importlib.import_module(dotted)
    except ImportError:
        module, _, name = dotted.rpartition('.')
        return getattr(importlib.import_module(module), name)
