"""
Tutorials: runnable documentation of the tools and of the Python API.

Every code step is executed on the example files installed with Moonstone (moonstone.examples) when a tutorial
is opened, and its real result is shown under the code. Nothing here is a pasted output, so a tutorial cannot
describe a behaviour the code no longer has: tests/ui/test_tutorials.py runs them all.
"""
import ast
import contextlib
import io
import warnings
from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class Step:
    """One step of a tutorial: Markdown text, then optional Python code whose last expression is shown."""
    text: str
    code: str = ''


@dataclass(frozen=True)
class Tutorial:
    """A tutorial. ``tool`` is the catalog key of the tool whose page also shows it."""
    key: str
    title: str
    group: str
    summary: str
    steps: tuple
    tool: Optional[str] = None


def run_tutorial(tutorial: Tutorial) -> list:
    """
    Run the code of every step, in order and in one shared namespace.

    :return: one (value of the last expression or None, text printed) per step.
    :raise RuntimeError: a step failed; the message names the tutorial and the step.
    """
    import matplotlib
    matplotlib.use('Agg')    # plots are saved or returned, never opened in a window
    namespace, results = {'__name__': f'tutorial_{tutorial.key}'}, []
    for number, step in enumerate(tutorial.steps, start=1):
        value, printed = None, io.StringIO()
        if step.code:
            try:
                tree = ast.parse(step.code)
                last = tree.body.pop() if isinstance(tree.body[-1], ast.Expr) else None
                with contextlib.redirect_stdout(printed), warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    exec(compile(tree, f'<{tutorial.key} step {number}>', 'exec'), namespace)
                    if last is not None:
                        value = eval(compile(ast.Expression(last.value), f'<{tutorial.key} step {number}>', 'eval'),
                                     namespace)
            except Exception as error:
                raise RuntimeError(f'Tutorial {tutorial.key!r}, step {number} failed: '
                                   f'{type(error).__name__}: {error}') from error
        results.append((value, printed.getvalue()))
    return results


# The tables used by most tutorials.
SETUP = '''\
from moonstone import examples
from moonstone.analysis.permanova import read_taxonomy_counts
from moonstone.parsers import MetadataParser

counts = read_taxonomy_counts(examples.path('species_counts.tsv'))
metadata = MetadataParser(examples.path('metadata.csv'), sep=',', index_col='sample').dataframe
'''
ONE_SAMPLE_PER_PATIENT = '''\
# each patient was sampled twice; the tests below need independent samples, so keep one time point
inclusion = metadata[metadata['Type'] == 'inclusion']
'''
EXAMPLE_NOTE = ('The example files are installed with Moonstone and all their values are simulated: 30 patients, '
                'each sampled at inclusion and at delivery, 12 of them with bacterial vaginosis (`bv`).')

TUTORIALS = {tutorial.key: tutorial for tutorial in [
    # ------------------------------------------------------------------------------------------------- Start here
    Tutorial(
        'first_analysis', 'Your first analysis, start to finish', 'Start here',
        'From two files to a result, in the interface and in Python. Read this one after a long break.',
        (
            Step('Moonstone starts from two files: a **count table** (taxa or genes in rows, samples in columns) '
                 f'and a **metadata** file (one row per sample). {EXAMPLE_NOTE} `examples.path(name)` gives the '
                 'location of an example file.\n\nHere is the count table: the taxonomy is the row index, and each '
                 'column is a sample.',
                 SETUP + 'counts'),
            Step('And the metadata. The first column of the file holds the sample IDs; the other columns are the '
                 'variables you can test.',
                 'metadata'),
            Step('**Check that the two files go together** before anything else. Samples are matched by ID, as '
                 'text: `01` and `1` are different samples, and nothing is renamed for you. An empty list from '
                 '`check_file_format` means the file is in the expected format.\n\nIn the interface this is the '
                 '*Check my files* page, and every tool page runs the same checks as you fill it in.',
                 '''\
from moonstone.validation import check_file_format, match_samples

print('format findings:', check_file_format(examples.path('species_counts.tsv')))
matching = match_samples(counts.columns, metadata.index)
print(len(matching.matched), 'samples are in both files')
'''),
            Step('**Choose the tool from the question.**\n\n'
                 '| Question | Tool |\n|---|---|\n'
                 '| Which taxa are worth keeping? | Read threshold |\n'
                 '| Do groups differ in overall composition? | PERMANOVA |\n'
                 '| Can the profile predict a variable? | Classification |\n'
                 '| Which features rise and fall together? | Correlation |\n'
                 '| Which features follow a variable, after adjusting for others? | Association |\n\n'
                 'Each has its own page, with its own tutorial.'),
            Step('**One sample per subject.** PERMANOVA, correlation and association assume independent samples. '
                 'Here every patient has two, so we keep the inclusion samples.',
                 ONE_SAMPLE_PER_PATIENT + "inclusion['disease'].value_counts()"),
            Step('**Ask the question.** Does the community differ between `bv` and `healthy`? `R2` is the share of '
                 'the variation between samples that the variable explains; `p_value_corrected` is the p-value after '
                 'correction for every test of the run. (999 permutations keep this tutorial fast; the default is '
                 '9,999.)',
                 '''\
from moonstone.analysis.permanova import run_permanova

results = run_permanova(counts, inclusion, columns=['disease'], permutations=999)
results[['level', 'variable', 'n_samples', 'R2', 'pseudo_F', 'p_value', 'p_value_corrected']]
'''),
            Step('**The same in the interface.** Open *PERMANOVA*, press *Use the example files*, then *Run*. The '
                 'page shows the command it runs, which you can copy into a terminal or a cluster job:\n\n'
                 '```bash\nmoonstone species_counts.tsv metadata.csv results --permanova disease '
                 '--select Type=inclusion --group_column patient\n```\n\n'
                 'Every run is written to a new directory that also holds `command.txt` (the exact command, the '
                 'Moonstone version and the time) and `run.log`. *Previous runs* reopens any of them.'),
        )),
    # ------------------------------------------------------------------------------------------------- Input files
    Tutorial(
        'input_tables', 'What the input files must look like', 'Input files',
        'The count-table formats, which tool reads which, and the metadata file.',
        (
            Step('**Taxonomy table**: one row per taxon, the taxonomy in columns named `kingdom` ... `species`, then '
                 'one column per sample. Tab- or comma-separated. Read by PERMANOVA, correlation and association. '
                 'This is the file as written:',
                 '''\
import pandas as pd
from moonstone import examples

pd.read_csv(examples.path('species_counts.tsv'), sep='\\t').iloc[:, :10]
'''),
            Step('`read_taxonomy_counts` is the reader those tools use. The taxonomy becomes the row index, which '
                 'lets Moonstone sum the counts to any level. Sample columns must be numeric with no blank cell '
                 '(write 0 when a taxon was not found), and must be counts, not relative abundances.',
                 '''\
from moonstone.analysis.permanova import read_taxonomy_counts

read_taxonomy_counts(examples.path('species_counts.tsv'))
'''),
            Step('**Plain feature table**: the first column names the features (here resistance genes), then one '
                 'column per sample. Correlation and association accept it too: `read_feature_table` reads either '
                 'kind, and returns the table with notes for the methods paragraph.',
                 '''\
from moonstone.analysis.compositional_pipeline import read_feature_table

genes, notes = read_feature_table(examples.path('resistance_genes.csv'))
genes
'''),
            Step('**Classification** reads a plain feature table, comma-separated only, and turns it so that '
                 'samples are rows:',
                 '''\
from moonstone.parsers import handleCounts

handleCounts.Inputs(examples.path('resistance_genes.csv')).opencounts()
'''),
            Step('**Read threshold** needs a comma-separated table with `species` and `NCBI_taxonomy_ID` columns '
                 '(the other taxonomy columns are optional and ignored). Rows with the same species and ID are '
                 'summed.',
                 '''\
from moonstone.filtering.read_threshold_analysis import load_species_table

load_species_table(examples.path('taxa_counts.csv'))
'''),
            Step('**Metadata**: comma-separated, one row per sample, sample IDs in the first column, written '
                 'exactly as in the header of the count table. Text columns are groups; numeric columns are '
                 'numbers (a numeric column of group codes such as 1/2/3 must be declared with `--categorical`).',
                 '''\
from moonstone.parsers import MetadataParser

metadata = MetadataParser(examples.path('metadata.csv'), sep=',', index_col='sample').dataframe
metadata.dtypes.rename('type').to_frame()
'''),
            Step('**In every file**: UTF-8 text, a dot as decimal mark, ISO dates (yyyy-mm-dd). Files saved by '
                 'Excel in a French locale differ on all of these: see *Checking and repairing a file saved by '
                 'Excel*.'),
        )),
    Tutorial(
        'check_and_repair', 'Checking and repairing a file saved by Excel', 'Input files',
        'check_file_format, convert_to_utf8, convert_locale and match_samples on a French-locale file.',
        (
            Step('A typical export from Excel in a French locale: cp1252 encoding, semicolons, decimal commas, '
                 'spaces as thousands separators, day-first dates. We write one to a temporary folder.',
                 '''\
import os
import tempfile

folder = tempfile.mkdtemp()
source = os.path.join(folder, 'export_excel.csv')
text = ('sample;âge;lectures;date\\n'
        'S001;25;12 500,5;03/02/2024\\n'
        'S002;31;8 730,25;17/02/2024\\n'
        'S003;28;950,0;21/02/2024\\n')
with open(source, 'w', encoding='cp1252') as handle:
    handle.write(text)
print(text)
'''),
            Step('`check_file_format` reads the file and reports what it finds, with the first offending row and '
                 'value. It never changes anything. An `error` means the file will load wrongly or not at all.',
                 '''\
import pandas as pd
from moonstone.validation import check_file_format

pd.DataFrame(check_file_format(source))[['level', 'code', 'column', 'row', 'example', 'message']]
'''),
            Step('**Step 1: the encoding.** `convert_to_utf8` writes a new UTF-8 file. The source encoding must be '
                 'given: cp1252, latin-1 and Mac Roman cannot be told apart reliably, so it is never guessed. '
                 'The result lists every line that changed.',
                 '''\
from moonstone.validation import convert_to_utf8

as_utf8 = os.path.join(folder, 'export_utf8.csv')
result = convert_to_utf8(source, as_utf8, 'cp1252')
result.changes_frame()
'''),
            Step('**Step 2: numbers, dates and separator.** `convert_locale` rewrites the file with English '
                 'conventions, from the conventions *you declare*. A column is converted only when every value '
                 'fits; otherwise it is left unchanged and reported in `result.notes`. Every changed value is '
                 'listed, so the conversion can be audited.',
                 '''\
from moonstone.validation import convert_locale

clean = os.path.join(folder, 'export_clean.csv')
result = convert_locale(as_utf8, clean, sep=';', decimal=',', date_order='dmy')
result.changes_frame()
'''),
            Step('The new file, and its check: an empty list means nothing is left to report. Neither function '
                 'overwrites a file: asking for an output that exists raises `FileExistsError`.',
                 '''\
with open(clean, encoding='utf-8') as handle:
    print(handle.read())
print('findings:', check_file_format(clean))
'''),
            Step('**Sample IDs.** `match_samples` pairs the count table and the metadata by ID, exactly as written. '
                 'IDs that differ only by case or spaces are *reported*, never matched: correct the files if they '
                 'are the same sample.',
                 '''\
from moonstone.validation import match_samples

matching = match_samples(['S001', 'S002', 'S31', 'S040'], ['S001', 's31', 'S002 ', 'S050'])
pd.DataFrame(matching.findings())[['level', 'code', 'message']]
'''),
        )),
    Tutorial(
        'parsers', 'Reading files with the parsers', 'Input files',
        'MetadataParser, GeneCountsParser and SunbeamKraken2Parser; the other parsers work the same way.',
        (
            Step('Every parser works the same way: give it a path, read the table from its `dataframe` property. '
                 'The parsers default to tab-separated files; pass `sep` for another separator.\n\n'
                 '`MetadataParser` indexes the table by the sample column.',
                 '''\
from moonstone import examples
from moonstone.parsers import MetadataParser

parser = MetadataParser(examples.path('metadata.csv'), sep=',', index_col='sample')
parser.dataframe
'''),
            Step('`get_stats` describes each column: its type, how many values and how many different ones. A '
                 'quick way to spot a column that is not what you expect.',
                 '''\
import pandas as pd

pd.DataFrame(parser.get_stats())[['col_name', 'python_col_type', 'n_values', 'n_uniq_values']]
'''),
            Step('Metadata typed by hand is rarely consistent. `cleaning_operations` applies transformations while '
                 'reading, as `{column: [(operation, options)]}`. Available operations are `to_slug` (lower case, '
                 'no spaces or accents) and `remove_trailing_spaces`.',
                 '''\
cleaned = MetadataParser(examples.path('metadata.csv'), sep=',', index_col='sample',
                         cleaning_operations={'patient': [('to_slug', {})]})
cleaned.dataframe['patient'].head()
'''),
            Step('`GeneCountsParser` reads a plain matrix: one row per gene, one column per sample.',
                 '''\
from moonstone.parsers import GeneCountsParser

GeneCountsParser(examples.path('resistance_genes.csv'), sep=',').dataframe
'''),
            Step('`SunbeamKraken2Parser` reads a Kraken2 table merged by the Sunbeam pipeline, where the lineage is '
                 'one text column (`k__Bacteria; p__Firmicutes; ...`). It splits the lineage into the taxonomy '
                 'index used everywhere else in Moonstone, and keeps the NCBI identifier as a column.',
                 '''\
from moonstone.parsers import SunbeamKraken2Parser

kraken2 = SunbeamKraken2Parser(examples.path('kraken2_counts.tsv')).dataframe
kraken2.drop(columns='NCBI_taxonomy_ID')
'''),
            Step('**Other parsers** follow the same pattern; no example file in their format is shipped, so they '
                 'are not run here:\n\n'
                 '- `Metaphlan2Parser`, `Metaphlan3Parser`: merged MetaPhlAn tables (`analysis_type` names the '
                 'MetaPhlAn output type).\n'
                 '- `MetaphlanReadStatsParser`: per-sample MetaPhlAn profiles made with `-t rel_ab_w_read_stats`, '
                 'joined into a table of estimated read counts.\n'
                 '- `Qiime2Parser`: a QIIME 2 table exported as CSV.\n'
                 '- `Picrust2PathwaysParser`: PICRUSt2 predicted pathway abundances (pathways by samples).\n'
                 '- `YAMLBasedMetadataParser`: a metadata parser whose cleaning operations are written in a YAML '
                 'file.\n\nTheir arguments are on the *Python API* page.'),
        )),
    # ------------------------------------------------------------------------------------------ Prepare the table
    Tutorial(
        'filtering', 'Filtering taxa, genes and samples', 'Prepare the table',
        'The filtering classes: by counts, by name, by mean, on a taxonomy level, and on metadata columns.',
        (
            Step('Every filter works the same way: `Filter(table, ...).filtered_df`. Filters are chained by giving '
                 'the result of one to the next.\n\n`NoCountsFiltering` removes rows and columns that hold only '
                 'zeros. We add an empty gene to see it go.',
                 SETUP + '''\
from moonstone.filtering import NoCountsFiltering
from moonstone.parsers import GeneCountsParser

genes = GeneCountsParser(examples.path('resistance_genes.csv'), sep=',').dataframe
with_empty_row = genes.copy()
with_empty_row.loc['mecA'] = 0.0
print('before:', list(with_empty_row.index))
print('after: ', list(NoCountsFiltering(with_empty_row).filtered_df.index))
'''),
            Step('`NamesFiltering` keeps (or, with `keep=False`, removes) rows by name; `axis=1` works on columns.',
                 '''\
from moonstone.filtering import NamesFiltering

without_two_samples = NamesFiltering(genes, ['S001', 'S002'], axis=1, keep=False).filtered_df
NamesFiltering(without_two_samples, ['tetM', 'ermB'], axis=0, keep=True).filtered_df
'''),
            Step('`MeanFiltering` removes items whose mean count over the samples is below a threshold. '
                 '`report_data` says what was removed. (Without a threshold, give `percentage_to_keep` and the '
                 'threshold is computed for you.)',
                 '''\
from moonstone.filtering import MeanFiltering

mean_filter = MeanFiltering(counts, threshold=100)
print(mean_filter.filtered_df.index.get_level_values('species').tolist())
mean_filter.report_data
'''),
            Step('On a taxonomy table, `TaxonomyNamesFiltering` selects by name at any level ...',
                 '''\
from moonstone.filtering import TaxonomyNamesFiltering

TaxonomyNamesFiltering(counts, ['Lactobacillus'], level='genus', keep=True).filtered_df
'''),
            Step('... and `TaxonomyMeanFiltering` keeps the taxa whose mean at a chosen level reaches a value. To '
                 'choose that value from the data instead of by hand, see the *Read threshold* tutorial.',
                 '''\
from moonstone.filtering import TaxonomyMeanFiltering

TaxonomyMeanFiltering(counts, 500, level='species').filtered_df
'''),
            Step('Two filters are useful on **metadata**. `NumberOfDifferentValuesFiltering` keeps the columns '
                 '(`axis=1`) with a number of different values in a range: here the variables with 2 or 3 groups. '
                 '`NaNPercentageFiltering` removes the columns with too many missing values.',
                 '''\
from moonstone.filtering.basics_filtering import NaNPercentageFiltering, NumberOfDifferentValuesFiltering

complete_enough = NaNPercentageFiltering(metadata, percentage=20, axis=1).filtered_df
NumberOfDifferentValuesFiltering(complete_enough, min=2, max=3, axis=1).filtered_df
'''),
        )),
    Tutorial(
        'read_threshold', 'Choosing a filtering threshold from the data', 'Prepare the table',
        'The read-threshold method: keep the taxa that are both prevalent and abundant, without a hand-picked cut.',
        (
            Step('**The question.** Which taxa are real members of these communities, and which are noise? A fixed '
                 'cut ("mean above 10 reads") is arbitrary. This method scores each taxon as '
                 '`prevalence x log10(mean reads + 1)` and keeps the taxa above the **knee** of the sorted scores: '
                 'the point where the score starts rising faster than its average rate. A taxon must be both '
                 'prevalent and abundant to pass.\n\n**The input**: a comma-separated table with `species` and '
                 '`NCBI_taxonomy_ID` columns and one column per sample, holding read counts (raw or normalized), '
                 'not relative abundances.',
                 '''\
import pandas as pd
from moonstone import examples

pd.read_csv(examples.path('taxa_counts.csv')).iloc[:, 6:12]
'''),
            Step('`analyze_normalized_reads` returns the statistics of the cut and the table of taxa kept. Here we '
                 'also save the plots and the table to a temporary folder, as the interface does in the run '
                 'directory.',
                 '''\
import os
import tempfile
from moonstone.filtering.read_threshold_analysis import analyze_normalized_reads

folder = tempfile.mkdtemp()
stats, kept = analyze_normalized_reads(examples.path('taxa_counts.csv'), save_plots=True, save_filtered=True,
                                       output_dir=folder)
kept
'''),
            Step('**Reading the statistics.** `species_removed` and `reads_pct_removed` say what the cut costs: a '
                 'good cut removes many taxa and few reads. `min_reads_retained` and `min_prevalence_retained` '
                 'describe the weakest taxon still kept.',
                 "pd.Series(stats).rename('value').to_frame()"),
            Step('**The knee plot** shows the sorted scores and where the cut falls. Look at it before trusting '
                 'the threshold: with few taxa, or scores that rise evenly, the knee is not well defined. (The '
                 'example has only eight species, far fewer than a real table.)',
                 "os.path.join(folder, 'combined_score_knee.png')"),
            Step('**Which taxa were kept**, on prevalence against abundance:',
                 "os.path.join(folder, 'combined_score_retained.png')"),
            Step('**From the command line or the interface**\n\n```bash\npython -m '
                 'moonstone.filtering.read_threshold_analysis taxa_counts.csv --save-filtered --save-plots '
                 '--output-dir results/\n```\n\nPERMANOVA applies this same filter at each taxonomic level unless '
                 '`--no_read_threshold` is given. The knee has been checked at species level only, so other '
                 'levels carry a warning in PERMANOVA results.'),
        ), tool='read_threshold'),
    Tutorial(
        'normalization', 'Normalizing counts', 'Prepare the table',
        'TotalCountsNormalization, GeometricMeanNormalization and random selection (rarefaction).',
        (
            Step('Samples are not sequenced to the same depth, so raw counts are not comparable between samples. '
                 'Every normalization works the same way: `Normalization(table).normalized_df`.\n\nThe depth of '
                 'the example samples:',
                 SETUP + "counts.sum().describe().rename('reads per sample').to_frame()"),
            Step('`TotalCountsNormalization` scales each sample by its total, so that all samples end with the '
                 'same total. `scaling_factors` holds the factor applied to each sample.',
                 '''\
import pandas as pd
from moonstone.normalization.counts.total_counts import TotalCountsNormalization

total = TotalCountsNormalization(counts)
pd.DataFrame({'reads before': counts.sum(), 'scaling factor': total.scaling_factors,
              'reads after': total.normalized_df.sum()}).head()
'''),
            Step('`GeometricMeanNormalization` is the median-of-ratios normalization of DESeq2: each sample is '
                 'compared with a reference made of the geometric mean of every taxon, and scaled by its median '
                 'ratio. It is less affected than total counts by one dominant taxon. Totals are *not* equal '
                 'afterwards.',
                 '''\
from moonstone.normalization.counts.geometric_mean import GeometricMeanNormalization

geometric = GeometricMeanNormalization(counts)
pd.DataFrame({'reads before': counts.sum(), 'scaling factor': geometric.scaling_factors,
              'reads after': geometric.normalized_df.sum().round()}).head()
'''),
            Step('**Random selection (rarefaction)** draws the same number of reads from every sample, at random, '
                 'in proportion to the counts. `TaxonomyRandomSelection` does it on a taxonomy table; by default '
                 'the number drawn is the depth of the shallowest sample. The seed makes it repeatable.',
                 '''\
from moonstone.normalization.counts.random_selection import TaxonomyRandomSelection

rarefied = TaxonomyRandomSelection(counts, random_seed=2935).normalized_df
pd.DataFrame({'reads before': counts.sum(), 'reads after': rarefied.sum()}).head()
'''),
            Step('`RandomSelection` is the same on a plain table (whole-number counts), with an explicit number of '
                 'reads to draw.',
                 '''\
from moonstone.normalization.counts.random_selection import RandomSelection

species = counts.droplevel(['kingdom', 'phylum', 'class', 'order', 'family', 'genus'])
RandomSelection(species, threshold=10000).normalized_df.sum().head().rename('reads after').to_frame()
'''),
        )),
    Tutorial(
        'transformers', 'Compositional transforms (CLR)', 'Prepare the table',
        'multiplicative_replacement, clr and rclr: what correlation and association do before testing.',
        (
            Step('Relative abundances are **compositional**: the values of a sample share a fixed total, so one '
                 'taxon going up pushes the others down. Correlations computed on them are partly artefacts. The '
                 'centred log-ratio (CLR) expresses each value relative to the geometric mean of its sample '
                 'instead. Logarithms need non-zero values, so zeros are dealt with first.\n\nThe resistance-gene '
                 'table has zeros:',
                 '''\
from moonstone import examples
from moonstone.parsers import GeneCountsParser

genes = GeneCountsParser(examples.path('resistance_genes.csv'), sep=',').dataframe
genes.iloc[:, :6]
'''),
            Step('`multiplicative_replacement` returns proportions (each sample sums to 1) in which every zero is '
                 'replaced by 65% of the smallest non-zero proportion of its sample, the other values being '
                 'scaled down to keep the total.',
                 '''\
from moonstone.transformers.compositional import multiplicative_replacement

multiplicative_replacement(genes).iloc[:, :6]
'''),
            Step('`clr` does both steps: zero replacement, then the centred log-ratio. A positive value means '
                 '"above the geometric mean of this sample". Each sample sums to zero. This is the transform that '
                 '`--correlation` and `--associate` apply.',
                 '''\
from moonstone.transformers.compositional import clr

clr(genes).iloc[:, :6]
'''),
            Step('`rclr` (robust CLR) does not replace zeros: it computes the CLR on the non-zero values of each '
                 'sample and leaves zeros missing (`NaN`), so what you do next must handle missing values.',
                 '''\
from moonstone.transformers.compositional import rclr

rclr(genes).iloc[:, :6]
'''),
        )),
    # ---------------------------------------------------------------------------------------------------- Describe
    Tutorial(
        'plots', 'Plotting counts and metadata', 'Describe',
        'PlotCountsStats, PlotTaxonomyCounts, PlotMetadataStats and the category-flow plots.',
        (
            Step('Plot classes take a table and offer one method per plot. Each method returns a plotly figure. In '
                 'a script, `show=True` (the default) opens it in the browser and `output_file="plot.html"` saves '
                 'it; here `show=False` hands the figure to the page.\n\n`PlotCountsStats`: the reads of each '
                 'sample ...',
                 SETUP + '''\
from moonstone.plot import PlotCountsStats

PlotCountsStats(counts, items_name='species').plot_sample_totals(binned=False, show=False)
'''),
            Step('... and the distribution of the mean reads of the items.',
                 "PlotCountsStats(counts, items_name='species').plot_mean_distribution(show=False)"),
            Step('`PlotTaxonomyCounts` works on a taxonomy table. The composition of each sample, for the most '
                 'abundant taxa:',
                 '''\
from moonstone.plot import PlotTaxonomyCounts

taxonomy_plots = PlotTaxonomyCounts(counts)
taxonomy_plots.plot_sample_composition_most_abundant_taxa(show=False)
'''),
            Step('The most abundant taxa, at any level (`taxa_level`) ...',
                 "taxonomy_plots.plot_most_abundant_taxa(taxa_level='genus', show=False)"),
            Step('... and the most prevalent ones.',
                 "taxonomy_plots.plot_most_prevalent_taxa(show=False)"),
            Step('`PlotMetadataStats` describes the cohort: a numeric column ...',
                 '''\
from moonstone.plot import PlotMetadataStats

metadata_plots = PlotMetadataStats(metadata)
metadata_plots.plot_age(show=False)
'''),
            Step('... or the groups of any column.',
                 "metadata_plots.plot_category_distribution('disease', show=False)"),
            Step('`plot_category_flows` shows how the samples are distributed across several columns at once: '
                 'which groups go together. `kind` is `"sankey"` or `"parcats"`.',
                 '''\
from moonstone.plot import plot_category_flows

plot_category_flows(metadata, ['Type', 'disease', 'antibiotics'], kind='parcats', color_by='disease', show=False)
'''),
            Step('`plot_taxonomy_flows` does the same down the taxonomy: how the reads split from phylum to genus '
                 'to species.',
                 '''\
from moonstone.plot import plot_taxonomy_flows

plot_taxonomy_flows(counts, ['phylum', 'genus', 'species'], show=False)
'''),
            Step('**Changing the look.** Every plot method takes `plotting_options`, a dictionary passed to plotly: '
                 '`traces` for the marks, `xaxes` / `yaxes` for the axes, `layout` for the rest.',
                 '''\
PlotCountsStats(counts).plot_sample_totals(
    binned=False, show=False,
    plotting_options={'traces': {'marker_color': 'teal'}, 'xaxes': {'tickangle': 90}})
'''),
        )),
    Tutorial(
        'alpha_diversity', 'Alpha diversity: how diverse is each sample?', 'Describe',
        'ShannonIndex, SimpsonInverseIndex and Chao1Index, and comparing them between groups.',
        (
            Step('Alpha diversity summarises one sample by one number. Every index works the same way: '
                 '`Index(table).diversity_indexes`. The indices are computed on the table as given, and depend on '
                 'sequencing depth: normalize or rarefy first (see *Normalizing counts*).\n\n`ShannonIndex` grows '
                 'with the number of taxa and with how evenly the reads are spread among them.',
                 SETUP + ONE_SAMPLE_PER_PATIENT + '''\
from moonstone.analysis.diversity import ShannonIndex
from moonstone.normalization.counts.random_selection import TaxonomyRandomSelection

rarefied = TaxonomyRandomSelection(counts[inclusion.index], random_seed=2935).normalized_df
shannon = ShannonIndex(rarefied)
shannon.diversity_indexes.head().to_frame()
'''),
            Step('`SimpsonInverseIndex` is driven by the dominant taxa (1 means a single taxon holds all the '
                 'reads); `Chao1Index` estimates the number of taxa, including those missed, from the taxa seen '
                 'once or twice. With eight species present in almost every sample, Chao1 says little here.',
                 '''\
import pandas as pd
from moonstone.analysis.diversity import Chao1Index, SimpsonInverseIndex

pd.DataFrame({'shannon': shannon.diversity_indexes,
              'simpson_inverse': SimpsonInverseIndex(rarefied).diversity_indexes,
              'chao1': Chao1Index(rarefied).diversity_indexes}).head()
'''),
            Step('**Comparing groups.** `analyse_groups` joins the index with the metadata, tests the groups '
                 'against each other and draws them. Pass `show=False, show_pval=False` in a script that should '
                 'not open the browser. It returns a dictionary; `fig` is the plot ...',
                 '''\
result = shannon.analyse_groups(inclusion, 'disease', show=False, show_pval=False)
result['fig']
'''),
            Step('... `pval` the p-value of each pair of groups, and `param` the test and correction that were '
                 'used, to quote in your methods.',
                 '''\
print(result['param'])
result['pval']
'''),
            Step('`FaithsPhylogeneticDiversity` works the same way but also needs a phylogenetic tree of the taxa '
                 '(a scikit-bio `TreeNode`); no tree is shipped with the examples, so it is not run here.'),
        )),
    Tutorial(
        'beta_diversity', 'Beta diversity: how different are the samples?', 'Describe',
        'BrayCurtis distances, the PCoA plot, and comparing distances within and between groups.',
        (
            Step('Beta diversity compares samples two by two. `BrayCurtis(table).beta_diversity_df` is the matrix '
                 'of distances: 0 for identical communities, 1 for no taxon in common. Like alpha diversity, it '
                 'is computed on the table as given: normalize first.',
                 SETUP + ONE_SAMPLE_PER_PATIENT + '''\
from moonstone.analysis.diversity import BrayCurtis
from moonstone.normalization.counts.total_counts import TotalCountsNormalization

normalized = TotalCountsNormalization(counts[inclusion.index]).normalized_df
bray_curtis = BrayCurtis(normalized)
bray_curtis.beta_diversity_df.iloc[:6, :6]
'''),
            Step('**Principal coordinates analysis (PCoA)** places the samples on a plane so that their distances '
                 'on the plot approximate the Bray-Curtis distances. Colour by a metadata column to see whether '
                 'groups separate, and whether one group is more spread out than the other.',
                 "bray_curtis.visualize_pcoa(inclusion, 'disease', show=False)"),
            Step('`analyse_groups` compares the distances within each group and between groups.',
                 '''\
result = bray_curtis.analyse_groups(inclusion, 'disease', show=False, show_pval=False)
result['fig']
'''),
            Step('**The test** for "do the groups differ as a whole" is PERMANOVA, available on every '
                 'beta-diversity object. It also runs the dispersion check (PERMDISP). See the *PERMANOVA* '
                 'tutorial for how to read it.',
                 "bray_curtis.permanova(inclusion, 'disease', permutations=999).to_frame()"),
            Step('`Jaccard` (in `moonstone.analysis.diversity.beta`) uses presence and absence only. '
                 '`WeightedUniFrac` and `UnweightedUniFrac` take the relatedness of the taxa into account and '
                 'need a phylogenetic tree; no tree is shipped with the examples, so they are not run here. All '
                 'have the same methods as `BrayCurtis`.'),
        )),
    # ----------------------------------------------------------------------------------------------------- Analyse
    Tutorial(
        'permanova', 'PERMANOVA: do groups differ in overall composition?', 'Analyse',
        'The screen, the model adjusted for covariates, the dispersion check, and how to read the table.',
        (
            Step('**The question.** Do samples from different groups (here `bv` and `healthy`) differ in their '
                 'community as a whole? PERMANOVA compares the Bray-Curtis distances between groups with the '
                 'distances within groups (a pseudo-F statistic), and gets its p-value by shuffling the group '
                 'labels many times.\n\n**The input**: a taxonomy count table and the metadata. Samples must be '
                 'independent: one per subject.',
                 SETUP + ONE_SAMPLE_PER_PATIENT + "inclusion.groupby('disease').size().rename('samples').to_frame()"),
            Step('**The screen** tests each variable alone. For each taxonomic level, counts are summed to that '
                 'level, taxa below the read-threshold knee are removed (a filter that never looks at the groups), '
                 'then PERMANOVA and PERMDISP are run.\n\n'
                 '- `R2`: share of the variation between samples explained by the variable.\n'
                 '- `p_value_corrected`: corrected over every test of the run. Use this one.\n'
                 '- `p_value_se`: Monte Carlo error of the p-value; more permutations make it smaller.\n'
                 '- `warnings`: read them before the p-values.\n\n'
                 '(999 permutations keep the tutorial fast; the default is 9,999.)',
                 '''\
from moonstone.analysis.permanova import run_permanova

screen = run_permanova(counts, inclusion, columns=['disease', 'antibiotics'], permutations=999)
screen[['level', 'variable', 'n_samples', 'n_groups', 'pseudo_F', 'R2', 'p_value', 'p_value_se',
        'p_value_corrected', 'warnings']]
'''),
            Step('**The dispersion check.** A significant PERMANOVA means the groups differ in their centre, in '
                 'their spread, or both. PERMDISP tests the spread alone. If `permdisp_p_value_corrected` is also '
                 'small, part or all of the result may come from one group being more variable, not from a '
                 'different composition.',
                 "screen[['variable', 'p_value_corrected', 'permdisp_F', 'permdisp_p_value', "
                 "'permdisp_p_value_corrected']]"),
            Step('Look at the PCoA before concluding: are the groups in different places, or is one simply more '
                 'spread out?',
                 '''\
from moonstone.analysis.diversity import BrayCurtis

BrayCurtis(counts[inclusion.index]).visualize_pcoa(inclusion, 'disease', show=False)
'''),
            Step('**Adjusting for covariates.** To test `disease` after age and antibiotic use, add covariates: the '
                 'variables and covariates form one model, and every term is tested after all the others (as '
                 '`vegan::adonis2` with `by="margin"`). `strata` keeps the shuffles within blocks, here the '
                 'sequencing run; it is not a term of the model.\n\nText columns are groups; numeric columns are '
                 'one slope (`term_type`). A sample with a missing value in any term is left out (`n_missing`).',
                 '''\
adjusted = run_permanova(counts, inclusion, columns=['disease'], covariates=['age', 'antibiotics'],
                         strata='batch', permutations=999)
adjusted[['variable', 'term_type', 'df', 'n_samples', 'n_missing', 'R2', 'pseudo_F', 'p_value',
          'p_value_corrected']]
'''),
            Step('**From the command line or the interface**\n\n```bash\n# the screen\nmoonstone species_counts.tsv '
                 'metadata.csv results --permanova disease antibiotics \\\n    --select Type=inclusion '
                 '--group_column patient\n\n# adjusted for covariates, shuffled within sequencing runs\nmoonstone '
                 'species_counts.tsv metadata.csv results --permanova disease \\\n    --covariates age antibiotics '
                 '--strata batch --select Type=inclusion --group_column patient\n```\n\n'
                 '`--select` keeps one time point, and `--group_column` makes the run stop if a subject still has '
                 'more than one sample. `--permanova` with no variable tests every metadata column. '
                 '`--taxa_level all` runs every level. The table goes to `permanova.csv`.\n\n'
                 '**Limits**: no interaction between terms, and no repeated measures.'),
        ), tool='permanova'),
    Tutorial(
        'classification', 'Classification: can the profile predict a variable?', 'Analyse',
        'Cross-validated SVM and random forest: what the score means and which features carry it.',
        (
            Step('**The question.** Can the counts of a sample predict a metadata variable, for example disease '
                 'status? A model is trained on part of the samples and scored on samples it has not seen '
                 '(cross-validation, 5 folds by default). Only scores on unseen samples are reported.\n\n'
                 '**The input**: a comma-separated feature table (features in rows) and the metadata. These are the '
                 'readers the command line uses; samples become rows.',
                 '''\
import tempfile
from moonstone import examples
from moonstone.parsers import handleCounts, handleMetadata

folder = tempfile.mkdtemp()    # the results are written here
counts = handleCounts.Inputs(examples.path('resistance_genes.csv')).opencounts()
metadata = handleMetadata.Inputs(examples.path('metadata.csv')).openmeta(folder)
counts.head()
'''),
            Step('**Run it.** Here the patients have two samples each. `group_column` keeps both samples of a '
                 'patient in the same fold; otherwise the model can recognise the patient instead of the disease, '
                 'and the score is too good.\n\nThe summary printed is what the command line prints. **Balanced '
                 'accuracy** is the average of the per-group hit rates: its chance level is 1 / number of groups '
                 'whatever the group sizes.',
                 '''\
from moonstone.analysis.classification import Classification

svm = Classification(counts, metadata, 'disease', folder, model='svm', group_column='patient')
results = svm.run(roc=True)
results['folds']
'''),
            Step('**Which features carry the prediction?** Permutation importance shuffles one feature among the '
                 'test samples and records the loss of balanced accuracy. A large loss means the model depends on '
                 'that feature. In the simulation, `tetM` follows *Gardnerella vaginalis*, which is higher with '
                 'bacterial vaginosis.',
                 '''\
import os
import pandas as pd

pd.read_csv(os.path.join(folder, 'disease_svm_permutation_importance.csv'), index_col=0)
'''),
            Step('**Another model.** A random forest is scored on the same folds (same seed), so the two can be '
                 'compared.',
                 '''\
forest = Classification(counts, metadata, 'disease', folder, model='random_forest', group_column='patient')
forest.run()['folds'][['balanced_accuracy', 'accuracy', 'auc']]
'''),
            Step('**The files written**, named `VARIABLE_MODEL_...`. The plots are PDF files; the interface shows '
                 'them on the results of a run.',
                 "pd.Series(sorted(os.listdir(folder))).rename('file').to_frame()"),
            Step('**From the command line or the interface**\n\n```bash\nmoonstone resistance_genes.csv metadata.csv '
                 'results --classify disease \\\n    --model svm random_forest --roc --group_column patient\n```\n\n'
                 '**Things that make a score look better than it is**, and what Moonstone does about them:\n\n'
                 '- *Several samples per subject*: give `--group_column`.\n'
                 '- *Choosing the best of several kernels*: `--explore_kernels` scores are exploratory; '
                 '`--kernel auto` chooses inside each training part and gives a fair score.\n'
                 '- *Selecting features on all samples*: `--feature_selection` selects inside each training part '
                 '(nested cross-validation).\n'
                 '- *Unequal groups*: read balanced accuracy, not accuracy.\n\n'
                 'A group with a single sample cannot be cross-validated and stops the run.'),
        ), tool='classify'),
    Tutorial(
        'correlation', 'Correlation: which features rise and fall together?', 'Analyse',
        'Compositional correlation with permutation p-values, within one table or across two.',
        (
            Step('**The question.** Which taxa (or genes) go up and down together across samples? The result is a '
                 'network: one coefficient per pair of features.\n\n**Why not correlate the counts directly?** '
                 'Abundances are compositional, so plain correlations are partly artefacts. Moonstone works on '
                 'CLR-transformed values (see *Compositional transforms*) and gets p-values by permutation, '
                 'because the CLR itself makes features slightly negatively correlated.\n\n'
                 '**The input**: a feature table and the metadata (used to select samples). The steps below are '
                 'the ones the command line applies, in this order: read, select samples, filter features.',
                 '''\
from moonstone import examples
from moonstone.analysis.compositional_pipeline import filter_features, read_feature_table
from moonstone.parsers import MetadataParser

species, notes = read_feature_table(examples.path('species_counts.tsv'), 'species')
metadata = MetadataParser(examples.path('metadata.csv'), sep=',', index_col='sample').dataframe
''' + ONE_SAMPLE_PER_PATIENT + '''\
species, note = filter_features(species[inclusion.index], min_prevalence=0.1)
print(note)
species
'''),
            Step('**Run it.** One row per pair of features:\n\n'
                 '- `r`: the coefficient (Spearman on CLR values by default).\n'
                 '- `p_value`: from the permutations, which shuffle every feature across samples and transform the '
                 'table again.\n'
                 '- `q_value`: corrected over all pairs (Benjamini-Hochberg). Use this one.\n'
                 '- `prevalence_a`, `prevalence_b`: share of samples where each feature is present.',
                 '''\
from moonstone.analysis.correlation import feature_correlation

edges = feature_correlation(species, permutations=999, seed=0)
edges.sort_values('p_value')
'''),
            Step('`filter_edges` keeps the significant pairs, `edges_to_matrix` turns the table into a matrix, and '
                 '`to_gephi` renames the columns for import in Gephi. The command line writes all three.',
                 '''\
from moonstone.analysis.correlation import edges_to_matrix, filter_edges

print(len(filter_edges(edges, alpha=0.05)), 'significant pairs')
edges_to_matrix(edges).round(2)
'''),
            Step('**Across two tables.** With a second table, for example resistance genes, only the pairs across '
                 'the two tables are tested, and each table is transformed on its own. In the simulation `tetM` '
                 'follows *Gardnerella vaginalis*.',
                 '''\
genes, _ = read_feature_table(examples.path('resistance_genes.csv'))
across = feature_correlation(species, genes[inclusion.index], permutations=999, seed=0)
across.sort_values('p_value').head()
'''),
            Step('**From the command line or the interface**\n\n```bash\n# network of species, at inclusion only\n'
                 'moonstone species_counts.tsv metadata.csv results --correlation \\\n    --select Type=inclusion '
                 '--group_column patient\n\n# species against resistance genes\nmoonstone species_counts.tsv '
                 'metadata.csv results --correlation \\\n    --second_countfile resistance_genes.csv --select '
                 'Type=inclusion\n```\n\n'
                 'Outputs: `correlation_edges.csv` (the table above), `correlation_matrix.csv`, '
                 '`correlation_gephi_edges.csv` and `methods.txt`, a methods paragraph describing exactly what was '
                 'run.\n\n**Good to know**: at least 100 permutations are required (999 by default); '
                 '`--correlation_method` also offers `pearson`, `kendall`, `rho` (proportionality) and `sparcc` '
                 '(which needs read counts); one sample per subject.'),
        ), tool='correlation'),
    Tutorial(
        'association', 'Association: which features follow a variable?', 'Analyse',
        'One linear model per feature on CLR values, adjusted for covariates, with ANCOM-BC as a check.',
        (
            Step('**The question.** Which taxa are linked to a clinical variable, once other variables are '
                 'accounted for? This is a regression, not a correlation: one linear model per feature on CLR '
                 'values (the linear model of MaAsLin2, without random effect), with the covariates you choose.\n\n'
                 '**The input**: a feature table and the metadata, one sample per subject.',
                 '''\
from moonstone import examples
from moonstone.analysis.compositional_pipeline import filter_features, read_feature_table
from moonstone.parsers import MetadataParser

species, notes = read_feature_table(examples.path('species_counts.tsv'), 'species')
metadata = MetadataParser(examples.path('metadata.csv'), sep=',', index_col='sample').dataframe
''' + ONE_SAMPLE_PER_PATIENT + '''\
species, note = filter_features(species[inclusion.index], min_prevalence=0.1)
inclusion[['disease', 'age', 'antibiotics']].head()
'''),
            Step('**Run it.** One row per feature and tested term:\n\n'
                 '- `level` and `reference`: for a group variable, the group compared and the group it is compared '
                 'with. Set the reference yourself, otherwise it is the first value in sorted order.\n'
                 '- `coef`: the difference on the CLR scale (per unit for a numeric variable), with its confidence '
                 'interval `ci_low`, `ci_high`.\n'
                 '- `q_value`: corrected over all tests. Use this one.\n'
                 '- `n_samples`, `n_dropped`: samples used, and left out for a missing value.\n\n'
                 'Covariates are adjusted for but not reported.',
                 '''\
from moonstone.analysis.association import feature_association

results = feature_association(species, inclusion, ['disease'], covariates=['age', 'antibiotics'],
                              reference={'disease': 'healthy'})
results[['feature', 'variable', 'level', 'reference', 'coef', 'ci_low', 'ci_high', 'p_value', 'q_value',
         'n_samples']].sort_values('p_value')
'''),
            Step('**A trap of the CLR scale.** Each CLR value is relative to the geometric mean of its sample. '
                 'When one dominant taxon drops sharply (here *Lactobacillus crispatus* with bacterial vaginosis), '
                 'the others shift on the CLR scale even if their own abundance did not change. In this '
                 'simulation *Escherichia coli* and *Lactobacillus iners* were given **no** change, yet the '
                 'model can report them.\n\n`ancombc_check` runs ANCOM-BC, which corrects for this bias, on the '
                 'same model. Its log-fold changes are not on the CLR scale: compare directions and significance, '
                 'not values. **Associations found by both methods are the robust ones.**',
                 '''\
from moonstone.analysis.association import ancombc_check

check = ancombc_check(species, inclusion, ['disease'], covariates=['age', 'antibiotics'],
                      reference={'disease': 'healthy'})
both = results.merge(check[['feature', 'variable', 'level', 'q_value', 'log_fold_change']],
                     on=['feature', 'variable', 'level'], suffixes=('', '_ancombc'))
both[['feature', 'coef', 'q_value', 'log_fold_change', 'q_value_ancombc']].sort_values('q_value')
'''),
            Step('**From the command line or the interface**\n\n```bash\nmoonstone species_counts.tsv metadata.csv '
                 'results --associate disease \\\n    --covariates age antibiotics --reference disease=healthy '
                 '--ancombc \\\n    --select Type=inclusion --group_column patient\n```\n\n'
                 'Outputs: `association.csv` (with the ANCOM-BC columns when `--ancombc` is given), `ancombc.csv` '
                 'and `methods.txt`.\n\n**Good to know**: a numeric column of group codes (1/2/3) is treated as a '
                 'number unless listed after `--categorical`; `--standardize` gives coefficients per standard '
                 'deviation instead of per unit; features present in fewer than 10% of the samples are left out '
                 '(`--min_prevalence`).'),
        ), tool='associate'),
    Tutorial(
        'differential', 'Differential analysis with simple tests', 'Analyse',
        'DifferentialAnalysis: one classical test per taxon between the groups of a variable.',
        (
            Step('`DifferentialAnalysis` runs one classical test per taxon between the groups of a metadata '
                 'variable, then corrects the p-values. It does not transform the counts and does not adjust for '
                 'other variables: for that, use *Association*. Give it comparable values (normalized counts).\n\n'
                 'Both tables must have the **samples as columns**, in the same order, so the metadata is turned '
                 '(`.T`).',
                 SETUP + ONE_SAMPLE_PER_PATIENT + '''\
from moonstone.analysis import DifferentialAnalysis
from moonstone.normalization.counts.total_counts import TotalCountsNormalization

normalized = TotalCountsNormalization(counts[inclusion.index]).normalized_df
species = normalized.droplevel(['kingdom', 'phylum', 'class', 'order', 'family', 'genus'])
analysis = DifferentialAnalysis(inclusion[['disease', 'batch']].T, species)
species.iloc[:, :5]
'''),
            Step('**Two groups**: `"dichotomic_features"` with `t_test` or `wilcoxon_rank_test`. The last argument '
                 'is the correction (`fdr_bh`, `bonferroni` ...).',
                 "analysis.differential_analysis_by_feature(['disease'], 'dichotomic_features', "
                 "'wilcoxon_rank_test', 'fdr_bh')"),
            Step('**More than two groups**: `"multiple_features"` with `kruskal_test` or `one_way_anova`. Here '
                 '`batch`, the sequencing run, which was simulated with no effect.',
                 "analysis.differential_analysis_by_feature(['batch'], 'multiple_features', 'kruskal_test', "
                 "'fdr_bh')"),
        )),
]}
