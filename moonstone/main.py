import argparse
import os
import sys
import logging

from moonstone.parsers import (handleCounts, handleMetadata)
from moonstone.analysis import stats
from moonstone.analysis.association import ancombc_check, feature_association
from moonstone.analysis import compositional_pipeline as pipeline
from moonstone.analysis.correlation import METHODS as CORRELATION_METHODS
from moonstone.analysis.correlation import edges_to_matrix, feature_correlation, filter_edges, to_gephi
from moonstone.analysis.classification import Classification, KERNELS, MODELS
from moonstone.analysis.diversity.beta import DEFAULT_PERMUTATIONS, MISSING_WARNING_FRACTION
from moonstone.analysis.permanova import TAXA_LEVELS, read_taxonomy_counts, run_permanova
from moonstone.filtering.mean_filtering import MeanFiltering
from moonstone.validation.sample_matching import match_samples


""" Moonstone is designed as a means of enabling machine learning analysis of metagenomic sequencing data.
The program englobes the following functionalities:
    > Import of count data (currently limited to SHAMAN-style counts output)
    > Import of clinical metadata (in development)
    > A basic statistical analysis of the imported data
    > Classification of a metadata variable from the counts
    > PERMANOVA with a dispersion check
    > Compositional correlation between features, and feature-metadata association with covariates"""


class RemovedFlag(argparse.Action):
    """A removed flag: fail and name its replacement, if there is one."""

    def __init__(self, option_strings, dest, replacement=None, **kwargs):
        self.replacement = replacement
        super().__init__(option_strings, dest, help=argparse.SUPPRESS, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        if self.replacement is None:
            parser.error(f'{option_string} has been removed.')
        parser.error(f'{option_string} has been removed. Use: {self.replacement.format(variable=values)}')


def parse_arguments():
    """
    Parser for commandline arguments:
    * Countfile and Metadata are required
    * Optional parameters include filtering, and the implemented analysis
    """
    parser = argparse.ArgumentParser(prog='moonstone', description='Microbiota Analysis Scripts for Machine LEarning',
                                     usage='%(prog)s countfile, metadata file, output directory [options]')
    parser.add_argument('countfile', type=str, help="Normalized count file input")
    parser.add_argument('metadata', type=str, help="Clinical data input file")
    parser.add_argument('outdir', nargs='?', type=str, help="Output files directory", default='output')
    parser.add_argument('-os', '--skip_prompt', help='Suppress output directory exists prompt.', action='store_true')
    parser.add_argument('-f', '--filtering', metavar='min_mean_read', type=float,
                        help="Minimum mean reads per variable: use a float >0")
    parser.add_argument('-c', '--classify', metavar='variable', type=str,
                        help='Cross-validated classification of a metadata variable from the counts')
    parser.add_argument('--model', nargs='+', choices=MODELS, default=['svm'],
                        help='Classifier(s) for --classify: Support Vector Machine and/or random forest '
                             '(default: svm)')
    parser.add_argument('--kernel', choices=KERNELS + ('auto',), default='linear',
                        help='Support Vector Machine kernel; auto picks one by nested cross-validation '
                             '(default: linear)')
    parser.add_argument('--explore_kernels', action='store_true',
                        help='With --classify, compare Support Vector Machine kernels (exploratory scores)')
    parser.add_argument('--folds', type=int, default=5, help='Number of cross-validation folds (default: 5)')
    parser.add_argument('--seed', type=int, default=0,
                        help='Random seed for folds and models, and for PERMANOVA and correlation permutations '
                             '(default: 0)')
    parser.add_argument('--group_column', metavar='column', type=str,
                        help='Metadata column naming the subject of each sample; keeps a subject in one fold. With '
                             '--correlation or --associate, the run stops if a subject has more than one sample')
    parser.add_argument('--roc', action='store_true',
                        help='With --classify, plot Receiver Operating Characteristic curves')
    parser.add_argument('--feature_selection', action='store_true',
                        help='With --classify, select features by nested recursive feature elimination')
    parser.add_argument('--permanova', metavar='variable', nargs='*',
                        help='PERMANOVA with a dispersion check (PERMDISP) for the given metadata variables, or all '
                             'variables if none are given. The count file must then be a normalized table with '
                             'taxonomy columns (kingdom ... species). Place this flag after the file arguments')
    parser.add_argument('--taxa_level', choices=TAXA_LEVELS + ['all'], default='species',
                        help='With --permanova, taxonomic level to sum the counts to, or all (default: species)')
    parser.add_argument('--no_read_threshold', action='store_true',
                        help='With --permanova, do not filter taxa by the read-threshold knee at each level')
    parser.add_argument('--permutations', type=int,
                        help=f'Number of permutations for --permanova (default: {DEFAULT_PERMUTATIONS}) and '
                             f'--correlation (default: {pipeline.DEFAULT_CORRELATION_PERMUTATIONS})')
    parser.add_argument('--threads', type=int, default=0,
                        help='With --permanova, number of tests run in parallel, one per processor core; 0 uses '
                             'every available core. Results do not depend on it (default: 0)')
    compositional = parser.add_argument_group(
        'compositional correlation and association',
        'Count files have one row per feature and one column per sample: with taxonomy columns (kingdom ... '
        'species), counts are summed to --taxa_level; otherwise the first column names the features. Features are '
        'filtered, then CLR-transformed. A methods paragraph is written to methods.txt')
    compositional.add_argument('--correlation', action='store_true',
                               help='Correlation between every pair of features, with permutation p-values')
    compositional.add_argument('--correlation_method', choices=CORRELATION_METHODS, default='spearman',
                               help='Correlation method (default: spearman, on CLR values)')
    compositional.add_argument('--second_countfile', metavar='file', type=str,
                               help='With --correlation, a second table (for example resistance genes): only pairs '
                                    'across the two tables are tested, each table transformed on its own')
    compositional.add_argument('--associate', metavar='variable', nargs='+',
                               help='Association of every feature with these metadata variables (one linear model '
                                    'per feature on CLR values). Place this flag after the file arguments')
    compositional.add_argument('--covariates', metavar='variable', nargs='+', default=[],
                               help='With --associate, metadata variables to adjust for')
    compositional.add_argument('--categorical', metavar='variable', nargs='+', default=[],
                               help='With --associate, numeric variables to treat as groups (for example 1/2/3 codes)')
    compositional.add_argument('--reference', metavar='VARIABLE=LEVEL', action='append', default=[],
                               help='With --associate, reference level of a categorical variable; repeatable')
    compositional.add_argument('--standardize', action='store_true',
                               help='With --associate, scale continuous variables to unit standard deviation')
    compositional.add_argument('--ancombc', action='store_true',
                               help='With --associate, also run ANCOM-BC on the same model as a sensitivity check')
    compositional.add_argument('--select', metavar='COLUMN=VALUE', action='append', default=[],
                               help='Keep only the samples with this metadata value, for example Type=inclusion; '
                                    'repeatable, all must match')
    compositional.add_argument('--min_prevalence', type=float, default=pipeline.DEFAULT_MIN_PREVALENCE,
                               help='Keep features present in at least this fraction of the selected samples '
                                    f'(default: {pipeline.DEFAULT_MIN_PREVALENCE})')
    compositional.add_argument('--top', type=int,
                               help='Then keep only the N features with the highest mean relative abundance')
    compositional.add_argument('--correction', choices=['fdr_bh', 'bonferroni', 'none'], default='fdr_bh',
                               help='Multiple-testing correction (default: fdr_bh, Benjamini-Hochberg)')
    compositional.add_argument('--alpha', type=float, default=0.05,
                               help='Significance threshold on corrected p-values (default: 0.05)')
    for flags, replacement in [(['-s', '--svm'], '--classify {variable}'),
                               (['-sr', '--svm_roc'], '--classify {variable} --roc'),
                               (['-sc', '--svm_classifier'], '--classify {variable} --feature_selection'),
                               (['-rf', '--random_forest'], '--classify {variable} --model random_forest')]:
        parser.add_argument(*flags, action=RemovedFlag, replacement=replacement, metavar='variable')
    parser.add_argument('-p', '--pca_plot', action=RemovedFlag, nargs=0)
    parser.add_argument('-k', '--k_means', action=RemovedFlag, metavar='nb_clusters')
    try:
        args = parser.parse_args()
        check_analysis_combinations(parser, args)
        return args
    except SystemExit as exit_error:
        # --help exits with 0; argument errors exit with 1.
        sys.exit(0 if exit_error.code == 0 else 1)


def check_analysis_combinations(parser, args):
    """Stop on flags that cannot be combined, before any file is read."""
    compositional = args.correlation or args.associate is not None
    for flag, used in [('--permanova', args.permanova is not None), ('--correlation or --associate', compositional)]:
        if not used:
            continue
        others = [other for other, other_used in [('-f/--filtering', args.filtering), ('-c/--classify', args.classify)]
                  if other_used]
        if others:
            parser.error(f"{flag} reads a count table with one row per feature and cannot be combined with "
                         f"{', '.join(others)}. Run them separately.")
    if args.permanova is not None and compositional:
        parser.error("--permanova cannot be combined with --correlation or --associate. Run them separately.")
    if args.permutations is not None and args.permutations < 1:
        parser.error("--permutations must be at least 1.")
    if not 0 < args.alpha < 1:
        parser.error("--alpha must be between 0 and 1.")
    needs = [('--second_countfile', args.second_countfile, '--correlation', args.correlation)]
    needs += [(flag, used, '--associate', args.associate is not None) for flag, used in [
        ('--covariates', args.covariates), ('--categorical', args.categorical), ('--reference', args.reference),
        ('--standardize', args.standardize), ('--ancombc', args.ancombc)]]
    needs += [(flag, used, '--correlation or --associate', compositional) for flag, used in [
        ('--select', args.select), ('--top', args.top is not None)]]
    for flag, used, required, present in needs:
        if used and not present:
            parser.error(f"{flag} is only used with {required}.")
    if args.correlation_method != 'spearman' and not args.correlation:
        parser.error("--correlation_method is only used with --correlation.")


def handle_output_directory(outdir_path, suppress_outdir_prompt):
    """
    Check if the output directory exists.
    Unless suppressed with -os argument, user will need to confirm overwriting to continue
    :return: path to output directory
    """
    yes_answers = ['Y', 'YES', 'MAKE IT SO', 'OUI']
    if os.path.isdir(outdir_path):
        if suppress_outdir_prompt:
            print('Output directory %s already exists, but away we go anyway!' % outdir_path)
        else:
            answer = input('Output directory %s already exists! Continue anyway [y/N]?' % outdir_path).upper()
            if answer in yes_answers:
                print('Okay. Files in %s will be overwritten.' % outdir_path)
            else:
                raise SystemExit("Not overwriting...aborting.")

    # Create directory and report success.
    if not os.path.isdir(outdir_path):
        try:
            os.mkdir(outdir_path)
        except OSError:
            print('There was a problem creating %s' % outdir_path)
        else:
            print('Output directory %s successfully created.' % outdir_path)
    return outdir_path


def get_count_file(countfile_path, outdir_path):
    """
    The following section concerns opening and reporting on the count file.
    :return: dataframe of counts
    """
    print("\nOpening {} which should contain normalized counts\n".format(countfile_path))
    count_read = handleCounts.Inputs(countfile_path)
    count_df = count_read.opencounts()
    num_samples, num_otus = count_df.shape
    print(f'Found {num_samples} samples and {num_otus} OTUs\n')

    run_statistics = stats.Descriptive(count_df, outdir_path)
    run_statistics.matrix_stats('starting_variables.csv')
    check_sparse = stats.Density(count_df)
    if check_sparse.is_sparse():
        print("Count Table is sparse with %2.3f%s non-zero values." % (check_sparse.percent_non_zeros(), "%"))
    else:
        print("Count Table has %2.3f%s non-zero values." % (check_sparse.percent_non_zeros(), "%"))
    return count_df


def filter_count_df(count_df, min_read_mean, outdir_path):
    """
    :return: filtered count dataframe
    """
    check_sparse = stats.Density(count_df)
    # MeanFiltering expects one row per taxon; the count table here has one row per sample.
    filtering = MeanFiltering(count_df.T, threshold=min_read_mean)
    count_df = filtering.filtered_df.T
    filtered_stats = stats.Descriptive(count_df, outdir_path)
    filtered_stats.matrix_stats('filtered_variables.csv')
    count_df.to_csv(path_or_buf=outdir_path+'/'+'filteredCountFile.csv')
    check_f_sparse = stats.Density(count_df)
    if check_f_sparse.is_sparse():
        print("Filtered Count Table is sparse with %2.3f%s non-zero values.\n" %
              (check_f_sparse.percent_non_zeros(), "%"))
    elif check_sparse.is_sparse():
        print("Filtered Count Table is no longer sparse: %2.3f%s non-zero values.\n" %
              (check_f_sparse.percent_non_zeros(), "%"))
    else:
        print("Filtered Count Table now has %2.3f%s non-zero reads.\n" % (check_f_sparse.percent_non_zeros(), "%"))
    return count_df


def get_metadata_file(metadata_path, count_df, outdir):
    """
    This step is to read in the variables file.
    This could be rendered optional depending the the analyses selected
    It also deal with reporting on the metadata file
    :return: dataframe of metadata
    """
    print("\nOpening {} which is expected to contain clinical metadata.".format(metadata_path))
    metadata_read = handleMetadata.Inputs(metadata_path)
    metadata_df = metadata_read.openmeta(outdir)
    report_sample_matching(count_df.index, metadata_df.index)
    return metadata_df


def report_sample_matching(count_samples, metadata_samples):
    """
    Match the samples of the count table and the metadata by ID, print what does not match, and stop
    if the files cannot be paired (no ID in common, or an ID used twice).
    """
    matching = match_samples(count_samples, metadata_samples)
    findings = matching.findings()
    for finding in findings:
        if finding.level != 'error':
            print(f'{finding.level.capitalize()}: {finding.message}.')
    errors = [finding.message for finding in findings if finding.level == 'error']
    if errors:
        raise SystemExit('Count and metadata samples cannot be matched by ID: ' + '; '.join(errors) + '.')
    print(f'{len(matching.matched)} samples in both the count and metadata files.\n')
    return matching


def permanova(args, outdir, logger):
    """
    Run PERMANOVA from a normalized count table with taxonomy and write permanova.csv to the output directory.
    """
    print(f"\nOpening {args.countfile} which should contain normalized counts with taxonomy columns\n")
    counts_df = read_taxonomy_counts(args.countfile)
    print(f'Found {counts_df.shape[1]} samples and {counts_df.shape[0]} taxa')
    metadata_df = handleMetadata.Inputs(args.metadata).openmeta(outdir)
    report_sample_matching(counts_df.columns, metadata_df.index)
    try:
        results = run_permanova(
            counts_df, metadata_df, columns=args.permanova or None, taxa_level=args.taxa_level,
            read_threshold=not args.no_read_threshold, permutations=args.permutations or DEFAULT_PERMUTATIONS, seed=args.seed,
            n_jobs=args.threads,
        )
    except ValueError as error:
        logger.error(str(error))
        raise SystemExit(f'PERMANOVA failed: {error}')
    output_file = os.path.join(outdir, 'permanova.csv')
    results.to_csv(output_file, index=False)
    for (level, warnings) in results[['level', 'warnings']].drop_duplicates().itertuples(index=False):
        if warnings:
            print(f'Warning ({level}): {warnings}')
    missing = results.groupby('variable', sort=False)[['n_samples', 'n_missing']].max()
    for variable, row in missing[missing['n_missing'] > 0].iterrows():
        share = row.n_missing / (row.n_missing + row.n_samples)
        check = ' - mostly missing, check and clean this variable' if share > MISSING_WARNING_FRACTION else ''
        print(f'{variable}: {int(row.n_missing)} sample(s) with a missing value left out{check}')
    significant = results[results['p_value_corrected'] < 0.05]
    print(f'\n{significant.shape[0]} of {results["p_value"].notna().sum()} tests significant after '
          f'false discovery rate correction (p < 0.05).')
    for row in significant.itertuples(index=False):
        dispersion = ' (groups also differ in spread)' if row.permdisp_p_value_corrected < 0.05 else ''
        print(f'  {row.level} / {row.variable}: R2 = {row.R2:.3f}, corrected p = {row.p_value_corrected:.3g}'
              f'{dispersion}')
    print(f'Results written to {output_file}')


def read_compositional_inputs(args, outdir):
    """
    Read the metadata and count table(s), select samples and filter features, in the order fixed by
    :mod:`moonstone.analysis.compositional_pipeline`.

    Returns:
        The metadata of the selected samples, the selection notes, and one (name, table, notes, filter_note) tuple
        per count table.
    """
    metadata = handleMetadata.Inputs(args.metadata).openmeta(outdir)
    paths = [('count table', args.countfile)]
    if args.correlation and args.second_countfile:
        paths = [('first count table', args.countfile), ('second count table', args.second_countfile)]
    raw_tables = []
    for name, path in paths:
        print(f"\nOpening {path}, one row per feature and one column per sample")
        table, notes = pipeline.read_feature_table(path, args.taxa_level)
        print(f'Found {table.shape[1]} samples and {table.shape[0]} features')
        report_sample_matching(table.columns, metadata.index)
        raw_tables.append((name, table, notes))
    metadata, selection_notes = pipeline.select_samples(
        metadata, [pipeline.parse_selection(text) for text in args.select])
    for note in selection_notes:
        print(f'Selected {note}')
    selected = set(metadata.index)
    samples = [s for s in raw_tables[0][1].columns if s in selected]
    for _, table, _ in raw_tables[1:]:
        samples = [s for s in samples if s in set(table.columns)]
    if len(samples) < 3:
        raise ValueError(f"Only {len(samples)} sample(s) are in every count table and the selected metadata.")
    print(f'{len(samples)} samples analysed.')
    tables = []
    for name, table, notes in raw_tables:
        filtered, filter_note = pipeline.filter_features(table[samples], args.min_prevalence, args.top)
        print(f'{name.capitalize()}: {filter_note}')
        tables.append((name, filtered, notes, filter_note))
    return metadata.loc[samples], selection_notes, tables


def correlation(args, outdir, metadata, selection_notes, tables):
    """Feature correlation: write the edge table, the significant edges for Gephi and the r matrix."""
    correction = None if args.correction == 'none' else args.correction
    permutations = args.permutations or pipeline.DEFAULT_CORRELATION_PERMUTATIONS
    edges = feature_correlation(
        tables[0][1], tables[1][1] if len(tables) > 1 else None, method=args.correlation_method,
        correction_method=correction, permutations=permutations, seed=args.seed,
        metadata=metadata if args.group_column else None, group_column=args.group_column,
    )
    edges.to_csv(os.path.join(outdir, 'correlation_edges.csv'), index=False)
    edges_to_matrix(edges).to_csv(os.path.join(outdir, 'correlation_matrix.csv'))
    significant = filter_edges(edges, alpha=args.alpha)
    to_gephi(significant).to_csv(os.path.join(outdir, 'correlation_gephi_edges.csv'), index=False)
    positive = int((significant['r'] > 0).sum())
    print(f'\n{len(significant)} of {edges["p_value"].notna().sum()} pairs significant '
          f'({"q" if correction else "p"} < {args.alpha}): {positive} positive, {len(significant) - positive} '
          'negative.')
    for row in significant.sort_values('p_value').head(10).itertuples(index=False):
        print(f'  {row.feature_a} - {row.feature_b}: r = {row.r:.3f}, q = {row.q_value:.3g}')
    print(f'Results written to {outdir}/correlation_edges.csv, correlation_matrix.csv and '
          'correlation_gephi_edges.csv')
    return pipeline.correlation_methods(
        args.correlation_method, len(metadata), [(name, notes, note) for name, _, notes, note in tables],
        selection_notes, permutations, args.seed, correction, args.alpha, args.group_column,
    )


def association(args, outdir, metadata, selection_notes, tables):
    """Feature-metadata association: write association.csv (and ancombc.csv)."""
    correction = None if args.correction == 'none' else args.correction
    reference = dict(pipeline.parse_selection(text) for text in args.reference)
    options = dict(covariates=args.covariates, categorical=args.categorical, reference=reference,
                   standardize=args.standardize, correction_method=correction, group_column=args.group_column)
    name, table, notes, filter_note = tables[0]
    results = feature_association(table, metadata, args.associate, **options)
    if args.ancombc:
        check = ancombc_check(table, metadata, args.associate, **options)
        check.to_csv(os.path.join(outdir, 'ancombc.csv'), index=False)
        keys = ['feature', 'variable', 'level']
        merged = check[keys + ['q_value', 'log_fold_change']].rename(
            columns={'q_value': 'ancombc_q_value', 'log_fold_change': 'ancombc_log_fold_change'})
        results = results.merge(merged, on=keys, how='left')
    results.to_csv(os.path.join(outdir, 'association.csv'), index=False)
    column = 'q_value' if correction else 'p_value'
    significant = results[results[column] < args.alpha]
    print(f'\n{len(significant)} of {results["p_value"].notna().sum()} tests significant '
          f'({"q" if correction else "p"} < {args.alpha}), on {results["n_samples"].iloc[0]} samples'
          + (f' ({results["n_dropped"].iloc[0]} left out for missing metadata).' if results["n_dropped"].iloc[0]
             else '.'))
    if args.ancombc:
        both = significant[significant['ancombc_q_value' if correction else column] < args.alpha]
        print(f'{len(both)} of them also significant with ANCOM-BC.')
    for row in significant.sort_values('p_value').head(10).itertuples(index=False):
        level = f' {row.level} vs {row.reference}' if isinstance(row.level, str) else ''
        print(f'  {row.feature} / {row.variable}{level}: coef = {row.coef:.3f} '
              f'[{row.ci_low:.3f}, {row.ci_high:.3f}], q = {row.q_value:.3g}')
    print(f'Results written to {outdir}/association.csv' + (' and ancombc.csv' if args.ancombc else ''))
    return pipeline.association_methods(
        args.associate, args.covariates, int(results['n_samples'].iloc[0]), int(results['n_dropped'].iloc[0]),
        ("count table", notes, filter_note), selection_notes, correction, args.alpha, args.standardize, args.group_column,
        args.ancombc, results.dropna(subset=['reference']).groupby('variable', sort=False)['reference'].first().to_dict(),
    )


def compositional(args, outdir, logger):
    """Run --correlation and/or --associate, then write methods.txt."""
    try:
        metadata, selection_notes, tables = read_compositional_inputs(args, outdir)
        paragraphs = []
        if args.correlation:
            paragraphs.append(correlation(args, outdir, metadata, selection_notes, tables))
        if args.associate is not None:
            paragraphs.append(association(args, outdir, metadata, selection_notes, tables))
    except ValueError as error:
        logger.error(str(error))
        raise SystemExit(f'Analysis failed: {error}')
    methods_file = os.path.join(outdir, 'methods.txt')
    with open(methods_file, 'w') as handle:
        handle.write('\n\n'.join(paragraphs) + '\n')
    print(f'Methods paragraph written to {methods_file}')


def classify(args, count_df, metadata_df, outdir, logger):
    """Cross-validated classification of args.classify, for each model."""
    if args.explore_kernels:
        exploration = Classification(count_df, metadata_df, args.classify, outdir, model='svm',
                                     n_folds=args.folds, seed=args.seed, group_column=args.group_column)
        try:
            exploration.explore_kernels()
        except ValueError as error:
            logger.error(str(error))
            raise SystemExit(f'Kernel exploration failed: {error}')
    for model in args.model:
        classification = Classification(count_df, metadata_df, args.classify, outdir, model=model,
                                        kernel=args.kernel, n_folds=args.folds, seed=args.seed,
                                        group_column=args.group_column)
        try:
            classification.run(roc=args.roc, feature_selection=args.feature_selection)
        except ValueError as error:
            logger.error(str(error))
            raise SystemExit(f'Classification failed: {error}')


def run():
    args = parse_arguments()
    outdir = handle_output_directory(args.outdir, args.skip_prompt)

    logging.basicConfig(level=logging.DEBUG, filename=outdir+'/moonstone.log', filemode='w',
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logger = logging.getLogger('moonstone_main')
    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.WARNING)
    error_log = logging.StreamHandler()
    error_log.setLevel(logging.WARNING)
    logger.addHandler(error_log)

    if args.permanova is not None:
        permanova(args, outdir, logger)
        return
    if args.correlation or args.associate is not None:
        compositional(args, outdir, logger)
        return

    # Parse input files
    logger.info('Starting handleCounts module, get_count_file function to import count data.')
    count_df = get_count_file(args.countfile, outdir)  # df
    if args.filtering:
        count_df = filter_count_df(count_df, args.filtering, outdir)
    metadata_df = get_metadata_file(args.metadata, count_df, args.outdir)  # dm

    if args.classify:
        classify(args, count_df, metadata_df, outdir, logger)


if __name__ == "__main__":
    run()
