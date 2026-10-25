import argparse
import os
import sys
import logging

from moonstone.parsers import (handleCounts, handleMetadata)
from moonstone.analysis import stats
from moonstone.analysis.classification import Classification, KERNELS, MODELS
from moonstone.analysis.diversity.beta import DEFAULT_PERMUTATIONS, MISSING_WARNING_FRACTION
from moonstone.analysis.permanova import TAXA_LEVELS, read_taxonomy_counts, run_permanova
from moonstone.filtering.mean_filtering import MeanFiltering


""" Moonstone is designed as a means of enabling machine learning analysis of metagenomic sequencing data.
The program englobes the following functionalities:
    > Import of count data (currently limited to SHAMAN-style counts output)
    > Import of clinical metadata (in development)
    > A basic statistical analysis of the imported data
    > Classification of a metadata variable from the counts
    > PERMANOVA with a dispersion check"""


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
                        help='Random seed for folds and models, and for PERMANOVA permutations (default: 0)')
    parser.add_argument('--group_column', metavar='column', type=str,
                        help='Metadata column naming the subject of each sample; keeps a subject in one fold')
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
    parser.add_argument('--permutations', type=int, default=DEFAULT_PERMUTATIONS,
                        help=f'With --permanova, number of permutations (default: {DEFAULT_PERMUTATIONS})')
    parser.add_argument('--threads', type=int, default=0,
                        help='With --permanova, number of tests run in parallel, one per processor core; 0 uses '
                             'every available core. Results do not depend on it (default: 0)')
    for flags, replacement in [(['-s', '--svm'], '--classify {variable}'),
                               (['-sr', '--svm_roc'], '--classify {variable} --roc'),
                               (['-sc', '--svm_classifier'], '--classify {variable} --feature_selection'),
                               (['-rf', '--random_forest'], '--classify {variable} --model random_forest')]:
        parser.add_argument(*flags, action=RemovedFlag, replacement=replacement, metavar='variable')
    parser.add_argument('-p', '--pca_plot', action=RemovedFlag, nargs=0)
    parser.add_argument('-k', '--k_means', action=RemovedFlag, metavar='nb_clusters')
    try:
        args = parser.parse_args()
        if args.permanova is not None:
            others = [flag for flag, used in [('-f/--filtering', args.filtering), ('-c/--classify', args.classify)]
                      if used]
            if others:
                parser.error(f"--permanova reads a count table with taxonomy and cannot be combined with "
                             f"{', '.join(others)}. Run them separately.")
        return args
    except SystemExit as exit_error:
        # --help exits with 0; argument errors exit with 1.
        sys.exit(0 if exit_error.code == 0 else 1)


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
    if count_df.shape[0] != metadata_df.shape[0]:
        print('Count samples: {} does not equal sample number in metadata file {}'.format(count_df.shape[0],
                                                                                          metadata_df.shape[0]))
        raise SystemExit("Count and Metadata sample numbers don't match...disaster!")
    else:
        print(f'Nice! {count_df.shape[0]} samples in count AND metadata files.\n')
    return metadata_df


def permanova(args, outdir, logger):
    """
    Run PERMANOVA from a normalized count table with taxonomy and write permanova.csv to the output directory.
    """
    print(f"\nOpening {args.countfile} which should contain normalized counts with taxonomy columns\n")
    counts_df = read_taxonomy_counts(args.countfile)
    print(f'Found {counts_df.shape[1]} samples and {counts_df.shape[0]} taxa')
    metadata_df = handleMetadata.Inputs(args.metadata).openmeta(outdir)
    n_common = len(set(counts_df.columns) & set(metadata_df.index))
    print(f'{n_common} samples in both the count and metadata files.\n')
    try:
        results = run_permanova(
            counts_df, metadata_df, columns=args.permanova or None, taxa_level=args.taxa_level,
            read_threshold=not args.no_read_threshold, permutations=args.permutations, seed=args.seed,
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

    # Parse input files
    logger.info('Starting handleCounts module, get_count_file function to import count data.')
    count_df = get_count_file(args.countfile, outdir)  # df
    if args.filtering:
        count_df = filter_count_df(count_df, args.filtering, outdir)
    metadata_df = get_metadata_file(args.metadata, count_df, args.outdir)  # dm

    # Run different analysis
    if args.classify:
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


if __name__ == "__main__":
    run()
