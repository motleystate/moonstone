import argparse
import os
import sys
import logging

from moonstone.parsers import (handleCounts, handleMetadata)
from moonstone.analysis import (classify, clustering, randomForest, stats)
from moonstone.filtering.mean_filtering import MeanFiltering


""" Moonstone is designed as a means of enabling machine learning analysis of metagenomic sequencing data.
The program englobes the following functionalities:
    > Import of count data (currently limited to SHAMAN-style counts output)
    > Import of clinical metadata (in development)
    > A basic statistical analysis of the imported data
    > Unsupervised clustering algorithms: PCA (other is development)
    > Supervised clustering (in development)
    > Classifier and biomarker discovery (in development)"""


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
    parser.add_argument('-p', '--pca_plot', help='Generates PCA plot of data', action='store_true')
    parser.add_argument('-k', '--k_means', metavar='nb_clusters', type=int,
                        help="Runs K-Means clustering. Provide number of clusters")
    parser.add_argument('-s', '--svm', metavar='variable', type=str, help='Run SVMachine Classification')
    parser.add_argument('-sr', '--svm_roc', metavar='variable', type=str, help='SVM Classifier with ROC Analysis')
    parser.add_argument('-sc', '--svm_classifier', metavar='variable', type=str,
                        help='SVM and output Classifier Components')
    parser.add_argument('-rf', '--random_forest', metavar='variable', type=str,
                        help="Random Forest Analysis, using supplied variable")
    try:
        return parser.parse_args()
    except SystemExit:
        sys.exit(1)


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
    filtering = MeanFiltering(count_df, threshold=min_read_mean)
    count_df = filtering.filtered_df
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

    # Parse input files
    logger.info('Starting handleCounts module, get_count_file function to import count data.')
    count_df = get_count_file(args.countfile, outdir)  # df
    if args.filtering:
        count_df = filter_count_df(count_df, args.filtering, outdir)
    metadata_df = get_metadata_file(args.metadata, count_df, args.outdir)  # dm

    # Run different analysis
    if args.pca_plot:
        pca = clustering.Unsupervised(count_df, metadata_df, outdir)
        pca.pca()

    if args.k_means:
        kmeans = clustering.Unsupervised(count_df, metadata_df, outdir)
        kmeans.kmeans('metaData_withKClusters.csv', n_clusters=args.k_means)

    if args.svm:
        svm = classify.ML(count_df, metadata_df, outdir)
        svm.svm(variable=args.svm)

    if args.svm_roc:
        roc = classify.ML(count_df, metadata_df, outdir)
        roc.roc_analysis(variable=args.svm_roc)

    if args.svm_classifier:
        scomponents = classify.ML(count_df, metadata_df, outdir)
        scomponents.feature_analysis(variable=args.svm_classifier)

    if args.random_forest:
        forest = randomForest.RandomForest(count_df, metadata_df, outdir, variable=args.random_forest)
        forest.forest('rf_Allfeatures.csv')


if __name__ == "__main__":
    run()
