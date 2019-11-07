import argparse
import handleCounts
import handleMetadata
import clustering
import stats
import filtering
import classify
from moonstone.analysis import randomForest

""" Moonstone is designed as a means of enabling machine learning analysis of metagenomic sequencing data.
The program englobes the following functionalities:
    > Import of count data (currently limited to SHAMAN-style counts output)
    > Import of clinical metadata (in development)
    > A basic statistical analysis of the imported data 
    > Unsupervised clustering algorithms: PCA (other is development)
    > Supervised clustering (in development)
    > Classifier and biomarker discovery (in development)"""

# Parser for commandline arguments:
# Countfile and Metadata are required
# Optional parameters include filtering, and the implemented analysis

parser = argparse.ArgumentParser(description='Microbiota Analysis Scripts for Machine LEarning')
parser.add_argument("countfile", type=str, help="Normalized count file input")
parser.add_argument("metadata", type=str, help="Clinical data input file")
parser.add_argument("-f", metavar='filtering', type=float, help="Minimum mean reads per variable: use a float >0")
parser.add_argument('-p', help='Generates PCA plot of data', action='store_true')
parser.add_argument('-k', metavar='clusters', help="Runs K-Means clustering. Provide number of clusters", type=int)
parser.add_argument('-s', metavar='variable', type=str, help='Run SVMachine Classification')
parser.add_argument('-sr', metavar='variable', type=str, help='SVM Classifier with ROC Analysis')
parser.add_argument('-sc', metavar='variable', type=str, help='SVM and output Classifier Components')
parser.add_argument('-rf', metavar='variable', type=str, help="Random Forest Analysis, using supplied variable")
args = parser.parse_args()

# The following section concerns opening and reporting on the count file
print("Opening {} which should contain normalized counts\n".format(args.countfile))
count_read = handleCounts.Inputs(args.countfile)
df = count_read.opencounts()
num_samples, num_otus = df.shape
print(f'Found {num_samples} samples and {num_otus} OTUs\n')
run_statistics = stats.Descriptive(df, 'starting_variables.csv')
run_statistics.matrix_stats()
check_sparse = stats.Density(df)
if check_sparse.is_sparse():
    print("Count Table is sparse with %2.3f%s non-zero values." % (check_sparse.percent_non_zeros(), "%"))
else:
    print("Count Table has %2.3f%s non-zero values." % (check_sparse.percent_non_zeros(), "%"))

# The second step is to read in the variables file.
# This could be rendered optional depending the the analyses selected
metadata_read = handleMetadata.Inputs(args.metadata)
dm = metadata_read.openmeta()

# The following section concerns opening and reporting on the metadata file
print("\nOpening {} which is expected to contain clinical metadata.".format(args.metadata))
meta_read = handleMetadata.Inputs(args.metadata)
dfc = meta_read.openmeta()
if df.shape[0] != dfc.shape[0]:
    print('Count samples: {} does not equal sample number in metadata file {}'.format(df.shape[0], dfc.shape[0]))
    raise SystemExit("Count and Metadata sample numbers don't match...disaster!")
else:
    print(f'Nice! {df.shape[0]} samples in count AND metadata files.\n')

if args.f:
    filtered = filtering.Filtering(df, args.f)
    df = filtered.by_mean()
    filtered_stats = stats.Descriptive(df, 'filtered_variables.csv')
    filtered_stats.matrix_stats()
    df.to_csv(path_or_buf='filteredCountFile.csv')
    check_f_sparse = stats.Density(df)
    if check_f_sparse.is_sparse():
        print("Filtered Count Table is sparse with %2.3f%s non-zero values.\n" %
              (check_f_sparse.percent_non_zeros(), "%"))
    elif check_sparse.is_sparse():
        print("Filtered Count Table is no longer sparse: %2.3f%s non-zero values.\n"
              % (check_f_sparse.percent_non_zeros(), "%"))
    else:
        print("Filtered Count Table now has %2.3f%s non-zero reads.\n" % (check_f_sparse.percent_non_zeros(), "%"))


if args.p:
    pca = clustering.Unsupervised(df, dm)
    pca.pca()

if args.k:
    kmeans = clustering.Unsupervised(df, dm)
    kmeans.kmeans(n_clusters=args.k)

if args.s:
    svm = classify.SVM(df, dm)
    svm.analyze(variable=args.s)

if args.sr:
    roc = classify.SVM(df, dm)
    roc.roc_analysis(variable=args.sr)

if args.sc:
    scomponents = classify.SVM(df, dm)
    scomponents.feature_analysis(variable=args.sc)

if args.rf:
    forest = randomForest.RandomForest(df, dm, variable=args.rf)
    forest.forest()
