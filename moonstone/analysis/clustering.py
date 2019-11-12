from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px


class Unsupervised(object):
    def __init__(self, count_matrix, metadata, outdir):
        self.count_matrix = count_matrix
        self.metadata = metadata
        self.outdir = outdir

    def pca(self):
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
        from mpl_toolkits.mplot3d import Axes3D

        df = self.count_matrix
        x = StandardScaler().fit_transform(df)
        clf = PCA(n_components=3)
        principal_components = clf.fit_transform(x)
        pc1_value, pc2_value, pc3_value = clf.explained_variance_ratio_
        pc1_label = '{:4} : {:5.2f}%'.format('PC1', pc1_value * 100)
        pc2_label = '{:4} : {:5.2f}%'.format('PC2', pc2_value * 100)
        pc3_label = '{:4} : {:5.2f}%'.format('PC3', pc3_value * 100)
        principal_df = pd.DataFrame(data=principal_components,
                                    columns=['pc1', 'pc2', 'pc3'])

        fig = plt.figure(figsize=(8, 8))
        ax = Axes3D(fig)
        ax.set_xlabel(pc1_label)
        ax.set_ylabel(pc2_label)
        ax.set_zlabel(pc3_label)
        ax.scatter3D(principal_df.loc[:, 'pc1'], principal_df.loc[:, 'pc2'], principal_df.loc[:, 'pc3'])
        plt.show()

    def kmeans(self, filename, n_clusters=2):
        from sklearn.cluster import KMeans
        from sklearn import preprocessing

        print('\nRunning K-Means analysis on {} samples and {} variables, setting {} clusters'
              .format(self.count_matrix.shape[0], self.count_matrix.shape[1], n_clusters))
        km = KMeans(init='k-means++', n_clusters=n_clusters, max_iter=1500)
        x = self.count_matrix
        x_normalized = preprocessing.maxabs_scale(x)
        x_standardized = preprocessing.StandardScaler().fit_transform(x)
        x_robust = preprocessing.RobustScaler().fit_transform(x)

        # Generate the data frame with the sample numbers and their cluster
        # In order of cluster members
        data_sets = (x, x_normalized, x_standardized, x_robust)
        set_names = ('Raw', 'Norm', 'Stand', 'Robust_S')
        # Setup the clinical data to which to add the cluster data.
        d_final = self.metadata
        for normalization, name in zip(data_sets, set_names):
            # Train and then fit the data
            y_clusters = km.fit_predict(normalization)

            # Do a count of the clusters and report for each type of normalization
            c = dict(Counter(y_clusters))
            print("Sample distribution over the %s clusters for %s data pretreatment: %s" % (n_clusters, name, c))

            # This section will add cluster labels for each sample in the metadata file
            # The first step is to import the metadata file. The Index is the sample name.
            samples = self.count_matrix.index
            dc = pd.DataFrame(zip(samples, y_clusters), columns=('Sample', name))

            # Count the number of samples in each cluster, and sort from largest to smallest cluster.
            # The grouping value becomes the index.
            cluster_summary = dc.groupby(name).count().sort_values(by="Sample", ascending=False)

            # Add a column from the index, which is now in order of cluster size.
            # This is useful as being independent of the number of clusters set
            # The original clusters are still the index
            cluster_summary['rename'] = np.sort(np.array(cluster_summary.index))

            # Build a dictionary of replacements to make, original clusters for new clusters
            # Add 1 to start counting clusters at 1 instead of 0
            rename_dictionary = dict(zip(cluster_summary.index, cluster_summary['rename']+1))
            dc[name].replace(to_replace=rename_dictionary, inplace=True)
            # Append the cluster information for each sample to the clinical parameters
            d_final = pd.merge(d_final, dc[name], left_index=True, right_index=True)
            d_final.index.name = 'sample'
        output_file = self.outdir+'/'+filename
        print("\n\tAppending KMeans clustering information to original Metadata file. Saving as %s" % output_file)
        d_final.to_csv(path_or_buf=output_file)
        print(d_final.head())

        data_clusters = px.data.gapminder()  # noqa
        fig = px  # noqa

    def gmm(self):
        from sklearn.mixture import GaussianMixture
        from sklearn.decomposition import PCA
        # x_train = self.count_matrix
        # clf = GaussianMixture(n_components=2, covariance_type='full')
        # clf.fit(x_train)
        pca = PCA(0.99)
        data = pca.fit_transform(self.count_matrix)
        n_components = np.arange(10, 500, 20)
        models = [GaussianMixture(n, covariance_type='full', random_state=0)
                  for n in n_components]
        aics = [model.fit(data).aic(data) for model in models]
        plt.plot(n_components, aics)
        plt.show()

        # print(gmm.means_)
        # print('\n')
        # print(gmm.covariances_)
