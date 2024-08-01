from sklearn.cluster import MiniBatchKMeans
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
import sklearn.metrics as metrics
from scipy.stats import ttest_ind
import pandas as pd
import numpy as np
import seaborn as sns


class CellNeighborhood(object):
    def __init__(self,cells, Xpos, YPos, ROI, CellType, outpath="CellNeighborhoodOutput", n_neighbors=7, method='Distancecutoff', method_param=50):
        self.name = outpath
        self.cells = cells
        self._set_column_names(Xpos, YPos, ROI, CellType)
        self.neighbor_nums = n_neighbors
        assert (method == 'Windowcutoff' or method == 'Distancecutoff')
        self.method = method
        self.method_param = method_param
        self.save_path = 'Distance' + str(self.method_param) + '_'

    def _set_column_names(self, X_position, Y_position, ROI, CellType):
        self.X_position = X_position
        self.Y_position = Y_position
        self.ROI = ROI
        self.CellType = CellType
        self.keep_cols = [X_position, Y_position, ROI, CellType]
        self.save_path = ''
        self.method = ''
        self.neighbor_nums = 0
        self.neighborhood_name = ''  # This is used to track whether neighborhoods are constructed or not


    def get_windows(self, job):
        idx, tissue_name, indices = job
        tissue = self.tissue_group.get_group(tissue_name)
        to_fit = tissue.loc[indices][[self.X_position, self.Y_position]].values
        fit = NearestNeighbors(n_neighbors=self.method_param).fit(tissue[[self.X_position, self.Y_position]].values)
        m = fit.kneighbors(to_fit)
        m = m[0], m[1]

        # sort_neighbors
        args = m[0].argsort(axis=1)
        add = np.arange(m[1].shape[0]) * m[1].shape[1]
        sorted_indices = m[1].flatten()[args + add[:, None]]
        neighbors = tissue.index.values[sorted_indices]
        return neighbors.astype(np.int32)  # returns k neighbor indices for each cell in the sample

    def get_neighbors_distance(self, job):
        idx, tissue_name, indices = job
        tissue = self.tissue_group.get_group(tissue_name)
        to_fit = tissue.loc[indices][[self.X_position, self.Y_position]].values

        point_tree = cKDTree(to_fit)
        neighbors = []
        indices_set = point_tree.query_ball_point(to_fit, self.method_param).flatten()
        for indices in indices_set:
            indexlist = [x for x in indices]
            neighbors.append(tissue.index.values[indexlist])
        return neighbors

    def identifyNeighborhoods(self):
        df = pd.concat([self.cells, pd.get_dummies(self.cells[self.CellType])], axis=1)  # converts cell types to indicators

        self.sum_cols = df[self.CellType].unique()  # sum_cols in the unique phenotypes
        self.values = df[self.sum_cols].values  # values in those indicator cell types

        self.tissue_group = df[[self.X_position, self.Y_position, self.ROI]].groupby(self.ROI)
        exps = list(df[self.ROI].unique())
        tissue_chunks = [(exps.index(t), t, a) for t, indices in self.tissue_group.groups.items() for a in np.array_split(indices, 1)]   ##### debug here

        if self.method == 'Windowcutoff':
            tissues = [self.get_windows(job) for job in tissue_chunks]
        if self.method == 'Distancecutoff':
            tissues = [self.get_neighbors_distance(job) for job in tissue_chunks]

        self.tissues=tissues
        self.tissue_chunks=tissue_chunks
        out_dict = {}

        for neighbors, job in zip(tissues, tissue_chunks):
            chunk = np.arange(len(neighbors))  # indices
            tissue_name = job[1]
            indices = job[2]
            if self.method == 'Windowcutoff':
                window = self.values[neighbors[chunk].flatten()].reshape(len(chunk), self.method_param,len(self.sum_cols)).sum(axis=1)
                out_dict[(tissue_name, self.method_param)] = (window.astype(np.float16), indices)

            if self.method == 'Distancecutoff':
                window = np.zeros((len(neighbors), len(self.sum_cols)))
                for i in range(len(neighbors)):
                   
                    window[i] += self.values[neighbors[i]].sum(axis=0)  # add counts of its neighbors celltypes together as its value (The total counts of its neighbor for each point)
                out_dict[(tissue_name, self.method_param)] = (window.astype(np.float16), indices)

        win = [pd.DataFrame(out_dict[(exp, self.method_param)][0], index=out_dict[(exp, self.method_param)][1].astype(int),columns=self.sum_cols) for exp in exps]
        allwindow = pd.concat(win, axis=0)  # all the groups concatenated
        allwindow = allwindow.loc[df.index.values]
        allwindow = pd.concat([df[self.keep_cols], allwindow], axis=1)


        # Now Perform the clustering
        self.km = MiniBatchKMeans(n_clusters=self.neighbor_nums, random_state=0)  # DBSCAN, KMedoid, Mean Shift
 
        df_c = allwindow[self.sum_cols]
        df_c_norm = df_c.div(df_c.sum(axis=1), axis=0)
        
        labelskm = self.km.fit_predict(df_c_norm.values)
      
        self.k_centroids = self.km.cluster_centers_

        self.neighborhood_name = "Neighborhood" + str(self.neighbor_nums)
        self.cells[self.neighborhood_name] = labelskm
        self.cells[self.neighborhood_name] = self.cells[self.neighborhood_name].astype('category')

        return(pd.concat([self.cells, df_c_norm], axis=1))



    def get_neighborhoods(self):
        if self.neighborhood_name == '':
            self.identifyNeighborhoods()
        return self.cells[
            self.neighborhood_name]  # returns a Pandas series containing neihborhood no for each cell type

    def scale_mat(self, matrix, range_max=5, range_min = -5):
        new_mat = matrix
        xmin = matrix.min()
        xmax = matrix.max()
        convert_to_range = lambda x, xmin, xmax, range_max, range_min: (range_max - range_min) * (x - xmin) / (xmax - xmin) + range_min
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                new_mat[i][j] = convert_to_range(matrix[i][j], xmin, xmax, range_max, range_min)
        return new_mat

    def save_clustermap(self,scale=True):
        self.clustermap(save=True)
        


    def get_neighborhood_lfc(self, ratio):
        if self.neighborhood_name == '':
            self.identifyNeighborhoods()

        niche_clusters = self.k_centroids
        tissue_avgs = self.values.mean(axis=0)

        fc = np.log2(
            ((niche_clusters + ratio*tissue_avgs) / (niche_clusters + ratio*tissue_avgs).sum(axis=1, keepdims=True)) / tissue_avgs)

        fc = self.scale_mat(fc, range_max=5, range_min = -5)
        fc = pd.DataFrame(fc, columns=self.sum_cols)
        
        print(fc.shape)
        return (fc)
 
    def clustermap(self, save=False):
        if self.neighborhood_name == '':
            self.identifyNeighborhoods()

        niche_clusters = self.k_centroids
        tissue_avgs = self.values.mean(axis=0)

        fc = np.log2(
            ((niche_clusters + 0.001*tissue_avgs) / (niche_clusters + 0.001*tissue_avgs).sum(axis=1, keepdims=True)) / tissue_avgs)
        fc = pd.DataFrame(fc, columns=self.sum_cols)
        hT = fc.T
        vmax=np.max(np.abs(hT.to_numpy()))
        s = sns.clustermap(hT.T,annot=True, cmap="bwr",row_cluster=False,col_cluster=False) #vmin=-vmax, vmax=vmax)
                           

        if save:
            s.savefig(self.save_path + self.name + 'ClusterMap' + str(self.method_param) + self.neighborhood_name + '.png')
            fc.to_csv(self.save_path + self.name + str(self.method_param) + 'ClusterMapDataframe.csv')
        else:
            return s, fc

    def save_neighborhoods(self):

        self.cells.to_csv(self.save_path + self.name + '_with_Neighborhood.csv',index=False)

    def set_group_plot_colums(self, patient, group):
        self.Group = group
        self.patient = patient

    def get_Groupwise_Stripplot(self):
        groupslist = self.cells[self.Group].unique()
        if len(groupslist) > 2:
            raise ValueError("More than 2 groups is not supported currently.")
        if len(groupslist) < 2:
            raise ValueError("At least 2 groups are required for group based analysis.")

        fc = self.cells.groupby([self.Patient, self.Group]).apply(
            lambda x: x[self.neighborhood_name].value_counts(sort=False, normalize=True))
        fc.columns = range(self.neighoborhood_nums)
        melt = pd.melt(fc.reset_index(), id_vars=[self.Patient, self.Group])
        melt = melt.rename(columns={'variable': 'neighborhood', 'value': 'frequency of neighborhood'})
        f, ax = plt.subplots(dpi=200, figsize=(10, 5))
        sns.stripplot(data=melt, hue=self.Group, dodge=True, alpha=.2, x='neighborhood',
                      y='frequency of neighborhood')
        sns.pointplot(data=melt, scatter_kws={'marker': 'd'}, hue=self.Group, dodge=.5, join=False,
                      x='neighborhood', y='frequency of neighborhood')
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[:2], labels[:2], title=self.Group, handletextpad=0, columnspacing=1, loc="upper right",
                  ncol=3, frameon=True)

        plt.savefig('GroupDiffWindow' + self.neighborhood_name + '.png')

        f = open('GroupStatWindow' + self.neighborhood_name + '.txt', 'w')
        for i in range(self.neighoborhood_nums):
            n2 = melt[melt['neighborhood'] == i]
            f.write(str(i) + ' ' + str(ttest_ind(n2[n2[self.Group] == groupslist[0]]['frequency of neighborhood'],n2[n2[self.Group] == groupslist[1]]['frequency of neighborhood'])[1]) + '\n')
        f.close()
