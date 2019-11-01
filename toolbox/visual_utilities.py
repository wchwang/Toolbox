# Created by woochanghwang at 2019-03-27

import seaborn as sns
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

class MidpointNormalize(mpl.colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=0, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        normalized_min = max(0, 1 / 2 * (1 - abs((self.midpoint - self.vmin) / (self.midpoint - self.vmax))))
        normalized_max = min(1, 1 / 2 * (1 + abs((self.vmax - self.midpoint) / (self.midpoint - self.vmin))))
        normalized_mid = 0.5
        x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
        return np.ma.masked_array(np.interp(value, x, y))

def read_network():
    print("test")

def set_clustermap_side_colors(side_series, color_set):
    '''
    Set row/col colors
    :param side_series: row_data from dataframe
    :param color_set: color set from colot palettes
    :return: row_colors
    '''
    rgb_colors = sns.color_palette(color_set, len(side_series.unique()))  # http://seaborn.pydata.org/tutorial/color_palettes.html

    tissue_type_color = dict(zip(side_series.unique(), rgb_colors))
    return side_series.map(tissue_type_color), tissue_type_color

def set_denrogram_row_color_legend(cm, tissue_type_color, legend_label, legend_loc):

    for label in legend_label.unique():
        cm.ax_row_dendrogram.bar(1,0, color=tissue_type_color[label],
                                label=label, linewidth=0)

    cm.ax_row_dendrogram.legend(bbox_to_anchor=(2.5, 1), loc=legend_loc)  # bbox_to_anchor: located out box
    # cm.ax_row_dendrogram.set_position([.15, .2, .03, .45])


def set_side_colormap_adjust(cm,row_ratio=1, col_ratio=1):

    row = cm.ax_row_dendrogram.get_position()
    row_c = cm.ax_row_colors.get_position()
    cm.ax_row_dendrogram.set_position([row.x0 + row_c.width * (1 - row_ratio), row.y0, row.width, row.height])
    cm.ax_row_colors.set_position([row_c.x0 + row_c.width * (1 - row_ratio), row_c.y0, row_c.width * row_ratio, row_c.height])

    # col = cm.ax_col_dendrogram.get_position()
    # col_c = cm.ax_col_colors.get_position()
    # cm.ax_col_dendrogram.set_position([col.x0 + col_c.width * (1 - col_ratio), row.y0, row.width, row.height])
    # cm.ax_col_colors.set_position(
    #     [col_c.x0 + col_c.width * (1 - col_ratio), col_c.y0, col_c.width * col_ratio, col_c.height])


def pca_plot(pca,pca_result,data_df):
    ######################
    ## Draw PCA Result
    ######################

    plt.figure(figsize=(8, 8))
    plt.scatter(pca_result[:, 0], pca_result[:, 1])

    for idx, symbol in enumerate(data_df.index):
        print(idx, symbol)

        # Add the text label
        labelpad = 0.01  # Adjust this based on your dataset

        plt.text(pca_result[idx, 0] + labelpad, pca_result[idx, 1] + labelpad, symbol, fontsize=7)

        # Mark the labeled observations with a star marker
        plt.scatter(pca_result[idx, 0], pca_result[idx, 1], s=100)

    # Add the axis labels
    plt.xlabel('PC 1 (%.2f%%)' % (pca.explained_variance_ratio_[0] * 100))
    plt.ylabel('PC 2 (%.2f%%)' % (pca.explained_variance_ratio_[1] * 100))

    # Done
    plt.show()

def draw_venn_2group(group1, group2, group_labels, save_addr,title='Comparision of significant gene set(2)'):
    from matplotlib_venn import venn2

    venn2([set(group1), set(group2)], set_labels=group_labels)

    plt.title(title)
    plt.savefig(save_addr)
    plt.show()

def draw_venn_3group(group1, group2, group3,group_labels,save_addr,title='Comparision of significant 3 gene set'):
    '''
    venn3([set(lst1), set(lst2), set(lst3)], set_labels = ('Drug A responsive genes', 'Drug B responsive genes', 'Drug C responsive genes'))
$ plt.title('Network of drug responsive genes: The great update\n')

    '''

    venn3([set(group1), set(group2), set(group3)], set_labels=group_labels)

    plt.title(title)
    plt.savefig(save_addr)
    plt.show()
# def cluster_dataframe(data_result_df, save_addr):
#
#     sub_subtypes= data_result_df.loc["disease_sub_subtype"]
#
#     print(depmap_lof_data_df)
#
#
#     ################################
#     ## Cluster
#     ################################
#
#
#     rgb_colors = sns.color_palette("hls",
#                                    len(sub_subtypes.unique()))  # http://seaborn.pydata.org/tutorial/color_palettes.html
#
#     # ccle_tissue_unique = ccle_tissue.unique()
#     # tmp = ccle_tissue_unique[0]
#     # ccle_tissue_unique[0] = ccle_tissue_unique[1]
#     # ccle_tissue_unique[1]=tmp
#     #
#     # print(ccle_tissue_unique)
#     #
#     sub_type_color = dict(zip(sub_subtypes.unique(), rgb_colors))
#     #
#     # print(tissue_type_color)
#     # print(ccle_tissue)
#     # print(data_result_df['Tissue'])
#     col_colors = sub_subtypes.map(sub_type_color)
#
#     # print(row_colors)
#     # cmap = sns.cubehelix_palette(start=0.3, rot= -.5, dark=0, light=1 , as_cmap=True, reverse=True) #color map
#     # cmap = sns.cubehelix_palette(start=0.3, rot= -.5, dark=0.8, light=0 , as_cmap=True)
#     # cmap = sns.color_palette("RdBu_r", 7)
#     # cmap = sns.color_palette("Reds")
#     # cm = sns.clustermap(depmap_lof_data_df, metric="euclidean", cmap=cmap, robust=True, method="average",
#     #                    row_colors=row_colors)  # Average is best
#     # cm = sns.clustermap(depmap_lof_data_df, metric="correlation", cmap=cmap, robust=True, method="average",row_colors=row_colors)  # Average is best
#
#     # depmap_lof_data_df.to_csv("test_2.csv")
#     cm = sns.clustermap(depmap_lof_data_df, metric="correlation", method="average",robust=True, row_cluster=False,col_colors=col_colors)  # Average is best, no y-label
#
#
#     # ########################################
#     # ## Adjust Clustermap, size, colour, label, rotate
#     # ##########################################
#     #
#     ## col dendrogram legend
#     for label in sub_subtypes.unique():
#         cm.ax_col_dendrogram.bar(1,0, color=sub_type_color[label],
#                                 label=label, linewidth=0)
#     cm.ax_col_dendrogram.legend(loc="upper right", ncol=1)
#     # cm.ax_col_dendrogram.legend(bbox_to_anchor=(2.5, 1), loc="upper left")  # bbox_to_anchor: located out box
#     # cm.ax_row_dendrogram.set_position([.15, .2, .03, .45])
#     #
#     # ## label
#     # # cm.cax.set_position([.5, .5, .03, .45])
#     cm.cax.set_position([.15, .2, .03, .45])
#     # plt.setp(cm.ax_heatmap.get_yticklabels(), rotation=0)  # For y axis
#     # plt.setp(cm.ax_heatmap.get_xticklabels(), rotation=90)  # For x axis
#     plt.tight_layout()
#     #
#     # ## Heatmap size adjust
#     # hm = cm.ax_heatmap.get_position()
#     # plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels(), fontsize=6)
#     # cm.ax_heatmap.set_position([hm.x0, hm.y0, hm.width * 0.25, hm.height])
#     # col = cm.ax_col_dendrogram.get_position()
#     #
#     # row = cm.ax_row_dendrogram.get_position()
#     # row_c = cm.ax_row_colors.get_position()
#     #
#     # cm.ax_col_dendrogram.set_position([col.x0, col.y0, col.width * 0.25, col.height * 0.25])
#     #
#     # # cm.ax_row_dendrogram.set_position([row.x0-row_c.width * 0.5, row.y0,row.width, row.height])
#     # # cm.ax_row_colors.set_position([row_c.x0-row_c.width * 0.5,row_c.y0, row_c.width * 0.5, row_c.height])
#     #
#     # cm.ax_row_dendrogram.set_position([row.x0+row_c.width * (1-0.3), row.y0,row.width, row.height])
#     # cm.ax_row_colors.set_position([row_c.x0+row_c.width * (1-0.3),row_c.y0, row_c.width * 0.3, row_c.height])
#     #
#     # #############
#     # ## https://github.com/mwaskom/seaborn/pull/1393
#     # #ClusterGrid layouts
#
#
#     ################
#     cm.savefig(save_addr+'.pdf')
#
#
#     plt.show()
#
#     # print(cm.dendrogram_col.reordered_ind)
#     # print(cm.dendrogram_row.reordered_ind)
#     #
#     # # print(data_result_df.loc[[233]])
#     # ccle_cluster_df = data_result_df.reindex(index = cm.dendrogram_row.reordered_ind, columns = cm.dendrogram_col.reordered_ind)
#     # ccle_cluster_df.to_csv(save_addr + "cluster.csv")
