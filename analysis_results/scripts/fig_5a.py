"""Plot figure 5a: dendogram"""
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hierarchy
import scipy.spatial.distance as distance
from scipy.cluster.hierarchy import dendrogram

import matplotlib as mplt
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.font_manager as font_manager
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.ticker import LogLocator, LogFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import axes3d

from copy import deepcopy
from collections import Counter
import cmocean
import cmasher
import seaborn
import random
import os


# set the font family style
mplt.rcParams['font.family'] = 'Myriad Pro'  # change to a font you have installed on your computer - checkout Google Fonts for free fonts available for download

# set some initial paths

# path to the directory where this script lives
thisdir = os.path.abspath('')

# path to the main directory of the repository
maindir = os.path.split(os.path.split(thisdir)[0])[0]

# path to the analysis_results subdirectory
analysisdir = os.path.split(thisdir)[0]

# path to the data subdirectory
datadir = os.path.join(os.path.split(os.path.split(thisdir)[0])[0], 'data')

# path to the figures subdirectory within analysis_results
figdir = os.path.join(analysisdir, 'figures')

location_file = os.path.join(analysisdir, 'location_country_names.csv')
locations_df = pd.read_csv(location_file, delimiter = ',')


def get_country_name(df,location):
    """
    Return country name of the location

    Args:
        df (pandas Dataframe) : dataframe containing locations and corresponding countries
        location (str)        : name of the location

    Returns:
        Name of the country the location is in.
    """
    d = df[df.location == location]
    return d.country.values[0]


def read_in_distance(measure, num_agebrackets):
    """ Read in the calculated distances between overall contact matrices.
        
        Args:
            measure (str)         : name of the distance measure used
            num_agebrackets (int) : the number of age brackets for the contact matrices

        Returns:
            A pandas dataframe.
    """
    file_path = os.path.join(analysisdir, 'distance_measures', measure + '_distance_M_all_locations_' + '%i' % num_agebrackets + '_new_formulation_matrices.csv')
    df = pd.read_csv(file_path, delimiter = ',')
    return df


def distance_matrix(locations, measure, distance_df, num_agebrackets):
    """
    Create distance matrix, square form matrix, and mappings from location id to index and back.

    Args:
        locations (list): list of locations
        measure (str): name of the distance measure used
        distance_df (pandas dataframe): distance dataframe
        num_agebrackets (int): the number of age brackets for the contact matrices

    Returns:
        A distance matrix, square form distance matrix, dictionary mapping location to index, dictionary mapping index to location

    """
    locations_list = deepcopy(locations)
    N = len(locations_list)
    M = np.zeros((N,N))

    locations_to_index = dict((location, n) for n, location in enumerate(locations_list))
    index_to_locations = dict((n, location) for n, location in enumerate(locations_list))

    for i in range(N):
        li = index_to_locations[i]
        for j in range(N):
            lj = index_to_locations[j]

            if M[i][j] == 0:

                dist = float(distance_df[(distance_df.location_1 == li) & (distance_df.location_2 == lj)].distance.values[0])

                if np.isnan(dist):
                    print(li,lj,dist)
                M[i][j] = dist
                M[j][i] = dist
        M[i][i] = 0

    squareM = distance.squareform(M)
    
    return M, squareM, locations_to_index, index_to_locations


def cluster_matrix(distance_df, M, squareM, locations_to_index, index_to_locations, method, distance_threshold, criterion):
    """
    Create distance matrix, square form matrix, and mappings from location id to index and back.

    Args:
        distance_df (pandas dataframe): distance dataframe
        M (np.ndarray): distance matrix
        squareM (np.ndarray): scipy.distance.squareform of distance matrix M
        locations_to_index (dict): dictionary mapping location to index
        index_to_locations (dict): dictionary mapping index to location
        method (str): linkage method for scipy.cluster.hierarchy
        distance_threshold (float): distance threshold for cutting the dendogram and creating clusters
        criterion (str): clustering criterion type

    Returns:
        A clustered matrix, list of clusters, clustered locations, list of clustered locations, cluster labels, cluster ids, cluster sizes, cluster size by location, clustering from fcluster

    """
    clustering = hierarchy.linkage(squareM, method=method) # cluster the locations together into one giant group step by step starting with the pair of locations most similar to each other
    clusters = hierarchy.fcluster(clustering, t=distance_threshold, criterion=criterion)
    C = Counter(clusters)
    N = len(locations_to_index)

    clustered_locations = [] # group locations together in a list if they're in the same cluster
    cluster_ids = set(clusters) # what are the cluster ids in the network? Get the set of them
    cluster_labels = {} # a dictionary that maps each locationid to the cluster id it belongs to
    cluster_sizes = [] # push back the size for each cluster c
    cluster_size_by_locations = {} # give the cluster size for the cluster each location n belongs to 


    for c in sorted(cluster_ids):
        clustered_locations.append([]) # create a list for each cluster where we will house the actual location ids within
        cluster_sizes.append(C[c]) # add the size for each cluster

    for location_index in range(N):
        c = clusters[location_index] # cluster location_index belongs to
        clustered_locations[c-1].append(location_index) # add the location index to the cluster list
        location_name = index_to_locations[location_index]
        cluster_labels[location_name] = c
        cluster_size_by_locations[location_index] = C[c]

    clustered_locations_list = [] # an ordered list of locations with their original index 
    for cindex in range(len(clustered_locations)):
        clustered_locations_list = clustered_locations_list + clustered_locations[cindex]

    clustered_Matrix = np.zeros((N,N))
    for i in range(N):
        location_index_i = clustered_locations_list[i]
        location_i = index_to_locations[location_index_i]
        for j in range(N):
            location_index_j = clustered_locations_list[j]
            location_j = index_to_locations[location_index_j]

            if clustered_Matrix[i][j] == 0:
                dist = float(distance_df[(distance_df.location_1 == location_i) & (distance_df.location_2 == location_j)].distance.values[0])
                clustered_Matrix[i][j] = dist
                clustered_Matrix[j][i] = dist

    return clustered_Matrix, clusters, clustered_locations, clustered_locations_list, cluster_labels, cluster_ids, cluster_sizes, cluster_size_by_locations, clustering


def plot_clustered_matrix(locations, measure, df, num_agebrackets, cmap=cmocean.cm.matter, save=False, show=False):
    """
    Plot clustered matrix of locations.
        
    Args:
        locations (list)                                        : list of locations
        measure (str)                                           : name of the distance measure used
        distance_df (pandas dataframe)                          : distance dataframe
        num_agebrackets (int)                                   : the number of age brackets for the contact matrices
        cmap (str or matplotlib.colors.LinearSegmentedColormap) : name of the colormap to use or the actual colormap
        save (bool)                                             : If True, save the figure
        show (bool)                                             : If True, show the figure

    Returns:
        Matplotlib figure.

    """
    if isinstance(cmap, str):
        cmap = mplt.cm.get_cmap(cmap)

    criterion = 'distance'
    method = 'single'
    distance_threshold = 1
    color_threshold = 1

    M, squareM, locations_to_index, index_to_locations = distance_matrix(locations, measure, df, num_agebrackets)

    clustered_Matrix, clusters, clustered_locations, clustered_locations_list, cluster_labels, cluster_ids, cluster_sizes, cluster_size_by_locations, clustering = cluster_matrix(df, M, squareM, locations_to_index, index_to_locations, method, distance_threshold, criterion)

    clustering = hierarchy.linkage(M, method=method)

    # figure and axis details
    fig = plt.figure(figsize=(20, 20))

    # dendogram 1
    R = dendrogram(clustering, 
        leaf_rotation = 90, above_threshold_color = 'k', color_threshold = color_threshold, orientation = 'left',
        distance_sort = 'descending',
        no_plot = True,
        leaf_font_size = 0,
    )

    # margins to split up the figure into space for the matrix, the dendogram, labels, and the legend
    small_margin = 0.01
    medium_margin = 0.03
    large_margin = 0.05
    tree_height = 0.4

    matrix_height = 1 - 2 * small_margin - large_margin - tree_height
    colorscale_width_percent = 0.03
    colorscale_padding_percent = 0.03

    colorscale_width = matrix_height * colorscale_width_percent
    colorscale_padding = matrix_height * colorscale_width_percent

    matrix_width = matrix_height + colorscale_width + colorscale_padding
    tree_width = matrix_height
    tree_bottom = 1 - small_margin - tree_height


    legend_padding = 0.05
    legend_width = colorscale_width
    legend_height = matrix_height

    axdendo = fig.add_axes([large_margin,tree_bottom,tree_width,tree_height])
    axmatrix = fig.add_axes([large_margin,large_margin,matrix_width,matrix_height])
    axright = fig.add_axes([large_margin + matrix_width + 3 * small_margin, large_margin, legend_width, legend_height])
    axlegend = fig.add_axes([large_margin + tree_width, tree_bottom, colorscale_width + colorscale_padding + 8 * legend_width, tree_height])
    axlegend.set_axis_off()

    axbottom = fig.add_axes([large_margin, 2 * small_margin, matrix_height, legend_width])

    # dendogram 2
    R2 = dendrogram(clustering, 
        leaf_rotation = 0, above_threshold_color = 'k', color_threshold = color_threshold,
        distance_sort = 'descending',
        orientation = 'top',
        no_plot = False,
        leaf_font_size = 0,
        ax = axdendo
        )

    # turn off the spines
    for a in ['right','top','bottom']:
        axdendo.spines[a].set_visible(False)
    axdendo.yaxis.set_ticks_position('left')

    # dendogram 3
    R3 = dendrogram(clustering,
        leaf_rotation = 0, above_threshold_color = 'k', color_threshold = distance_threshold,
        distance_sort = 'descending',
        orientation = 'top',
        no_plot = True
        )

    # remove axis ticks
    axdendo.set_xticks([])
    axdendo.tick_params(axis = 'y', labelsize = 26, labelrotation = 90)
    idx = R['leaves']
    idx2 = R3['leaves']
    labels = []
    labels_to_show = []
    labels_to_show_ticks = []

    locations_to_label = [
                          'Daman_and_Diu', 'Yamalo_Nenets_Autonomous_Okrug', 'Chukotka', 'Nunavut', 'Northwest_Territories', 'Yukon', 
                          'District_of_Columbia', 'Australian_Capital_Territory', 'Puerto_Rico', 'Hawaii',
                          'Tuva', 'Meghalaya', 'Arunachal_Pradesh', 'Ingushetia', 'Chechnya', 'Andaman_and_Nicobar_Islands',
                          'Mpumalanga', 'KwaZulu-Natal', 'Dagestan'
                          ]

    for ii in range(len(idx)):

        i = idx[ii]
        name = index_to_locations[i]

        labels.append(name.replace('_',' ').replace('-ken','').replace('-to','').replace('-fu','').replace('-',' ').replace('Distict_of_Columbia','Distict_of\nColumbia').replace('Australian_Capital_Territory','Australian\nCapital\nTerritory').replace('Andaman_and_Nicobar_Islands','Andaman & Nicobar\nIslands'))
        if name in locations_to_label:

            labels_to_show.append(name.replace('_',' ').replace('-ken','').replace('-to','').replace('-fu','').replace('-',' ').replace('Autonomous','Auton.').replace('Distict_of_Columbia','Distict_of\nColumbia').replace('Australian_Capital_Territory','Australian\nCapital\nTerritory').replace('Andaman_and_Nicobar_Islands','Andaman & Nicobar\nIslands'))
            labels_to_show_ticks.append(ii)

    D = deepcopy(M)
    D = D[idx,:]
    D = D[:,idx2]

    flipD = D.copy()

    im = axmatrix.matshow(flipD, aspect = 'auto', origin = 'upper', interpolation = 'nearest', cmap = cmap)

    axmatrix.set_xticks([]) # turn off xticks
    axmatrix.set_yticks([])

    divider = make_axes_locatable(axmatrix)
    cax = divider.append_axes('right',size = '3%', pad = '3%')
    cbar = fig.colorbar(im, cax = cax)
    cbar.ax.set_yticklabels([int(cl) for cl in cbar.ax.get_yticks()], rotation = 90, fontsize = 24)

    ### Add a colour bar indicating the country for each location around the left and bottom sides ###
    country_colour_dic = {}

    country_colour_dic['Australia'] = '#0000ff'
    country_colour_dic['United_States'] = '#00ace7'
    country_colour_dic['China'] = '#fcc200'
    country_colour_dic['Canada'] = '#369f5c'
    country_colour_dic['Europe'] = '#941cca'
    country_colour_dic['India'] = 'darkorange'
    country_colour_dic['Israel'] = '#dddddd'
    country_colour_dic['Japan'] = '#000098'
    country_colour_dic['Russia'] = '#dc142b'
    country_colour_dic['South_Africa'] = '#b5d93c'

    countries = ['Australia', 'Canada', 'China', 'Europe', 'India', 'Israel', 'Japan', 'South_Africa', 'Russia', 'United_States']
    country_indices = dict((c, n) for n, c in enumerate(countries))

    ### you must include all countries otherwise it will reset colorbar so that the lowest country index in your list maps to the lowest possible country ###
    ### if Australia is not included but the color is still there in the list, then country index 1 (Canada) maps to Australia's colour, country index 2 (China) maps to Canada's colour, etc. 

    cmap_colourbar = mplt.colors.ListedColormap([ country_colour_dic[country] for n,country in enumerate(countries) ], N = len(countries))
    dcountries = np.zeros((len(locations),1))
    yticks = []
    yticklabels = []

    yticklabel_adjustments = dict.fromkeys(locations_to_label, 0)
    yticklabel_adjustments['Daman_and_Diu'] = 1
    yticklabel_adjustments['Yamalo_Nenets_Autonomous_Okrug'] = 11
    yticklabel_adjustments['Chukotka'] = 24
    yticklabel_adjustments['Nunavut'] = -9
    yticklabel_adjustments['Northwest_Territories'] = 0
    yticklabel_adjustments['Yukon'] = 10
    yticklabel_adjustments['Tuva'] = -42
    yticklabel_adjustments['District_of_Columbia'] = -22
    yticklabel_adjustments['Australian_Capital_Territory'] = -1
    yticklabel_adjustments['Puerto_Rico'] = 15
    yticklabel_adjustments['Hawaii'] = 25
    yticklabel_adjustments['Arunachal_Pradesh'] = -34
    yticklabel_adjustments['Meghalaya'] = -24
    yticklabel_adjustments['Chechnya'] = -15
    yticklabel_adjustments['Ingushetia'] = -5
    yticklabel_adjustments['Dagestan'] = 2
    yticklabel_adjustments['Andaman_and_Nicobar_Islands'] = 5
    yticklabel_adjustments['Mpumalanga'] = 24.5
    yticklabel_adjustments['KwaZulu-Natal'] = 15.5

    for ii in range(len(idx)):
        i = idx[ii]
        name = index_to_locations[i]
        country = get_country_name(locations_df,name)
        dcountries[ii] = country_indices[country]
        if name in locations_to_label:
            ytick_pos = ii
            ytick_pos += yticklabel_adjustments[name]
            yticks.append(ytick_pos)

            name = name.replace('_',' ').replace('-ken','').replace('-to','').replace('-fu','').replace('-',' ').replace('Autonomous','Auton.').replace('District of Columbia','District of\nColumbia').replace('Australian Capital Territory','Australian        \nCapital     \nTerritory').replace('Andaman and Nicobar Islands','Andaman & Nicobar        \nIslands').replace('Nenets ','Nenets      \n')
            if 'Islands' not in name:
                name += '        '
            elif 'Islands' in name:
                name = name.rstrip('        ')

            yticklabels.append(name)

    im_right = axright.matshow(dcountries, aspect = 'auto', origin = 'upper', cmap = cmap_colourbar)

    axright.set_xticks([])
    axright.tick_params(axis='y', right=True, labelright=True, left=False, labelleft=False, length=0)
    axright.set_yticks(yticks)
    axright.set_yticklabels(yticklabels, rotation=180, fontsize=22, fontweight='heavy', fontstyle='oblique')

    im_bottom = axbottom.matshow(dcountries.T, aspect = 'auto', origin = 'upper', cmap = cmap_colourbar)
    axbottom.set_xticks([])
    axbottom.set_yticks([])

    if show:
        plt.show()

    if save:
        file_name = 'fig_5a.pdf'
        fig_path = os.path.join(figdir, file_name)
        fig.savefig(fig_path, format='pdf')


if __name__ == '__main__':

    measure = 'canberra'
    num_agebrackets = 85

    df = read_in_distance(measure, num_agebrackets)
    locations = list(set(df.location_1.values))

    criterion = 'distance'
    method = 'average'
    distance_threshold = 1
    color_threshold = 1
    save = True  # save to pdf in figures folder
    show = False  # show on screen

    plot_clustered_matrix(locations, measure, df, num_agebrackets, cmap=cmocean.cm.matter, save=save, show=show)







