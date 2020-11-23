"""Plot survey and synthetic matrices for France, Japan, and Shanghai, China as shown in figure 3."""
import numpy as np
import pandas as pd
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
import matplotlib.backends.backend_pdf
import cmocean
import cmasher
import seaborn
import copy
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


def read_contact_matrix(location, country, level, setting, num_agebrackets=85):
    """
    Read in the contact for each setting.

    Args:
        location (str)        : name of the location
        country (str)         : name of the country
        level (str)           : name of level (country or subnational)
        setting (str)         : name of the contact setting
        num_agebrackets (int) : the number of age brackets for the matrix

    Returns:
        np.ndarray: A numpy matrix of contact.
    """
    setting_type, setting_suffix = 'F', 'setting'
    if setting == 'overall':
        setting_type, setting_suffix = 'M', 'contact_matrix'

    if level == 'country':
        file_name = country + '_' + level + '_level_' + setting_type + '_' + setting + '_' + setting_suffix + '_' + '%i' % num_agebrackets + '.csv'
    else:
        file_name = country + '_' + level + '_' + location + '_' + setting_type + '_' + setting + '_' + setting_suffix + '_' + '%i' % num_agebrackets + '.csv'
    file_path = os.path.join(datadir, 'contact_matrices', file_name)
    M = np.loadtxt(file_path, delimiter=',')
    return M


def read_validation_contact_matrix(location, matrix_type):
    """
    Read in either the survey or synthetic contact matrix used in the validation analysis of figure 3b.


    Args:
        location (str)    : name of the location
        matrix_type (str) : the type of contact matrix, either 'Model' for the overall contact matrix from a linear combination of the synthetic setting contact matrices, or 'Survey' for the empirical survey matrix.

    Returns:
        np.ndarray: A numpy matrix of contact.
    """
    file_path = os.path.join(analysisdir, 'survey_model_matrix_comparisons', matrix_type + '_' + location + '.csv')
    delimiter = ' '
    M = np.loadtxt(file_path, delimiter=delimiter)
    return M


def plot_matrix(location, matrix_type, cmap=cmocean.cm.deep_r, save=False, show=False):
    """
    Plot the matrices from figure 3b.

    Args:
        location (str): name of the location
        matrix_type (str) : the type of contact matrix, either 'Model' for the overall contact matrix from a linear combination of the synthetic setting contact matrices, or 'Survey' for the empirical survey matrix.
        cmap (str or matplotlib.colors.LinearSegmentedColormap) : name of the colormap to use
        save (bool)                                             : If True, save the figure
        show (bool)                                             : If True, show the figure
    
    Returns:
        Matplotlib figure.
    """

    if isinstance(cmap, str):
        cmap = mplt.cm.get_cmap(cmap)

    fontsizes = {'colorbar': 30, 'colorbarlabels': 22, 'title': 44, 'ylabel': 28, 'xlabel': 28, 'xticks': 24, 'yticks': 24}

    # defaults without colorbar label
    left = 0.155
    right = 0.935
    top = 0.90
    bottom = 0.12

    # move margins in a little
    left -= 0.00
    right -= 0.07
    top -= 0.03
    bottom += 0.03

    fig = plt.figure(figsize=(8, 8))
    fig.subplots_adjust(left=left, right=right, top=top, bottom=bottom)
    ax = fig.add_subplot(111)

    matrix = read_validation_contact_matrix(location, matrix_type)

    num_agebrackets = len(matrix)
    for a in range(num_agebrackets):
        matrix[a, :] = matrix[a, :]/np.sum(matrix[a, :])
    matrix = matrix/np.sum(matrix)


    min_CM = 1e-1
    max_CM = 1e-4
    im = ax.imshow(matrix.T, origin='lower', interpolation='nearest', cmap=cmap, norm=LogNorm(vmin=min_CM, vmax=max_CM))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='4%', pad=0.15)
    cbar = fig.colorbar(im, cax=cax)
    cbar.ax.tick_params(labelsize=fontsizes['colorbarlabels'])
    cbar.ax.set_ylabel('Frequency of contacts', fontsize=fontsizes['ylabel'])

    title_dic = {'Model': 'Synthetic', 'Survey': 'Survey'}
    ax.set_title(title_dic[matrix_type] + ' Matrix', fontsize=fontsizes['title'])
    ax.set_ylabel('Age of contact', fontsize=fontsizes['ylabel'])
    ax.set_xlabel('Age', fontsize=fontsizes['xlabel'])


    if num_agebrackets == 15:
        age_brackets = [str(5*i) + '-' + str(5 * (i+1)-1) for i in range(14)] + ['70+']

    elif num_agebrackets == 17:
        age_brackets = ['0-2', '3-6', '7-9'] + [str(5*i) + '-' + str(5 * (i+1) - 1) for i in range(2, 15)] + ['75+']

    ax.set_xticks(np.arange(len(age_brackets)))
    ax.set_xticklabels(age_brackets, rotation=60)
    ax.set_yticks(np.arange(len(age_brackets)))
    ax.set_yticklabels(age_brackets)
    ax.tick_params(labelsize=fontsizes['xticks'])

    if show:
        plt.show()

    if save:

        file_name = 'fig_3b_' + matrix_type + '_' + location + '.pdf'
        fig_path = os.path.join(figdir, file_name)
        fig.savefig(fig_path, format='pdf')


if __name__ == '__main__':
    
    # location = 'France'
    # location = 'Japan'
    location = 'Shanghai'

    # matrix_type = 'Model'
    matrix_type = 'Survey'

    save = True
    show = False

    plot_matrix(location, matrix_type, save=save, show=show)



