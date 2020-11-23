"""Plot figure 4: multiple overall contact matrices."""
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


def plot_multi_matrix(locations, countries, levels, setting, num_agebrackets=85, cmap=cmocean.cm.deep_r, min_CM=1e-4, max_CM=1e-1, save=False, show=False, stacked=False, rows=1, cols=3):
    """
    Plot a single contact matrix.

    Args:
        locations (list)                                        : list of the locations
        countries (list)                                        : list of the countries corresponding to the locations
        levels (list)                                           : list of levels (country or subnational) corresponding to the locations
        setting (str)                                           : name of the contact setting
        num_agebrackets (int)                                   : the number of age brackets for the matrix
        cmap (str or matplotlib.colors.LinearSegmentedColormap) : name of the colormap to use or the actual colormap
        min_CM (float)                                          : minimum value for the colorscale
        max_CM (float)                                          : maximum value for the colorscale
        save (bool)                                             : If True, save the figure
        show (bool)                                             : If True, show the figure
        stacked (bool)                                          : If True, arrange the subplots into multiple rows

    Returns:
        Matplotlib figure.
    """

    if isinstance(cmap, str):
        cmap = mplt.cm.get_cmap(cmap)

    fontsizes = {'colorbar': 30, 'colorbarlabels': 24, 'title': 44, 'ylabel': 28, 'xlabel': 28, 'xticks': 24, 'yticks': 24}

    n = len(locations)  # how many places to plot

    width = 9
    height = 8

    if stacked:
        if n < rows * cols:
            raise ValueError('Please increase either the number of rows or columns to fit all of the subplots in the figure.')
        width = height
        fig, ax = plt.subplots(rows, cols, figsize=(cols * width, rows * height))
    else:
        fig, ax = plt.subplots(1, n, figsize=(n * width, height))

    if stacked:
        ax_list = ax.copy()
        ax = []
        for r in range(rows):
            for c in range(cols):
                ax.append(ax_list[r][c])

    # defaults without colorbar label
    if stacked:
        left = 0.06
        right = 0.92
        top = 0.92
        bottom = 0.08
        wspace = 0.35
        hspace = 0.35

    else:
        left = 0.05
        right = 0.94
        top = 0.90
        bottom = 0.115
        wspace = 0.35
        hspace = 0.5

    fig.subplots_adjust(left=left, right=right, top=top, bottom=bottom, hspace=hspace, wspace=wspace)

    matrices = []
    ims = []
    for i, location in enumerate(locations):
        country = countries[i]
        level = levels[i]
        matrix = read_contact_matrix(location, country, level, setting, num_agebrackets)
        matrices.append(matrix)

        im = ax[i].imshow(matrix.T, origin='lower', interpolation='nearest', cmap=cmap, norm=LogNorm(vmin=min_CM, vmax=max_CM))
        divider = make_axes_locatable(ax[i])
        cax = divider.append_axes('right', size='4%', pad=0.15)
        cbar = fig.colorbar(im, cax=cax)
        cbar.ax.tick_params(labelsize=fontsizes['colorbarlabels'])

        if setting == 'overall': cbar.set_label('$M_{ij}$', fontsize = fontsizes['colorbar'])
        if setting == 'household': cbar.set_label('$H_{ij}$', fontsize = fontsizes['colorbar'])
        if setting == 'school': cbar.set_label('$S_{ij}$', fontsize = fontsizes['colorbar'])
        if setting == 'work': cbar.set_label('$W_{ij}$', fontsize = fontsizes['colorbar'])
        if setting == 'community': cbar.set_label('$R_{ij}$', fontsize = fontsizes['colorbar'])

        title = location.replace('_', ' ').replace('-ken', '').replace('-to', '').replace('-fu', '').replace('-', ' ')
        ax[i].set_title(title, fontsize=fontsizes['title'])
        ax[i].set_ylabel('Age of contact', fontsize=fontsizes['ylabel'])
        ax[i].set_xlabel('Age', fontsize=fontsizes['xlabel'])
        ax[i].set_xticks(np.arange(0, 81, 10))
        ax[i].set_yticks(np.arange(0, 81, 10))
        ax[i].tick_params(labelsize=fontsizes['xticks'])

    if show:
        plt.show()

    if save:
        setting_type, setting_suffix = 'F', 'setting'
        if setting == 'overall': 
            setting_type, setting_suffix = 'M', 'contact_matrix'

        if stacked:
            file_name = 'multi_stacked_' + setting + '_' + setting_suffix + '_' + '%i' % num_agebrackets + '.pdf'
        else:
            file_name = 'multi_' + setting + '_' + setting_suffix + '_' + '%i' % num_agebrackets + '.pdf'
        fig_path = os.path.join(figdir, file_name)
        fig.savefig(fig_path, format='pdf')


if __name__ == '__main__':

    # Examples of use

    setting = 'overall'
    num_agebrackets = 85

    cmap = cmocean.cm.deep_r

    save = True
    show = False
    stacked = False
    rows, cols = 1, 3

    # these are colorbar limits that work best for the separate social layers
    min_CM, max_CM = 1e-4, 1e-1
    # the overall combined matrix has different units and thus different colorbar limits 
    if setting == 'overall':
        min_CM, max_CM = 1e-2, 1e1

    # 3 locations in 1 row
    locations = ['Beijing', 'New_York', 'Maharashtra']
    countries = ['China', 'United_States', 'India']
    levels = ['subnational', 'subnational', 'subnational']

    plot_multi_matrix(locations, countries, levels, setting, num_agebrackets, cmap=cmap, min_CM=min_CM, max_CM=max_CM, save=save, show=show, rows=rows, cols=cols, stacked=stacked)


    # 4 locations in 2 rows and 2 columns (stacked)
    locations = ['Beijing', 'New_York', 'Maharashtra', 'Meghalaya']
    countries = ['China', 'United_States', 'India', 'India']
    levels = ['subnational', 'subnational', 'subnational', 'subnational']

    stacked = True
    rows, cols = 2, 2

    plot_multi_matrix(locations, countries, levels, setting, num_agebrackets, cmap=cmap, min_CM=min_CM, max_CM=max_CM, save=save, show=show, rows=rows, cols=cols, stacked=stacked)

