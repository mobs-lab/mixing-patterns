"""Plot inferred contact matrices by setting and overall contact matrices."""
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


def plot_single_matrix(location, country, level, setting, num_agebrackets=85, cmap=cmocean.cm.deep_r, min_CM=1e-4, max_CM=1e-1, save=False, show=False):
    """
    Plot a single contact matrix.

    Args:
        location (str)                                          : name of the location
        country (str)                                           : name of the country
        level (str)                                             : name of level (country or subnational)
        setting (str)                                           : name of the contact setting
        num_agebrackets (int)                                   : the number of age brackets for the matrix
        cmap (str or matplotlib.colors.LinearSegmentedColormap) : name of the colormap to use
        min_CM (float)                                          : minimum value for the colorscale
        max_CM (float)                                          : maximum value for the colorscale
        save (bool)                                             : If True, save the figure
        show (bool)                                             : If True, show the figure

    Returns:
        Matplotlib figure.
    """

    if isinstance(cmap, str):
        cmap = mplt.cm.get_cmap(cmap)

    fontsizes = {'colorbar': 30, 'colorbarlabels': 22, 'title': 44, 'ylabel': 28, 'xlabel': 28, 'xticks': 24, 'yticks': 24}

    # defaults without colorbar label
    left = 0.12
    right = 0.925
    top = 0.90
    bottom = 0.095

    # move margins in a little
    left -= 0.00
    right -= 0.07
    top -= 0.03
    bottom += 0.03

    fig = plt.figure(figsize=(8, 8))
    fig.subplots_adjust(left=left, right=right, top=top, bottom=bottom)
    ax = fig.add_subplot(111)

    matrix = read_contact_matrix(location, country, level, setting, num_agebrackets)

    im = ax.imshow(matrix.T, origin='lower', interpolation='nearest', cmap=cmap, norm=LogNorm(vmin=min_CM, vmax=max_CM))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='4%', pad=0.15)
    cbar = fig.colorbar(im, cax=cax)
    cbar.ax.tick_params(labelsize=fontsizes['colorbarlabels'])

    if setting == 'overall': cbar.set_label('$M_{ij}$', fontsize = fontsizes['colorbar'])
    if setting == 'household': cbar.set_label('$H_{ij}$', fontsize = fontsizes['colorbar'])
    if setting == 'school': cbar.set_label('$S_{ij}$', fontsize = fontsizes['colorbar'])
    if setting == 'work': cbar.set_label('$W_{ij}$', fontsize = fontsizes['colorbar'])
    if setting == 'community': cbar.set_label('$R_{ij}$', fontsize = fontsizes['colorbar'])

    title = location.replace('_', ' ').replace('-ken', '').replace('-to', '').replace('-fu', '').replace('-', ' ')
    ax.set_title(title, fontsize=fontsizes['title'])
    ax.set_ylabel('Age of contact', fontsize=fontsizes['ylabel'])
    ax.set_xlabel('Age', fontsize=fontsizes['xlabel'])
    ax.set_xticks(np.arange(0, 81, 10))
    ax.set_yticks(np.arange(0, 81, 10))
    ax.tick_params(labelsize=fontsizes['xticks'])

    if show:
        plt.show()

    if save:
        setting_type, setting_suffix = 'F', 'setting'
        if setting == 'overall':
            setting_type, setting_suffix = 'M', 'contact_matrix'

        file_name = country + '_' + level + '_' + location + '_' + setting_type + '_' + setting + '_' + setting_suffix + '_' + '%i' % num_agebrackets + '.pdf'
        fig_path = os.path.join(figdir, file_name)
        fig.savefig(fig_path, format='pdf')


if __name__ == '__main__':

    # Examples of use

    # plotting a country
    location = 'China'
    country = 'China'
    level = 'country'

    # # plotting a subnational location
    # location = 'Beijing'
    # country = 'China'
    # level = 'subnational'

    settings = ['household', 'school', 'work', 'community', 'overall']

    num_agebrackets = 85
    cmap = cmocean.cm.deep_r  # define the colormap you want to use

    save = True  # saves in the figures folder
    show = False  # show to screen

    for setting in settings:

        # these are colorbar limits that work best for the separate social layers
        min_CM, max_CM = 1e-4, 1e-1
        # the overall combined matrix has different units and thus different colorbar limits 
        if setting == 'overall':
            min_CM, max_CM = 1e-2, 1e1

        # plot the matrix
        plot_single_matrix(location, country, level, setting, num_agebrackets, cmap=cmap, min_CM=min_CM, max_CM=max_CM, save=save, show=show)







