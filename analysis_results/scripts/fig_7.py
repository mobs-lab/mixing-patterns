"""Plot figure 7b: 3 scatter plots and 3 maps"""
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hierarchy
import scipy.spatial.distance as distance
from scipy.cluster.hierarchy import dendrogram

import geopandas
from descartes import PolygonPatch
import shapely.geometry as sgeom
from cartopy.feature import ShapelyFeature
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature

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
import cmasher as cmr
import seaborn as sns
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
        str: Name of the country the location is in.
    """
    d = df[df.location == location]
    return d.country.values[0]


def get_locations_by_country(df,country):
    """
    Return locations in the country

    Args:
        df (pandas Dataframe) : dataframe containing locations and corresponding countries
        country (str)         : name of the country

    Returns:
        list: Available locations in the country.
    """
    locations = list(df[df.country == country].location.values)
    return locations


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
        np.ndarray: A numpy square matrix of contact.
    """
    setting_type, setting_suffix = 'F', 'setting'
    if setting == 'overall':
        setting_type, setting_suffix = 'M', 'contact_matrix'

    if country == 'Europe':
        country = location
        level = 'country'

    if level == 'country':
        file_name = country + '_' + level + '_level_' + setting_type + '_' + setting + '_' + setting_suffix + '_' + '%i' % num_agebrackets + '.csv'
    else:
        file_name = country + '_' + level + '_' + location + '_' + setting_type + '_' + setting + '_' + setting_suffix + '_' + '%i' % num_agebrackets + '.csv'
    file_path = os.path.join(datadir, 'contact_matrices', file_name)
    M = np.loadtxt(file_path, delimiter=',')
    return M


def get_attack_rates_df(reference_location, reference_scenario, R0_star, susceptibility, gamma_inverse, num_agebrackets):
    """
    Get attack rates dataframe.

    Args:
        reference_location (str) : name of reference location or locations
        reference_scenario (str) : specific reference scenario
        R0_star (float)          : the basic reproduction number
        susceptibility (float)   : susceptibility of adults to under 18
        gamma_inverse (float)    : the mean recovery period
        num_agebrackets (int)    : the number of age brackets for the matrix
    
    Returns:
        Pandas dataframe of attack rates by location names
    """

    file_path = os.path.join(analysisdir, 'reference_location_polymod_and_tomsk', 'mcmc_beta_and_dropfactor_scenario_' + reference_scenario + '_country_R0',
                             'all_locations_attack_rates_by_age_reference_scenario_' + reference_scenario + '_country_R0_' + '%.2f' % R0_star  + '_susceptibility_' + '%.2f' % susceptibility + '_gamma_inverse_' + '%.1f' % gamma_inverse + '_' + str(num_agebrackets) + '.csv')
    return pd.read_csv(file_path)


def get_attack_rates_country_R0_df(reference_location, reference_scenario, R0_star, susceptibility, gamma_inverse, num_agebrackets):
    """
    Get attack rates dataframe.

    Args:
        reference_location (str) : name of reference location or locations
        reference_scenario (str) : specific reference scenario
        R0_star (float)          : basic reproduction number for the country matrix
        susceptibility (float)   : susceptibility of adults to under 18
        gamma_inverse (float)    : the mean recovery period
        num_agebrackets (int)    : the number of age brackets for the matrix
    
    Returns:
        Pandas dataframe of attack rates by location names
    """

    file_path = os.path.join(analysisdir, 'reference_location_polymod_and_tomsk_with_country_matrix', 'mcmc_beta_and_dropfactor_scenario_' + reference_scenario + '_country_R0',
                             'all_locations_attack_rates_by_age_reference_scenario_' + reference_scenario + '_country_R0_' + '%.2f' % R0_star + '_susceptibility_' + '%.2f' % susceptibility + '_gamma_inverse_' + '%.1f' % gamma_inverse + '_' + str(num_agebrackets) + '.csv')
    return pd.read_csv(file_path)


def get_attack_rate(df, location):
    """
    Get the attack rate for the location.

    Args:
        df (pd.DataFrame) : pandas DataFrame with attack rates for this scenario
        location (str)    : name of the location

    Returns:
        float: Attack rate as a fraction from SIR model, values between 0 and 1. 
    """
    return df.loc[df['location'] == location]['artotal'].values[0]


def rescale_x(x,min_x,max_x):
    """
    Rescale or nornalize value from 0 to 1 with defined limits
    
    Args:
        x (float)     : value to be rescaled
        min_x (float) : lower limit of range
        max_x (float) : upper limit of range
    
    Returns:
        float: Value between 0 and 1 representing a normalized value of x.
    """
    return (x - min_x)/(max_x - min_x)


def get_eigenvalue(matrix):
    """
    Get the real component of the leading eigenvalue of a square matrix.

    Args:
        matrix (np.ndarray): square matrix

    Returns:
        float: Real component of the leading eigenvalue of the matrix.
    """
    eigenvalue = max(np.linalg.eigvals(matrix)).real
    return eigenvalue


def get_beta(R0, gamma_inverse, matrix):
    """
    Get the transmissibilty, beta, for an SIR compartmental model given the basic reproduction number, the mean recovery period, and the age specific contact matrix.
    
    Args:
        R0 (float)            : the basic reproduction number
        gamma_inverse (float) : the mean recovery period
        matrix (np.ndarray)   : the contact matrix

    Returns:
        float: The transmissibility, beta, for an SIR compartmental model with age specific contact patterns.
    """
    gamma = float(1)/gamma_inverse
    eigenvalue = get_eigenvalue(matrix)
    return R0 * gamma / eigenvalue


def get_R0(beta, gamma_inverse, matrix):
    """
    Get the basic reproduction number, R0, for an SIR compartmental model given the basic reproduction number, the mean recovery period, and the age specific contact matrix.
    
    Args:
        beta (float)          : the transmissibility
        gamma_inverse (float) : the mean recovery period
        matrix (np.ndarray)   : the contact matrix

    Returns:
        float: The basic reproduction number R0 for an SIR compartmental model with age specific contact patterns.
    """
    gamma = 1./gamma_inverse
    eigenvalue = get_eigenvalue(matrix)
    return beta * eigenvalue / gamma


def get_factor_vector(susceptibility_drop_factor, susceptibility_drop_age, num_agebrackets):
    """
    Get a vector with the age susceptibility drop factor.

    Args:
        susceptibility_drop_factor (float) : susceptibility drop factor
        susceptibility_drop_age (int)      : age at which susceptibility drops
        num_agebrackets (int)              : the number of age brackets

    Returns:
        np.ndarray: A numpy array that gives the relative susceptibility drop factor for ages past the susceptibility_drop_age.
    """
    factor_vector = np.ones((num_agebrackets,1))
    factor_vector[susceptibility_drop_age:] *= susceptibility_drop_factor
    return factor_vector


def get_age_effective_contact_matrix(contact_matrix, susceptibility_drop_factor, susceptibility_drop_age):
    """
    Get an effective age specific contact matrix with an age dependent susceptibility drop factor.

    Args:
        contact_matrix (np.ndarray)        : the contact matrix
        susceptibility_drop_factor (float) : susceptibility drop factor
        susceptibility_drop_age (int)      : age at which susceptibility drops

    Returns:
        np.ndarray: A numpy square matrix that gives the effective contact matrix given an age dependent susceptibility drop factor.
    """
    factor_vector = get_factor_vector(susceptibility_drop_factor,susceptibility_drop_age,len(contact_matrix))
    effective_matrix = contact_matrix * factor_vector
    return effective_matrix


def get_R0_with_factor_vector(beta, factor_vector, gamma_inverse, contact_matrix):
    """
    Get the basic reproduction number R0 for an SIR model with an age dependent susceptibility drop factor.

    Args:
        beta (float)                       : the transmissibility
        factor_vector (np.ndarray)         : array of relative susceptibility factor by age
        gamma_inverse (float)              : the mean recovery period
        contact_matrix (np.ndarray)        : the contact matrix

    Returns:
        float: The basic reproduction number for an SIR model with an age dependent susceptibility drop factor.
    """
    gamma = 1./gamma_inverse
    effective_matrix = contact_matrix.T * factor_vector
    eigenvalue = get_eigenvalue( effective_matrix )
    R0 = beta * eigenvalue / gamma
    return R0


def get_beta_with_beta_vector(factor_vector, R0, gamma_inverse, contact_matrix):
    """
    Get the transmissibility beta for an SIR model with an age dependent susceptibility drop factor.

    Args:
        R0 (float)                  : the basic reproduction number
        factor_vector (np.ndarray)  : array of relative susceptibility factor by age
        gamma_inverse (float)       : the mean recovery period
        contact_matrix (np.ndarray) : the contact matrix

    Returns:
        float: The transmissibility beta for an SIR model with an age dependent susceptibility drop factor.
    """
    gamma = 1./gamma_inverse
    effective_matrix = contact_matrix.T * factor_vector
    eigenvalue = get_eigenvalue( effective_matrix )
    beta = R0 * gamma / eigenvalue
    return beta


def get_shapefile_record_name_df():
    """
    Returns a pandas DataFrame of shapefile record names and the associated locations.
    """
    df = pd.read_csv(os.path.join(analysisdir, 'shapefile_record_names_dictionary.csv'))
    return df


def get_shapefile_record_name_from_location(shapefile_record_names_df, location):
    """
    Get the shapefile record name for the location given.

    Args:
        shapefile_record_names_df (pandas DataFrame) : dataframe of shapefile record names and the associated locations
        location (str)                               : name of the location

    Returns:
        str: Record name for the location as listed in the shapefiles used.
    """
    return shapefile_record_names_df.loc[shapefile_record_names_df['location'] == location].record_name.values[0]


def get_location_from_shapefile_record_name(shapefile_record_names_df, record_name):
    """
    Get location for the shapefile record name given.

    Args:
        shapefile_record_names_df (pandas DataFrame) : dataframe of shapefile record names and the associated locations
        location (str)                               : name of the location

    Returns:
        str: Location for the shapefile record name as listed in the shapefiles used.
    """
    return shapefile_record_names_df.loc[shapefile_record_names_df['record_name'] == record_name].location.values[0]


def map_location_from_shapefile_record_name(shapefile_record, country, country_categories, shapefile_record_names_df):
    """
    Map the shapefile record name to the location name depending on how the shapefile records are organized.

    Args:
        shapefile_record (cartopy.io.shapereader.FionaRecord) : the shapefile FionaRecord
        country (str)                                         : name of the country
        country_categories (dict)                             : dictionary of countries separated into different categories depending on how the shapefile FionaRecord is organized
        shapefile_record_names_df (pandas DataFrame)          : dataframe of shapefile record names and the associated locations

    Returns:
        str: Location name for the shapefile record. 
    """

    if country in country_categories[1]:
        state_name = shapefile_record.attributes['STATE_NAME'].replace(' ', '_')

    elif country in country_categories[2]:
        state_name = shapefile_record.attributes['NAME_1'].replace(' And ', ' and ').replace(' ', '_')

        if state_name in shapefile_record_names_df.record_name.values:
            state_name = get_location_from_shapefile_record_name(shapefile_record_names_df, state_name)

    elif country in country_categories[3]:
        state_name = shapefile_record.attributes['NAME_1']
        state_name = state_name.replace(' ', '_')

        if state_name in shapefile_record_names_df.record_name.values:
            state_name = get_location_from_shapefile_record_name(shapefile_record_names_df, state_name)

    elif country in country_categories[4]:
        state_name = shapefile_record.attributes['NAME_ENGLI']

        if state_name in shapefile_record_names_df.record_name.values:
            state_name = get_location_from_shapefile_record_name(shapefile_record_names_df, state_name)

    elif country in country_categories[5]:
        state_name = shapefile_record.attributes['NAME_0']

    return state_name


def define_map_extents():
    """ Define map extents for the different countries."""
    countries = ['Australia','Canada','China','India','Europe','Israel','Russia','South_Africa','Japan','United_States']

    llcrnrlat_dic = dict.fromkeys(countries,-90)
    llcrnrlon_dic = dict.fromkeys(countries,-180)
    urcrnrlat_dic = dict.fromkeys(countries,90)
    urcrnrlon_dic = dict.fromkeys(countries,180)

    llcrnrlat_dic['China'] = 16
    llcrnrlat_dic['India'] = 6
    llcrnrlat_dic['United_States'] = 17
    llcrnrlat_dic['Japan'] = 30
    llcrnrlat_dic['Russia'] = 40
    llcrnrlat_dic['Australia'] = -45 
    llcrnrlat_dic['Canada'] = 40
    llcrnrlat_dic['Europe'] = 31
    llcrnrlat_dic['Israel'] = 29
    llcrnrlat_dic['South_Africa'] = -36


    llcrnrlon_dic['China'] = 73
    llcrnrlon_dic['India'] = 65
    llcrnrlon_dic['United_States'] = -140
    llcrnrlon_dic['Japan'] = 128
    llcrnrlon_dic['Australia'] = 112
    llcrnrlon_dic['Canada'] = -142
    llcrnrlon_dic['Europe'] = -12
    llcrnrlon_dic['Israel'] = 33
    llcrnrlon_dic['South_Africa'] = 15


    urcrnrlat_dic['China'] = 49
    urcrnrlat_dic['India'] = 37
    urcrnrlat_dic['United_States'] = 50
    urcrnrlat_dic['Japan'] = 47
    urcrnrlat_dic['Russia'] = 90
    urcrnrlat_dic['Australia'] = -6.5
    urcrnrlat_dic['Canada'] = 84
    urcrnrlat_dic['Europe'] = 71
    urcrnrlat_dic['Israel'] = 35
    urcrnrlat_dic['South_Africa'] = -21

    urcrnrlon_dic['China'] = 139
    urcrnrlon_dic['India'] = 100
    urcrnrlon_dic['United_States'] = -66
    urcrnrlon_dic['Japan'] = 150
    urcrnrlon_dic['Australia'] = 155.5
    urcrnrlon_dic['Canada'] = -50
    urcrnrlon_dic['Europe'] = 36
    urcrnrlon_dic['Israel'] = 36
    urcrnrlon_dic['South_Africa'] = 33

    return llcrnrlon_dic, llcrnrlat_dic, urcrnrlon_dic, urcrnrlat_dic


def plot_fig(countries, reference_scenario, R0_star, gamma_inverse, beta, susceptibility_drop_factor, num_agebrackets):
    """
    Plot figure of the difference in attack rate and R0 for each country and the states within it in one large figure with multiple panels.
    
    Args:
        countries (list)                   : list of countries
        reference_scenario (str)           : label of the reference scenario
        R0_star (float)                    : the basic reproduction number
        gamma_inverse (float)              : the mean recovery period
        beta (float)                       : the transmissibility
        susceptibility_drop_factor (float) : susceptibility of adults to those under 18
        num_agebrackets (int)              : the number of age brackets

    Returns:
        Matplotlib figure.
    """
    shapefile_record_names_df = get_shapefile_record_name_df()

    susceptibility_drop_age = 18

    country_colour_dic = {}
    country_colour_dic['Australia'] = '#0000ff'
    country_colour_dic['Canada'] = '#2ab207'
    country_colour_dic['China'] = '#fcc200'
    country_colour_dic['Europe'] = '#941cca'
    country_colour_dic['India'] = 'darkorange'
    country_colour_dic['Israel'] = '#9b9b9b'
    country_colour_dic['Japan'] = '#000098'
    country_colour_dic['Russia'] = '#dc142b'
    country_colour_dic['South_Africa'] = '#b5d93c'
    country_colour_dic['United_States'] = '#00ace7'

    country_codes = {'United_States': 'USA', 'China': 'CHN', 'Australia': 'AUS', 'Canada': 'CAN', 'India': 'IND',
                      'South_Africa': 'ZAF', 'Japan': 'JPN', 'Russia': 'RUS', 'Europe': 'EUR', 'Israel': 'ISR'}

    shapefiledir = os.path.join('/Users', 'dina', 'Dropbox (MOBS)', 'shapefiles')  # change this before pushing to repo
    shapefiles = {}

    country_categories = {}
    country_categories[1] = ['United_States']
    country_categories[2] = ['China', 'Australia', 'Canada', 'India']
    country_categories[3] = ['South_Africa', 'Japan', 'Russia']
    country_categories[4] = ['Europe']
    country_categories[5] = ['Israel']

    locations_not_included = set(['Dadra_and_Nagar_Haveli', 'Lakshadweep', 'Chandigarh', 'Ashmore_and_Cartier_Islands','Coral_Sea_Islands','Jervis_Bay_Territory',
                              'Albania','Armenia','Azerbaijan','Belarus','Belgium','Bosnia and Herzegovina','Croatia','Faroe Islands',
                             'Iceland','Kosovo','Macedonia','Moldova','Montenegro','Northern Cyprus','Poland','Russia','Serbia',
                              'Turkey','Ukraine','Russian_Federation', 'Ryazan', "Ryazan'"])

    for country in countries:
        if country in country_categories[1]:
            shpf = os.path.join(shapefiledir, 'USA_Qian', 'states_PRI.shp')

        elif country in country_categories[2]:
            shpf = os.path.join(shapefiledir, country_codes[country] + '_adm_shp', country_codes[country] + '_adm1.shp')

        elif country in country_categories[3]:
            shpf = os.path.join(shapefiledir, 'gadm36_' + country_codes[country] + '_shp', 'gadm36_' + country_codes[country] + '_1.shp')

        elif country in country_categories[4]:
            shpf = os.path.join(shapefiledir, country, country + '.shp')

        elif country in country_categories[5]:
            shpf = os.path.join(shapefiledir, 'gadm36_' + country_codes[country] + '_shp', 'gadm36_' + country_codes[country] + '_1.shp')
        shapefiles[country] = shpf

    attack_rates_df = get_attack_rates_df(reference_location, reference_scenario, R0_star, susceptibility_drop_factor, gamma_inverse, num_agebrackets)
    attack_rates_country_R0_df = get_attack_rates_country_R0_df(reference_location, reference_scenario, R0_star, susceptibility_drop_factor, gamma_inverse, num_agebrackets)

    width = 9.
    height = 9
    marker_size = 50

    fontsizes = {'legend': 20, 'title': 48, 'xlabel': 28, 'ylabel': 28, 'xticks': 24, 'yticks': 24}
    for f in fontsizes:
        fontsizes[f] *= 1.25

    left = 0.06
    right = 0.97
    bottom = 0.04
    top = 0.92
    hspace = 0.30
    wspace = 0.30

    cmap = cmocean.cm.matter_r

    if susceptibility_drop_factor == 1:
        max_ar_percent_diff = 15
        min_ar_percent_diff = -15
    elif susceptibility_drop_factor == 0.63:
        max_ar_percent_diff = 20
        min_ar_percent_diff = -20
    else:
        max_ar_percent_diff = 25
        min_ar_percent_diff = -25

    fig, ax = plt.subplots(2, len(countries), figsize=(width * len(countries), height * 2))
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top, hspace=hspace, wspace=wspace)
    ax = ax.reshape((ax.size,))
    leg = []

    ax2 = fig.add_axes([-0.2, -0.2, 0.1, 0.1])  # won't show up
    L = 200
    cbrange = np.zeros((8, L + 1))
    for i in range(L + 1):
        cbrange[:, i] = (i-L/2)/(L/2) * max_ar_percent_diff

    im = ax2.imshow(cbrange, cmap=cmap)

    llcrnrlon_dic, llcrnrlat_dic, urcrnrlon_dic, urcrnrlat_dic = define_map_extents()

    for c, country in enumerate(countries):
        locations = get_locations_by_country(locations_df, country)

        if country in locations and country != 'Israel':
            locations.remove(country)

        if country == 'India':
            locations.remove('Dadra_and_Nagar_Haveli')
            locations.remove('Chandigarh')
            locations.remove('Lakshadweep')

        if country == 'Russia':
            locations.remove('Russian_Federation')

        R0_list = []
        R0_country_list = []
        attack_rates_list = []
        attack_rates_country_list = []

        ar_dic = {}
        country_ar_dic = {}
        ar_percent_diff_dic = {}

        matrix_location = country
        if country == 'Europe':
            matrix_location = 'United-Kingdom'
        elif country == 'Russia':
            matrix_location = 'Russian_Federation'

        reference_matrix = read_contact_matrix(matrix_location, country, 'country', 'overall')
        factor_vector = get_factor_vector(susceptibility_drop_factor,susceptibility_drop_age,num_agebrackets)

        reference_beta = get_beta_with_beta_vector(factor_vector,R0_star,gamma_inverse,reference_matrix)

        for i, location in enumerate(locations):
            matrix = read_contact_matrix(location, country, 'subnational', 'overall')
            R0 = get_R0_with_factor_vector(reference_beta, factor_vector, gamma_inverse, matrix)
            R0_list.append(R0)
            R0_country_list.append(R0_star)

            ar = get_attack_rate(attack_rates_df, location) * 100
            car = get_attack_rate(attack_rates_country_R0_df, location) * 100

            attack_rates_list.append(ar)
            attack_rates_country_list.append(car)

            ar_dic[location] = ar
            country_ar_dic[location] = car
            ar_percent_diff_dic[location] = (car - ar)/car * 100

        ax[c].scatter(R0_list, attack_rates_list, color = country_colour_dic[country], s=marker_size, edgecolor = country_colour_dic[country], label = 'State matrix results')
        ax[c].scatter(R0_country_list, attack_rates_country_list, color='k', s=marker_size, edgecolor='k', label = 'Country matrix results')

        for i in range(len(R0_list)):
            ax[c].plot([R0_list[i], R0_country_list[i]], [attack_rates_list[i], attack_rates_country_list[i]], color = country_colour_dic[country], lw=2)

        if susceptibility_drop_factor == 1:
            xlim_min, xlim_max = 1.35, 1.65
            ax[c].set_xlim(xlim_min, xlim_max)
            ax[c].set_xticks(np.arange(xlim_min, xlim_max + 0.01, 0.05))

            ylim_min, ylim_max = 40, 65

            ax[c].set_ylim(ylim_min, ylim_max)
            ax[c].set_yticks(np.arange(ylim_min, ylim_max + 1, 5))

        elif round(susceptibility_drop_factor, 2) == 0.63:

            xlim_min, xlim_max = 1.35, 1.65
            ax[c].set_xlim(xlim_min, xlim_max)
            ax[c].set_xticks(np.arange(xlim_min, xlim_max + 0.01, 0.05))

            ylim_min, ylim_max = 15, 55

            ax[c].set_ylim(ylim_min, ylim_max)
            ax[c].set_yticks(np.arange(ylim_min, ylim_max + 1, 5))

        ax[c].set_title(country.replace('_', ' ').replace('-', ' '), fontsize = fontsizes['title'])
        ax[c].set_xlabel(r'$R_0$', fontsize = fontsizes['xlabel'])
        ax[c].set_ylabel('Attack Rates (%)', fontsize = fontsizes['ylabel'])
        ax[c].tick_params(labelsize=fontsizes['xticks'])
        leg.append(ax[c].legend(loc=2, fontsize=fontsizes['legend']))
        leg[c].draw_frame(False)
        if c < len(countries)-1:
            nn = 20
            x = np.ones(nn) * 1.68
            y = np.linspace(7, 68.5, nn)
            ax[c].plot(x, y, color = 'k', ls = (0, (13.5, 4)), lw = 1, clip_on=False)

        ax[c + len(countries)].axis('off')
        ax[c + len(countries)] = fig.add_subplot(2, len(countries), c + len(countries) + 1, projection=ccrs.PlateCarree(), frameon=False)
        ax[c + len(countries)].set_xticks([])
        ax[c + len(countries)].set_yticks([])

        extent = [llcrnrlon_dic[country], urcrnrlon_dic[country], llcrnrlat_dic[country], urcrnrlat_dic[country]]
        reader = shpreader.Reader(shapefiles[country])

        for shapefile_record in reader.records():
            location = map_location_from_shapefile_record_name(shapefile_record, country, country_categories, shapefile_record_names_df)
            if country == 'Europe' and location == 'Georgia':
                pass

            try:
                if country == 'Japan':
                    if location in shapefile_record_names_df.location:
                        x = ar_percent_diff_dic[location]
                    else:
                        try:
                            x = ar_percent_diff_dic[location + '-to']
                        except:
                            try:
                                x = ar_percent_diff_dic[location + '-ken']
                            except:
                                x = ar_percent_diff_dic[location + '-fu']
                elif country == 'Russia':
                    if location in shapefile_record_names_df.location:
                        x = ar_percent_diff_dic[location]
                    else:
                        try:
                            x = ar_percent_diff_dic[location + '_Oblast']
                        except:
                            x = ar_percent_diff_dic[location]

                else:
                    x = ar_percent_diff_dic[location]
            except:
                continue

            x = rescale_x(x, min_ar_percent_diff, max_ar_percent_diff)

            ax[c + len(countries)].add_geometries([shapefile_record.geometry], ccrs.PlateCarree(), facecolor=cmap(x), edgecolor='k', lw=0.2, clip_on=False)

        cax, kw = mplt.colorbar.make_axes(ax[c + len(countries)],location='bottom',pad=0.05, aspect=25, shrink=0.9)
        cbar = fig.colorbar(im, cax=cax, **kw)
        cbar.ax.tick_params(labelsize=24)
        cbar.set_label('% Variation of Attack Rate', fontsize=28)

        ax[c + len(countries)].set_extent(extent, crs=ccrs.PlateCarree())
        ax[c + len(countries)].tick_params(colors='white')


    # fig_path = os.path.join(figdir, 'fig_7.pdf')  # use at your own peril - pdfs are gigantic for this figure
    # fig.savefig(fig_path, format='pdf')

    fig_path = os.path.join(figdir, 'fig_7.png')
    fig.savefig(fig_path, format='png', dpi=150)


if __name__ == '__main__':
    
    beta = 0.04752
    susceptibility_drop_factor = 1.
    gamma_inverse = 2.6
    R0_star = 1.5
    num_agebrackets = 85

    reference_location = 'polymod_and_tomsk'
    reference_scenario = 'all_locations'

    countries = ['China', 'United_States', 'India']

    plot_fig(countries, reference_scenario, R0_star, gamma_inverse, beta, susceptibility_drop_factor, num_agebrackets)





