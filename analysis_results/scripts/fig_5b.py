"""Plot figure 5b: map"""
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
mplt.rcParams['font.family'] = 'Myriad Pro'  # change to a font you have installed on your computer - checkout Google Fonts for free fonts available for download, otherwise this defaults to the system default

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
    Get the country name of the location.

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


def read_in_distance(measure, num_agebrackets):
    """ Read in the calculated distances between overall contact matrices.
        
        Args:
            measure (str)         : name of the distance measure used
            num_agebrackets (int) : the number of age brackets for the contact matrices

        Returns:
            A pandas DataFrame of distance measures between two locations.
    """
    file_path = os.path.join(analysisdir, 'distance_measures', measure + '_distance_M_all_locations_' + '%i' % num_agebrackets + '_new_formulation_matrices.csv')
    df = pd.read_csv(file_path, delimiter = ',')
    return df


def get_distance(distance_df, location_1, location_2):
    """
    Return distance value between two locations.

    Args:
        distance_df (pandas DataFrame) : distance dataframe
        location_1 (str)               : name of the first location
        location_2 (str)               : name of the second location

    Returns:
        float: Distance between two locations.
    """
    dist = float(distance_df[(distance_df.location_1 == location_1) & (distance_df.location_2 == location_2)].distance.values[0])
    if location_1 == location_2:
        dist = 0
    return dist


# def get_shapefile_state_to_ids(country):
#     """
#     Return dictionary mapping state to id
    
#     Args:
#         country (str): name of the country

#     Returns:
#         dict: A dictionary mapping state names to shapefile ID numbers.
#     """
#     file_path = os.path.join(datadir, 'shapefiles', country + '_shapefile_state_ids.csv')
#     df = pd.read_csv(file_path, delimiter = ',')
#     state_to_ids = dict(zip(df.location,df.location_id))
#     return state_to_ids


# def get_shapefile_ids_to_state(country):
#     """
#     Return dictionary mapping id to state
    
#     Args:
#         country (str): name of the country

#     Returns:
#         dict: A dictionary mapping shapefile ID numbers to state names.
#     """
#     file_path = os.path.join(datadir, 'shapefiles', country + '_shapefile_state_ids.csv')
#     df = pd.read_csv(file_path, delimiter = ',')
#     ids_to_state = dict(zip(df.location_id,df.location))


def rescale_x(x,min_x,max_x):
    """
    Rescale or nornalize value from 0 to 1 with defined limits.
    
    Args:
        x (float)     : value to be rescaled
        min_x (float) : lower limit of range
        max_x (float) : upper limit of range
    
    Returns:
        float: A value between 0 and 1 representing the normalized value of x given defined limits.
    """
    return (x - min_x)/(max_x - min_x)


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


def plot_map(countries, reference_location, measure, shapefile_record_names_df, num_agebrackets=85):
    """
    Plot map of states colored by distance from reference location.

    Args:
        countries (list)                             : list of countries
        reference_location (str)                     : name of reference location
        measure (str)                                : name of distance measure used
        shapefile_record_names_df (pandas DataFrame) : dataframe of shapefile record names and the associated locations
        num_agebrackets (int)                        : number of age brackets for the matrices used to calculate the distances

    Returns:
        Matplotlib figure.

    """
    fig = plt.figure(figsize=(10, 5))
    fig.subplots_adjust(bottom=0.08, top =0.99, left=0.04, right=0.96)
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    ax.coastlines(resolution='110m', color='k', lw=0.4)

    country_codes = {'United_States': 'USA', 'China': 'CHN', 'Australia': 'AUS', 'Canada': 'CAN', 'India': 'IND',
                      'South_Africa': 'ZAF', 'Japan': 'JPN', 'Russia': 'RUS', 'Europe': 'EUR', 'Israel': 'ISR'}

    shpfilename = shpreader.natural_earth(resolution='110m', category='cultural', name='admin_0_countries')
    reader = shpreader.Reader(shpfilename)
    ne_countries = reader.records()

    for nec in ne_countries:
        ax.add_geometries([nec.geometry],
         crs=ccrs.PlateCarree(),
         facecolor='#cccccc',
         edgecolor='k', lw=0.4)

    distance_df = read_in_distance(measure, num_agebrackets)

    shapefiledir = os.path.join('/Users', 'dina', 'Dropbox (MOBS)', 'shapefiles')  # change this before pushing to repo

    shapefiles = {}

    country_categories = {}
    country_categories[1] = []
    country_categories[2] = ['China', 'Australia', 'Canada', 'India', 'United_States']
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

    p = []
    distance_list = []
    min_x = 0
    max_x = 2500
    reader = shpreader.Reader(shpf)

    cmap = cmocean.cm.matter

    distance_dic = {}

    for country in countries:
        locations = set(get_locations_by_country(locations_df, country))
        if country != 'Israel' and country in locations:
            locations.remove(country)
        locations = locations - locations_not_included

        for location in sorted(locations):
            dist = get_distance(distance_df, reference_location, location)
            distance_dic[(country, location)] = dist

    for country in countries:
        reader = shpreader.Reader(shapefiles[country])

        for shapefile_record in reader.records():

            location = map_location_from_shapefile_record_name(shapefile_record, country, country_categories, shapefile_record_names_df)
            if country == 'Europe' and location == 'Georgia':
                pass

            poly = shapefile_record.geometry

            try:
                if country == 'Japan':
                    if location in shapefile_record_names_df.location:
                        x = distance_dic[(country, location)]
                    else:
                        try:
                            x = distance_dic[(country, location + '-to')]
                        except:
                            try:
                                x = distance_dic[(country, location + '-ken')]
                            except:
                                x = distance_dic[(country, location + '-fu')]
                elif country == 'Russia':
                    if location in shapefile_record_names_df.location:
                        x = distance_dic[(country, location)]
                    else:
                        try:
                            x = distance_dic[(country, location + '_Oblast')]
                        except:
                            x = distance_dic[(country, location)]

                else:
                    x = distance_dic[(country, location)]

                x = rescale_x(x,min_x,max_x)

                if poly.geom_type == 'MultiPolygon':
                    for pol in poly:
                        p.append(PolygonPatch(pol))
                        distance_list.append(x)
                else:
                    p.append(PolygonPatch(poly))
                    distance_list.append(x)

            except:
                continue

            ax.add_geometries([shapefile_record.geometry], ccrs.PlateCarree(), facecolor=cmap(x), edgecolor='k', lw=0.2)

    ax.set_extent([-180, 180, -57, 110], crs=ccrs.PlateCarree())

    ax2 = fig.add_axes([-0.08, -0.14, 0.84, 0.01])  # don't show the axis itself, just need the colorbar
    L = 200
    cbrange = np.zeros((8, L + 1))
    for i in range(L+1):
        cbrange[:, i] = float(i)/(L+0) * max_x
    im = ax2.imshow(cbrange, cmap=cmap)

    cax, kw = mplt.colorbar.make_axes(ax,location='bottom',pad=0.05, aspect=33)
    cbar = fig.colorbar(im, cax=cax, **kw)
    cbar.ax.tick_params(labelsize=20)
    cbar.set_label(measure.replace('_distance_M', ' ').title() + ' Distance from ' + reference_location.replace('_', ' ').replace('-', ' '), fontsize=24)

    ax.set_frame_on(False)

    # fig_path = os.path.join(figdir, 'fig_5b.pdf')
    # fig.savefig(fig_path, format='pdf')  # use at your own peril - pdfs are gigantic for this figure

    fig_path = os.path.join(figdir, 'fig_5b.png')
    fig.savefig(fig_path, format='png', dpi=450)


if __name__ == '__main__':

    measure = 'canberra'
    shapefile_record_names_df = get_shapefile_record_name_df()
    num_agebrackets = 85

    countries = [
    'Australia', 
    'Canada', 
    'China', 
    'Europe',
    'India', 
    'Israel', 
    'Japan',
    'Russia',
    'South_Africa',
    'United_States'
    ]
    reference_location = 'New_York'
    plot_map(countries, reference_location, measure, shapefile_record_names_df)

