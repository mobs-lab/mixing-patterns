"""
Fill in the table listing year of source data. Currently, only for the age
distribution data. This information is also reported in our article

Mistry, D., Litvinova, M., Pastore y Piontti, A. et al. Inferring
high-resolution human mixing patterns for disease modeling. Nat Commun 12, 323
(2021). https://doi.org/10.1038/s41467-020-20544-y

"""

import numpy as np
import pandas as pd
import copy
import os


# set some initial paths

# path to the directory where this script lives
thisdir = os.path.abspath('')

# path to the main directory of the repository
maindir = os.path.split(os.path.split(thisdir)[0])[0]

# path to the analysis_results subdirectory
analysisdir = os.path.split(thisdir)[0]

# path to the data subdirectory
datadir = os.path.join(os.path.split(os.path.split(thisdir)[0])[0], 'data')

location_file = os.path.join(analysisdir, 'location_country_names.csv')
locations_df = pd.read_csv(location_file, delimiter=',')


def get_country_name(df, location):
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


def get_locations_by_country(df, country):
    """
    Return locations in the country

    Args:
        df (pandas Dataframe) : dataframe containing locations and corresponding countries
        country (str)         : name of the country

    Returns:
        str: Name of the country the location is in.
    """
    locations = list(df[df.country == country].location.values)
    return locations


def update_age_source_year(df, locations, num_agebrackets, source_year):
    """
    Update the data for year of the source data in the dataframe.
    """
    for li, location in enumerate(locations):
        try:
            country = get_country_name(locations_df, location)
        except:
            country = location
        if country == 'Europe':
            country = location
            level = 'country'
        if country == location:
            level = 'country'
        else:
            level = 'subnational'

        if level == 'country':
            file_name = f"{country}_{level}_level_age_distribution_{num_agebrackets:.0f}.csv"
        else:
            file_name = f"{country}_{level}_{location}_age_distribution_{num_agebrackets:.0f}.csv"
        file_path = os.path.join(datadir, 'age_distributions', file_name)
        df.loc[df['file_name'] == file_name, 'year_of_source_data'] = source_year

    return df


if __name__ == '__main__':

    # read in the dataframe for source year
    df_path = os.path.join(datadir, 'age_distribution_source_year.csv')
    df = pd.read_csv(df_path, delimiter=',')
    df.sort_values(by='file_name', inplace=True)

    # update source year for all locations

    country = 'Australia'
    locations = get_locations_by_country(locations_df, country)
    locations.append(country)
    source_year = 2011
    for num_agebrackets in [18, 85]:
        df = update_age_source_year(df, locations, num_agebrackets, source_year)

    country = 'Canada'
    locations = get_locations_by_country(locations_df, country)
    locations.append(country)
    source_year = 2011
    for num_agebrackets in [18, 85]:
        df = update_age_source_year(df, locations, num_agebrackets, source_year)

    country = 'China'
    locations = get_locations_by_country(locations_df, country)
    locations.append(country)
    source_year = 2010
    for num_agebrackets in [18, 85]:
        df = update_age_source_year(df, locations, num_agebrackets, source_year)

    country = 'Europe'
    locations = get_locations_by_country(locations_df, country)
    source_year = 2011
    for num_agebrackets in [18, 85]:
        df = update_age_source_year(df, locations, num_agebrackets, source_year)

    country = 'India'
    locations = get_locations_by_country(locations_df, country)
    locations.append(country)
    source_year = 2011
    for num_agebrackets in [18, 85]:
        df = update_age_source_year(df, locations, num_agebrackets, source_year)

    country = 'Israel'
    locations = get_locations_by_country(locations_df, country)
    locations.append(country)
    source_year = 2008
    for num_agebrackets in [18, 85]:
        df = update_age_source_year(df, locations, num_agebrackets, source_year)

    country = 'Japan'
    locations = get_locations_by_country(locations_df, country)
    locations.append(country)
    source_year = 2010
    for num_agebrackets in [18, 85]:
        df = update_age_source_year(df, locations, num_agebrackets, source_year)

    country = 'Russia'
    locations = get_locations_by_country(locations_df, country)
    source_year = 2010
    for num_agebrackets in [18, 85]:
        df = update_age_source_year(df, locations, num_agebrackets, source_year)

    country = 'South_Africa'
    locations = get_locations_by_country(locations_df, country)
    locations.append(country)
    source_year = 2011
    for num_agebrackets in [18, 85]:
        df = update_age_source_year(df, locations, num_agebrackets, source_year)

    country = 'United_States'
    locations = get_locations_by_country(locations_df, country)
    locations.append(country)
    source_year = 2010
    for num_agebrackets in [18, 85]:
        df = update_age_source_year(df, locations, num_agebrackets, source_year)

    # write the updated dataframe to file - only need to run this once
    write_flag = True
    if write_flag:
        df.to_csv(df_path, index=False)
