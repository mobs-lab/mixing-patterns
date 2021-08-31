"""
How to re-aggregate the contact matrices to different age group or bracket
mappings.

The original contact matrices we offer in this repo and report in

Mistry, D., Litvinova, M., Pastore y Piontti, A. et al. Inferring
high-resolution human mixing patterns for disease modeling. Nat Commun 12, 323
(2021). https://doi.org/10.1038/s41467-020-20544-y

are aggregated to 85 age brackets. Here, we show how to work with these matrices
to reconstruct the symmetric matrices, aggregate to a different grouping of ages
and calculate the new asymmetric contact matrices. For more details on the methods
shown here, we refer you to the Methods section of our paper.

Note: In order to use the new matrices with a different age grouping, you'll need
an age distribution that is also mapped in the same grouping. Here, we only provide
the age distribution for places mapped to the original 85 brackets and the ages
aggregated to 18 brackets as an example. The methods below should show you how to
write the matrices and age distribution to different age groupings.

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

setting_fractions_file = os.path.join(analysisdir, 'summary_fractions_by_location.csv')
setting_fractions_df = pd.read_csv(setting_fractions_file)
setting_codes = ['H', 'S', 'W', 'R']


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
        A numpy matrix of contact.
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


def get_ages(location, country, level, num_agebrackets=85):
    """
    Get the age count for the synthetic population of the location.

    Args:
        location (str)        : name of the location
        country (str)         : name of the country
        level (str)           : name of level (country or subnational)
        num_agebrackets (int) : the number of age brackets

    Returns:
        dict: A dictionary of the age count.
    """

    if country == 'Europe':
        country = location
        level = 'country'

    if level == 'country':
        file_name = country + '_' + level + '_level_age_distribution_' + '%i' % num_agebrackets + '.csv'
    else:
        file_name = country + '_' + level + '_' + location + '_age_distribution_' + '%i' % num_agebrackets + '.csv'
    file_path = os.path.join(datadir, 'age_distributions', file_name)
    df = pd.read_csv(file_path, delimiter=',', header=None)
    df.columns = ['age', 'age_count']
    ages = dict(zip(df.age.values.astype(int), df.age_count.values))
    return ages


def get_available_age_brackets_and_mapping():
    """
    Create a mapping of the number of age brackets to the brackets and for each
    set of brackets, a mapping of age to bracket or group.
    """
    brackets_dic = {}
    dict_of_age_by_brackets = {}

    for num_agebrackets in [85, 20, 18, 15, 12]:
        brackets = []
        if num_agebrackets == 85:
            for i in range(84):
                brackets.append([i])
            brackets.append(np.arange(84, 101))
        # # num_agebracket = 20 only works if you have an age distribution and
        # # matrices that go in more detail than age 84+, so if you wanted 
        # # brackets of 85-89, 90-94, 95-100+, etc. it would be hard unless 
        # # you have those matrices (which we don't because of the European 
        # # matrices for Hungary and Portugal)

        # if num_agebrackets == 20:
        #     for i in range(19):
        #         brackets.append(np.arange(5 * i, 5 * (i + 1)))
        #     brackets.append(np.arange(95, 101))
        if num_agebrackets == 18:
            for i in range(16):
                brackets.append(np.arange(5 * i, 5 * (i + 1)))
            brackets.append(np.arange(80, 84))
            brackets.append(np.arange(84, 101))
        if num_agebrackets == 15:
            for i in range(14):
                brackets.append(np.arange(5 * i, 5 * (i + 1)))
            brackets.append(np.arange(70, 101))
        if num_agebrackets == 12:
            for i in range(11):
                brackets.append(np.arange(5 * i, 5 * (i + 1)))
            brackets.append(np.arange(55, 101))

        age_by_brackets_dic = dict.fromkeys(np.arange(101), 0)
        for n, b in enumerate(brackets):
            for a in b:
                age_by_brackets_dic[a] = n

        brackets_dic[num_agebrackets] = brackets
        dict_of_age_by_brackets[num_agebrackets] = age_by_brackets_dic

    return brackets_dic, dict_of_age_by_brackets


def get_age_brackets(available_age_brackets, age_by_brackets_mapping, num_agebrackets):
    """
    Return the age brackets and mapping of age to bracket for a specific mapping indexed by the number of brackets.
    """
    return available_age_brackets[num_agebrackets], age_by_brackets_mapping[num_agebrackets]


def get_aggregate_ages(ages, age_by_brackets_dic):
    num_agebrackets = len(set(age_by_brackets_dic.values()))
    aggregate_ages = dict.fromkeys(np.arange(num_agebrackets), 0)
    for a in ages:
        b = age_by_brackets_dic[a]
        aggregate_ages[b] += ages[a]
    return aggregate_ages


def get_aggregate_matrix(M, age_by_brackets_dic):
    N = len(M)
    num_agebrackets = len(set(age_by_brackets_dic.values()))
    M_agg = np.zeros((num_agebrackets, num_agebrackets))
    for i in range(N):
        bi = age_by_brackets_dic[i]
        for j in range(N):
            bj = age_by_brackets_dic[j]
            M_agg[bi][bj] += M[i][j]
    return M_agg


def get_asymmetric_matrix(symmetric_matrix, aggregate_ages):
    M = copy.deepcopy(symmetric_matrix)
    for a in aggregate_ages:
        M[a, :] = M[a, :] / float(aggregate_ages[a])
    return M


def recalculate_symmetric_matrix(asymmetric_matrix, ages):
    symmetric_matrix = asymmetric_matrix.copy()
    for a in ages:
        symmetric_matrix[a, :] = symmetric_matrix[a, :] * ages[a]
    return symmetric_matrix


def get_aggregate_symmetric_community_matrices(aggregate_ages):
    num_agebrackets = len(aggregate_ages)
    M = np.ones((num_agebrackets, num_agebrackets))
    for a in aggregate_ages:
        M[a, :] *= float(aggregate_ages[a])
        M[:, a] *= float(aggregate_ages[a])
    return M


def get_new_aggregate_matrices(location, num_agebrackets, available_age_brackets, age_by_brackets_mapping):
    """
    For a location, recombine the contact matrices to a new mapping of age
    brackets or groups. The last age bracket must include all ages from 84 on
    up.

    """
    country = get_country_name(locations_df, location)
    if location == 'Israel':
        ages = get_ages(location, country, 'country')
    else:
        ages = get_ages(location, country, 'subnational')
    aggregate_ages = get_aggregate_ages(ages, age_by_brackets_dic)

    new_matrices = dict()

    for layer in ['household', 'school', 'work']:
        if location == 'Israel':
            matrix = read_contact_matrix(location, country, 'country', layer)
        else:
            matrix = read_contact_matrix(location, country, 'subnational', layer)
        symmetric_matrix = recalculate_symmetric_matrix(matrix, ages)
        assert np.allclose(symmetric_matrix, symmetric_matrix.T), f"Check failed for {location}, {layer}"
        aggregate_symmetric_matrix = get_aggregate_matrix(symmetric_matrix, age_by_brackets_dic)
        aggregate_matrix = get_asymmetric_matrix(aggregate_symmetric_matrix, aggregate_ages)
        new_matrices[layer] = aggregate_matrix

    if location == 'Israel':
        og_matrix = read_contact_matrix(location, country, 'country', 'community')
    else:
        og_matrix = read_contact_matrix(location, country, 'subnational', 'community')
    # community matrix
    aggregate_symmetric_matrix = get_aggregate_symmetric_community_matrices(aggregate_ages)
    new_matrices['community'] = get_asymmetric_matrix(aggregate_symmetric_matrix, aggregate_ages)
    for a in aggregate_ages:
        new_matrices['community'][a, :] = new_matrices['community'][a, :] / new_matrices['community'][a, :].sum()

    return new_matrices


def write_new_aggregated_matrices(location, new_matrices, datadir):
    num_agebrackets = len(new_matrices[list(new_matrices.keys())[0]])
    country = get_country_name(locations_df, location)

    if country == 'Europe':
        country = location
        level = 'country'
    if country == location:
        level = 'country'
    else:
        level = 'subnational'

    for layer in ['household', 'school', 'work', 'community']:
        setting = layer
        setting_type, setting_suffix = 'F', 'setting'
        if setting == 'overall':
            setting_type, setting_suffix = 'M', 'contact_matrix'

        if level == 'country':
            file_name = f"{country}_{level}_level_{setting_type}_{setting}_{setting_suffix}_{num_agebrackets:.0f}.csv"
        else:
            file_name = f"{country}_{level}_{location}_{setting_type}_{setting}_{setting_suffix}_{num_agebrackets:.0f}.csv"
        file_path = os.path.join(datadir, 'contact_matrices', file_name)
        print(location, layer, new_matrices[setting].sum() / num_agebrackets)
        np.savetxt(file_path, new_matrices[setting], fmt='%.16f', delimiter=',')
    return


def write_new_aggregated_ages(location, aggregate_ages, datadir):
    num_agebrackets = len(aggregate_ages)
    country = get_country_name(locations_df, location)

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
    # print(file_path)
    f = open(file_path, 'w+')
    for a in range(len(aggregate_ages)):
        f.write(f"{a:.16f},{aggregate_ages[a]:.16f}\n")
    return


if __name__ == '__main__':

    available_age_brackets, age_by_brackets_mapping = get_available_age_brackets_and_mapping()

    num_agebrackets = 18
    age_brackets = available_age_brackets[num_agebrackets]
    age_by_brackets_dic = age_by_brackets_mapping[num_agebrackets]

    countries = ['Australia', 'Canada', 'China', 'Europe', 'India', 'Israel', 'Japan', 'Russia', 'South_Africa', 'United_States']
    locations = []

    for country in countries:
        locations += get_locations_by_country(locations_df, country)
        if country in locations and country != 'Israel':
            locations.remove(country)

        if country == 'India':
            locations.remove('Dadra_and_Nagar_Haveli')
            locations.remove('Chandigarh')
            locations.remove('Lakshadweep')

    if 'Russian_Federation' in locations:
        locations.remove('Russian_Federation')
    # locations = locations[0:1]

    write_flag = True
    # write_flag = False

    for li, location in enumerate(locations):
        country = get_country_name(locations_df, location)
        if location == 'Israel':
            ages = get_ages(location, country, 'country')
        else:
            ages = get_ages(location, country, 'subnational')
        aggregate_ages = get_aggregate_ages(ages, age_by_brackets_dic)

        new_matrices = get_new_aggregate_matrices(location, num_agebrackets, available_age_brackets, age_by_brackets_mapping)
        print(location, country)
        if write_flag:
            write_new_aggregated_matrices(location, new_matrices, datadir)
            write_new_aggregated_ages(location, aggregate_ages, datadir)
