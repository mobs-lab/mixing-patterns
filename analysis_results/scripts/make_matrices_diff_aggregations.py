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
write the matrices and age distribution to different age groupings. Also provided
are age mappings for 15 and 12 age brackets.

"""
from glob import glob as ls
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

location_file = os.path.join(datadir, 'location_to_country.csv')
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
    if country == location:
        level = 'country'
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
    if country == location:
        level = 'country'

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

    for num_agebrackets in [85, 18, 15, 12]:
        brackets = []
        if num_agebrackets == 85:
            for i in range(84):
                brackets.append([i])
            brackets.append(np.arange(84, 101))
        # # num_agebracket = 20 only works if you have an age distribution and
        # # matrices that go in more detail than age 84+, so if you wanted
        # # brackets of 85-89, 90-94, 95-100+, etc. it would be hard unless
        # # you have those matrices (which we don't because of the European
        # # matrices)

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

    Args:
        available_age_brackets (dict)  : a dictionary of available mappings
        age_by_brackets_mapping (dict) : a dictionary of mappings for single years of age to the bracket or bin they map to depending on the number of age brackets defined
        num_agebrackets (int)          : number of age brackets desired

    Returns:
        The age brackets and mapping of age to bracket for a specific mapping
        indexed by the number of brackets.
    """
    return available_age_brackets[num_agebrackets], age_by_brackets_mapping[num_agebrackets]


def get_aggregate_ages(ages, age_by_brackets_dic):
    """
    Return an aggregated age count distribution.

    Args:
        ages (dict)                : original age count distribution
        age_by_brackets_dic (dict) : a dictionary mapping single years of age to the bracket or bin they map to

    Returns:
        An aggregated age count distribution from the mapping age_by_brackets_dic.
    """
    num_agebrackets = len(set(age_by_brackets_dic.values()))
    aggregate_ages = dict.fromkeys(np.arange(num_agebrackets), 0)
    for a in ages:
        b = age_by_brackets_dic[a]
        aggregate_ages[b] += ages[a]
    return aggregate_ages


def get_aggregate_matrix(M, age_by_brackets_dic):
    """
    Return an aggregated age mixing contact matrix.

    Args:
        M (np.array): original age mixing contact matrix
        age_by_brackets_dic (dict) : a dictionary mapping single years of age to the bracket or bin they map to

    Returns:
        np.array: An aggregated age mixing contact matrix using the age_by_brackets_dic mapping.
    """
    N = len(M)
    num_agebrackets = len(set(age_by_brackets_dic.values()))
    assert N >= num_agebrackets, 'The number of age brackets to aggregate to the matrix to is greater than the original matrix. We cannot aggregate the contact matrix with this mapping.'
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
    # each individual can interact with others in their own age group, but we shouldn't count self interactions
    for a in aggregate_ages:
        M[a, a] -= aggregate_ages[a]
    return M


def get_new_aggregate_matrices(location, num_agebrackets, available_age_brackets, age_by_brackets_mapping, verbose=False):
    """
    For a location, recombine the contact matrices to a new mapping of age
    brackets or groups. The last age bracket must include all ages from 84 on
    up.

    """
    country = get_country_name(locations_df, location)
    if location == 'Israel' or location == country:
        ages = get_ages(location, country, 'country')
    else:
        ages = get_ages(location, country, 'subnational')
    aggregate_ages = get_aggregate_ages(ages, age_by_brackets_dic)

    new_matrices = dict()
    for layer in ['household', 'school', 'work']:
        if location == 'Israel' or location == country:
            matrix = read_contact_matrix(location, country, 'country', layer)
        else:
            matrix = read_contact_matrix(location, country, 'subnational', layer)
        symmetric_matrix = recalculate_symmetric_matrix(matrix, ages)
        assert np.allclose(symmetric_matrix, symmetric_matrix.T), f"Check failed for {location}, {layer}"
        aggregate_symmetric_matrix = get_aggregate_matrix(symmetric_matrix, age_by_brackets_dic)
        aggregate_matrix = get_asymmetric_matrix(aggregate_symmetric_matrix, aggregate_ages)
        new_matrices[layer] = aggregate_matrix

    if location == 'Israel' or location == country:
        og_matrix = read_contact_matrix(location, country, 'country', 'community')
    else:
        og_matrix = read_contact_matrix(location, country, 'subnational', 'community')
    # community matrix
    aggregate_symmetric_matrix = get_aggregate_symmetric_community_matrices(aggregate_ages)
    new_matrices['community'] = get_asymmetric_matrix(aggregate_symmetric_matrix, aggregate_ages)
    for a in aggregate_ages:
        new_matrices['community'][a, :] = new_matrices['community'][a, :] / new_matrices['community'][a, :].sum()

    if verbose:
        print(location)
        for layer in ['household', 'school', 'work', 'community']:
            print(layer, new_matrices[layer].sum() / num_agebrackets)
    return new_matrices


def write_new_aggregated_matrices(location, new_matrices, datadir, overwrite=False):
    num_agebrackets = len(new_matrices[list(new_matrices.keys())[0]])
    country = get_country_name(locations_df, location)

    if country == 'Europe':
        country = location
        level = 'country'
    if country == location:
        level = 'country'
    else:
        level = 'subnational'

    for layer in new_matrices.keys():
        setting = layer
        setting_type, setting_suffix = 'F', 'setting'
        if setting == 'overall':
            setting_type, setting_suffix = 'M', 'contact_matrix'

        if level == 'country':
            file_name = f"{country}_{level}_level_{setting_type}_{setting}_{setting_suffix}_{num_agebrackets:.0f}.csv"
        else:
            file_name = f"{country}_{level}_{location}_{setting_type}_{setting}_{setting_suffix}_{num_agebrackets:.0f}.csv"
        file_path = os.path.join(datadir, 'contact_matrices', file_name)
        e = os.path.exists(file_path)
        if e:
            print(f'{file_path} already exists.')
            if overwrite:
                print('Overwriting the file.')
                np.savetxt(file_path, new_matrices[setting], fmt='%.16f', delimiter=',')
            else:
                print('Not overwriting the existing file.')
        else:
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
    f = open(file_path, 'w+')
    for a in range(len(aggregate_ages)):
        f.write(f"{a:.16f},{aggregate_ages[a]:.16f}\n")
    f.close()
    return


def check_matrix_files_delimiter(num_agebrackets=85):
    matrix_files = ls(os.path.join(datadir, 'contact_matrices', f'*{num_agebrackets:.0f}*'))
    matrix_files = sorted(matrix_files)
    s = False
    for file_path in matrix_files:
        try:
            M = np.loadtxt(file_path, delimiter=',')
        except:
            print(f'{file_path} has a different delimiter than the expected ",".')
            s = True
    if s:
        print('Found at least one contact matrix using a different delimiter than the expected ",".')


def check_age_distribution_count(locations, num_agebrackets=18):
    for location in locations:
        country = get_country_name(locations_df, location)
        level = 'subnational'
        ages_85 = get_ages(location, country, level, num_agebrackets=85)
        ages = get_ages(location, country, level, num_agebrackets=num_agebrackets)
        # print(location, country)
        assert sum(ages.values()) == sum(ages_85.values()), f'Sum of aggregated ages ({num_agebrackets} age brackets) does not match original age count (85 age brackets).'


def combine_matrices(location, num_agebrackets=85, weights=None, verbose=False):
    country = get_country_name(locations_df, location)
    level = 'subnational'

    matrices = dict()
    matrices['overall'] = np.zeros((num_agebrackets, num_agebrackets))

    # use default weights if none provided - see paper for details on how these were derived
    if weights is None:
        weights = dict()
        weights['household'] = 4.1103
        weights['school'] = 11.4069
        weights['work'] = 8.0746
        weights['community'] = 2.7868

    elif not isinstance(weights, dict):
        raise ValueError('weights must be a dictionary.')

    for l in ['household', 'school', 'work', 'community']:
        matrices[l] = read_contact_matrix(location, country, level, l, num_agebrackets)
        matrices['overall'] += matrices[l] * weights[l]

    if verbose:
        print(location, country, matrices['overall'].sum() / num_agebrackets)

    for l in ['household', 'school', 'work', 'community']:
        matrices.pop(l)

    return matrices


if __name__ == '__main__':

    # get all of the locations possible from the countries/regions we provide data for
    countries = ['Australia', 'Canada', 'China', 'Europe', 'India', 'Israel', 'Japan', 'Russia', 'South_Africa', 'United_States']
    locations = []

    for country in countries:
        locations += get_locations_by_country(locations_df, country)
        if country not in locations and country not in ['Russia', 'Europe']:
            locations.append(country)

        # a few territories in India we did not model
        if country == 'India':
            locations.remove('Dadra_and_Nagar_Haveli')
            locations.remove('Chandigarh')
            locations.remove('Lakshadweep')

    # find out what age brackets and mappings are available -- if you want to use different ones than available, look at the details of this method to add them
    available_age_brackets, age_by_brackets_mapping = get_available_age_brackets_and_mapping()

    # let's create contact matrices aggregated to 18 age brackets: 5 year bins, where the last bin is 84 years old and up
    num_agebrackets = 18
    age_brackets = available_age_brackets[num_agebrackets]
    age_by_brackets_dic = age_by_brackets_mapping[num_agebrackets]

    # first a check: all the original matrix files for 85 age brackets have ',' as the delimiter
    check_matrix_files_delimiter(85)

    # next: write the setting matrices aggregated to num_agebrackets
    write_flag = False  # set to True to save matrices to disk
    if write_flag:
        for location in locations:
            country = get_country_name(locations_df, location)
            if location == 'Israel' or location == country:
                ages = get_ages(location, country, 'country')
            else:
                ages = get_ages(location, country, 'subnational')
            aggregate_ages = get_aggregate_ages(ages, age_by_brackets_dic)

            new_matrices = get_new_aggregate_matrices(location, num_agebrackets, available_age_brackets, age_by_brackets_mapping)
            write_new_aggregated_matrices(location, new_matrices, datadir, overwrite)
            write_new_aggregated_ages(location, aggregate_ages, datadir)

    # check those matrices use the correct delimiter of ','
    check_matrix_files_delimiter(num_agebrackets=num_agebrackets)  # this will only work if you've already run the code with write_flag = True before or at the same time
    # check the sum of the aggregated age distribution count for num_agebrackets matches the original count
    check_age_distribution_count(locations, num_agebrackets=num_agebrackets)

    # write the overall matrix for the new aggregation
    write_flag = False  # set to True to save matrices to disk
    if write_flag:
        for location in locations:
            matrix = combine_matrices(location, num_agebrackets=num_agebrackets)  # just returns a new overall matrix (weighted linear combination) for num_agebrackets
            write_new_aggregated_matrices(location, matrix, datadir)
