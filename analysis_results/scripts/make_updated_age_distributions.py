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

popsizes_file = os.path.join(datadir, 'population_size.csv')
popsizes_df = pd.read_csv(popsizes_file, delimiter=',')


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


def get_popsize(popsizes_df, location):
    """
    Return the population size for the location (according to the year of the original data - either 2010 or 2011).

    Args:
        popsizes_df (pandas Dataframe) : dataframe containing population sizes for each location, the year of the count, and the source
        location (str) : name of the location

    Returns:
        int: An integer of the estimated population size for the location.
    """
    return popsizes_df.loc[popsizes_df['location'] == location]['population_size'].values[0]


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


def get_rescaled_ages(location, popsize):
    country = get_country_name(locations_df, location)
    level = 'subnational'
    num_agebrackets = 85
    ages = get_ages(location, country, level, num_agebrackets)
    age_dist = ages.copy()

    for a in ages:
        age_dist[a] = ages[a] / sum(ages.values())

    original_popsize = sum(ages.values())
    ratio = popsize / original_popsize
    age_range = np.arange(len(ages))

    rescaled_ages = dict()
    for a in ages:
        rescaled_ages[a] = np.round(ages[a] * ratio)

    diff = int(popsize - sum(rescaled_ages.values()))

    if diff > 0:
        for n in range(diff):
            a = np.random.choice(age_range, p=[age_dist[a] for a in age_range])
            rescaled_ages[a] += 1

    elif diff < 0:
        for n in range(-diff):
            a = np.random.choice(age_range, p=[age_dist[a] for a in age_range])
            while rescaled_ages[a] <= 0:
                a = np.random.choice(age_range, p=[age_dist[a] for a in age_range])
            rescaled_ages[a] -= 1

    assert sum(rescaled_ages.values()) == popsize, f'New age distribution does not match the size of popsize that was desired. The difference is {sum(rescaled_ages.values()) - popsize} (sum(new age distribution) - popsize).'
    print(f'{location}, {country}, {ratio:0.3f}, {original_popsize}, {popsize}, {sum(rescaled_ages.values()) - popsize}')
    return rescaled_ages


def write_rescaled_ages(location, rescaled_ages, datadir, overwrite=False):
    num_agebrackets = len(rescaled_ages)
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
    file_path = os.path.join(datadir, 'population_rescaled_age_distributions', file_name)
    e = os.path.exists(file_path)
    if e:
        print(f'{file_path} already exists.')
        if overwrite:
            print('Overwriting the file.')
            f = open(file_path, 'w+')
            for a in rescaled_ages:
                f.write(f"{a:.16f},{rescaled_ages[a]:.16f}\n")
            f.close()
        else:
            print('Not overwriting the existing file.')
    else:
        f = open(file_path, 'w+')
        for a in rescaled_ages:
            f.write(f"{a:.16f},{rescaled_ages[a]:.16f}\n")
        f.close()
    return


def write_rescaled_aggregated_ages(location, aggregate_ages, datadir, overwrite=False):
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
    file_path = os.path.join(datadir, 'population_rescaled_age_distributions', file_name)
    e = os.path.exists(file_path)
    if e:
        print(f'{file_path} already exists.')
        if overwrite:
            print('Overwriting the file.')
            f = open(file_path, 'w+')
            for a in range(len(aggregate_ages)):
                f.write(f"{a:.16f},{aggregate_ages[a]:.16f}\n")
            f.close()
        else:
            print('Not overwriting the existing file.')
    else:
        f = open(file_path, 'w+')
        for a in range(len(aggregate_ages)):
            f.write(f"{a:.16f},{aggregate_ages[a]:.16f}\n")
        f.close()
    return



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

    # create population size rescaled ages with sizes
    write_flag = False  # set to True to save matrices to disk
    for location in locations:
        popsize = get_popsize(popsizes_df, location)
        rescaled_ages = get_rescaled_ages(location, popsize)

        aggregate_ages = get_aggregate_ages(rescaled_ages, age_by_brackets_dic)
        assert sum(rescaled_ages.values()) == sum(aggregate_ages.values())
        # write to disk
        if write_flag:
            write_rescaled_ages(location, rescaled_ages, datadir)
            write_rescaled_aggregated_ages(location, aggregate_ages, datadir)
