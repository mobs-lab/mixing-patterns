"""Plot figure 6: attack rates vs R0."""
import numpy as np
import pandas as pd
from scipy.stats import linregress

import matplotlib as mplt
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.font_manager as font_manager

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

location_file = os.path.join(analysisdir, 'location_country_names.csv')
locations_df = pd.read_csv(location_file, delimiter = ',')

setting_fractions_file = os.path.join(analysisdir, 'summary_fractions_by_location.csv')
setting_fractions_df = pd.read_csv(setting_fractions_file)
setting_codes = ['H','S','W','R']


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
        str: Name of the country the location is in.
    """
    locations = list(df[df.country == country].location.values)
    return locations


def get_fractions(setting_fractions_df, location):
    """
    Get the fraction of people with contacts in each setting of household (H), school (S), or work (W).

    Args:
        setting_fractions_df (pandas DataFrame) : a dataframe
        location (str)                          : name of the location

    Returns:
        dict: A dictionary of fractions of people with contacts in each setting for the location.
    """
    fractions = dict.fromkeys(['H','S','W','R'],1.)
    d = setting_fractions_df[setting_fractions_df.location == location]
    fractions['H'] = d.NhN.values[0]
    fractions['S'] = d.NsN.values[0]
    fractions['W'] = d.NwN.values[0]
    return fractions


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


def get_average_age(ages):
    """
    Get the average age from a dictionary of age counts.

    Args:
        ages (dict): dictionary of age counts

    Return:
        float: The average age given the age count.
    """
    average_age = 0
    total_population = sum(ages.values())
    for a in ages:
        average_age += ages[a] * a
    average_age = average_age/total_population
    return average_age


def get_school_age_distribution(location):
    """
    Get the age count of people active in the school setting.

    Args:
        location (str): name of the location

    Returns:
        dict: Age count of people active in the school setting.
    """
    ages = {}
    file_path = os.path.join(analysisdir, 'schools_age_distributions',
                'schools_age_distributions_' + location + '.dat')
    df = pd.read_csv(file_path, delimiter = ',')
    ages = dict(zip(df.age.values, df.setting_count.values))
    return ages


def get_percent_in_school(ages, school_ages):
    """
    Get the percent of people in school.

    Args: 
        ages (dict)        : age count
        school_ages (dict) : school age count

    Returns:
        float: The percent of people in the school setting.
    """
    total_in_school = np.sum([v for v in school_ages.values()], dtype = float)
    total_population = np.sum([v for v in ages.values()], dtype = float)
    return total_in_school/total_population * 100


def get_attack_rates_df(reference_location, reference_scenario, beta, susceptibility_drop_factor, gamma_inverse, num_agebrackets):
    """
    Get attack rates dataframe for an SIR compartmental model with age specific contact patterns.

    Args:
        reference_location (str)           : name of reference location or locations
        reference_scenario (str)           : specific reference scenario
        beta (float)                       : the transmissibilty
        susceptibility_drop_factor (float) : susceptibility of adults to those under 18
        gamma_inverse (float)              : the mean recovery period
        num_agebrackets (int)              : the number of age brackets for the matrix
    
    Returns:
        Pandas dataframe of attack rates by location
    """

    file_path = os.path.join(analysisdir, 'reference_location_' + reference_location, 'mcmc_beta_and_dropfactor_scenario_' + reference_scenario,
                             'all_locations_attack_rates_by_age_reference_scenario_' + reference_scenario + '_beta_' + '%.2f' % beta + '_susceptibility_' + '%.2f' % susceptibility + '_gamma_inverse_' + '%.1f' % gamma_inverse + '_' + str(num_agebrackets) + '.csv')
    return pd.read_csv(file_path)


def get_attack_rate(df, location):
    """
    Get the total attack rate for the location from the dataframe for an SIR compartmental model with age specific contact patterns.

    Args:
        df (pd.DataFrame) : a dataframe of attack rates by location
        location (str)    : name of the location

    Returns:
        float: The total attack rate for the location as a fraction from an SIR compartmental model with age specific contact patterns. Values between 0 and 1.
    """
    return df.loc[df['location'] == location]['artotal'].values[0]


def get_homogeneous_attack_rate_df(gamma_inverse, num_agebrackets):
    """
    Get a dataframe with the attack rate for an SIR model with the homogeneous mixing assumption for different basic reproduction, R0, values with a given average recovery period.

    Args:
        gamma_inverse (float): the mean recovery period
        num_agebrackets (int): the number of age brackets for the matrix

    Returns:
        Pandas dataframe of attack rates by R0 value
    """
    file_path = os.path.join(analysisdir, 'homogeneous_sir_attack_rates', 'attack_rates_SIR_homogeneous_mixing.csv')
    df = pd.read_csv(file_path)
    return df


def get_homogeneous_attack_rate(df, R0_star):
    """
    Get the attack rate for an SIR model with the homogeneous mixing assumption for a given basic reproduction, R0_star, value.

    Args:
        df (pd.DataFrame) : a dataframe of attack rates by the basic reproduction number
        R0_star (float)   : the basic reproduction number

    Returns:
        float: The total attack rate as a fraction from an SIR model with homogeneous mixing assumptions and specified basic reproduction number R0_star. Values between 0 and 1.
    """
    return df.loc[df['R0'] == R0_star]['attack_rate'].values[0]


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


def get_R0(beta, gamma_inverse, matrix):
    """
    Get the basic reproduction number, R0, for an SIR compartmental model given the basic reproduction number, the mean recovery period, and the contact matrix.
    
    Args:
        beta (float)          : the transmissibility beta
        gamma_inverse (float) : the mean recovery period
        matrix (np.ndarray)   : the contact matrix

    Returns:
        float: The basic reproduction number R0 for an SIR compartmental model.
    """

    gamma = 1./gamma_inverse
    eigenvalue = get_eigenvalue(matrix)
    return beta * eigenvalue / gamma


def get_beta(R0, gamma_inverse, matrix):
    """
    Get the transmissibility from an SIR model with age specific contact patterns and basic reproduction number.

    Args:
        R0_star (float)       : the basic reproduction number
        gamma_inverse (float) : the mean recovery period
        matrix (np.ndarray)   : the age specific contact matrix

    Returns:
        float: The transmissibility from an SIR model with age specific contact patterns.
    """
    gamma = float(1)/gamma_inverse
    eigenvalue = get_eigenvalue(matrix)
    return R0 * gamma / eigenvalue


def linear_function(x,m,b):
    """
    Get the y value of a linear function given the x value.

    Args:
        m (float): the slope
        b (float): the intercept

    Returns:
        The expected y value from a linear function at some specified x value.
    """
    return m*x + b


def plot_fig(countries, reference_location, reference_scenario, beta, susceptibility_drop_factor, gamma_inverse, num_agebrackets):
    """
    Plot the attack rates from an SIR model with age specific contact patterns vs the basic reproduction, the average age, 
    and the percent of the population with contacts in the school layer for subnational locations.

    Args:
        countries (list)                   : list of countries
        reference_location (str)           : the reference location (or set of locations) the transmissibilty is calibrated to
        reference_scenario (str)           : label of the reference scenario
        beta (float)                       : the transmissibility
        susceptibility_drop_factor (float) : susceptibility of adults to those under 18
        gamma_inverse (float)              : the mean recovery period
        num_agebrackets (int)              : the number of age brackets

    Returns:
        Matplotlib figure.
    """
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

    hdf = get_homogeneous_attack_rate_df(gamma_inverse, num_agebrackets)
    df = get_attack_rates_df(reference_location, reference_scenario, beta, susceptibility_drop_factor, gamma_inverse, num_agebrackets)

    width = 16
    height = 5
    left = 0.06
    right = 0.865
    bottom = 0.16
    top = 0.88
    wspace = 0.32

    fig, ax = plt.subplots(1, 3, figsize=(width, height))
    fig.subplots_adjust(left = left, right = right, top = top, bottom = bottom, wspace = wspace)


    leg_left = right + 0.01
    leg_right = 0.985
    leg_bottom = bottom
    leg_top = top
    leg_width = leg_right - leg_left
    leg_height = leg_top - leg_bottom

    fontsize = 20

    axleg = fig.add_axes([leg_left, leg_bottom, leg_width, leg_height])
    axleg.axis('off')

    attack_rates_list = []
    R0_list = []
    average_age_list = []
    percent_in_school_list = []
    color_list = []

    beta_drop_age = 18

    for n, location in enumerate(locations):
        country = get_country_name(locations_df, location)
        if location == 'Israel':
            matrix = read_contact_matrix(location, country, 'country', 'overall')
            ages = get_ages(location, country, 'country')
        else:
            matrix = read_contact_matrix(location, country, 'subnational', 'overall')
            ages = get_ages(location, country, 'subnational')

        R0 = get_R0(beta, gamma_inverse, matrix)
        R0_list.append(R0)

        average_age = get_average_age(ages)
        average_age_list.append(average_age)

        if country != 'Europe':
            school_ages = get_school_age_distribution(location)
            percent_in_school = get_percent_in_school(ages, school_ages)

        else:
            fractions = get_fractions(setting_fractions_df, location)
            percent_in_school = fractions['S'] * 100

        percent_in_school_list.append(percent_in_school)

        ar = get_attack_rate(df, location) * 100
        attack_rates_list.append(ar)
        color = country_colour_dic[country]
        color_list.append(color)

    homogeneous_attack_rates_list = hdf.attack_rate.values * 100
    homogeneous_R0_list = hdf.R0.values

    size = 12

    leg = []

    ax[0].scatter(R0_list, attack_rates_list, marker='o', color=color_list, s=size)
    ax[0].plot(homogeneous_R0_list, homogeneous_attack_rates_list, color='k', lw=1.5, label='Homogenous \nmixing model')
    leg.append(ax[0].legend(loc=2, fontsize=17))
    ax[0].set_xlim(1.5, 2.0)
    ax[0].set_xticks(np.arange(1.5, 2.01, 0.1))
    ax[0].set_xlabel(r'$R_0$', fontsize = fontsize)
    ax[0].text(1.39, 50, 'a', fontsize=fontsize+24, fontstyle='oblique')


    m,b,r,p,std_err = linregress(average_age_list,attack_rates_list)
    y_theory = np.array([linear_function(a,m,b) for a in average_age_list])
    ax[1].plot(average_age_list, y_theory, color = 'k', lw = 1.5)
    ax[1].scatter(average_age_list, attack_rates_list, marker='o', color=color_list, s=size)
    ax[1].text(45,73, r'$\rho$ = ' + '%.2f' % r , fontsize = 18, verticalalignment = 'center', horizontalalignment = 'center', color = 'k')
    ax[1].set_xlim(20, 55)
    ax[1].set_xticks(np.arange(20, 51, 10))
    leg.append(ax[1].legend(loc = 7, fontsize = 16, ncol = 1))
    ax[1].set_xlabel('Average Age', fontsize = fontsize)
    ax[1].text(13.5, 50, 'b', fontsize=fontsize+24, fontstyle='oblique')


    m,b,r,p,std_err = linregress(percent_in_school_list,attack_rates_list)
    y_theory = np.array([linear_function(a,m,b) for a in percent_in_school_list])
    ax[2].plot(percent_in_school_list, y_theory, color = 'k', lw = 1.5)
    ax[2].scatter(percent_in_school_list, attack_rates_list, marker='o', color=color_list, s=size)
    ax[2].text(20,73, r'$\rho$ = ' + '%.2f' % r , fontsize = 18, verticalalignment = 'center', horizontalalignment = 'center', color = 'k')
    ax[2].set_xlim(10, 45)
    ax[2].set_xticks(np.arange(10, 50, 10))
    ax[2].set_xlabel('% In Educational Institutions', fontsize = fontsize)
    ax[2].text(3.5, 50, 'c', fontsize=fontsize+24, fontstyle='oblique')


    for country in ['Australia','Canada','China','Europe','India','Israel','Japan','Russia','South_Africa','United_States']:
        ax[2].scatter(0,0, color = country_colour_dic[country], s = size * 2, label = country.replace('-',' ').replace('_',' '))
    leg.append(ax[2].legend(loc = 7, fontsize = 17, ncol = 1, bbox_to_anchor = (0.6,0.4,1,0.2)))

    for i in range(len(ax)):
        ax[i].set_ylabel('Attack Rate (%)', fontsize = fontsize)
        ax[i].set_yticks(np.arange(55, 81, 5))
        ax[i].set_ylim(54, 80)
        ax[i].tick_params(labelsize=fontsize-2)
        ax[i].tick_params(axis='x', which='minor', bottom=False)
        ax[i].tick_params(axis='y', which='minor', left=False)
        leg[i].draw_frame(False)

    plt.minorticks_off()

    fig_path = os.path.join(figdir, 'fig_6.pdf')
    fig.savefig(fig_path, format='pdf')



if __name__ == '__main__':
    
    reference_location = 'polymod_and_tomsk'
    reference_scenario = 'all_locations'

    beta = 0.04752
    susceptibility = 1.0
    gamma_inverse = 2.6
    num_agebrackets = 85

    countries = ['Australia', 'Canada', 'China', 'Europe', 'India', 'Israel', 'Japan', 'Russia', 'South_Africa', 'United_States']

    plot_fig(countries, reference_location, reference_scenario, beta, susceptibility, gamma_inverse, num_agebrackets)









