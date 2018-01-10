from datetime import datetime, timedelta

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import src.periodogram_calculation as pc
from src import amplitude_phase, g_factor, clustering

CIRCADIAN_TIME = HOURS_PER_DAY = 24
MINUTES_PER_DAY = 1440

global window_size

number_of_days = 3  # TODO need to calculate this value from the user's input of relevant days

# if false, all the '-' values in the data will translate as 0, otherwise as None and will be out of calculations
ignore_nan_values = False

# The relevant columns to parse
# cols_labels = 2  # col 'C'
# time_labels = 3  # col 'D'
# data_col = 5  # col 'F'


def xls_to_csv(path, sheet, csv_file):
    """
    converting the xls file to csv
    :param path: the path to the xls file
    :param sheet: name of sheet to convert
    :param csv_file: name of the output csv file
    :return:
    """
    print("path: ", path)
    data_xls = pd.read_excel(path, sheet, index_col=None)
    data_xls.to_csv(csv_file, encoding='utf-8')


def parse_input(input_file, types_names, wells_names_by_type, num_of_groups, start_time, total_days, sampling_intervals,
                cols_labels, time_labels, data_col, none_val_to_num): ### TODO ###
    """

    :param input_file: input filename
    :return:
    """

    selected_cols = [cols_labels, time_labels, data_col]
    temp_labels = 'type', 'time', 'data'
    # TODO find how to generalize the header value
    global experiment_title
    experiment_title = pd.read_csv(input_file, header=None, nrows=1)[data_col][0]

    all_data = pd.read_csv(input_file, header=3, names=temp_labels, usecols=selected_cols, low_memory=False)
    print("all_data: ", all_data)

    return {types_names[i]: create_data_table(all_data, wells_names_by_type[i], start_time, total_days,
                                              sampling_intervals, none_val_to_num) for i in range(num_of_groups)}


def time_range(start, end, delta):
    """
    Generates values of time range between start and end with delta interval
    :param start:
    :param end:
    :param delta:
    :return:
    """
    curr = start
    while curr < end:
        yield curr
        curr += delta


def time_range_list(start, days, interval):
    """
    create a list of time using time_range function from start_time up to total_days days with interval of
    sampling_intervals
    :param start:
    :param days:
    :param interval:
    :return:
    """
    return [result.time().isoformat() for result in time_range(start, start.replace() + timedelta(days=days),
                                                               timedelta(minutes=interval))]


def circadian_time_list(interval, num_of_days):
    """

    :param interval:
    :param num_of_days:
    :return:
    """
    num_of_values_per_day = MINUTES_PER_DAY / interval
    x = float(60/interval)
    return [int((i/x))+float((i % x)/x) for i in range(0, int(num_of_values_per_day*num_of_days)+1)]


def create_data_table(data, selected_wells, start_time, total_days, sampling_intervals, unknown_val): ### TODO ###
    """
    creates data table for specific type of sample group (for example to WT sample)
    :param data:
    :param selected_wells:
    :return:
    """

    list_time = time_range_list(start_time, total_days, sampling_intervals)

    # return all rows that relate to the wt wells.
    data_rows = [[row for index, row in data.iterrows() if row['type'] == well] for well in selected_wells]

    # create list of data values lists
    data_values = [[float(data[2]) if data[2] != '-' else unknown_val for data in inner_list] for inner_list in
                   data_rows]

    # dictionary of well name and its data values
    data_dict = dict(zip(selected_wells, data_values))

    # Create DataFrame for each wt and mutant data
    df = pd.DataFrame.from_dict(data_dict)

    # Set DataFrame labels
    df.index = list_time

    return df


def calc_statistic_values(df, axis_to_calc):
    """

    :param df:
    :param axis_to_calc: 1 - calculating value for every row, 0 calculating value for every col
    :return:
    """

    # calc mean values of each row
    mean_value = df.mean(axis=axis_to_calc, skipna=True)

    # calc std values of each row
    std_value = df.std(axis=axis_to_calc, skipna=True)

    # calc standard error values of each row
    se_value = df.sem(axis=axis_to_calc, skipna=True)

    return mean_value, std_value, se_value


def calculate_activity_patterns(data_table, method, start_time, total_days, sampling_intervals, ignore_part_beg,
                                ignore_part_end, samples_per_hour): ### TODO ###
    """
    calculates period by fourier (0), period by chi-square (1), amplitude and phase (2), g factor (3)
    the calculation is made for all larvas of one group type.
    the return value is a list of values of each larva
    :param data_table: data table of one group
    :param method:
    :return:
    """

    values = []
    columns_headers = data_table.columns.values.tolist()

    # runs on every larva in the group and calculates the chosen calculation
    for i, header in enumerate(columns_headers):
        samples_data = data_table[header].tolist()

        # removes the ignored hours - i.e the irrelevant rows
        relevant_data = np.array(samples_data[ignore_part_beg: len(samples_data) - ignore_part_end])
        if method == "fourier":
            value = pc.calculate_periodogram(relevant_data, i, "fourier", sampling_intervals)
        elif method == "chi_square":
            value = pc.calculate_periodogram(relevant_data, i, "chi_square", sampling_intervals)
        elif method == "amplitude_phase":
            times = time_range_list(start_time, total_days, sampling_intervals)
            circadian_time = circadian_time_list(sampling_intervals, total_days)
            value = amplitude_phase.amplitude_phase(samples_data, sampling_intervals, window_size, samples_per_hour,
                                                    times, circadian_time)
        elif method == "g_factor":
            value = g_factor.g_factor_calculation(relevant_data, number_of_days, samples_per_hour)
            # TODO add a call for the function g_factor_significant()
        else:
            value = None
        values.append(value)

    return values


# def smooth_group_data(data_table, window=20):
#     """
#     Runs over all larvas in one group and smooths the data using moving average
#     :param data_table:
#     :param window:
#     :return:
#     """
#     columns_headers = data_table.columns.values.tolist()
#     smoothed_data_array = []
#     # runs on every larva in the group and smooth the data
#     for header in columns_headers:
#         smoothed_data_array.append(amplitude_phase.moving_average(data_table[header].tolist(), window))
#
#     smoothed_data_array = pd.DataFrame(smoothed_data_array).T
#     # TODO: add the correct time list as the size of data is smaller
#     return smoothed_data_array


def graph_data(sampling_intervals, total_days, types_names, DF, samples_per_hour): ### TODO ###
    """
    Create a graph of the average data of each group
    (contains all days of experiment)
    :return:
    """
    ct = circadian_time_list(sampling_intervals, total_days)
    ct = [float("%.2f" % elem) for elem in ct]

    plt.clf()  # clears the entire current figure
    fig, a = plt.subplots(frameon=False)

    for name in types_names:
        # calc mean values of each row in the table
        data_mean = DF[name].mean(axis=1, skipna=True).values
        plt.plot(data_mean)
    # the space between ticks needed in order to have a tick every 12 hours
    frequency = int(samples_per_hour * CIRCADIAN_TIME / 2)
    plt.xticks(np.arange(0, len(ct)+1, frequency), ct[::frequency])

    # major ticks every 12 hours, minor ticks every 1 hour
    ax = plt.gca()
    minor_ticks = np.arange(0, len(ct)+1, samples_per_hour)
    ax.set_xticks(minor_ticks, minor=True)

    plt.title("Average data per experimental group by time")
    plt.xlabel("Time (hrs)")
    plt.ylabel(experiment_title)
    plt.legend(types_names, loc='upper left')
    plt.savefig("data_average_graph")
    return fig


def groups_calculation(method, types_names, DF, num_of_groups, start_time, total_days, sampling_intervals,
                       ignore_part_beg, ignore_part_end, samples_per_hour): ### TODO ###

    groups_calc_values = {}

    for name in types_names:
        values = calculate_activity_patterns(DF[name], method, start_time, total_days, sampling_intervals,
                                             ignore_part_beg, ignore_part_end, samples_per_hour)
        groups_calc_values[name] = values

    return groups_calc_values


def period_results_graph(results, group_names, method):

    fig = plt.figure(frameon=False)
    mean_value, std_value, se_value = calc_statistic_values(pd.DataFrame.from_dict(results), 0)  # TODO 0 for calc on cols
    print("mean_value: ", mean_value)
    mean_value = list(mean_value)
    se_value = list(se_value)
    print("mean_value: ", mean_value)
    print("se_value: ", se_value)

    width = 0.5       # the width of the bars: can also be len(x) sequence
    y_pos = np.arange(len(group_names))
    plt.bar(y_pos, mean_value, align="center", width=width, yerr=se_value)
    plt.ylabel(method)
    plt.title(method)
    plt.xticks(y_pos, group_names)

    plt.savefig(method + " results graph")
    return fig


def g_factor_results_graph(data_to_pres, group_names, action):
    fig = plt.figure(frameon=False)

    data_as_list = [value for value in data_to_pres.values()]
    bp = plt.boxplot(data_as_list, zorder=0, showcaps=False, showbox=False, labels=group_names,
                     medianprops={"color": "b", "linewidth": 2})

    for i in range(len(group_names)):
        y = data_as_list[i]
        x = np.random.normal(1+i, 0, size=len(y))
        plt.plot(x, y, 'k.', alpha=1)
    plt.ylabel(action)
    plt.title(action)
    plt.savefig("g factor results graph")
    return fig


def amp_phase_results_graph(results, group_names, method, title):

    fig = plt.figure(figsize=(8, 4.8), dpi=100, frameon=False)

    amplitude_results = {key: results[key]["amplitude"] for key in results}
    phase_results = {key: results[key]["phase c.t."] for key in results}

    amp_mean, amp_std, amp_se = calc_statistic_values(pd.DataFrame.from_dict(amplitude_results), 0)
    print("amp_mean: ", amp_mean)
    amp_mean = list(amp_mean)
    amp_se = list(amp_se)
    print("amp_mean: ", amp_mean)
    print("amp_se: ", amp_se)

    phase_mean, phase_std, phase_se = calc_statistic_values(pd.DataFrame.from_dict(phase_results), 0)
    print("phase_mean: ", phase_mean)
    phase_mean = list(phase_mean)
    phase_se = list(phase_se)
    print("phase_mean: ", phase_mean)
    print("phase_se: ", phase_se)

    width = 0.5       # the width of the bars: can also be len(x) sequence
    y_pos = np.arange(len(group_names))

    plt.subplot(121)
    plt.bar(y_pos, amp_mean, align="center", width=width, yerr=amp_se)
    plt.ylabel("Average amplitude")
    # plt.title("Average amplitude and standard error per experimental group")
    plt.xticks(y_pos, group_names)

    plt.subplot(122)
    plt.bar(y_pos, phase_mean, align="center", width=width, yerr=phase_se)
    plt.ylabel("Average phase")
    # plt.title("Average phase and standard error per experimental group")
    plt.xticks(y_pos, group_names)

    plt.savefig("amplitude & phase results graph")
    return fig


def f(t):
    return np.exp(-t) * np.cos(2*np.pi*t)


def manage_statistic_tests(groups, method, group_names):
    """
    Run statistic test for every two groups
    :param groups: dictionary with groups results values. key - group name, value - list of group values
    :param method:
    :param group_names:
    :return:
    """
    print("method ", method)
    num_of_groups = len(group_names)
    p_values = {}
    for i in range(num_of_groups):
        for j in range(i+1, num_of_groups):
            statistic, p_value = statistic_tests(groups[group_names[i]], groups[group_names[j]], method)
            p_values[(group_names[i], group_names[j])] = p_value
    print("p_value dict  :", p_values)
    return p_values


def statistic_tests(group1, group2, method):
    """
    take a statistic test between the two groups according to the method used
    :param group1:
    :param group2:
    :param method:
    :return:
    """

    statistic = p_value = 0
    if method == "g_factor":
        statistic, p_value = g_factor.g_factor_significant(group1, group2)
        print("          g factor statistic, p_value ", statistic, p_value)
    else:
        statistic, p_value = stats.ttest_ind(group1, group2)

    return statistic, p_value

