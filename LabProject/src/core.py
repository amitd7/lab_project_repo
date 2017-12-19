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

    print("all_larvas_values: ", values)
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
    fig, a = plt.subplots()

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

    # manage_statistic_tests(method, groups_calc_values, groups_amps_values, groups_phase_values, num_of_groups,
    # types_names)

    return groups_calc_values


def plot_results_graph(results, group_names, method):
    fig, ax = plt.subplots()

    mean_value, std_value, se_value = calc_statistic_values(pd.DataFrame.from_dict(results), 0)  # TODO 0 for calc on cols
    print("mean_value: ", mean_value)
    print("se_value: ", type(se_value))
    mean_value = list(mean_value)
    se_value = list(se_value)
    print("mean_value: ", mean_value)
    print("se_value: ", se_value)
    ind = np.arange(len(group_names))    # the x locations for the groups
    width = 0.35       # the width of the bars: can also be len(x) sequence

    rects = ax.bar(ind, mean_value, width, yerr=se_value)

    plt.ylabel(method)
    plt.title("---------")
    plt.xticks(ind, group_names)
    # plt.legend((p1[0], p2[0]), ('Men', 'Women'))
    autolabel(rects, ax)
    plt.show()


def autolabel(rects, ax):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height, "%.2f" % int(height), ha='center', va='bottom')



# def manage_statistic_tests(method, groups_values, groups_amps_values, groups_phase_values, num_of_groups,
# types_names): ### TODO ###
#
#     if num_of_groups > 2:
#         statistic_anova_test(groups_values, groups_amps_values, groups_phase_values)
#     else:
#         type1 = types_names[0]
#         type2 = types_names[1]
#         if method == "amplitude_phase":
#             statistic, p_value = statistic_tests(groups_amps_values[type1][0], groups_amps_values[type2][0], method)
#             print("statistic amp: ", statistic, "p_value amp: ", p_value)
#             statistic, p_value = statistic_tests(groups_phase_values[type1][0], groups_phase_values[type2][0], method)
#             print("statistic phase: ", statistic, "p_value phase: ", p_value)
#         else:
#             statistic, p_value = statistic_tests(groups_values[type1][0], groups_values[type2][0], method)
#             print("statistic: ", statistic, "p_value: ", p_value)


# def statistic_anova_test(groups_values, groups_amps_values, groups_phase_values):
#     stats.f_oneway()
#     # TODO see how to send all the lists in the same time
#     statistic = p_value = 0
#
#     return statistic, p_value


# def statistic_tests(group1, group2, method):
#     """
#     take a statistic test between the two groups according to the method used
#     :param group1:
#     :param group2:
#     :param method:
#     :return:
#     """
#
#     statistic = p_value = 0
#     if method == "fourier" or method == "chi_square" or method == "amplitude_phase":
#         statistic, p_value = stats.ttest_ind(group1, group2)
#     elif method == "g_factor":
#         statistic, p_value = g_factor.g_factor_significant(group1, group2)
#
#     return statistic, p_value

