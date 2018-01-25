from datetime import datetime, timedelta
from scipy import stats

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


import src.periodogram_calculation as pc
from src import amplitude_phase, g_factor, clustering

CIRCADIAN_TIME = HOURS_PER_DAY = 24
MINUTES_PER_HOUR = 60
MINUTES_PER_DAY = MINUTES_PER_HOUR * HOURS_PER_DAY

global window_size

number_of_days = 3  # TODO need to calculate this value from the user's input of relevant days


def xls_to_csv(path, sheet, csv_file):
    """
    converting the xls file to csv
    :param path: the path to the xls file
    :param sheet: name of sheet to convert
    :param csv_file: name of the output csv file
    :return:
    """
    data_xls = pd.read_excel(path, sheet, index_col=None)
    data_xls.to_csv(csv_file, encoding='utf-8')


def parse_input(input_file, types_names, wells_names_by_type, num_of_groups, start_time, sampling_intervals,
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

    return {types_names[i]: create_data_table(all_data, wells_names_by_type[i], start_time, sampling_intervals,
                                              none_val_to_num) for i in range(num_of_groups)}


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


def circadian_time_list(interval, ct_zero, rec_start_time):
    """

    :param interval:
    :param ct_zero:
    :return:
    """
    # samples_per_hour = float(60/interval)
    # diff_times = (rec_start_time[0] - ct_zero[0]) + ((rec_start_time[1] - ct_zero[1])/MINUTES_PER_HOUR)
    # elems_to_ct_zero = diff_times * samples_per_hour
    #
    # # ct_zero_formatted = ct_zero[0] + (ct_zero[1]/60)
    # print("ct_zero  circadian time list  ", ct_zero)
    # cz_zero_index = (ct_zero[0]*samples_per_hour) + ((ct_zero[1]/MINUTES_PER_HOUR) * samples_per_hour)
    # print("cz_zero_index:  ", cz_zero_index)
    # ct = [int((i/samples_per_hour))+float((i % samples_per_hour)/samples_per_hour) for i in
    #       # range(0, num_of_values_per_larva + int(cz_zero_index) + 1)]
    #       range(0, num_of_values_per_larva - elems_to_ct_zero + 1)]
    #
    # before_ct_zero = [i for i in g_factor.frange(0-diff_times, elems_to_ct_zero, interval/MINUTES_PER_HOUR)]
    # return [float("%.2f" % elem) for elem in ct]

    samples_per_hour = float(60/interval)
    diff_times = abs(ct_zero[0] - rec_start_time[0]) + ((ct_zero[1] - rec_start_time[1])/MINUTES_PER_HOUR)
    # print("ctl - diff_times: ", diff_times)
    diff_to_add = 0
    if diff_times > 0:
        diff_to_add = (HOURS_PER_DAY - diff_times) * samples_per_hour
        # print("ctl - diff_to_add: ", diff_to_add)

    elems_to_ct_zero = diff_times * samples_per_hour
    # print("ctl - elems_to_ct_zero: ", elems_to_ct_zero)
    ct = [int((i/samples_per_hour))+float((i % samples_per_hour)/samples_per_hour) for i in
          # range(0, int(num_of_values_per_larva - elems_to_ct_zero) + 1)]
          range(0, int(num_of_values_per_larva + diff_to_add) + 1)]
    # print("ctl - ct: ", ct)
    # print("ctl - len(ct): ", len(ct))

    # time_bin_dec = interval/MINUTES_PER_HOUR
    # before_ct_zero = [float("%.2f" % i) for i in g_factor.frange(-diff_times, 0, time_bin_dec)]
    # print("ctl - before_ct_zero: ", before_ct_zero)
    # return before_ct_zero[:int(elems_to_ct_zero)] + [float("%.2f" % elem) for elem in ct]
    return [float("%.2f" % elem) for elem in ct], diff_to_add


def create_data_table(data, selected_wells, start_time, sampling_intervals, unknown_val): ### TODO ###
    """
    creates data table for specific type of sample group (for example to WT sample)
    :param data:
    :param selected_wells:
    :return:
    """
    # return all rows that relate to the wt wells.
    data_rows = [[row for index, row in data.iterrows() if row['type'] == well] for well in selected_wells]

    # create list of data values lists
    data_values = [[float(data[2]) if data[2] != '-' else unknown_val for data in inner_list] for inner_list in
                   data_rows]

    global num_of_values_per_larva, list_time

    num_of_values_per_larva = len(data_values[0])
    total_days = (num_of_values_per_larva / (60 / sampling_intervals)) / 24

    list_time = time_range_list(start_time, total_days, sampling_intervals)

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


def calculate_activity_patterns(data_table, method, start_time, circadian_time, sampling_intervals, ignore_part_beg,
                                ignore_part_end, samples_per_hour, diff_time): ### TODO ###
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
        relevant_data = np.array(samples_data[ignore_part_beg: ignore_part_end])
        if method == "fourier":
            value = pc.calculate_periodogram(relevant_data, i, "fourier", sampling_intervals)
        elif method == "chi_square":
            value = pc.calculate_periodogram(relevant_data, i, "chi_square", sampling_intervals)
        elif method == "amplitude_phase":
            value = amplitude_phase.amplitude_phase(samples_data, sampling_intervals, window_size, samples_per_hour,
                                                    list_time, circadian_time, diff_time)
        elif method == "g_factor":
            value = g_factor.g_factor_calculation(relevant_data, number_of_days, samples_per_hour)
        else:
            value = None
        values.append(value)

    return values


def graph_data(sampling_intervals, c_times, diff_to_add, types_names, DF, samples_per_hour): ### TODO ###
    """
    Create a graph of the average data of each group
    (contains all days of experiment)
    :return:
    """
    print("c_times: ", c_times)
    print("len c_times: ", len(c_times))
    plt.clf()  # clears the entire current figure
    fig, a = plt.subplots(figsize=(7, 4.8), dpi=100, frameon=False)

    for name in types_names:
        # calc mean values of each row in the table
        data_mean = DF[name].mean(axis=1, skipna=True).values

        none_list_to_add = [None]*int(diff_to_add)
        concatenated_data = np.concatenate((none_list_to_add, data_mean))

        plt.plot(concatenated_data)
    # the space between ticks needed in order to have a tick every 12 hours
    frequency = int(samples_per_hour * CIRCADIAN_TIME / 2)
    plt.xticks(np.arange(0, len(c_times)+1, frequency), c_times[::frequency])

    # major ticks every 12 hours, minor ticks every 1 hour
    ax = plt.gca()
    minor_ticks = np.arange(0, len(c_times)+1, samples_per_hour)
    ax.set_xticks(minor_ticks, minor=True)

    plt.title("Average data per experimental group by time")
    plt.xlabel("Time (hrs)")
    plt.ylabel(experiment_title)
    plt.legend(types_names, loc='upper left')
    plt.savefig("data_average_graph")
    return fig


def groups_calculation(method, types_names, DF, num_of_groups, start_time, c_times, sampling_intervals,
                       ignore_part_beg, ignore_part_end, samples_per_hour, diff_time): ### TODO ###

    groups_calc_values = {}

    for name in types_names:
        values = calculate_activity_patterns(DF[name], method, start_time, c_times, sampling_intervals,
                                             ignore_part_beg, ignore_part_end, samples_per_hour, diff_time)
        groups_calc_values[name] = values

    return groups_calc_values


def period_results_graph(results, group_names, method):

    fig = plt.figure(frameon=False)
    mean_value, std_value, se_value = calc_statistic_values(pd.DataFrame.from_dict(results), 0)  # TODO 0 for calc on cols

    mean_value = [mean_value[name] for name in group_names]
    se_value = [se_value[name] for name in group_names]

    width = 0.5       # the width of the bars: can also be len(x) sequence
    y_pos = np.arange(len(group_names))
    rects = plt.bar(y_pos, mean_value, align="center", width=width, yerr=se_value)

    ax = plt.gca()

    # Attach a text label above each bar displaying its height
    for i, rect in enumerate(rects):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 0.5*height, "%.2f" % height, ha="center", va="bottom")
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height, "%.2f" % se_value[i], ha="center", va="bottom")

    plt.ylabel("Period (hours)")
    plt.title("Period by groups")
    plt.xticks(y_pos, group_names)

    plt.savefig(method + " results graph")
    return fig


def g_factor_results_graph(data_to_pres, group_names, action):
    fig = plt.figure(frameon=False)

    data_as_list = [value for value in data_to_pres.values()]
    plt.boxplot(data_as_list, zorder=0, showcaps=False, showbox=False, labels=group_names,
                medianprops={"color": "b", "linewidth": 2})

    patch = mpatches.Patch(color='blue', label="Median")
    plt.legend(handles=[patch])

    for i in range(len(group_names)):
        y = data_as_list[i]
        x = np.random.normal(1+i, 0, size=len(y))
        plt.plot(x, y, 'k.', markersize=5, alpha=1)
        # plt.figtext(i/2 + 1/4, medians[i] + 0.2, medians[i],  horizontalalignment='center',  color='k', weight='semibold')

    plt.ylabel("G factor")
    plt.title("G factor by group")
    plt.savefig("g factor results graph")
    return fig


def amp_phase_results_graph(results, group_names, method, title):

    fig = plt.figure(figsize=(8, 4.8), dpi=100, frameon=False)
    if method == "average_amplitude":
        fig.suptitle("Average amplitude & phase by group")
    else:
        fig.suptitle("Amplitude & phase by group")

    amplitude_results = {key: results[key]["amplitude"] for key in results}
    phase_results = {key: results[key]["phase c.t."] for key in results}

    amp_mean, amp_std, amp_se = calc_statistic_values(pd.DataFrame.from_dict(amplitude_results), 0)
    amp_mean = [amp_mean[name] for name in group_names]
    amp_se = [amp_se[name] for name in group_names]

    phase_mean, phase_std, phase_se = calc_statistic_values(pd.DataFrame.from_dict(phase_results), 0)
    phase_mean = [phase_mean[name] for name in group_names]
    phase_se = [phase_se[name] for name in group_names]

    width = 0.5       # the width of the bars: can also be len(x) sequence
    y_pos = np.arange(len(group_names))

    plt.subplot(121)
    rects1 = plt.bar(y_pos, amp_mean, align="center", width=width, yerr=amp_se)
    plt.ylabel("Average amplitude")
    plt.xticks(y_pos, group_names)

    ax = plt.gca()
    print("ax1: ", ax)

    # Attach a text label above each bar displaying its height
    for i, rect in enumerate(rects1):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 0.5*height, "%.2f" % height, ha="center", va="bottom")
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height, "%.2f" % amp_se[i], ha="center", va="bottom")

    plt.subplot(122)
    rects2 = plt.bar(y_pos, phase_mean, align="center", width=width, yerr=phase_se)
    plt.ylabel("Average phase")
    plt.xticks(y_pos, group_names)

    ax = plt.gca()
    print("ax2: ", ax)
    # Attach a text label above each bar displaying its height
    for i, rect in enumerate(rects2):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 0.5*height, "%.2f" % height, ha="center", va="bottom")
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height, "%.2f" % phase_se[i], ha="center", va="bottom")

    plt.savefig(method + " results graph")
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

    if method == "g_factor":
        statistic, p_value = g_factor.g_factor_significant(group1, group2)
        print("          g factor statistic, p_value ", statistic, p_value)
    else:
        statistic, p_value = stats.ttest_ind(group1, group2)

    return statistic, p_value

