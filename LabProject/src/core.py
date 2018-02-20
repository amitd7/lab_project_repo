from datetime import timedelta
from scipy import stats

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xlrd


import src.periodogram_calculation as pc
from src import amplitude_phase, g_factor

CIRCADIAN_TIME = HOURS_PER_DAY = 24
MINUTES_PER_HOUR = 60
MINUTES_PER_DAY = MINUTES_PER_HOUR * HOURS_PER_DAY

global window_size


def xls_to_csv(path, sheet, csv_file):
    """
    converting the xls file to csv
    :param path: the path to the xls file
    :param sheet: name of sheet to convert
    :param csv_file: name of the output csv file
    :return:
    """
    try:
        data_xls = pd.read_excel(path, sheet, index_col=None)
        data_xls.to_csv(csv_file, encoding='utf-8')
    except xlrd.XLRDError:
        print("Error Catch!")
        return "XLRDError"


def parse_input(input_file, types_names, wells_names_by_type, num_of_groups, start_time, sampling_intervals,
                excel_data_row, cols_labels, data_col, none_val_to_num):
    """

    :param input_file:
    :param types_names:
    :param wells_names_by_type:
    :param num_of_groups:
    :param start_time:
    :param sampling_intervals:
    :param excel_data_row:
    :param cols_labels:
    :param data_col:
    :param none_val_to_num:
    :return:
    """
    selected_cols = [cols_labels, data_col]
    temp_labels = ["type", "data"]

    global experiment_title
    try:
        experiment_title = pd.read_csv(input_file, header=None, nrows=1)[data_col][0]
        # experiment_title = pd.read_excel(input_file, header=None, nrows=1)[data_col][0]
    except FileNotFoundError:
        return "FileNotFoundError"
    except KeyError:
        return "UnexpectedError"

    try:
        all_data = pd.read_csv(input_file, header=None, skip_blank_lines=True, skiprows=(excel_data_row-1),
                               names=temp_labels, usecols=selected_cols, low_memory=False)
        # all_data = pd.read_excel(input_file, header=3, parse_cols=selected_cols)
    except FileNotFoundError:
        return "FileNotFoundError"

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
    :param rec_start_time:
    :return:
    """
    samples_per_hour = float(60/interval)
    diff_times = ct_zero[0] - rec_start_time[0] + ((ct_zero[1] - rec_start_time[1])/MINUTES_PER_HOUR)

    if diff_times > 0:
        diff_to_add = (HOURS_PER_DAY - diff_times) * samples_per_hour
    else:
        diff_to_add = int(abs(diff_times) * samples_per_hour)

    ct = [int((i/samples_per_hour))+float((i % samples_per_hour)/samples_per_hour) for i in
          range(0, int(num_of_values_per_larva + diff_to_add) + 1)]

    return [float("%.2f" % elem) for elem in ct], diff_to_add


def create_data_table(data, selected_wells, start_time, sampling_intervals, unknown_val):
    """
    creates data table for specific type of sample group (for example to WT sample)
    :param data:
    :param selected_wells:
    :param start_time:
    :param sampling_intervals:
    :param unknown_val:
    :return:
    """
    # return all rows that relate to the wt wells.
    data_rows = [[row for index, row in data.iterrows() if row["type"] == well] for well in selected_wells]

    try:
        # create list of data values lists
        data_values = [[float(data[1]) if data[1] != "-" else unknown_val for data in inner_list] for inner_list in
                       data_rows]
    except ValueError:
        return "ValueError"

    global num_of_values_per_larva, list_time

    if data_values:
        num_of_values_per_larva = len(data_values[0])
        total_days = (num_of_values_per_larva / (MINUTES_PER_HOUR / sampling_intervals)) / HOURS_PER_DAY

        list_time = time_range_list(start_time, total_days, sampling_intervals)

        # dictionary of well name and its data values
        data_dict = dict(zip(selected_wells, data_values))

        # Create DataFrame for each wt and mutant data
        df = pd.DataFrame.from_dict(data_dict)
        # Set DataFrame labels
        df.index = list_time
    else:
        df = pd.DataFrame()

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


def calculate_activity_patterns(data_table, method, circadian_time, sampling_intervals, ignore_part_beg,
                                ignore_part_end, samples_per_hour, diff_time, path_to_save, days_for_calc):
    """
    calculates period by fourier, period by chi-square, amplitude and phase, g factor
    the calculation is made for all larvas of one group type.
    the return value is a list of values of each larva
    :param data_table: data table of one group
    :param method:
    :param circadian_time:
    :param sampling_intervals:
    :param ignore_part_beg:
    :param ignore_part_end:
    :param samples_per_hour:
    :param diff_time:
    :param path_to_save:
    :param days_for_calc:
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
            value = pc.calculate_periodogram(relevant_data, header, "fourier", sampling_intervals, path_to_save)
        elif method == "chi_square":
            value = pc.calculate_periodogram(relevant_data, header, "chi_square", sampling_intervals, path_to_save)
        elif method == "amplitude_phase":
            value = amplitude_phase.amplitude_phase(samples_data, window_size, samples_per_hour,
                                                    list_time, circadian_time, diff_time)
        elif method == "g_factor":
            value = g_factor.g_factor_calculation(relevant_data, days_for_calc, samples_per_hour)
        else:
            value = None
        values.append(value)

    return values


def graph_data(c_times, diff_to_add, types_names, df, samples_per_hour, path_to_save, title, filename):
    """
    Create a graph of the average data of each group
    (contains all days of experiment)
    :param c_times:
    :param diff_to_add:
    :param types_names:
    :param df:
    :param samples_per_hour:
    :param path_to_save:
    :param title:
    :param filename:
    :return:
    """
    plt.clf()  # clears the entire current figure
    fig, a = plt.subplots(figsize=(7, 4.8), dpi=100, frameon=False)

    for name in types_names:
        # calc mean values of each row in the table
        data_mean = df[name].mean(axis=1, skipna=True).values

        none_list_to_add = [None]*int(diff_to_add)
        concatenated_data = np.concatenate((none_list_to_add, data_mean))

        data_se = df[name].sem(axis=1, skipna=True).values
        concatenated_data_se = np.concatenate((none_list_to_add, data_se))
        # plt.plot(concatenated_data)
        plt.errorbar(list(range(len(concatenated_data))), concatenated_data, concatenated_data_se, elinewidth=0.3)
    # the space between ticks needed in order to have a tick every 12 hours
    frequency = int(samples_per_hour * CIRCADIAN_TIME / 2)
    plt.xticks(np.arange(0, len(c_times)+1, frequency), c_times[::frequency])

    # major ticks every 12 hours, minor ticks every 1 hour
    ax = plt.gca()
    minor_ticks = np.arange(0, len(c_times)+1, samples_per_hour)
    ax.set_xticks(minor_ticks, minor=True)
    plt.title(title)
    plt.xlabel("Time (hrs)")
    plt.ylabel(experiment_title)
    try:
        plt.legend(types_names, loc="best")
    except IndexError:
        return "IndexError"
    plt.savefig(path_to_save + filename)
    return fig


def groups_calculation(method, types_names, df, c_times, sampling_intervals, ignore_part_beg, ignore_part_end,
                       samples_per_hour, diff_time, path, days_for_calc):
    """

    :param method:
    :param types_names:
    :param df:
    :param c_times:
    :param sampling_intervals:
    :param ignore_part_beg:
    :param ignore_part_end:
    :param samples_per_hour:
    :param diff_time:
    :param path:
    :param days_for_calc:
    :return:
    """

    groups_calc_values = {}

    for name in types_names:
        values = calculate_activity_patterns(df[name], method, c_times, sampling_intervals, ignore_part_beg,
                                             ignore_part_end, samples_per_hour, diff_time, path, days_for_calc)
        groups_calc_values[name] = values

    return groups_calc_values


def period_results_graph(results, group_names, method, path):
    """

    :param results:
    :param group_names:
    :param method:
    :param path:
    :return:
    """

    fig = plt.figure(frameon=False)

    dd = pd.DataFrame.from_dict(results)

    mean_value, std_value, se_value = calc_statistic_values(dd, 0)

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
        ax.text(rect.get_x() + rect.get_width()/2., height+se_value[i]+0.05, "%.2f" % se_value[i], ha="center",
                va="bottom")

    plt.ylabel("Period (hours)")
    plt.title("Period by groups (method used: %s)" % method)
    plt.xticks(y_pos, group_names)
    plt.savefig(path + method + " results graph")
    return fig


def g_factor_results_graph(data_to_pres, group_names, _, path):
    """

    :param data_to_pres:
    :param group_names:
    :param _:
    :param path:
    :return:
    """
    fig = plt.figure(frameon=False)
    # data_as_list = [value for value in data_to_pres.values()]
    data_as_list = [[x for x in lst if x is not None] for lst in data_to_pres.values()]
    plt.boxplot(data_as_list, showcaps=False, showbox=False, sym="", whiskerprops={"linestyle": ""}, labels=group_names,
                medianprops={"color": "r", "linewidth": 2})

    patch = mpatches.Patch(color="red", label="Median")
    plt.legend(handles=[patch])

    for i in range(len(group_names)):
        y = data_as_list[i]
        x = np.random.normal(1+i, 0, size=len(y))
        plt.plot(x, y, ".", markersize=5)

    plt.ylabel("G-factor")
    plt.title("G-factor by group")
    plt.savefig(path + "g-factor results graph")
    return fig


def amp_phase_results_graph(results, group_names, method, path):
    """

    :param results:
    :param group_names:
    :param method:
    :param path:
    :return:
    """

    fig = plt.figure(figsize=(8, 4.8), dpi=100, frameon=False)
    if method == "average_amplitude":
        fig.suptitle("Average amplitude & phase by group")
    else:
        fig.suptitle("Amplitude & phase by group")

    amplitude_results = {key: results[key]["amplitude"] for key in results}
    phase_results = {key: results[key]["phase CT"] for key in results}

    amplitude_results = pad_to_match_lengths(amplitude_results)
    phase_results = pad_to_match_lengths(phase_results)

    amp_mean, amp_std, amp_se = calc_statistic_values(pd.DataFrame.from_dict(amplitude_results), 0)
    amp_mean = [amp_mean[name] for name in group_names]
    amp_se = [amp_se[name] for name in group_names]

    phase_mean, phase_std, phase_se = calc_statistic_values(pd.DataFrame.from_dict(phase_results), 0)
    phase_mean = [phase_mean[name] for name in group_names]
    phase_se = [phase_se[name] for name in group_names]

    width = 0.5       # the width of the bars: can also be len(x) sequence
    y_pos = np.arange(len(group_names))

    ax = plt.subplot(121)
    rects1 = plt.bar(y_pos, amp_mean, align="center", width=width, yerr=amp_se)
    plt.ylabel("Average amplitude")
    plt.xticks(y_pos, group_names)

    # ax = plt.gca()

    # Attach a text label above each bar displaying its height
    for i, rect in enumerate(rects1):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 0.5*height, "%.2f" % height, ha="center", va="bottom")
        ax.text(rect.get_x() + rect.get_width()/2., float("%.2f" % height)+amp_se[i]+0.05, "%.2f" % amp_se[i],
                ha="center", va="bottom")

    ax = plt.subplot(122)
    rects2 = plt.bar(y_pos, phase_mean, align="center", width=width, yerr=phase_se)
    plt.ylabel("Average phase")
    plt.xticks(y_pos, group_names)

    # ax = plt.gca()
    # Attach a text label above each bar displaying its height
    for i, rect in enumerate(rects2):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 0.5*height, "%.2f" % height, ha="center", va="bottom")
        ax.text(rect.get_x() + rect.get_width()/2., float("%.2f" % height)+phase_se[i]+0.05, "%.2f" % phase_se[i],
                ha="center", va="bottom")

    plt.savefig(path + method + " results graph")
    return fig


def manage_statistic_tests(groups, method, group_names):
    """
    Run statistic test for every two groups
    :param groups: dictionary with groups results values. key - group name, value - list of group values
    :param method:
    :param group_names:
    :return:
    """
    num_of_groups = len(group_names)
    p_values = {}
    for i in range(num_of_groups):
        for j in range(i+1, num_of_groups):
            statistic, p_value = statistic_tests(groups[group_names[i]], groups[group_names[j]], method)
            p_values[(group_names[i], group_names[j])] = p_value
    return p_values


def statistic_tests(group1, group2, method):
    """
    take a statistic test between the two groups according to the method used
    :param group1:
    :param group2:
    :param method:
    :return:
    """
    group1 = [elem for elem in group1 if elem is not None]
    group2 = [elem for elem in group2 if elem is not None]

    if method == "g_factor":
        statistic, p_value = g_factor.g_factor_significant(group1, group2)
    else:

        statistic, p_value = stats.ttest_ind(group1, group2, nan_policy="omit")
    return statistic, p_value


def pad_to_match_lengths(dict_of_lists):
    max_len = 0
    for key in dict_of_lists.keys():
        lst_len = len(dict_of_lists[key])
        if lst_len > max_len:
            max_len = lst_len

    for key in dict_of_lists.keys():
        lst_len = len(dict_of_lists[key])
        if lst_len < max_len:
            dict_of_lists[key] += [None] * (max_len-lst_len)
    return dict_of_lists
