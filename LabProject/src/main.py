from datetime import datetime, timedelta

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import src.periodogram_calculation as pc
from src import amplitude_phase, g_factor, clustering

CIRCADIAN_TIME = HOURS_PER_DAY = 24
MINUTES_PER_DAY = 1440


plate_size = (6, 8)
num_of_groups = 2  # number of mutations and wt groups
types_names = ["WT", "MUT"]

sampling_intervals = 10  # in minutes

samples_per_hour = 6
number_of_days = 3
# select the relevant days
hours_to_ignore_beg = 24  # hours to ignore from the beginning
hours_to_ignore_end = 24  # hours to ignore from the end
ignore_part_beg = hours_to_ignore_beg * samples_per_hour  # number of samples to ignore from the beginning
ignore_part_end = hours_to_ignore_end * samples_per_hour  # number of samples to ignore from the end

# if false, all the '-' values in the data will translate as 0, otherwise as None and will be out of calculations
ignore_nan_values = False

total_days = 5  # TODO should be in days??

start_time = datetime(datetime.now().year, datetime.now().month, datetime.now().day, hour=9, minute=0)

total_num_of_values = int(MINUTES_PER_DAY / sampling_intervals * total_days)

# The relevant columns to parse
cols_labels = 2  # col 'C'
time_labels = 3  # col 'D'
data_col = 5  # col 'F'

# TODO generalize the lists
# wells_names_by_type is a list of lists. each inner list contains the wells names of each type
wells_names_by_type = [["A2", "A4", "A6", "A8", "B2", "B4", "B6", "B8", "C2", "C4", "C6", "C8", "D1", "D3", "D5",
                        "D7", "E1", "E3", "E5", "E7", "F1", "F3", "F5", "F7"],
                       ["A1", "A3", "A5", "A7", "B1", "B3", "B5", "B7", "C1", "C3", "C5", "C7", "D2", "D4", "D6",
                        "D8", "E2", "E4", "E6", "E8", "F2", "F4", "F6", "F8"]]
# wells_names_by_type = [["A2", "A4"],["A1", "A3"]]


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


def parse_input(input_file):
    """

    :param input_file: input filename
    :return:
    """

    selected_cols = [cols_labels, time_labels, data_col]
    temp_labels = 'type', 'time', 'data'
    # TODO find how to generalize the header value
    all_data = pd.read_csv(input_file, header=3, names=temp_labels, usecols=selected_cols, low_memory=False)

    return {types_names[i]: create_data_table(all_data, wells_names_by_type[i]) for i in range(0, num_of_groups)}


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
    return [result.time().isoformat() for result in time_range(start, start.replace() +
                                                    timedelta(days=days), timedelta(minutes=interval))]


def circadian_time_list(interval, num_of_days):
    """

    :param interval:
    :param num_of_days:
    :return:
    """
    num_of_values_per_day = MINUTES_PER_DAY / interval
    x = float(60/interval)
    return [int((i/x))+float((i % x)/x) for i in range(0, int(num_of_values_per_day*num_of_days)+1)]


def create_data_table(data, selected_wells):
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
    if not ignore_nan_values:
        unknown_val = 0
    else:
        unknown_val = None
    data_values = [[float(data[2]) if data[2] != '-' else unknown_val for data in inner_list] for inner_list in
                   data_rows]

    # dictionary of well name and its data values
    data_dict = dict(zip(selected_wells, data_values))

    # Create DataFrame for each wt and mutant data
    df = pd.DataFrame.from_dict(data_dict)

    # Set DataFrame labels
    df.index = list_time

    return df


def calc_statistic_values(df):
    """

    :param df:
    :return:
    """

    # calc mean values of each row
    mean_value = df.mean(axis=1, skipna=True)

    # calc std values of each row
    std_value = df.std(axis=1, skipna=True)

    # calc standard error values of each row
    se_value = df.sem(axis=1, skipna=True)

    # add values to table
    # df['mean'] = mean_value
    # df['std'] = std_value
    # df['s. error'] = se_value
    # TODO make sure id adds the values and they can be seen outside the function

    return mean_value, std_value, se_value


def calculate_activity_patterns(data_table, method):
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
        relevant_data = np.array(samples_data[ignore_part_beg - 1: len(samples_data) - ignore_part_end - 1])
        if method == "fourier":
            value = pc.calculate_periodogram(relevant_data, i, "fourier")
        elif method == "chi_square":
            value = pc.calculate_periodogram(relevant_data, i, "chi_square")
        elif method == "amplitude_phase":
            window_size = 20
            times = time_range_list(start_time, total_days, sampling_intervals)
            circadian_time = circadian_time_list(sampling_intervals, total_days)
            value = amplitude_phase.amplitude_phase(moving_average(samples_data, window_size), window_size,
                                                    samples_per_hour, times, circadian_time)
        elif method == "g_factor":
            value = g_factor.g_factor_calculation(relevant_data, number_of_days, samples_per_hour)
            # TODO add a call for the function g_factor_significant()
        elif method == "clustering":
            value = None
            print("clustering")  # TODO
        else:
            value = None
        values.append(value)

    print("all_larvas_values: ", values)
    return values


def smooth_group_data(data_table, window=20):
    """
    Runs over all larvas in one group and smooths the data using moving average
    :param data_table:
    :param window:
    :return:
    """
    columns_headers = data_table.columns.values.tolist()
    smoothed_data_array = []
    # runs on every larva in the group and smooth the data
    for header in columns_headers:
        smoothed_data_array.append(moving_average(data_table[header].tolist(), window))

    smoothed_data_array = pd.DataFrame(smoothed_data_array).T
    # TODO: add the correct time list as the size of data is smaller
    return smoothed_data_array


def moving_average(x, window=20):
    """
    Smooths the data in the column vector x using a moving average
    :param x: data list
    :param window: The default span for the moving average is 20
    :return: list of smoothed data
    """
    cum_sum = np.cumsum(np.insert(x, 0, 0))
    return (cum_sum[window:] - cum_sum[:-window]) / window


def graph_data():
    """
    Create a graph of the average data of each group
    (contains all days of experiment)
    :return:
    """

    ct = circadian_time_list(sampling_intervals, total_days)
    ct = [float("%.2f" % elem) for elem in ct]

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

    plt.legend(types_names, loc='upper left')
    plt.savefig("data_average_graph")
    plt.show()
    plt.clf()  # clears the entire current figure


def groups_calculation(method):

    # TODO make the function more efficient?
    groups_calc_values = {}
    groups_amps_values = {}
    groups_phase_values = {}
    for name in types_names:
        values = calculate_activity_patterns(DF[name], method)
        if method == "amplitude_phase":
            formatted_values = [(float("%.2f" % elem[0]), (elem[1][0], float("%.2f" % elem[1][1]))) for elem in values]

            amps = [float("%.2f" % elem[0]) for elem in values]
            phase_circ = [float("%.2f" % elem[1][1]) for elem in values]  # only circadian values

            groups_amps_values[name] = (amps, ("mean", float("%.2f" % np.mean(amps))), ("sem",
                                                                                        float("%.2f" % stats.sem(amps))))
            groups_phase_values[name] = (phase_circ, ("mean", float("%.2f" % np.mean(phase_circ))),
                                         ("sem", float("%.2f" % stats.sem(phase_circ))))

            groups_calc_values[name] = formatted_values
        else:
            formatted_values = [float("%.2f" % elem) for elem in values]
            groups_calc_values[name] = (formatted_values, ("mean", float("%.2f" % np.mean(formatted_values))),
                                        ("sem", float("%.2f" % stats.sem(formatted_values))))

    manage_statistic_tests(method, groups_calc_values, groups_amps_values, groups_phase_values)

    return groups_calc_values


def manage_statistic_tests(method, groups_values, groups_amps_values, groups_phase_values):

    if num_of_groups > 2:
        statistic_anova_test(groups_values, groups_amps_values, groups_phase_values)
    else:
        type1 = types_names[0]
        type2 = types_names[1]
        if method == "amplitude_phase":
            statistic, p_value = statistic_tests(groups_amps_values[type1][0], groups_amps_values[type2][0], method)
            print("statistic amp: ", statistic, "p_value amp: ", p_value)
            statistic, p_value = statistic_tests(groups_phase_values[type1][0], groups_phase_values[type2][0], method)
            print("statistic phase: ", statistic, "p_value phase: ", p_value)
        else:
            statistic, p_value = statistic_tests(groups_values[type1][0], groups_values[type2][0], method)
            print("statistic: ", statistic, "p_value: ", p_value)


def statistic_anova_test(groups_values, groups_amps_values, groups_phase_values):
    stats.f_oneway()
    # TODO see how to send all the lists in the same time
    statistic = p_value = 0

    return statistic, p_value


def statistic_tests(group1, group2, method):
    """
    take a statistic test between the two groups according to the method used
    :param group1:
    :param group2:
    :param method:
    :return:
    """

    statistic = p_value = 0
    if method == "fourier" or method == "chi_square" or method == "amplitude_phase":
        statistic, p_value = stats.ttest_ind(group1, group2)
    elif method == "g_factor":
        statistic, p_value = g_factor.g_factor_significant(group1, group2)

    return statistic, p_value


if __name__ == "__main__":

    # xls_to_csv("files/example3.xlsx", "Analysis", "files/csv_file.csv")
    DF = parse_input("C:/Users/Amit/PycharmProjects/lab_project_repo/LabProject/files/csv_file.csv")

    # TODO calc_statistic_values(tables_dict[types_names[i]])
    # graph_data()
    # calculate_activity_patterns(DF['WT'], "amplitude_phase")
    # print(groups_calculation("amplitude_phase"))

    ### run mean amplitude of one group ###
    # ct_list = circadian_time_list(sampling_intervals, total_days)
    # window_s = 20
    # ct_list = ct_list[int(window_s/2): - int(window_s/2)]
    #
    # data = smooth_group_data(DF["WT"], window_s)
    #
    # print(amplitude_phase.average_amplitude(float(12), float(36), data, ct_list))
    ####################################


    # print(clustering.cluster(data))

