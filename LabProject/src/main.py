from datetime import datetime, timedelta

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import src.periodogram_calculation as pc
from src import amplitude_phase, g_factor

CIRCADIAN_TIME = HOURS_PER_DAY = 24
MINUTES_PER_DAY = 1440


plate_size = (6, 8)
num_of_dif_types = 2  # number of mutations and wt
types_names = ["WT", "MUT"]

start_time = 9.00
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

start_time2 = datetime(datetime.now().year, datetime.now().month, datetime.now().day, hour=9, minute=0)

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

    return {types_names[i]: create_data_table(all_data, wells_names_by_type[i]) for i in range(0, num_of_dif_types)}


# def time_list_by_decimal():
#     """
#     Creates a list of time by decimal
#     :return:
#     """
#     list_time = []
#     list_time2 = []
#
#     num_of_values_per_day = MINUTES_PER_DAY / sampling_intervals
#     x = float(60/sampling_intervals)
#
#     for i in range(1, int(num_of_values_per_day*total_days)):
#         list_time.append(int(start_time+(i/x))+float((i % x)/x))
#         list_time2.append(float(start_time+((i % 60)/100)+int(i/60)))  # time for showing
#     return list_time, list_time2


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
    create a list of time using time_range function from start_time2 up to total_days days with interval of
    sampling_intervals
    :param start:
    :param days:
    :param interval:
    :return:
    """
    # TODO change the start_time2 name (check if still need the start_time)

    return [result.time().isoformat() for result in time_range(start, start.replace() +
                                                    timedelta(days=days), timedelta(minutes=interval))]


def create_data_table(data, selected_wells):
    """
    creates data table for specific type of sample group (for example to WT sample)
    :param data:
    :param selected_wells:
    :return:
    """

    list_time = time_range_list(start_time2, total_days, sampling_intervals)

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
    df['mean'] = mean_value
    df['std'] = std_value
    df['s. error'] = se_value
    # TODO make sure id adds the values and they can be seen outside the function


def graph_maker(df_dict, type_name, time_list):  # TODO check if I really need this function
    print("in graph")
    # print(df_dict.get(type_name[0])['s. error'].tolist())

    # first_dimension = 0
    time_list = [start_time] + time_list

    # l = [first_dimension] + df_dict.get(type_name[0])['s. error'].tolist()
    # print(l)

    arr1 = np.array(df_dict.get(type_name[0])['s. error'].tolist())
    # arr2 = np.array([first_dimension] + df_dict.get(type_name[1])['s. error'].tolist())

    print("time_list: ", time_list)
    print("y: ", arr1)

    # create list of the times in numbers
    plt.plot(time_list, arr1)
    plt.show()
    # plt.plot(arr2, time_list)

    # plt.plot(df.get("MUT")['s. error'], time_list)


def calculate_activity_patterns(data_table, method):
    """
    calculates period by fourier (0), period by chi-square (1), amplitude and phase (2), g factor (3)
    the calculation is made for all larvas of one group type.
    the return value is a list of values of each larva
    :param data_table:
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
            times = time_range_list(start_time2, total_days, sampling_intervals)
            value = amplitude_phase.amplitude_phase(moving_average(samples_data, window_size), window_size,
                                                    samples_per_hour, times)
        elif method == "g_factor":
            value = g_factor.g_factor_calculation(relevant_data, number_of_days, samples_per_hour)
            # TODO add a call for the function g_factor_significant()
        else:
            value = None
        values.append(value)

    print("all_larvas_values: ", values)


def smooth_data(data_dict):

    smooth_data_dict = {"WT": smooth_group_data(DF['WT'])}

    print("f  ", smooth_data_dict['WT'])


def smooth_group_data(data_table):
    """
    Runs over all larvas in one group and smooths the data using moving average
    :param data_table:
    :return:
    """
    columns_headers = data_table.columns.values.tolist()
    smoothed_data_array = []
    # runs on every larva in the group and smooth the data
    for header in columns_headers:
        smoothed_data_array.append(moving_average(data_table[header].tolist(), 20))
    # TODO add the columns names
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


if __name__ == "__main__":

    # xls_to_csv("files/example3.xlsx", "Analysis", "files/csv_file.csv")
    DF = parse_input("C:/Users/Amit/PycharmProjects/LabProject/files/csv_file.csv")

    # TODO calc_statistic_values(tables_dict[types_names[i]])

    calculate_activity_patterns(DF['WT'], "fourier")

    # g_factor_significant()
    # graph_maker(DF, "WT", time_list_by_decimal())

















