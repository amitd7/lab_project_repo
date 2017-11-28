from datetime import datetime, timedelta
from matplotlib import style
from tkinter import filedialog

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import string
import tkinter as tk

import src.periodogram_calculation as pc
from src import amplitude_phase, g_factor, clustering

CIRCADIAN_TIME = HOURS_PER_DAY = 24
MINUTES_PER_DAY = 1440

style.use("ggplot")

# plate_size = (6, 8)
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


def xls_to_csv(sheet, csv_file):
    """
    converting the xls file to csv
    :param path: the path to the xls file
    :param sheet: name of sheet to convert
    :param csv_file: name of the output csv file
    :return:
    """
    data_xls = pd.read_excel(input_file_path, sheet, index_col=None)
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

    plt.legend(types_names, loc='upper left')
    plt.savefig("data_average_graph")
    # plt.show()
    return fig


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


# GUI and GUI-related functions:

def select_input_file():
    global input_file_path
    input_file_path = filedialog.askopenfile()
    browse_entry.insert(0, input_file_path.name)


# def add_entry(master, text):
#
#     column, row = master.grid_size()
#
#     label = tk.Label(master, text=text)
#     label.grid(row=row, column=0, sticky=tk.E, padx=2)
#
#     entry = tk.Entry(master)
#     entry.grid(row=row, column=1, sticky=tk.E+tk.W)
#
#     return entry


def plate_table(master, rows, cols):
    letters = list(string.ascii_uppercase[:cols])

    for col in range(cols):
        label = tk.Label(master, bd=2, text=letters[col])
        label.grid(row=0, column=col+1, padx=2, pady=2)

    for row in range(rows):
        label = tk.Label(master, bd=2, text=str(row+1))
        label.grid(row=row+1, column=0, padx=2, pady=2)

    for r in range(1, rows+1):
        for c in range(1, cols+1):
            entry_2 = tk.Entry(master, bd=2, width=4)
            entry_2.grid(row=r, column=c, padx=2, pady=2)


def select_groups_screen(previous_root):
    previous_root.withdraw()
    root = tk.Tk()
    root.geometry("425x600")
    root.resizable(height=False, width=False)

    top_frame = tk.Frame(root)
    select_groups_label = tk.Label(master=top_frame, font=14,
                                   text="Please write in every cell the number of group it belongs to ", pady=20)
    names_label = tk.Label(master=top_frame, text=separate_names_to_show())

    select_groups_label.pack()
    names_label.pack()
    top_frame.pack(fill="x")

    table_frame = tk.Frame(root)
    plate_table(table_frame, int(plate_size_entry.get()), int(plate_size_entry_2.get()))
    table_frame.pack(expand=True)
    # table_frame.place(relx=0.5, rely=0.2)

    bottom_frame = tk.Frame(root)

    submit_button = tk.Button(master=bottom_frame, text="  Submit  ", bg='gainsboro',
                              command=lambda: submit_data(previous_root))
    submit_button.grid(ipadx=5, padx=5, pady=2)

    back_button = tk.Button(master=bottom_frame, text="  Back  ", bg='gainsboro',
                            command=lambda: show_preview_screen(root, previous_root))
    back_button.grid(ipadx=5, padx=5, pady=2)

    bottom_frame.pack(expand=True)

    # root.mainloop()


def separate_by_comma(text):
    return text.split(",")


def separate_names_to_show():

    names = separate_by_comma(types_names_entry.get())
    names_str = ""
    for i, name in enumerate(names):
        names_str += str(i) + " - " + name + "  \n"
    return names_str


def show_preview_screen(curr_root, previous_root):
    previous_root.update()
    previous_root.deiconify()
    curr_root.destroy()


def validate_user_input():
    print("validate")
    # make sure that number values are really numbers
    # TODO check that all is filled in second screen
    # TODO check that everything is legal


def submit_data(previous_root):

    # TODO check that all is filled
    # TODO check that everything is legal

    # xls_to_csv("Analysis", "files/csv_file.csv")

    # global full_data_table
    # TODO change to relative path
    # full_data_table = parse_input("C:/Users/Amit/PycharmProjects/lab_project_repo/LabProject/files/csv_file.csv")
    # TODO open a new screen with calculation options

    choose_calculation_screen(previous_root)

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
def choose_calculation_screen(previous_root):
    previous_root.destroy()  # TODO if hide, then need to make sure that when backing it opens the right screen
    root = tk.Tk()
    root.geometry("800x600")
    root.resizable(height=False, width=False)

    # create graph of all data

    top_frame = tk.Frame(root)

    fig = graph_data()

    canvas = FigureCanvasTkAgg(fig, master=top_frame)
    canvas.show()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    top_frame.pack()

    root.mainloop()


def enable_button(*args):
    a = brows_sv.get()
    b = plate_size_x_sv.get()
    c = plate_size_y_sv.get()
    d = num_of_groups_sv.get()
    e = types_names_sv.get()
    f = samples_per_hour_sv.get()
    g = sampling_intervals_sv.get()
    h = total_days_sv.get()
    i = start_time_hour_sv.get()
    j = start_time_minutes_sv.get()

    if a and b and c and d and e and f and g and h and i and j:
        next_button.config(state='normal')
    else:
        next_button.config(state='disabled')


def main():
    root = tk.Tk()

    root.geometry("425x600")
    root.resizable(height=False, width=False)
    root.wm_title("Hello")

    top_frame = tk.Frame(master=root)

    title = tk.Label(master=top_frame, text="Welcome to ... tool!", font=("Arial", 24, "bold"))
    title.pack(fill="x")
    subtitle = tk.Label(master=top_frame, text="Tool for analysis of Zebrafish circadian rhythms", font=("Arial", 15))
    subtitle.pack(fill="x")
    top_frame.pack(fill="x")

    middle_frame = tk.Frame(root)

    logo = tk.PhotoImage(file="sand-time.gif")
    logo_label = tk.Label(master=top_frame, image=logo, justify='center')
    logo_label.image = logo
    logo_label.pack(fill='x', pady=20)

    global browse_entry, plate_size_entry, plate_size_entry_2, num_of_groups_entry, types_names_entry, \
        samples_per_hour_entry, sampling_intervals_entry, total_days_entry, start_time_hour_entry, \
        start_time_minutes_entry, next_button

    global brows_sv, plate_size_x_sv, plate_size_y_sv, num_of_groups_sv, types_names_sv, samples_per_hour_sv, \
        sampling_intervals_sv, total_days_sv, start_time_hour_sv, start_time_minutes_sv

    # create a browse button and locate it
    browse_button = tk.Button(master=middle_frame, text="  Browse  ", bg='gainsboro',
                              command=lambda: select_input_file())
    browse_button.grid(ipadx=5, padx=5, pady=2, sticky=tk.E)

    # create StringVar objects to trace after user input (to then enable the next button)
    brows_sv = tk.StringVar()
    brows_sv.trace("w", enable_button)

    plate_size_x_sv = tk.StringVar()
    plate_size_x_sv.trace("w", enable_button)

    plate_size_y_sv = tk.StringVar()
    plate_size_y_sv.trace("w", enable_button)

    num_of_groups_sv = tk.StringVar()
    num_of_groups_sv.trace("w", enable_button)

    types_names_sv = tk.StringVar()
    types_names_sv.trace("w", enable_button)

    samples_per_hour_sv = tk.StringVar()
    samples_per_hour_sv.trace("w", enable_button)

    sampling_intervals_sv = tk.StringVar()
    sampling_intervals_sv.trace("w", enable_button)

    total_days_sv = tk.StringVar()
    total_days_sv.trace("w", enable_button)

    start_time_hour_sv = tk.StringVar()
    start_time_hour_sv.trace("w", enable_button)

    start_time_minutes_sv = tk.StringVar()
    start_time_minutes_sv.trace("w", enable_button)

    # create the widgets
    browse_entry = tk.Entry(master=middle_frame, bd=2, textvariable=brows_sv)

    plate_size_label = tk.Label(master=middle_frame, text="Plate size: ")
    plate_size_entry = tk.Entry(master=middle_frame, bd=2, width=5, textvariable=plate_size_x_sv)
    plate_size_label_2 = tk.Label(master=middle_frame, text=" x ")
    plate_size_entry_2 = tk.Entry(master=middle_frame, bd=2, width=5, textvariable=plate_size_y_sv)

    num_of_groups_label = tk.Label(master=middle_frame, text="Number of groups: ")
    num_of_groups_entry = tk.Entry(master=middle_frame, bd=2, textvariable=num_of_groups_sv)

    types_names_label = tk.Label(master=middle_frame, text="Groups names: ")
    types_names_entry = tk.Entry(master=middle_frame, bd=2, textvariable=types_names_sv)

    samples_per_hour_label = tk.Label(master=middle_frame, text="Samples per 1 hour: ")
    samples_per_hour_entry = tk.Entry(master=middle_frame, bd=2, textvariable=samples_per_hour_sv)

    sampling_intervals_label = tk.Label(master=middle_frame, text="Sampling intervals (in minutes): ")
    sampling_intervals_entry = tk.Entry(master=middle_frame, bd=2, textvariable=sampling_intervals_sv)

    # ignore_nan_values = False
    check_ignore_nan_iv = tk.IntVar()
    check_ignore_cb = tk.Checkbutton(master=middle_frame, text="Set \'-\' values as \'0\'",
                                     variable=check_ignore_nan_iv, onvalue=0, offvalue=1)

    total_days_label = tk.Label(master=middle_frame, text="Total Number of days of experiment: ")
    total_days_entry = tk.Entry(master=middle_frame, bd=2, textvariable=total_days_sv)

    start_time_label = tk.Label(master=middle_frame, text="Start time: ")
    start_time_hour_label = tk.Label(master=middle_frame, text="hour: ")
    start_time_minutes_label = tk.Label(master=middle_frame, text=" minutes:           ")

    start_time_hour_entry = tk.Entry(master=middle_frame, bd=2, width=5, textvariable=start_time_hour_sv)
    start_time_minutes_entry = tk.Entry(master=middle_frame, bd=2, width=5, textvariable=start_time_minutes_sv)

    # locate all widgets on screen
    browse_entry.grid(row=0, column=1, columnspan=2, sticky=tk.W, pady=2, ipadx=46)
    plate_size_label.grid(row=1, column=0, sticky=tk.E, pady=2)
    plate_size_entry.grid(row=1, column=1, sticky=tk.W, pady=2)
    plate_size_label_2.grid(row=1, column=1, pady=2)
    plate_size_entry_2.grid(row=1, column=1, sticky=tk.E, pady=2)
    num_of_groups_label.grid(row=2, column=0, sticky=tk.E, pady=2)
    num_of_groups_entry.grid(row=2, column=1, pady=2)

    types_names_label.grid(row=3, column=0, sticky=tk.E, pady=2)
    types_names_entry.grid(row=3, column=1, pady=2)

    start_time_label.grid(row=4, column=0, sticky=tk.E, pady=2)
    start_time_hour_label.grid(row=4, column=1, sticky=tk.W, pady=2)
    start_time_hour_entry.grid(row=4, column=1, pady=2)
    start_time_minutes_label.grid(row=4, column=2, sticky=tk.W, padx=2, pady=2)
    start_time_minutes_entry.grid(row=4, column=2, sticky=tk.E, pady=2)

    samples_per_hour_label.grid(row=5, column=0, sticky=tk.E, pady=2)
    samples_per_hour_entry.grid(row=5, column=1, pady=2)

    sampling_intervals_label.grid(row=6, column=0, sticky=tk.E, pady=2)
    sampling_intervals_entry.grid(row=6, column=1, pady=2)

    total_days_label.grid(row=7, column=0, sticky=tk.E, pady=2)
    total_days_entry.grid(row=7, column=1, pady=2)

    check_ignore_cb.select()
    check_ignore_cb.grid(row=8, column=0, columnspan=3)

    middle_frame.pack(fill="x")  # (expand=True)

    bottom_frame = tk.Frame(root)

    next_button = tk.Button(master=bottom_frame, text="  Next  ", font=14, bg='gainsboro', state='disabled',
                            command=lambda: select_groups_screen(root))
    next_button.grid(ipadx=10, pady=10)

    bottom_frame.pack(fill="x", side=tk.BOTTOM)
    bottom_frame.place(relx=0.5, rely=0.9, anchor=tk.CENTER)

    root.mainloop()


if __name__ == "__main__":

    # xls_to_csv("Analysis", "files/csv_file.csv")
    DF = parse_input("C:/Users/Amit/PycharmProjects/lab_project_repo/LabProject/files/csv_file.csv")

    # TODO calc_statistic_values(tables_dict[types_names[i]])
    # graph_data()
    # print(f)

    main()
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

