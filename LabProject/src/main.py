from datetime import datetime, timedelta
from matplotlib import style
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import filedialog
from tkinter import ttk

import ast
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# from scipy import stats
import string
import tkinter as tk

# import src.periodogram_calculation as pc
# from src import amplitude_phase, g_factor, clustering
from src import core

CIRCADIAN_TIME = HOURS_PER_DAY = 24
MINUTES_PER_HOUR = 60
MINUTES_PER_DAY = 1440

style.use("ggplot")
# style.use("seaborn-dark")


# GUI and GUI-related functions:


def select_input_file():
    """

    :return:
    """
    global input_file_path
    input_file_path = filedialog.askopenfile().name
    browse_entry.insert(0, input_file_path)


def plate_table(master, rows, cols):
    """
    Create a table of entries with the size and shape of the experiment plate so the user can mark every cell in which
    group it belongs to
    :param master:
    :param rows:
    :param cols:
    :return:
    """
    letters = list(string.ascii_uppercase[:cols])
    entries_table = [[None for x in range(cols)] for y in range(rows)]

    for col in range(cols):
        label = tk.Label(master, bd=2, text=str(col+1))
        label.grid(row=0, column=col+1, padx=2, pady=2)

    for row in range(rows):
        label = tk.Label(master, bd=2, text=letters[row])
        label.grid(row=row+1, column=0, padx=2, pady=2)

    for r in range(1, rows+1):
        for c in range(1, cols+1):
            entry = tk.Entry(master, bd=2, width=4)
            entries_table[r-1][c-1] = entry
            entry.grid(row=r, column=c, padx=2, pady=2)

    return entries_table


def select_groups_screen(previous_root):
    """
    second screen, associate every larva to it's group
    :param previous_root:
    :return:
    """
    previous_root.withdraw()
    root = tk.Tk()
    root.geometry("425x600")
    root.resizable(height=False, width=False)

    set_global_values()

    top_frame = tk.Frame(root)
    select_groups_label = tk.Label(master=top_frame, font=14,
                                   text="Please write in every cell the number of group it belongs to ", pady=20)
    names_label = tk.Label(master=top_frame, text=separate_names_to_show())

    select_groups_label.pack()
    names_label.pack()
    top_frame.pack(fill="x")

    table_frame = tk.Frame(root)
    global submit_button
    plate_table_values = plate_table(table_frame, int(plate_size_entry.get()), int(plate_size_entry_2.get()))
    table_frame.pack(expand=True)

    bottom_frame = tk.Frame(root)

    submit_button = tk.Button(master=bottom_frame, text="  Submit  ", bg='gainsboro',
                              command=lambda: submit_data(root, plate_table_values))
    submit_button.grid(ipadx=5, padx=5, pady=2)

    back_button = tk.Button(master=bottom_frame, text="  Back  ", bg='gainsboro',
                            command=lambda: show_previous_screen(root, previous_root))
    back_button.grid(ipadx=5, padx=5, pady=2)

    bottom_frame.pack(expand=True)

    # root.mainloop()


def is_groups_table_full(plate_values):
    """
    checks if all entries in the select groups screen are filled
    :param plate_values:
    :return:
    """
    row = len(plate_values)
    col = len(plate_values[0])

    for r in range(row):
        for c in range(col):
            if not plate_values[r][c].get():
                print("empty")
                return False
    return True


def separate_names_to_show():
    """
    create labels of groups and groups numbers to show for easier group selection
    :return:
    """
    names = types_names_entry.get().split(",")
    names_str = ""
    for i, name in enumerate(names):
        names_str += str(i) + " - " + name + "  \n"
    return names_str


def show_previous_screen(curr_root, previous_root):
    """
    show previous screen and destroy current screen
    :param curr_root:
    :param previous_root:
    :return:
    """
    previous_root.update()
    previous_root.deiconify()
    curr_root.destroy()


# def validate_user_input():
#     """
#
#     :return:
#     """
#     print("validate")
#     # make sure that number values are really numbers
#     # TODO check that all is filled in second screen
#     # TODO check that everything is legal


def submit_data(previous_root, plate_values):
    """
    after selecting the groups, parse the input file
    :param previous_root:
    :param plate_values:
    :return:
    """

    # TODO check that all is filled
    # TODO check that everything is legal
    core.xls_to_csv(input_file_path, "Analysis", "C:/Users/Amit/PycharmProjects/lab_project_repo/LabProject/files/csv_file.csv")
    input_file = "C:/Users/Amit/PycharmProjects/lab_project_repo/LabProject/files/csv_file.csv"  # TODO change to relative path

    global full_data_table
    global wells_names_by_type  # TODO should define in another place??
    if is_groups_table_full(plate_values):
        plate_table_vals = [[val.get() for val in lst] for lst in plate_values]
        wells_names_by_type = separate_by_group(plate_table_vals)

        full_data_table = core.parse_input(input_file, types_names, wells_names_by_type, num_of_groups, start_time,
                                           total_days, sampling_intervals, excel_well_labels, excel_time, excel_data,
                                           none_to_num)
        choose_calculation_screen(previous_root)


def group_table_to_names(table_vals):
    """
    Create a table of well names at the same size as the zebrafish plate
    :param table_vals: the table received from the user for choosing groups
    :return:
    """
    letters = list(string.ascii_uppercase)
    numbers = range(len(letters))
    row = len(table_vals)
    col = len(table_vals[0])
    return [[str(letters[i]) + str(numbers[j+1]) for j in range(col)]for i in range(row)]


def separate_by_group(table_vals):
    """
    Separate the wells by groups.
    Creates a nested list of the wells names of each experiment group
    :param table_vals:
    :return:
    """

    names_table = group_table_to_names(table_vals)

    row = len(table_vals)
    col = len(table_vals[0])

    wells_by_group = [[] for x in range(num_of_groups)]

    for i in range(row):
        for j in range(col):
            wells_by_group[int(table_vals[i][j])].append(str(names_table[i][j]))
    return wells_by_group


def choose_calculation_screen(previous_root):
    """
    create a screen with the calculation options for the user to choose
    :param previous_root:
    :return:
    """

    previous_root.withdraw()  # TODO if hide, then need to make sure that when backing it opens the right screen
    root = tk.Tk()
    root.geometry("640x640")
    root.resizable(height=False, width=False)

    top_frame = tk.Frame(root)
    # create graph of all data
    fig = core.graph_data(sampling_intervals, total_days, types_names, full_data_table, samples_per_hour)  # todo remark added function input
    canvas = FigureCanvasTkAgg(fig, master=top_frame)
    canvas.show()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    top_frame.pack()

    bottom_frame = tk.Frame(root)

    period_btn = tk.Button(master=bottom_frame, text="  Calculate Period  ", bg='gainsboro',
                           command=lambda: period_g_settings_popup(root, "period"))
    period_btn.grid(ipadx=5, padx=5, pady=3)

    g_factor_btn = tk.Button(master=bottom_frame, text="  Calculate G factor  ", bg='gainsboro',
                             command=lambda: period_g_settings_popup(root, "g_factor"))
    g_factor_btn.grid(ipadx=5, padx=5, pady=3)

    amplitude_btn = tk.Button(master=bottom_frame, text="  Calculate Amplitude & Phase  ", bg='gainsboro',
                              command=lambda: amplitude_settings_popup(root, "amplitude_phase"))
    amplitude_btn.grid(ipadx=5, padx=5, pady=3)

    # clustering_btn = tk.Button(master=bottom_frame, text="  Clustering  ", bg='gainsboro')
    # clustering_btn.grid(ipadx=5, padx=5, pady=3)

    export_data_btn = tk.Button(master=bottom_frame, text="  Export Raw Data  ", bg='gainsboro',
                                command=lambda: export_raw_data(full_data_table))
    export_data_btn.grid(ipadx=5, padx=5, pady=3)

    back_btn = tk.Button(master=bottom_frame, text="  Back  ", bg='gainsboro',
                         command=lambda: show_previous_screen(root, previous_root))
    back_btn.grid(ipadx=5, padx=5, pady=3)

    bottom_frame.pack()
    root.mainloop()


def period_g_settings_popup(previous_root, action):
    """

    :param previous_root:
    :param action:
    :return:
    """
    period_settings = tk.Toplevel(master=previous_root)
    period_settings.wm_title("Settings")

    period_settings.grab_set()  # set the focus on the popup only
    print("action: ", action)
    top_frame = tk.Frame(master=period_settings)

    global from_day_entry, for_days_entry, from_period_entry, to_period_entry

    value = tk.StringVar()
    box = ttk.Combobox(top_frame, textvariable=value, state='readonly')  # TODO can I remove the value variable?
    if action == "period":
        methods_label = tk.Label(master=top_frame, text="Choose method: ")
        methods_label.grid(row=0, column=0, sticky=tk.E, pady=3)

        methods = ["fourier", "chi_square"]

        box["values"] = methods
        box.current(0)
        box.grid(row=0, column=1, padx=2, pady=3, columnspan=4, sticky=tk.W)

        range_label = tk.Label(master=top_frame, text="Periods range (hours): ")
        range_label.grid(row=2, column=0, sticky=tk.E, pady=3)

        from_period_label = tk.Label(master=top_frame, text="from ")
        from_period_label.grid(row=2, column=1, pady=3, sticky=tk.W)

        from_period_entry = tk.Entry(master=top_frame, bd=3, width=5)
        from_period_entry.insert(0, 16)
        from_period_entry.grid(row=2, column=2, pady=3)

        to_period_label = tk.Label(master=top_frame, text=" to ")
        to_period_label.grid(row=2, column=3, pady=3)

        to_period_entry = tk.Entry(master=top_frame, bd=2, width=5)
        to_period_entry.insert(0, 32)
        to_period_entry.grid(row=2, column=4, pady=3)

    days_label = tk.Label(master=top_frame, text="Start calculate from ")
    days_label.grid(row=1, column=0, sticky=tk.E, pady=3)

    from_day_label = tk.Label(master=top_frame, text="hours for ")
    from_day_label.grid(row=1, column=2, pady=3, sticky=tk.W)

    from_day_entry = tk.Entry(master=top_frame, bd=2, width=5)
    from_day_entry.grid(row=1, column=1, pady=3)

    for_days_entry = tk.Entry(master=top_frame, bd=2, width=5)
    for_days_entry.grid(row=1, column=3, pady=3)

    to_day_label = tk.Label(master=top_frame, text=" days ")
    to_day_label.grid(row=1, column=4, pady=3)

    cancel_btn = tk.Button(master=top_frame, text="  Cancel  ", command=lambda: cancel_action(period_settings,
                                                                                              previous_root))
    cancel_btn.grid(row=4, column=0, columnspan=2, pady=10)

    calculate_btn = tk.Button(master=top_frame, text="  Calculate  ", command=lambda: calc_action(period_settings,
                                                                                                  action, box.get()))
    calculate_btn.grid(row=4, column=2, columnspan=2, pady=10)

    top_frame.pack()
    period_settings.mainloop()


def amplitude_settings_popup(previous_root, action):
    """

    :param previous_root:
    :param action:
    :return:
    """
    amplitude_settings = tk.Toplevel(master=previous_root)
    amplitude_settings.wm_title("Settings")

    amplitude_settings.grab_set()  # set the focus on the popup only
    print("action: ", action)
    top_frame = tk.Frame(master=amplitude_settings)

    global window_entry, amp_from_day_entry, avg_amp_from_time_entry, avg_amp_days_entry

    window_label = tk.Label(master=top_frame, text="Moving average window size: ")
    window_label.grid(row=0, column=0, pady=3, sticky=tk.E)

    window_entry = tk.Entry(master=top_frame, bd=3, width=5)
    window_entry.insert(0, 20)
    window_entry.grid(row=0, column=1, pady=3, sticky=tk.W)

    amp_day_label = tk.Label(master=top_frame, text="Calculate amplitude from ")
    amp_day_label.grid(row=1, column=0,  pady=3, sticky=tk.E)

    # amp_from_day_label = tk.Label(master=top_frame, text="from ")
    # amp_from_day_label.grid(row=1, column=1, pady=3, sticky=tk.W)

    amp_from_day_entry = tk.Entry(master=top_frame, bd=2, width=5)
    amp_from_day_entry.grid(row=1, column=1, pady=3)

    amp_to_day_label = tk.Label(master=top_frame, text="hours ")
    amp_to_day_label.grid(row=1, column=2, pady=3)

    calculate_btn = tk.Button(master=top_frame, text="  Calculate  ", command=lambda: calc_action(amplitude_settings,
                                                                                                  action, None))
    calculate_btn.grid(row=2, column=2, columnspan=3, pady=10)

    days_label = tk.Label(master=top_frame, text="For average amplitude calculate from ")
    days_label.grid(row=3, column=0, sticky=tk.E, pady=3)

    avg_amp_from_time_entry = tk.Entry(master=top_frame, bd=2, width=5)
    avg_amp_from_time_entry.grid(row=3, column=1, pady=3)

    from_day_label = tk.Label(master=top_frame, text="hours for ")
    from_day_label.grid(row=3, column=2, pady=3, sticky=tk.W)

    avg_amp_days_entry = tk.Entry(master=top_frame, bd=2, width=5)
    avg_amp_days_entry.grid(row=3, column=3, pady=2)

    to_day_label = tk.Label(master=top_frame, text=" days ")
    to_day_label.grid(row=3, column=4, pady=3)

    cancel_btn = tk.Button(master=top_frame, text="  Cancel  ", command=lambda: cancel_action(amplitude_settings,
                                                                                              previous_root))
    cancel_btn.grid(row=4, column=0, columnspan=2, pady=10)

    average_amp_calculate_btn = tk.Button(master=top_frame, text="  Calculate average amplitude   ",
                                          command=lambda: average_amp_calc_action(amplitude_settings))
    average_amp_calculate_btn.grid(row=4, column=2, columnspan=3, pady=10)

    top_frame.pack()
    amplitude_settings.mainloop()


def calc_action(previous_root, action, method_type):
    """

    :param previous_root:
    :param action:
    :param method_type:
    :return:
    """
    print("method:  ", method_type)

    ignore_part_beg = 0
    ignore_part_end = 0

    if action == "amplitude_phase":
        core.window_size = int(window_entry.get())
        core.amplitude_phase.day_start = int(amp_from_day_entry.get())
    else:
        hours_to_ignore_beg = int(from_day_entry.get())
        period_to_day = int(for_days_entry.get()) * HOURS_PER_DAY + hours_to_ignore_beg
        hours_to_ignore_end = (total_days * HOURS_PER_DAY) - period_to_day

        # select the relevant days
        ignore_part_beg = int(hours_to_ignore_beg * samples_per_hour)  # number of samples to ignore from the beginning
        ignore_part_end = int(hours_to_ignore_end * samples_per_hour)  # number of samples to ignore from the end
        if action == "period":
            action = method_type
            core.pc.from_period = int(from_period_entry.get()) * MINUTES_PER_HOUR
            core.pc.to_period = int(to_period_entry.get()) * MINUTES_PER_HOUR

    results_values = core.groups_calculation(action, types_names, full_data_table, num_of_groups, start_time,
                                             total_days, sampling_intervals, ignore_part_beg, ignore_part_end,
                                             samples_per_hour)
    print("results_values:   ", results_values)

    if action == "fourier" or action == "chi_square":
        # period_results_screen(previous_root, results_values, action, core.period_results_graph, "  T-test  ")
        results_screen(previous_root, results_values, action, "core.period_results_graph", "  T-test  ")
    elif action == "amplitude_phase":
        results_values = break_amp_phase_to_dict(results_values)
        amp_phase_results_screen(previous_root, results_values, action)
    elif action == "g_factor":
        # g_factor_results_screen(previous_root, results_values, action)
        results_screen(previous_root, results_values, action, "core.g_factor_results_graph", "  KS-test  ")


def average_amp_calc_action(previous_root):

    print("in average_amp_calc_action")
    c_times = core.circadian_time_list(sampling_intervals, total_days)
    average_amp = core.amplitude_phase.average_amplitude_wrap(int(avg_amp_from_time_entry.get()),
                                                              int(avg_amp_days_entry.get()), full_data_table,
                                                              c_times, types_names, int(window_entry.get()))
    print("average_amp: ", average_amp)

    results_screen(previous_root, average_amp, "average_amplitude", "core.period_results_graph", "  T-test  ")


def cancel_action(root, previous_root):
    """

    :param root:
    :param previous_root:
    :return:
    """
    root.destroy()
    previous_root.grab_release()


def eval_text_entry(entry_str):
    if not entry_str:
        return None
    else:
            return ast.literal_eval(entry_str)


def groups_stat_tests_call(entry, results_values, action):

    p_values = core.manage_statistic_tests(results_values, action, types_names)
    entry.configure(state=tk.NORMAL)
    entry.delete(0, tk.END)
    entry.insert(tk.INSERT, str(p_values))
    entry.configure(state="readonly")


def results_screen(previous_root, results_values, action, func_name, btn_name):

    # previous_root.destroy()
    root = tk.Tk()
    root.geometry("640x580")
    root.resizable(height=False, width=False)

    top_frame = tk.Frame(root)
    # create graph of results
    # fig = core.period_results_graph(results_values, types_names, action)
    fig = eval(func_name)(results_values, types_names, action)
    canvas = FigureCanvasTkAgg(fig, master=top_frame)
    canvas.show()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    top_frame.pack()

    bottom_frame = tk.Frame(root)

    test_score = tk.Entry(master=bottom_frame, width=30)
    test_score.grid(row=0, column=1, ipadx=5, padx=5, pady=3)
    test_score.configure(state="readonly")

    test_btn = tk.Button(master=bottom_frame, text=btn_name, bg='gainsboro',
                         command=lambda: groups_stat_tests_call(test_score, results_values, action))
    test_btn.grid(row=0, column=0, ipadx=5, padx=5, pady=3)

    export_btn = tk.Button(master=bottom_frame, text="  Export results  ", bg='gainsboro',
                           command=lambda: save_to_excel(results_values, action, eval_text_entry(test_score.get())))
    export_btn.grid(ipadx=5, padx=5, pady=3, columnspan=2)

    bottom_frame.pack()
    root.mainloop()


def amp_phase_results_screen(previous_root, results_values, action):

    # previous_root.destroy()
    root = tk.Tk()
    root.geometry("800x580")
    root.resizable(height=False, width=False)

    top_frame = tk.Frame(root)
    # create graph of results
    fig = core.amp_phase_results_graph(results_values, types_names, action, action)
    canvas = FigureCanvasTkAgg(fig, master=top_frame)
    canvas.show()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    top_frame.pack()

    bottom_frame = tk.Frame(root)

    amplitude_results = {key: results_values[key]["amplitude"] for key in results_values}
    phase_results = {key: results_values[key]["phase c.t."] for key in results_values}
    full_phase_results = results_values
    for key in full_phase_results.keys():
        del full_phase_results[key]["amplitude"]

    amp_ttest_score = tk.Entry(master=bottom_frame, width=30)
    amp_ttest_score.grid(row=0, column=1, ipadx=5, padx=5, pady=3)
    amp_ttest_score.configure(state="readonly")

    phase_ttest_score = tk.Entry(master=bottom_frame, width=30)
    phase_ttest_score.grid(row=0, column=3, ipadx=5, padx=5, pady=3)
    phase_ttest_score.configure(state="readonly")

    amp_ttest_btn = tk.Button(master=bottom_frame, text="  Amplitude T-test  ", bg='gainsboro',
                              command=lambda: groups_stat_tests_call(amp_ttest_score, amplitude_results, action))
    amp_ttest_btn.grid(row=0, column=0, ipadx=5, padx=5, pady=3)

    phase_ttest_btn = tk.Button(master=bottom_frame, text="  Phase T-test  ", bg='gainsboro',
                                command=lambda: groups_stat_tests_call(phase_ttest_score, phase_results, action))
    phase_ttest_btn.grid(row=0, column=2, ipadx=5, padx=5, pady=3)

    amp_export_btn = tk.Button(master=bottom_frame, text="  Export amplitude results  ", bg='gainsboro',
                               command=lambda: save_to_excel(amplitude_results, "amplitude",
                                                             eval_text_entry(amp_ttest_score.get())))
    amp_export_btn.grid(row=1, column=0, ipadx=5, padx=5, pady=3, columnspan=2)

    phase_export_btn = tk.Button(master=bottom_frame, text="  Export phase results  ", bg='gainsboro',
                                 command=lambda: save_to_excel(full_phase_results, "phase",
                                                               eval_text_entry(phase_ttest_score.get())))
    phase_export_btn.grid(row=1, column=2, ipadx=5, padx=5, pady=3, columnspan=2)

    bottom_frame.pack()
    root.mainloop()


def save_to_excel(data_to_write, filename, p_values_d):

    # data_to_write has to be a dictionary
    results_df_for_export = pd.DataFrame()
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(filename+".xlsx", engine="xlsxwriter")
    n = 0
    for i, name in enumerate(types_names):
        if filename == "phase":
            for key in data_to_write[name]:
                results_df_for_export[name+" "+key] = data_to_write[name][key]
                if n < len(results_df_for_export[name+" "+key]):
                    n = len(results_df_for_export[name+" "+key])
        else:
            results_df_for_export[name] = data_to_write[name]
            if n < len(results_df_for_export[name]):
                n = len(results_df_for_export[name])

    # Change the index numbers (rows names) to start from 1 and not from 0 as in default
    results_df_for_export.index = range(1, n + 1)

    mean_value, std_value, se_value = core.calc_statistic_values(results_df_for_export, 0)  # TODO 0 for calc on cols
    print("mean_value: ", mean_value)
    print("mean_value type: ", type(mean_value))

    mean_value.name = "average"
    std_value.name = "std"
    se_value.name = "se"
    results_df_for_export = results_df_for_export.append(mean_value)
    results_df_for_export = results_df_for_export.append(std_value)
    results_df_for_export = results_df_for_export.append(se_value)
    if p_values_d is not None:
        p_value_df = pd.DataFrame(list(p_values_d.items()))
        p_value_df.index = ["p value"] * len(p_values_d)

        results_df_for_export = results_df_for_export.append(p_value_df)

    print("df:  ", results_df_for_export)

    # Convert the dataframe to an XlsxWriter Excel object.
    results_df_for_export.to_excel(writer, sheet_name=filename)

    # Close the Pandas Excel writer and output the Excel file.
    writer.save()


# TODO full_data_table is global but is not recognized here
def export_raw_data(data_tables):
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter("raw_data.xlsx", engine="xlsxwriter")

    for name in types_names:
        mean_value, std_value, se_value = core.calc_statistic_values(data_tables[name], 1)
        # add values to table
        data_tables[name]["average"] = mean_value
        data_tables[name]["std"] = std_value
        data_tables[name]["s. error"] = se_value

        # Convert the dataframe to an XlsxWriter Excel object.
        data_tables[name].to_excel(writer, sheet_name=name)

    # Close the Pandas Excel writer and output the Excel file.
    writer.save()


def break_amp_phase_to_dict(data):
    """
    convert to list of triplets inside every key to dictionary
    :param data: as dictionary
    :return:
    """
    full_dict = {}
    for i, name in enumerate(types_names):
        group_dict = {"amplitude": [data[name][i][0] for i in range(len(data[name]))]}  # TODO how to turn it to dict literal?
        group_dict["phase"] = [data[name][i][1][0] for i in range(len(data[name]))]
        group_dict["phase c.t."] = [data[name][i][1][1] for i in range(len(data[name]))]
        full_dict[name] = group_dict  # pd.DataFrame.from_dict(group_dict)
    return full_dict


def enable_button(*args):
    """

    :param args:
    :return:
    """
    print("in")
    # a = brows_sv.get()
    # b = plate_size_x_sv.get()
    # c = plate_size_y_sv.get()
    #
    # if a and b and c and num_of_groups and types_names and samples_per_hour and sampling_intervals and total_days and \
    #         start_time_hour and start_time_minutes:
    #     next_button.config(state='normal')
    # else:
    #     next_button.config(state='disabled')


def separate_time(time):

    return time.split(":")


def time_bin_minutes(lst):
    h = 0
    m = 1
    s = 2

    # TODO is 00:10:30 an option or it will be just one option filled?
    if lst[s] == "00" and lst[h] == "00":
        return int(lst[m])
    elif lst[m] == "00" and lst[h] == "00":
        return int(lst[s]) / MINUTES_PER_HOUR  # TODO may return float value, will everything work with float?


def parse_excel_cols_names():

    e_well_labels = excel_well_labels_sv.get()
    e_time = excel_time_sv.get()
    e_data = excel_data_sv.get()

    if not e_well_labels.isdigit():
        e_well_labels = ord(e_well_labels.lower()) % 32
    else:
        e_well_labels = int(e_well_labels)
    if not e_time.isdigit():
        e_time = ord(e_time.lower()) % 32
    else:
        e_time = int(e_time)
    if not e_data.isdigit():
        e_data = ord(e_data.lower()) % 32
    else:
        e_data = int(e_data)

    # Indexing in excel columns start from zero
    e_well_labels -= 1
    e_time -= 1
    e_data -= 1

    return e_well_labels, e_time, e_data


def set_global_values():
    """

    :return:
    """

    global num_of_groups, types_names, samples_per_hour, sampling_intervals, total_days, start_time_hour, \
        start_time_minutes, start_time, excel_well_labels, excel_time, excel_data, none_to_num  #, total_num_of_values
    # data from user input
    num_of_groups = int(num_of_groups_sv.get())  # number of mutations and wt groups
    types_names = types_names_sv.get().split(",")

    sampling_intervals = int(time_bin_minutes(separate_time(sampling_intervals_sv.get())))  # in minutes
    samples_per_hour = int(MINUTES_PER_HOUR / sampling_intervals)
    total_days = int(total_days_sv.get())

    excel_well_labels, excel_time, excel_data = parse_excel_cols_names()
    none_to_num = int(check_ignore_sv.get())

    start_time_separated = separate_time(start_time_sv.get())
    start_time_hour = int(start_time_separated[0])
    start_time_minutes = int(start_time_separated[1])

    start_time_minutes += sampling_intervals
    if start_time_minutes == MINUTES_PER_HOUR:
        start_time_hour += 1
        start_time_minutes = 0
    if start_time_hour == HOURS_PER_DAY:
        start_time_hour = 0

    start_time = datetime(datetime.now().year, datetime.now().month, datetime.now().day, hour=int(start_time_hour),
                          minute=int(start_time_minutes))

    # total_num_of_values = int(MINUTES_PER_DAY / sampling_intervals * total_days)


def main():
    root = tk.Tk()

    root.geometry("425x600")
    root.resizable(height=False, width=False)
    root.wm_title("Hello")

    top_frame = tk.Frame(master=root)

    title = tk.Label(master=top_frame, text="Welcome to Amitool!", font=("Arial", 24, "bold"))
    title.pack(fill="x")
    subtitle = tk.Label(master=top_frame, text="Zebrafish Circadian Rhythms Analysis Tool", font=("Arial", 15))
    subtitle.pack(fill="x")
    top_frame.pack(fill="x")

    middle_frame = tk.Frame(root)

    logo = tk.PhotoImage(file="Yoav_lab.gif")
    logo_label = tk.Label(master=top_frame, image=logo, justify='center')
    logo_label.image = logo
    logo_label.pack(fill='x', pady=20)

    global browse_entry, plate_size_entry, plate_size_entry_2, num_of_groups_entry, types_names_entry, \
        sampling_intervals_entry, total_days_entry, start_time_entry, next_button  # TODO check if can remove this

    global brows_sv, plate_size_x_sv, plate_size_y_sv, num_of_groups_sv, types_names_sv, samples_per_hour_sv, \
        sampling_intervals_sv, total_days_sv, start_time_sv, check_ignore_sv, excel_well_labels_sv, excel_time_sv, \
        excel_data_sv

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

    start_time_sv = tk.StringVar()
    start_time_sv.trace("w", enable_button)

    check_ignore_sv = tk.StringVar()
    check_ignore_sv.trace("w", enable_button)

    excel_well_labels_sv = tk.StringVar()
    excel_well_labels_sv.trace("w", enable_button)
    excel_time_sv = tk.StringVar()
    excel_time_sv.trace("w", enable_button)
    excel_data_sv = tk.StringVar()
    excel_data_sv.trace("w", enable_button)

    # create the widgets
    browse_entry = tk.Entry(master=middle_frame, bd=2, textvariable=brows_sv)

    plate_size_label = tk.Label(master=middle_frame, text="Plate size ")
    plate_size_row_label = tk.Label(master=middle_frame, text="row ")
    plate_size_entry = tk.Entry(master=middle_frame, bd=2, width=5, textvariable=plate_size_x_sv)
    # plate_size_label_2 = tk.Label(master=middle_frame, text=" x ")
    plate_size_col_label = tk.Label(master=middle_frame, text="column ")
    plate_size_entry_2 = tk.Entry(master=middle_frame, bd=2, width=5, textvariable=plate_size_y_sv)

    num_of_groups_label = tk.Label(master=middle_frame, text="Number of groups ")
    num_of_groups_entry = tk.Entry(master=middle_frame, bd=2, textvariable=num_of_groups_sv)

    types_names_label = tk.Label(master=middle_frame, text="Groups names (name1,name2,...) ")
    types_names_entry = tk.Entry(master=middle_frame, bd=2, textvariable=types_names_sv)

    sampling_intervals_label = tk.Label(master=middle_frame, text="Time bin (hh:mm:ss) ")
    sampling_intervals_entry = tk.Entry(master=middle_frame, bd=2, textvariable=sampling_intervals_sv)

    # ignore_nan_values = False
    check_ignore_label = tk.Label(master=middle_frame, text="Set \'-\' values as ")
    check_ignore_entry = tk.Entry(master=middle_frame, bd=2, width=5, textvariable=check_ignore_sv)
    check_ignore_entry.insert(0, 0)

    total_days_label = tk.Label(master=middle_frame, text="Total Number of days of experiment ")
    total_days_entry = tk.Entry(master=middle_frame, bd=2, textvariable=total_days_sv)

    start_time_label = tk.Label(master=middle_frame, text="Start time (hh:mm) ")
    start_time_entry = tk.Entry(master=middle_frame, bd=2, textvariable=start_time_sv)

    excel_cols_label = tk.Label(master=middle_frame, text="Relevant Excel columns ")
    excel_well_label = tk.Label(master=middle_frame, text="Labels ")
    excel_time_label = tk.Label(master=middle_frame, text="Time ")
    excel_data_label = tk.Label(master=middle_frame, text="Data ")
    excel_well_entry = tk.Entry(master=middle_frame, bd=2, width=5, textvariable=excel_well_labels_sv)
    excel_time_entry = tk.Entry(master=middle_frame, bd=2, width=5, textvariable=excel_time_sv)
    excel_data_entry = tk.Entry(master=middle_frame, bd=2, width=5, textvariable=excel_data_sv)

    # locate all widgets on screen
    browse_entry.grid(row=0, column=1, columnspan=2, sticky=tk.W, pady=2, ipadx=46)
    plate_size_label.grid(row=1, column=0, sticky=tk.E, pady=2)
    plate_size_row_label.grid(row=1, column=1, sticky=tk.W, pady=2)
    plate_size_entry.grid(row=1, column=1, pady=2)
    # plate_size_label_2.grid(row=1, column=1, pady=2)
    plate_size_col_label.grid(row=1, column=1, sticky=tk.E, pady=2)
    plate_size_entry_2.grid(row=1, column=2, sticky=tk.W, pady=2)
    num_of_groups_label.grid(row=2, column=0, sticky=tk.E, pady=2)
    num_of_groups_entry.grid(row=2, column=1, pady=2)

    types_names_label.grid(row=3, column=0, sticky=tk.E, pady=2)
    types_names_entry.grid(row=3, column=1, pady=2)

    start_time_label.grid(row=4, column=0, sticky=tk.E, pady=2)
    start_time_entry.grid(row=4, column=1, pady=2)

    sampling_intervals_label.grid(row=6, column=0, sticky=tk.E, pady=2)
    sampling_intervals_entry.grid(row=6, column=1, pady=2)

    total_days_label.grid(row=7, column=0, sticky=tk.E, pady=2)
    total_days_entry.grid(row=7, column=1, pady=2)

    check_ignore_label.grid(row=8, column=0, sticky=tk.E, pady=2)
    check_ignore_entry.grid(row=8, column=1, pady=2)

    excel_cols_label.grid(row=9, column=0, sticky=tk.E, pady=2)
    excel_well_label.grid(row=9, column=1, sticky=tk.W, pady=2)
    excel_time_label.grid(row=10, column=1, sticky=tk.W, pady=2)
    excel_data_label.grid(row=11, column=1, sticky=tk.W, pady=2)
    excel_well_entry.grid(row=9, column=1, pady=2)
    excel_time_entry.grid(row=10, column=1,pady=2)
    excel_data_entry.grid(row=11, column=1, pady=2)

    middle_frame.pack(fill="x")  # (expand=True)

    bottom_frame = tk.Frame(root)

    next_button = tk.Button(master=bottom_frame, text="  Next  ", font=14, bg='gainsboro',# state='disabled',
                            command=lambda: select_groups_screen(root))
    next_button.grid(ipadx=10, pady=10)

    bottom_frame.pack(fill="x", side=tk.BOTTOM)
    bottom_frame.place(relx=0.5, rely=0.95, anchor=tk.CENTER)

    root.mainloop()


if __name__ == "__main__":

    # xls_to_csv(input_file_path, "Analysis", "files/csv_file.csv")
    # DF = parse_input("C:/Users/Amit/PycharmProjects/lab_project_repo/LabProject/files/csv_file.csv")


    main()



