from datetime import datetime, timedelta
from matplotlib import style
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import filedialog, ttk, messagebox

import ast
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
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
# style.use("bmh")

# GUI and GUI-related functions:


def select_input_file():
    """

    :return:
    """
    global input_file_path
    input_file_path = filedialog.askopenfile().name
    brows_sv.set(input_file_path)


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

    top_frame = tk.Frame(root)
    select_groups_label = tk.Label(master=top_frame, font=14,
                                   text="Please write in every cell the number of group it belongs to ", pady=20)
    names_label = tk.Label(master=top_frame, font=14, text=separate_names_to_show())

    select_groups_label.pack()
    names_label.pack()
    top_frame.pack(fill="x")

    table_frame = tk.Frame(root)
    global submit_button
    plate_table_values = plate_table(table_frame, int(plate_size_entry.get()), int(plate_size_entry_2.get()))
    table_frame.pack(expand=True)

    loading_frame = tk.Frame(root)

    wait_label = tk.Label(master=loading_frame, font=("Arial", 14, "bold"), foreground="red")
    wait_label.grid(pady=20)

    loading_frame.pack(expand=True)

    bottom_frame = tk.Frame(root)

    submit_button = tk.Button(master=bottom_frame, text="  Submit  ", bg="gainsboro",
                              command=lambda: submit_data(root, plate_table_values, wait_label))
    submit_button.grid(ipadx=5, padx=5, pady=2)

    back_button = tk.Button(master=bottom_frame, text="  Back  ", bg="gainsboro",
                            command=lambda: show_previous_screen(root, previous_root))
    back_button.grid(ipadx=5, padx=5, pady=2)

    bottom_frame.pack(expand=True)

    # root.mainloop()


def is_groups_table_valid(plate_values):
    """
    checks if all entries in the select groups screen are filled
    :param plate_values:
    :return:
    """
    row = len(plate_values)
    col = len(plate_values[0])

    for r in range(row):
        for c in range(col):
            cell_val = plate_values[r][c].get()
            if not cell_val:
                messagebox.showinfo("Error", "The table is not full")
                return False
            elif not cell_val.isdigit():
                messagebox.showinfo("Error", "Value in cell number (%d, %d) is not a number" % (r+1, c+1))
                return False
            elif int(cell_val) >= num_of_groups:
                messagebox.showinfo("Error", "Cell number (%d, %d) has an invalid input. \nValue needs to be a number "
                                             "between 0 to %d" % (r+1, c+1, num_of_groups-1))
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


def clicked_next(previous_root):

    group_names_pattern = "([^:|,]+)(,([^:|,]+))+"  # TODO make sure always work correctly
    start_time_pattern = "(?:[01]\d|2[0123]):(?:[012345]\d)"
    time_bin_pattern = "(?:[01]\d|2[0123]):(?:[012345]\d):(?:[012345]\d)"

    if not plate_size_x_sv.get().isdigit() or not plate_size_y_sv.get().isdigit() or not num_of_groups_sv.get().isdigit() \
       or not check_ignore_sv.get().isdigit():
        tk.messagebox.showinfo("Error", "Input is not a number")
        return
    if (int(plate_size_x_sv.get()) * int(plate_size_y_sv.get())) > 96:  # TODO is the calculation correct?
        tk.messagebox.showinfo("Error", "Plate size is not valid (??)")
        return
    if not re.match(group_names_pattern, types_names_sv.get(), flags=0):
        tk.messagebox.showinfo("Error", "Group names's input is not in the right format")
        return
    if not re.match(start_time_pattern, start_time_sv.get(), flags=0):
        tk.messagebox.showinfo("Error", "Start time has an invalid input")
        return
    if not re.match(time_bin_pattern, sampling_intervals_sv.get(), flags=0):
        tk.messagebox.showinfo("Error", "Time bin has an invalid input")
        return

    set_global_values()
    select_groups_screen(previous_root)


def submit_data(previous_root, plate_values, wait_label):
    """
    after selecting the groups, parse the input file
    :param previous_root:
    :param plate_values:
    :param wait_label:
    :return:
    """
    global full_data_table, c_times, diff_times_to_add
    table_input_validation = is_groups_table_valid(plate_values)

    wait_label.config(text="It might take a while...\nPlease be patient")
    wait_label.update_idletasks()

    if table_input_validation:
        core.xls_to_csv(input_file_path, "Analysis", "C:/Users/Amit/PycharmProjects/lab_project_repo/LabProject/files/csv_file.csv")
        input_file = "C:/Users/Amit/PycharmProjects/lab_project_repo/LabProject/files/csv_file.csv"  # TODO change to relative path

        plate_table_vals = [[val.get() for val in lst] for lst in plate_values]
        wells_names_by_type = separate_by_group(plate_table_vals)
        print("wells_names_by_type: ", wells_names_by_type)
        full_data_table = core.parse_input(input_file, types_names, wells_names_by_type, num_of_groups, start_time,
                                           sampling_intervals, excel_well_labels, excel_time, excel_data,
                                           none_to_num)

        wait_label.config(text="")

        # diff_times_to_add is the number of empty values added to complete total amount that divides by 24
        c_times, diff_times_to_add = core.circadian_time_list(sampling_intervals, ct_zero, rec_start_time_list)
        diff_times_to_add = int(diff_times_to_add)
        print("diff_times_to_add: ", diff_times_to_add)

        choose_calculation_screen(previous_root)
    else:
        return


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
    root.geometry("700x680")
    root.resizable(height=False, width=False)

    top_frame = tk.Frame(root)
    # diff_times_to_add is the number of empty values added to complete total amount that divides by 24
    # c_times, diff_times_to_add = core.circadian_time_list(sampling_intervals, ct_zero, rec_start_time_list)

    # create graph of all data
    fig = core.graph_data(sampling_intervals, c_times, diff_times_to_add, types_names, full_data_table, samples_per_hour)  # todo remark added function input
    canvas = FigureCanvasTkAgg(fig, master=top_frame)
    canvas.show()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    top_frame.pack()

    bottom_frame = tk.Frame(root)

    period_btn = tk.Button(master=bottom_frame, text="  Calculate Period  ", bg='gainsboro',
                           command=lambda: period_g_settings_popup(root, "period"))
    period_btn.grid(row=0, columnspan=3, ipadx=5, padx=5, pady=3)

    g_factor_btn = tk.Button(master=bottom_frame, text="  Calculate G factor  ", bg='gainsboro',
                             command=lambda: period_g_settings_popup(root, "g_factor"))
    g_factor_btn.grid(row=1, columnspan=3, ipadx=5, padx=5, pady=3)

    amplitude_btn = tk.Button(master=bottom_frame, text="  Calculate Amplitude & Phase  ", bg='gainsboro',
                              command=lambda: amplitude_settings_popup(root, "amplitude_phase"))
    amplitude_btn.grid(row=2, columnspan=3, ipadx=5, padx=5, pady=3)

    export_data_btn = tk.Button(master=bottom_frame, text="  Export Raw Data  ", bg='gainsboro',
                                command=lambda: export_raw_data(full_data_table, "raw data"))
    export_data_btn.grid(row=3, columnspan=3, ipadx=5, padx=5, pady=3)

    export_smoothed_data_btn = tk.Button(master=bottom_frame, text="  Export Smoothed Data  ", bg='gainsboro',
                                         command=lambda: export_raw_data(smooth_all_data(full_data_table,
                                                                         smoothed_window_entry.get()), "smoothed data"))
    export_smoothed_data_btn.grid(row=4, column=0, ipadx=5, padx=5, pady=3)

    smoothed_window_entry = tk.Entry(master=bottom_frame, bd=2, width=5)
    smoothed_window_entry.insert(0, 20)
    smoothed_window_entry.grid(row=4, column=1, ipadx=5, pady=3)

    window_after_text = tk.Label(master=bottom_frame, text="(time points window)")
    window_after_text.grid(row=4, column=2, sticky=tk.E, pady=3)

    back_btn = tk.Button(master=bottom_frame, text="  Back  ", bg='gainsboro',
                         command=lambda: show_previous_screen(root, previous_root))
    back_btn.grid(row=5, columnspan=3, ipadx=5, padx=5, pady=3)

    bottom_frame.pack()
    root.mainloop()


def period_g_factor_settings_input_check(previous_root, action, method_type):

    c_times_min = c_times[0]
    c_times_max = c_times[len(c_times)-1]

    if not from_day_entry.get() or not for_days_entry.get():
        messagebox.showinfo("Error", "One or more of the inputs is empty")
        return
    elif action == "period" and (not from_period_entry.get() or not to_period_entry.get()):
        messagebox.showinfo("Error", "One or more of the inputs is empty")
        return
    elif not from_day_entry.get().isdigit() or not for_days_entry.get().isdigit():
        messagebox.showinfo("Error", "One or more of the inputs is not a number")
        return
    elif action == "period" and (not from_period_entry.get().isdigit() or not to_period_entry.get().isdigit()):
        messagebox.showinfo("Error", "One or more of the inputs is empty")
        return
    elif int(from_day_entry.get()) < c_times_min or int(from_day_entry.get()) >= c_times_max:
        messagebox.showinfo("Error", "From day is not in range")
        return
    # elif int(for_days_entry.get()) < 1 or int(for_days_entry.get()) > total_days:  # TODO
    #     messagebox.showinfo("Error", "Number of days is not in range")
    #     return
    elif (int(from_day_entry.get()) + (int(for_days_entry.get()) * HOURS_PER_DAY)) > c_times_max:
        messagebox.showinfo("Error", "Number of days is not in range")
        return
    elif action == "period":
        if int(from_period_entry.get()) < 0 or int(to_period_entry.get()) < 0 or int(from_period_entry.get()) > \
                int(to_period_entry.get()):  # TODO can they have the same value??
            messagebox.showinfo("Error", "Period range is not valid")
            return
    print("period_g_factor_settings_input_check before calc action")
    calc_action(previous_root, action, method_type)


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

    days_label = tk.Label(master=top_frame, text="Start calculate from c.t. ")
    days_label.grid(row=1, column=0, sticky=tk.E, pady=3)

    from_day_label = tk.Label(master=top_frame, text=" for ")
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

    calculate_btn = tk.Button(master=top_frame, text="  Calculate  ", command=lambda:
                              period_g_factor_settings_input_check(period_settings, action, box.get()))
    calculate_btn.grid(row=4, column=2, columnspan=2, pady=10)

    top_frame.pack()
    period_settings.mainloop()


def amp_settings_input_check(previous_root, action, method_type):

    c_times_min = c_times[0]
    c_times_max = c_times[len(c_times)-1]

    if not window_entry.get() or not amp_from_day_entry.get():
        messagebox.showinfo("Error", "One or more of the inputs is empty")
        return
    elif not window_entry.get().isdigit() or not amp_from_day_entry.get().isdigit():
        messagebox.showinfo("Error", "One or more of the inputs is not a number")
        return
    elif int(amp_from_day_entry.get()) < c_times_min or int(amp_from_day_entry.get()) > (c_times_max - HOURS_PER_DAY):
        messagebox.showinfo("Error", "From hour is not in range")
        return
    elif int(window_entry.get()) <= 0 or int(window_entry.get()) > len(c_times):
        messagebox.showinfo("Error", "Window size is not in range")
        return

    calc_action(previous_root, action, method_type)


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

    window_sec_label = tk.Label(master=top_frame, text="(number of time points) ")
    window_sec_label.grid(row=0, column=2, pady=3)

    amp_day_label = tk.Label(master=top_frame, text="Calculate amplitude & phase from c.t. ")
    amp_day_label.grid(row=1, column=0,  pady=3, sticky=tk.E)

    amp_from_day_entry = tk.Entry(master=top_frame, bd=2, width=5)
    amp_from_day_entry.grid(row=1, column=1, pady=3)

    calculate_btn = tk.Button(master=top_frame, text="  Calculate  ", command=lambda:
                              amp_settings_input_check(amplitude_settings, action, None))
    calculate_btn.grid(row=2, column=2, columnspan=3, pady=10)

    days_label = tk.Label(master=top_frame, text="For average amplitude & phase calculate from c.t. ")
    days_label.grid(row=3, column=0, sticky=tk.E, pady=3)

    avg_amp_from_time_entry = tk.Entry(master=top_frame, bd=2, width=5)
    avg_amp_from_time_entry.grid(row=3, column=1, pady=3)

    from_day_label = tk.Label(master=top_frame, text=" for ")
    from_day_label.grid(row=3, column=2, pady=3, sticky=tk.W)

    avg_amp_days_entry = tk.Entry(master=top_frame, bd=2, width=5)
    avg_amp_days_entry.grid(row=3, column=3, pady=2)

    to_day_label = tk.Label(master=top_frame, text=" days ")
    to_day_label.grid(row=3, column=4, pady=3)

    cancel_btn = tk.Button(master=top_frame, text="  Cancel  ", command=lambda: cancel_action(amplitude_settings,
                                                                                              previous_root))
    cancel_btn.grid(row=4, column=0, columnspan=2, pady=10)

    average_amp_calculate_btn = tk.Button(master=top_frame, text="  Calculate   ",
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

    ignore_part_beg, samples_to_ignore_end = 0, 0

    if action == "amplitude_phase":
        core.window_size = int(window_entry.get())
        core.amplitude_phase.day_start = int(amp_from_day_entry.get())
    else:
        hours_to_ignore_beg = int(from_day_entry.get())  # in hours
        period_to_day = int(for_days_entry.get()) * HOURS_PER_DAY  # amount of hours excluding the hours to ignore from the beginning
        # hours_to_ignore_end = (total_days * HOURS_PER_DAY) - period_to_day

        # select the relevant days
        # remove the values that were added to circadian time list to complete 24 hours cycle from ct 0
        # number of samples to ignore from the beginning
        ignore_part_beg = int(hours_to_ignore_beg * samples_per_hour) - diff_times_to_add
        samples_to_ignore_end = int(period_to_day * samples_per_hour) + ignore_part_beg
        print("ignore_part_beg: ", ignore_part_beg)
        print("samples_to_ignore_end: ", samples_to_ignore_end)
        if action == "period":
            action = method_type
            core.pc.from_period = int(from_period_entry.get()) * MINUTES_PER_HOUR
            core.pc.to_period = int(to_period_entry.get()) * MINUTES_PER_HOUR

    results_values = core.groups_calculation(action, types_names, full_data_table, num_of_groups, start_time, c_times,
                                             sampling_intervals, ignore_part_beg, samples_to_ignore_end,
                                             samples_per_hour, diff_times_to_add)
    print("results_values:   ", results_values)

    if action == "fourier" or action == "chi_square":
        # period_results_screen(previous_root, results_values, action, core.period_results_graph, "  T-test  ")
        results_screen(previous_root, results_values, action, "core.period_results_graph", "  T-test  ")
    elif action == "amplitude_phase":
        results_values = break_amp_phase_to_dict(results_values)
        print("results_values in amplitude_phase: ", results_values)
        amp_phase_results_screen(previous_root, results_values, action)
    elif action == "g_factor":
        # g_factor_results_screen(previous_root, results_values, action)
        results_screen(previous_root, results_values, action, "core.g_factor_results_graph", "  KS-test  ")


def average_amp_calc_action(previous_root):

    c_times_min = c_times[0]
    c_times_max = c_times[len(c_times)-1]

    if not avg_amp_from_time_entry.get() or not avg_amp_days_entry.get():
        messagebox.showinfo("Error", "One or more of the inputs is empty")
        return
    elif not avg_amp_from_time_entry.get().isdigit() or not avg_amp_days_entry.get().isdigit():
        messagebox.showinfo("Error", "One or more of the inputs is not a number")
        return
    elif int(avg_amp_from_time_entry.get()) < c_times_min or int(avg_amp_from_time_entry.get()) >= c_times_max:
        messagebox.showinfo("Error", "From day is not in range")
        return
    elif int(avg_amp_days_entry.get()) < 1:  # or int(avg_amp_days_entry.get()) > total_days: TODO
        messagebox.showinfo("Error", "Number of days is not in range")
        return
    elif (int(avg_amp_from_time_entry.get()) + (int(avg_amp_days_entry.get()) * HOURS_PER_DAY)) > c_times_max:
        messagebox.showinfo("Error", "From day and number of days are not in range")
        return

    average_amp_phase_results = core.amplitude_phase.average_amplitude_wrap(int(avg_amp_from_time_entry.get()),
                                                              int(avg_amp_days_entry.get()), full_data_table,
                                                              c_times, types_names, int(window_entry.get()),
                                                              sampling_intervals, diff_times_to_add)
    print("average_amp: ", average_amp_phase_results)
    amp_phase_results_screen(previous_root, average_amp_phase_results, "average_amplitude")
    # results_screen(previous_root, average_amp, "average_amplitude", "core.period_results_graph", "  T-test  ")


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
    p_values = {key: float("%.4f" % p_values[key]) for key in p_values.keys()}
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

    if action == "average_amplitude":
        a_filename = "average amplitude"
        p_filename = "average phase"
    else:
        a_filename = "amplitude"
        p_filename = "phase"

    amp_export_btn = tk.Button(master=bottom_frame, text="  Export amplitude results  ", bg='gainsboro',
                               command=lambda: save_to_excel(amplitude_results, a_filename,
                                                             eval_text_entry(amp_ttest_score.get())))
    amp_export_btn.grid(row=1, column=0, ipadx=5, padx=5, pady=3, columnspan=2)

    phase_export_btn = tk.Button(master=bottom_frame, text="  Export phase results  ", bg='gainsboro',
                                 command=lambda: save_to_excel(full_phase_results, p_filename,
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
        if filename == "phase" or filename == "average phase":
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
    print("mean_value excel : ", mean_value)
    print("mean_value excel type: ", type(mean_value))

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
def export_raw_data(data_tables, filename):
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(filename+".xlsx", engine="xlsxwriter")

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


def smooth_group_data(full_data, name, window=20):

    columns_headers = full_data[name].columns.values.tolist()
    s_data = {}
    for header in columns_headers:
        samples_data = full_data[name][header].tolist()
        s_data[header] = core.amplitude_phase.smooth_data(samples_data, window)
    df = pd.DataFrame.from_dict(s_data)
    # Set DataFrame labels
    df.index = full_data[types_names[0]].index
    return df


def smooth_all_data(data, window):

    # if window.isdigit(): TODO
    window = int(window)
    return {types_names[i]: smooth_group_data(data, types_names[i], window) for i in range(num_of_groups)}


def break_amp_phase_to_dict(data):
    """
    convert to list of triplets inside every key to dictionary
    :param data: as dictionary
    :return:
    """
    full_dict = {}
    for i, name in enumerate(types_names):
        group_dict = {"amplitude": [data[name][i][0] for i in range(len(data[name]))],
                      "phase": [data[name][i][1][0] for i in range(len(data[name]))],
                      "phase c.t.": [data[name][i][1][1] for i in range(len(data[name]))]}
        full_dict[name] = group_dict  # pd.DataFrame.from_dict(group_dict)
    return full_dict


def enable_next_button(*args):
    """

    :param args:
    :return:
    """
    a = brows_sv.get()
    b = plate_size_x_sv.get()
    c = plate_size_y_sv.get()
    d = num_of_groups_sv.get()
    e = types_names_sv.get()
    f = start_time_sv.get()
    g = circadian_start_time_sv.get()
    h = sampling_intervals_sv.get()
    i = check_ignore_sv.get()
    j = excel_well_labels_sv.get()
    k = excel_time_sv.get()
    l = excel_data_sv.get()

    if a and b and c and d and e and f and g and h and i and j and k and l:
        next_button.config(state='normal')
    else:
        next_button.config(state='disabled')


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

    global num_of_groups, types_names, samples_per_hour, sampling_intervals, ct_zero, start_time, excel_well_labels, \
        excel_time, excel_data, none_to_num, rec_start_time_list
    # data from user input
    num_of_groups = int(num_of_groups_sv.get())  # number of mutations and wt groups
    types_names = types_names_sv.get().split(",")

    sampling_intervals = int(time_bin_minutes(separate_time(sampling_intervals_sv.get())))  # in minutes
    samples_per_hour = int(MINUTES_PER_HOUR / sampling_intervals)

    excel_well_labels, excel_time, excel_data = parse_excel_cols_names()
    none_to_num = int(check_ignore_sv.get())

    rec_start_time_list = start_time_sv.get().split(":")
    rec_start_time_list = [int(elem) for elem in rec_start_time_list]

    start_time_hour = int(rec_start_time_list[0])
    start_time_minutes = int(rec_start_time_list[1])

    start_time_minutes += sampling_intervals
    if start_time_minutes == MINUTES_PER_HOUR:
        start_time_hour += 1
        start_time_minutes = 0
    if start_time_hour == HOURS_PER_DAY:
        start_time_hour = 0

    ct_zero = circadian_start_time_sv.get().split(":")
    ct_zero = [int(elem) for elem in ct_zero]
    print("ct_zero in set globals  ", ct_zero)

    start_time = datetime(datetime.now().year, datetime.now().month, datetime.now().day, hour=int(start_time_hour),
                          minute=int(start_time_minutes))


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

    logo = tk.PhotoImage(file="Yoav_lab.gif")
    logo_label = tk.Label(master=top_frame, image=logo, justify='center')
    logo_label.image = logo
    logo_label.pack(fill='x', pady=20)

    top_frame.pack(fill="x")

    global browse_entry, plate_size_entry, plate_size_entry_2, num_of_groups_entry, types_names_entry, \
        sampling_intervals_entry, start_time_entry, next_button  # TODO check if can remove this

    global brows_sv, plate_size_x_sv, plate_size_y_sv, num_of_groups_sv, types_names_sv, \
        sampling_intervals_sv, start_time_sv, circadian_start_time_sv, check_ignore_sv, excel_well_labels_sv, \
        excel_time_sv, excel_data_sv

    middle_frame = tk.Frame(root)

    # create a browse button and locate it
    browse_button = tk.Button(master=middle_frame, text="  Browse  ", bg='gainsboro',
                              command=lambda: select_input_file())
    browse_button.grid(ipadx=5, padx=5, pady=2, sticky=tk.E)

    # create StringVar objects to trace after user input (to then enable the next button)
    brows_sv = tk.StringVar()
    plate_size_x_sv = tk.StringVar()
    plate_size_y_sv = tk.StringVar()
    num_of_groups_sv = tk.StringVar()
    types_names_sv = tk.StringVar()
    sampling_intervals_sv = tk.StringVar()
    start_time_sv = tk.StringVar()
    circadian_start_time_sv = tk.StringVar()
    check_ignore_sv = tk.StringVar()
    excel_well_labels_sv = tk.StringVar()
    excel_time_sv = tk.StringVar()
    excel_data_sv = tk.StringVar()

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

    start_time_label = tk.Label(master=middle_frame, text="Recording start time (hh:mm) ")
    start_time_entry = tk.Entry(master=middle_frame, bd=2, textvariable=start_time_sv)

    circadian_start_time_label = tk.Label(master=middle_frame, text="CT 0 (hh:mm) ")
    circadian_start_time_entry = tk.Entry(master=middle_frame, bd=2, textvariable=circadian_start_time_sv)

    sampling_intervals_label = tk.Label(master=middle_frame, text="Time bin (hh:mm:ss) ")
    sampling_intervals_entry = tk.Entry(master=middle_frame, bd=2, textvariable=sampling_intervals_sv)

    # ignore_nan_values = False
    check_ignore_label = tk.Label(master=middle_frame, text="Set \'-\' values as ")
    check_ignore_entry = tk.Entry(master=middle_frame, bd=2, width=5, textvariable=check_ignore_sv)
    check_ignore_entry.insert(0, 0)

    excel_cols_label = tk.Label(master=middle_frame, text="Relevant Excel columns ")
    excel_well_label = tk.Label(master=middle_frame, text="Arenas ")
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
    plate_size_col_label.grid(row=1, column=1, sticky=tk.E, pady=2)
    plate_size_entry_2.grid(row=1, column=2, sticky=tk.W, pady=2)
    num_of_groups_label.grid(row=2, column=0, sticky=tk.E, pady=2)
    num_of_groups_entry.grid(row=2, column=1, pady=2)

    types_names_label.grid(row=3, column=0, sticky=tk.E, pady=2)
    types_names_entry.grid(row=3, column=1, pady=2)

    start_time_label.grid(row=4, column=0, sticky=tk.E, pady=2)
    start_time_entry.grid(row=4, column=1, pady=2)

    circadian_start_time_label.grid(row=5, column=0, sticky=tk.E, pady=2)
    circadian_start_time_entry.grid(row=5, column=1, pady=2)

    sampling_intervals_label.grid(row=6, column=0, sticky=tk.E, pady=2)
    sampling_intervals_entry.grid(row=6, column=1, pady=2)

    check_ignore_label.grid(row=8, column=0, sticky=tk.E, pady=2)
    check_ignore_entry.grid(row=8, column=1, pady=2)

    excel_cols_label.grid(row=9, column=0, sticky=tk.E, pady=2)
    excel_well_label.grid(row=9, column=1, sticky=tk.W, pady=2)
    excel_time_label.grid(row=10, column=1, sticky=tk.W, pady=2)
    excel_data_label.grid(row=11, column=1, sticky=tk.W, pady=2)
    excel_well_entry.grid(row=9, column=1, pady=2)
    excel_time_entry.grid(row=10, column=1,pady=2)
    excel_data_entry.grid(row=11, column=1, pady=2)

    middle_frame.pack(fill="x")

    bottom_frame = tk.Frame(root)

    next_button = tk.Button(master=bottom_frame, text="  Next  ", font=14, bg='gainsboro', state='disabled',
                            command=lambda: clicked_next(root))
    next_button.grid(ipadx=10, pady=10)

    bottom_frame.pack(fill="x", side=tk.BOTTOM)
    bottom_frame.place(relx=0.5, rely=0.95, anchor=tk.CENTER)

    brows_sv.trace("w", enable_next_button)
    plate_size_x_sv.trace("w", enable_next_button)
    plate_size_y_sv.trace("w", enable_next_button)
    num_of_groups_sv.trace("w", enable_next_button)
    types_names_sv.trace("w", enable_next_button)
    sampling_intervals_sv.trace("w", enable_next_button)
    start_time_sv.trace("w", enable_next_button)
    circadian_start_time_sv.trace("w", enable_next_button)
    check_ignore_sv.trace("w", enable_next_button)
    excel_well_labels_sv.trace("w", enable_next_button)
    excel_time_sv.trace("w", enable_next_button)
    excel_data_sv.trace("w", enable_next_button)

    root.mainloop()










if __name__ == "__main__":

    # xls_to_csv(input_file_path, "Analysis", "files/csv_file.csv")
    # DF = parse_input("C:/Users/Amit/PycharmProjects/lab_project_repo/LabProject/files/csv_file.csv")

    # input_file = "C:/Users/Amit/PycharmProjects/lab_project_repo/LabProject/files/csv_file.csv"
    #
    # start_time = datetime(datetime.now().year, datetime.now().month, datetime.now().day, hour=int(9),
    #                       minute=int(0))
    # wells_names_by_type = [['A2', 'A4', 'A6', 'A8', 'B2', 'B4', 'B6', 'B8', 'C2', 'C4', 'C6', 'C8', 'D1', 'D3', 'D5',
    #                         'D7', 'E1', 'E3', 'E5', 'E7', 'F1', 'F3', 'F5', 'F7'], ['A1', 'A3', 'A5', 'A7', 'B1', 'B3',
    #                                                                                 'B5', 'B7', 'C1', 'C3', 'C5', 'C7',
    #                                                                                 'D2', 'D4', 'D6', 'D8', 'E2', 'E4',
    #                                                                                 'E6', 'E8', 'F2', 'F4', 'F6', 'F8']]
    # group_names = ["wt", "mut"]
    # full_data = core.parse_input(input_file, group_names, wells_names_by_type, 2, start_time, 10, 2, 3, 5, 0)
    # num_of_groups = 2
    #
    # print("smooth_all_data(): ", smooth_all_data())

    main()




