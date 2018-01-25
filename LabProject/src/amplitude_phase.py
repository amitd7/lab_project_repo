import matplotlib.pyplot as plt
import math
import numpy as np


HOURS_PER_DAY = 24
MINUTES_PER_DAY = 1440

global day_start  # initialized in the main file


def amplitude_phase(measurements, interval_duration, window_size, samples_per_hour, times, circadian_time, diff_time):

    day_index = circadian_time.index(day_start)
    samples_per_day = int(samples_per_hour * HOURS_PER_DAY)

    smoothed_data = moving_average(measurements, window_size)
    diff_time = int(diff_time)
    after_smooth_index = int(day_index - diff_time - math.ceil(window_size/2) + 1)
    smoothed_day_data = smoothed_data[after_smooth_index:after_smooth_index + samples_per_day]

    day_times = times[day_index - diff_time: day_index - diff_time + samples_per_day]
    circadian_time = circadian_time[1:]

    day_ct = circadian_time[day_index: day_index + samples_per_day]
    peak = max(smoothed_day_data)
    trough = min(smoothed_day_data)

    amplitude_value = float("%.2f" % ((peak - trough) / 2))
    phase_value = [day_times[np.argmax(smoothed_day_data)], float("%.2f" % day_ct[np.argmax(smoothed_day_data)])]
    return [amplitude_value, phase_value]


def average_amplitude(from_ct, num_of_days, larva_data, circadian_times, interval_dur, diff_time):
    # the data and the circadian time list have to be of the same size
    # the data is of one larva and is smoothed
    amplitudes, phases = 0, 0
    index = circadian_times.index(from_ct)
    samples_per_day = int(MINUTES_PER_DAY/interval_dur)
    diff_time = int(diff_time)

    for i in range(num_of_days):
        curr_data = larva_data[index - diff_time: index - diff_time + samples_per_day]  # data of one day
        curr_ct = circadian_times[index:index + samples_per_day]

        peak = max(curr_data)
        trough = min(curr_data)
        amplitudes += ((peak - trough) / 2)

        phases += curr_ct[np.argmax(curr_data)]
        index += samples_per_day

    return float("%.2f" % (amplitudes / num_of_days)), float("%.2f" % (phases / num_of_days))


def average_amplitude_wrap(from_ct, num_of_days, data_table, circadian_times, group_names, window, interval_dur,
                           diff_time):

    print("average_amplitude_wrap")
    avg_by_name = {}

    for name in group_names:
        columns_headers = data_table[name].columns.values.tolist()
        avg_val_dict = {"amplitude": [], "phase c.t.": []}

        for i, header in enumerate(columns_headers):
            samples_data = data_table[name][header].tolist()
            smoothed_data = smooth_data(samples_data, window)
            amp, phs = average_amplitude(from_ct, num_of_days, smoothed_data, circadian_times, interval_dur, diff_time)
            avg_val_dict["amplitude"].append(amp)
            avg_val_dict["phase c.t."].append(phs)

        avg_by_name[name] = avg_val_dict

    return avg_by_name


def moving_average(x, window=20):
    """
    Smooths the data in the column vector x using moving average
    :param x: data list
    :param window: The default span for the moving average is 20
    :return: list of smoothed data
    """
    cum_sum = np.cumsum(np.insert(x, 0, [0]))
    return (cum_sum[window:] - cum_sum[:-window]) / window


def smooth_data(data, window):
    index = math.ceil(window / 2)
    return list([None]*(index-1)) + list(moving_average(data, window)) + list([None]*(int(window - index)))
    # return list(np.zeros(index-1)) + list(moving_average(data, window)) + list(np.zeros(int(window - index)))
