import matplotlib.pyplot as plt
import math
import numpy as np


HOURS_PER_DAY = 24
MINUTES_PER_DAY = 1440

global day_start
# TODO need to get these numbers somehow
sampling_intervals = 10


def amplitude_phase(measurements, interval_duration, window_size, samples_per_hour, times, circadian_time):

    # remove the partial days data from the beginning and from the end
    # as the moving average removed half window_size from the beginning and from the end
    # part_to_slice = int(samples_per_hour * HOURS_PER_DAY - (window_size / 2))
    #
    # spliced_moving_average_data = smoothed_measurements[part_to_slice: -part_to_slice or None]
    #
    # # plt.plot(spliced_moving_average_data)
    #
    # n = len(spliced_moving_average_data)
    # num_of_days = int(n / (samples_per_hour * HOURS_PER_DAY))  # TODO if it not an int?? could it happen?

    # time2 = ["0", "12"] * num_of_days + ["0"]
    #step_value = int((n + 1) / (len(time2) - 1))
    # Set number and labels of ticks for x-axis
    # plt.xticks(np.arange(int(n+1), step=step_value), time2)

    # plt.show()
    a, p = amplitude_phase_calculation(measurements, samples_per_hour, times, circadian_time, window_size)
    # print("amplitude: ", a, "phase: ", p)

    return [a, p]


# def amplitude_phase_calculation(data, day_num_to_use, num_of_days_of_measurements, samples_per_hour, times, ct):
def amplitude_phase_calculation(data, samples_per_hour, times, ct, window):

    day_index = ct.index(day_start)

    samples_per_day = int(samples_per_hour * HOURS_PER_DAY)

    smoothed_data = moving_average(data, window)

    after_smooth_index = day_index - math.ceil(window/2) + 1
    smoothed_day_data = smoothed_data[after_smooth_index:after_smooth_index + samples_per_day]

    # time_last = times[len(times)-1]
    # time_to_add = times[times.index(time_last) + 1]
    # times = times[1:] + [time_to_add]  # TODO will work if there is only one day??????

    day_times = times[day_index: day_index + samples_per_day]

    ct = ct[1:]
    # print("ct  2:  ", ct)
    day_ct = ct[day_index: day_index + samples_per_day]

    peak = max(smoothed_day_data)
    trough = min(smoothed_day_data)

    amplitude_value = float("%.2f" % ((peak - trough) / 2))

    phase_value = [day_times[np.argmax(smoothed_day_data)], float("%.2f" % day_ct[np.argmax(smoothed_day_data)])]
    #
    #
    #
    #
    # # TODO day_num_to_use should not be under 1 and no more then the number of days. if so raise a message
    # if not ((day_num_to_use > 1) and (day_num_to_use < num_of_days_of_measurements)):
    #     print("Error")
    #     # TODO raise a message and try again
    #
    # measurements_per_day = int(samples_per_hour * HOURS_PER_DAY)
    #
    # measurements_to_slice = int((day_num_to_use - 1) * measurements_per_day)
    #
    # times = times[: num_of_days_of_measurements * measurements_per_day]
    # circadian_times = ct[: num_of_days_of_measurements * measurements_per_day]
    #
    # one_day_measurements = data[measurements_to_slice + 1: measurements_to_slice + measurements_per_day + 1]
    # one_day_times = times[measurements_to_slice: measurements_to_slice + measurements_per_day]
    # one_day_circadian_times = circadian_times[measurements_to_slice: measurements_to_slice + measurements_per_day]
    #
    # # plt.plot(one_day_measurements)
    # # plt.show()
    #
    # peak = max(one_day_measurements)
    # trough = min(one_day_measurements)
    #
    # amplitude_value = (peak - trough) / 2
    #
    # phase_value = (one_day_times[np.argmax(one_day_measurements)], one_day_circadian_times[np.argmax(one_day_measurements)])
    return amplitude_value, phase_value


def average_amplitude(from_ct, num_of_days, larva_data, circadian_times):
    # the data and the circadian time list have to be of the same size
    # the data is of one larva and is smoothed
    amplitudes = 0
    index = circadian_times.index(from_ct)
    samples_per_day = int(MINUTES_PER_DAY/sampling_intervals)

    for i in range(num_of_days):
        curr_data = larva_data[index: index + samples_per_day]

        peak = max(curr_data)
        trough = min(curr_data)
        amplitudes += ((peak - trough) / 2)

        index += samples_per_day

    return float("%.2f" % (amplitudes / num_of_days))


def average_amplitude_wrap(from_ct, num_of_days, data_table, circadian_times, group_names, window):

    print("average_amplitude_wrap")
    avg_amp_by_name = {}
    average_amp = []
    for name in group_names:
        columns_headers = data_table[name].columns.values.tolist()
        for i, header in enumerate(columns_headers):
            samples_data = data_table[name][header].tolist()
            smoothed_data = smooth_data(samples_data, window)
            # print("smoothed_data: ", smoothed_data)
            # print("length smoothed_data: ", len(smoothed_data))
            average_amp.append(average_amplitude(from_ct, num_of_days, smoothed_data, circadian_times))
        avg_amp_by_name[name] = average_amp
        average_amp = []
        print(name, " average amplitude:  ", avg_amp_by_name[name])
    return avg_amp_by_name


def moving_average(x, window=20):
    """
    Smooths the data in the column vector x using a moving average
    :param x: data list
    :param window: The default span for the moving average is 20
    :return: list of smoothed data
    """
    cum_sum = np.cumsum(np.insert(x, 0, [0]))
    return (cum_sum[window:] - cum_sum[:-window]) / window


def smooth_data(data, window):
    index = math.ceil(window / 2)
    return list(np.zeros(index-1)) + list(moving_average(data, window)) + list(np.zeros(int(window - index)))
