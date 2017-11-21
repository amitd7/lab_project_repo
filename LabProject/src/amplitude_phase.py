import matplotlib.pyplot as plt
import numpy as np


HOURS_PER_DAY = 24
MINUTES_PER_DAY = 1440

# TODO need to get these numbers somehow
sampling_intervals = 10
total_days = 5


def amplitude_phase(smoothed_measurements, window_size, samples_per_hour, times, circadian_time):

    # remove the partial days data from the beginning and from the end
    # as the moving average removed half window_size from the beginning and from the end
    part_to_slice = int(samples_per_hour * HOURS_PER_DAY - (window_size / 2))

    spliced_moving_average_data = smoothed_measurements[part_to_slice: -part_to_slice or None]

    # plt.plot(spliced_moving_average_data)

    n = len(spliced_moving_average_data)
    num_of_days = int(n / (samples_per_hour * HOURS_PER_DAY))  # TODO if it not an int?? could it happen?

    # time2 = ["0", "12"] * num_of_days + ["0"]
    #step_value = int((n + 1) / (len(time2) - 1))
    # Set number and labels of ticks for x-axis
    # plt.xticks(np.arange(int(n+1), step=step_value), time2)

    # plt.show()
    a, p = amplitude_phase_calculation(spliced_moving_average_data, 2, num_of_days, samples_per_hour, times,
                                       circadian_time)
    # print("amplitude: ", a, "phase: ", p)

    return a, p


def amplitude_phase_calculation(data, day_num_to_use, num_of_days_of_measurements, samples_per_hour, times, ct):

    # TODO day_num_to_use should not be under 1 and no more then the number of days. if so raise a message
    if not ((day_num_to_use > 1) and (day_num_to_use < num_of_days_of_measurements)):
        print("Error")
        # TODO raise a message and try again

    measurements_per_day = int(samples_per_hour * HOURS_PER_DAY)

    measurements_to_slice = int((day_num_to_use - 1) * measurements_per_day)

    times = times[: num_of_days_of_measurements * measurements_per_day]
    circadian_times = ct[: num_of_days_of_measurements * measurements_per_day]

    one_day_measurements = data[measurements_to_slice + 1: measurements_to_slice + measurements_per_day + 1]
    one_day_times = times[measurements_to_slice: measurements_to_slice + measurements_per_day]
    one_day_circadian_times = circadian_times[measurements_to_slice: measurements_to_slice + measurements_per_day]

    # plt.plot(one_day_measurements)
    # plt.show()

    peak = max(one_day_measurements)
    trough = min(one_day_measurements)

    amplitude_value = (peak - trough) / 2

    phase_value = (one_day_times[np.argmax(one_day_measurements)], one_day_circadian_times[np.argmax(one_day_measurements)])

    return amplitude_value, phase_value


def average_amplitude_one_organism(from_ct, to_ct, data, circadian_times):
    # the data and the circadian time list have to be of the same size
    # the data is of one larva and is smoothed
    amplitudes = 0
    num_of_days = int((to_ct - from_ct) / HOURS_PER_DAY)
    index = circadian_times.index(from_ct)
    samples_per_day = int(MINUTES_PER_DAY/sampling_intervals)

    for i in range(num_of_days):
        curr_data = data[index: index + samples_per_day]  # make sure always true

        peak = max(curr_data)
        trough = min(curr_data)
        amplitudes += (peak - trough) / 2

        index += samples_per_day

    return amplitudes / num_of_days


def average_amplitude(from_ct, to_ct, group_data, circadian_times):
    average_amp = []
    n = group_data.shape[1]
    for i in range(n):
        average_amp.append(average_amplitude_one_organism(from_ct, to_ct, group_data[i], circadian_times))
    average_amp = [float("%.2f" % elem) for elem in average_amp]
    return average_amp
