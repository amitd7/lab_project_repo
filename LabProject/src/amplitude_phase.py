import matplotlib.pyplot as plt
import numpy as np


HOURS_PER_DAY = 24


def amplitude_phase(smoothed_measurements, window_size, samples_per_hour, times):

    # remove the partial days data from the beginning and from the end
    # as the moving average removed half window_size from the beginning and from the end
    part_to_slice = int(samples_per_hour * HOURS_PER_DAY - (window_size / 2))

    spliced_moving_average_data = smoothed_measurements[part_to_slice: -part_to_slice or None]

    plt.plot(spliced_moving_average_data)

    n = len(spliced_moving_average_data)
    num_of_days = int(n / (samples_per_hour * HOURS_PER_DAY))  # TODO if it not an int?? could it happen?

    time2 = ["0", "12"] * num_of_days + ["0"]
    step_value = int((n + 1) / (len(time2) - 1))
    # Set number and labels of ticks for x-axis
    plt.xticks(np.arange(int(n+1), step=step_value), time2)

    plt.show()
    a, p = amplitude_phase_calculation(spliced_moving_average_data, 2, num_of_days, samples_per_hour, times)
    print("amplitude: ", a, "phase: ", p)

    return a, p


def amplitude_phase_calculation(data, day_num_to_use, num_of_days_of_measurements, samples_per_hour, times):

    # TODO day_num_to_use should not be under 1 and no more then the number of days. if so raise a message
    if not ((day_num_to_use > 1) and (day_num_to_use < num_of_days_of_measurements)):
        print("Error")
        # TODO raise a message and try again

    measurements_per_day = int(samples_per_hour * HOURS_PER_DAY)

    measurements_to_slice = int((day_num_to_use - 1) * measurements_per_day)

    times = times[: num_of_days_of_measurements * measurements_per_day]

    one_day_measurements = data[measurements_to_slice + 1: measurements_to_slice + measurements_per_day + 1]
    one_day_times = times[measurements_to_slice: measurements_to_slice + measurements_per_day]

    plt.plot(one_day_measurements)
    plt.show()

    peak = max(one_day_measurements)
    trough = min(one_day_measurements)

    amplitude_value = (peak - trough) / 2

    phase_value = one_day_times[np.argmax(one_day_measurements)]  # TODO translate to circadian time

    return amplitude_value, phase_value
