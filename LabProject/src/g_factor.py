import numpy as np
from scipy import stats


HOURS_PER_DAY = CIRCADIAN_TIME = 24


def g_factor_calculation(larva_measurements_vec, number_of_days, samples_per_hour):
    """
    Calculate g factor for one larva vector of measurements
    :param larva_measurements_vec: vector of the relevant data
    :param number_of_days:
    :param samples_per_hour:
    :return:
    """

    tic_time = 1 / samples_per_hour
    number_of_time_points = HOURS_PER_DAY * number_of_days * samples_per_hour
    real_is = 1 / CIRCADIAN_TIME
    freq_step = 1 / (HOURS_PER_DAY * number_of_days)
    f = [j for j in frange(freq_step, 1/(2*tic_time), freq_step)]

    location_circadian = f.index(real_is)
    larva_measurements_vec = [value for value in larva_measurements_vec if not np.isnan(value)]  # remove nan values
    y = np.fft.fft(larva_measurements_vec, number_of_time_points)
    py = y * y.conj() / number_of_time_points
    py = py.real  # [float(i) for i in py]

    # test_factor = py[location_circadian+1]
    return float("%.2f" % float(py[location_circadian+1] / np.sum(np.array(py[1: (len(f)+1)]))))  # py[location_circadian+1] / np.sum(np.array(py[1: (len(f)+1)]))


# same as python range but for floats too
def frange(start, stop, span):
    """
    same as python range but for floats
    :param start: Starting number of the sequence
    :param stop: Generate numbers up to, but not including this number
    :param span: Difference between each number in the sequence
    :return:
    """
    while start < stop:
        yield start
        start += span


# TODO how to make the calculation if I have more then two types (i.e more then 2 vectors)?
def g_factor_significant(g_factor_group_1, g_factor_group_2):
    """
    receives two lists of g factor calculated from two groups and runs Kolmogorov-Smirnov statistic test on every
     2 samples
    :param g_factor_group_1:
    :param g_factor_group_2:
    :return:
    """
    ks, p_value = stats.ks_2samp(g_factor_group_1, g_factor_group_2)
    return ks, p_value
