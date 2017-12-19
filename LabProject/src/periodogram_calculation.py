import numpy as np
import matplotlib.pyplot as plt
import math
from math import cos, pi, sin
from pprint import pprint


global from_period, to_period
# from_period = 960  # TODO if I to init the global variables, will they be updated from the user's input?
# to_period = 1920


# calculate the period of one larva by the chosen method (fourier or chi square)
def calculate_periodogram(measurements, i, method, interval):

    global interval_duration
    interval_duration = interval  # TODO make sure it is initizlized

    if method == "fourier":
        periods, periodogram_values, p_values = calculate_fourier_periodogram(measurements)
    else:
        periods, periodogram_values, p_values = calc_chi_square_periodogram(measurements)

    # incorporate calibration
    periods = inc_calib(np.array(periods))
    # convert units from minutes to hours
    periods = [x / 60 for x in periods]

    max_index = np.argmax(periodogram_values)
    significant_period = periods[max_index]
    mark_max_value = [significant_period] * len(periods)

    if method == "fourier":
        larva_period_value = float("%.2f" % significant_period)

        plt.figure()
        plt.plot(periods, periodogram_values, 'b', mark_max_value, periodogram_values, 'k:')
        plt.title("Periodogram (Fourier) - #%d" % (i+1))
        plt.xlabel("Period (h)")
        plt.ylabel("R^2")
        plt.savefig("Periodogram (Fourier) - #%d" % (i+1))
        plt.close()
        # plt.clf()
    else:
        # check if the max value (the period) is above the p value i.e. is significant
        if p_values[max_index] < np.amax(periodogram_values):
            larva_period_value = float("%.2f" % significant_period)
        else:
            larva_period_value = 0  # TODO: which number to insert when it is not significant 0 or nan??

        plt.figure()
        plt.plot(periods, periodogram_values, 'b', periods, p_values, 'r--', mark_max_value, periodogram_values, 'k:')
        plt.title("Periodogram (Chi-Square) - #%d" % (i+1))
        plt.xlabel("Period (h)")
        plt.ylabel("Qp")
        plt.savefig("Periodogram (Chi-Square) - #%d" % (i+1))
        plt.close()
        # plt.clf()

    return larva_period_value


def calculate_fourier_periodogram(measurements, p_value=0.05):
    n = len(measurements)
    sum_r2 = 0.0
    fp = int(from_period / interval_duration)
    tp = int(to_period / interval_duration)

    periods = list(range(fp, tp))
    periodogram_values = []

    for i in range(tp - fp):
        j = float(n) / (fp + i)
        x = r2(measurements, j)
        periodogram_values.append(x)
        sum_r2 += x

    threshold = sum_r2 * (1 - ((p_value / n) ** (1.0 / (n - 1))))
    p_values = [threshold] * (tp - fp)

    return periods, periodogram_values, p_values


def r2(measurements, j):
    a_j, b_j = 0.0, 0.0
    n = len(measurements)

    for i, measurement in enumerate(measurements):
        arg = 2 * pi * j * i / n
        if not math.isnan(measurement):
            a_j += measurement * cos(arg)
            b_j += measurement * sin(arg)

    a_j *= (2 / n)
    b_j *= (2 / n)

    return a_j ** 2 + b_j ** 2


def calc_chi_square_periodogram(measurements, p_value=0.05):
    M = 0
    n = len(measurements)
    fp = int(from_period / interval_duration)
    tp = int(to_period / interval_duration)

    periodogram_values = [None] * (tp - fp)
    period = [None] * (tp - fp)
    p_values = [None] * (tp - fp)

    for i in range(0, n):
        if not math.isnan(measurements[i]):
            M += measurements[i]
    M /= n

    mse = 0
    for i in range(0, n):
        if not math.isnan(measurements[i]):
            diff = measurements[i] - M
            mse += diff * diff

    mse /= n

    for P in range(fp, tp):
        Qp = 0
        K = n / P  # integer division
        for h in range(0, P):
            Mh = 0
            for k in range(0, int(K)):
                if not math.isnan(measurements[k]):
                    Mh += measurements[h + k * P]
            Mh /= K
            diff = Mh - M
            Qp += diff * diff

        Qp = Qp * K / mse
        period[P - fp] = P
        periodogram_values[P - fp] = float(Qp)
        pv = p_value / 10.0
        p_values[P - fp] = float(chi_square_cdf_inv(1 - pv, P-1))

    return period, periodogram_values, p_values


def chi_square_cdf_inv(x, k):

    arg = norm_cdf_inv(x)

    h = 1 - (2.0 * k * k) / (3.0 * k * k)
    p = 1.0 / k
    m = (h - 1) * (1 - 3 * h)

    s = arg * h * math.sqrt(2 * p) * (1 + 0.5 * m * p) + (1 + h * p * (h - 1 - 0.5 * (2 - h) * m * p))
    return k * math.pow(s, 1.0 / h)


def norm_cdf_inv(p):
    return math.sqrt(2) * erf_inv(2 * p - 1)


a = 8 * (math.pi - 3) / 3 / math.pi / (4 - math.pi)


def erf_inv(x):
    lg = math.log(1 - x * x)
    s = 2 / (math.pi * a) + lg / 2
    return x / np.abs(x) * math.sqrt(math.sqrt(s * s - lg / a) - s)


# incorporate calibration
def inc_calib(periods):
    factor = interval_duration
    for i in range(periods.size):
        periods[i] *= factor
    return periods
