import numpy as np
import matplotlib.pyplot as plt
import math
from math import cos, pi, sin


global from_period, to_period


# calculate the period of one larva by the chosen method (fourier or chi square)
def calculate_periodogram(measurements, header, method, interval, path):

    global interval_duration
    interval_duration = interval

    if method == "fourier":
        periods, periodogram_values, p_values = calculate_fourier_periodogram(measurements)
    else:
        periods, periodogram_values, p_values = calc_chi_square_periodogram(measurements)

    max_ind = peak_finder(periodogram_values, p_values)

    # incorporate calibration
    periods = inc_calib(np.array(periods))
    # convert units from minutes to hours
    periods = [x / 60 for x in periods]

    significant_period = periods[max_ind]

    larva_period_value = float("%.4f" % significant_period)
    if method == "fourier":

        plt.figure()
        plt.plot(periods, periodogram_values, "b")
        plt.axvline(x=significant_period, color="k")
        ax = plt.gca()
        (y_bottom, y_top) = ax.get_ylim()
        plt.text(significant_period, y_top+0.03, float("%.2f" % larva_period_value))
        plt.title("Periodogram (Fourier) - %s" % header, y=1.05)
        plt.xlabel("Period (h)")
        plt.ylabel("R^2")
        plt.savefig(path + "Periodogram (Fourier) - %s" % header)
        plt.close()
        # plt.clf()
    else:
        plt.figure()
        plt.plot(periods, periodogram_values, 'b', periods, p_values, "r--")
        plt.axvline(x=significant_period, color="k")
        ax = plt.gca()
        (y_bottom, y_top) = ax.get_ylim()
        plt.text(significant_period, y_top+0.03, float("%.2f" % larva_period_value))
        plt.title("Periodogram (Chi-Square) - %s" % header, y=1.05)
        plt.xlabel("Period (h)")
        plt.ylabel("Qp")
        plt.savefig(path + "Periodogram (Chi-Square) - %s" % header)
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

    return float(a_j ** 2 + b_j ** 2)


def calc_chi_square_periodogram(measurements, p_value=0.05):
    measurements_mean = 0
    n = len(measurements)
    fp = int(from_period / interval_duration)
    tp = int(to_period / interval_duration)

    periodogram_values = [None] * (tp - fp)
    period = [None] * (tp - fp)
    p_values = [None] * (tp - fp)

    for i in range(0, n):
        if not math.isnan(measurements[i]):
            measurements_mean += measurements[i]
    measurements_mean /= n

    mse = 0
    for i in range(0, n):
        if not math.isnan(measurements[i]):
            diff = measurements[i] - measurements_mean
            mse += (diff * diff)
    mse /= n

    for P in range(fp, tp):
        q_p = 0
        num_of_blocks = int(n / P)  # integer division
        for h in range(0, P):
            mean_of_block = 0
            for k in range(0, num_of_blocks):
                if not math.isnan(measurements[h + k * P]):
                    mean_of_block += measurements[h + k * P]
            mean_of_block /= num_of_blocks
            diff = mean_of_block - measurements_mean
            q_p += (diff * diff)

        q_p = q_p * num_of_blocks / mse
        period[P - fp] = P
        periodogram_values[P - fp] = float(q_p)
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


def peak_finder(values, p_values):
    # find peaks
    # use p-values as references
    relatives = list(values)

    # use p values as references
    for i in range(0, len(values)):
        relatives[i] -= p_values[i]
    return relatives.index(max(relatives))
