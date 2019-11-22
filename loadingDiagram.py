# python script to create loading diagram
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy import integrate

# import big ass list
# Ldistr = [......]

def genLoadDiagram(Ldistr):
    Ldistr[0] = 0
    Ldistr.append(0)
    LdistrArr = -1 * np.array(Ldistr)
    Ldistr = list(LdistrArr)

    # input variables
    L = 12.009  # half wing span
    Lmax = 2500
    npoints = len(Ldistr)-1
    dy = L / npoints

    # find p position
    Ppos = round(0.35 * npoints)


    def p(n, Ppos):
        if n == Ppos:
            return 20267
        else:
            return 0


    # calculate wing weight along y
    def W_w(y):
        return 391.2366 * (-0.215585 * y + 3.695654)


    W_w_distr = []
    for n in range(npoints + 1):
        W_w_distr.append(W_w((L / npoints) * n))
    W_w_distr[0] = 0
    W_w_distr[npoints] = 0

    halfWingWeight, error = sp.integrate.quad(W_w, 0, L)
    print(halfWingWeight)
    # calculate shear
    sumdistr = []
    for n in range(npoints + 1):
        sumdistr.append((W_w_distr[n] + Ldistr[n]))

    shear_internal_distr = []
    for n in range(npoints + 1):
        sum_integral = []
        for n_1 in range(n, npoints + 1):
            sum_integral.append((sumdistr[n_1] * dy + p(n_1, Ppos)))
        sum_contribution = sum(sum_integral)
        shear_internal_distr.append(sum_contribution)

    moment_internal_distr = []
    for n in range(npoints + 1):
        sum_integral = []
        for n_1 in range(n, npoints + 1):
            sum_integral.append((shear_internal_distr[n_1] * dy))
        sum_contribution = sum(sum_integral)
        moment_internal_distr.append(sum_contribution)

    # make graphs
    ytab = []
    for n in range(npoints + 1):
        ytab.append(n)

    plt.subplot(211)
    plt.hlines(0, 0, npoints)
    plt.plot(ytab, Ldistr, label="Lift")
    plt.plot(ytab, W_w_distr, label="Weight")
    plt.legend(loc="upper left")
    plt.subplot(212)
    plt.hlines(0, 0, npoints)
    plt.plot(ytab, shear_internal_distr, label="Shear")
    plt.plot(ytab, moment_internal_distr, label="Moment")
    plt.legend(loc="upper left")
    plt.show()

def getMomentDistr(Ldistr):
    Ldistr[0] = 0
    Ldistr.append(0)
    LdistrArr = -1 * np.array(Ldistr)
    Ldistr = list(LdistrArr)

    # input variables
    L = 12.009  # half wing span
    Lmax = 2500
    npoints = len(Ldistr) - 1
    dy = L / npoints

    # find p position
    Ppos = round(0.35 * npoints)

    def p(n, Ppos):
        if n == Ppos:
            return 20267
        else:
            return 0


    # calculate wing weight along y
    def W_w(y):
        return 391.2366 * (-0.215585 * y + 3.695654)

    W_w_distr = []
    for n in range(npoints + 1):
        W_w_distr.append(W_w((L / npoints) * n))
    W_w_distr[0] = 0
    W_w_distr[npoints] = 0

    halfWingWeight, error = sp.integrate.quad(W_w, 0, L)
    print(halfWingWeight)
    # calculate shear
    sumdistr = []
    for n in range(npoints + 1):
        sumdistr.append((W_w_distr[n] + Ldistr[n]))

    shear_internal_distr = []
    for n in range(npoints + 1):
        sum_integral = []
        for n_1 in range(n, npoints + 1):
            sum_integral.append((sumdistr[n_1] * dy + p(n_1, Ppos)))
        sum_contribution = sum(sum_integral)
        shear_internal_distr.append(sum_contribution)

    moment_internal_distr = []
    for n in range(npoints + 1):
        sum_integral = []
        for n_1 in range(n, npoints + 1):
            sum_integral.append((shear_internal_distr[n_1] * dy))
        sum_contribution = sum(sum_integral)
        moment_internal_distr.append(sum_contribution)

    # make graphs
    ytab = []
    for n in range(npoints + 1):
        ytab.append(n)

    return moment_internal_distr