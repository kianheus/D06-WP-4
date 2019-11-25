import numpy as np
import scipy as sp
from scipy import integrate

# HELLO DAAN, HERE ARE THE MOMENT AND TORSION FUNCTIONS
# USE getMomentDistr(Ldistr) aAND getTorsionDistribution(Ldistr, M_distr, rho, V, T)

# get internal moment distribution
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

# get internal torsion distribution
def getTorsionDistribution(Ldistr, M_distr, rho, V, T):
    # Define inputs
    def c(y):
        c = -0.21558573 * y + 3.6956254
        return c

    b_half = 12.009

    W_e = (2066 + (872.57 / 2)) * 9.81
    z_e = 0.7149
    x_e = 0.4661 + 0.15 * c(0.35 * b_half)
    y_e_ratio = 0.35

    # Along span
    npoints = len(Ldistr)
    y_e_point = round(npoints * y_e_ratio)

    c_distr = []
    for n in range(npoints):
        c_distr.append(c((b_half / npoints) * n))

    # Starting at wingtip
    torsion_distr = []
    for n in range(npoints):
        if abs(n - y_e_point) < 0.020833 * npoints:
            torsion_distr.append(
                (T * z_e - W_e * x_e) / (0.020833 * 2 * b_half) + Ldistr[n] * 0.15 * c_distr[n] + M_distr[n])
        else:
            torsion_distr.append(Ldistr[n] * 0.15 * c_distr[n] + M_distr[n])

    torsion_distr_new = []
    for n in range(npoints):
        if n > y_e_point:
            total_torsion = Ldistr[n] * 0.15 * c_distr[n] + M_distr[n]
            torsion_distr_new.append(total_torsion)
        else:
            total_torsion = Ldistr[n] * 0.15 * c_distr[n] + M_distr[n]
            torsion_distr_new.append(total_torsion)

    # Integrating
    dy = b_half / npoints

    torsion_int_distr = []

    def engine(n, T, z_e, W_e, x_e, y_e_point, npoints, b_half):
        if abs(n - y_e_point) < 0.020833 * npoints:
            return (T * z_e - W_e * x_e) / (0.020833 * 2 * b_half)
        else:
            return 0

    for n in range(npoints):
        sum_integral = []
        for n_1 in range(n, npoints):
            sum_integral.append(
                torsion_distr[n_1] * dy + engine(n_1, T, z_e, W_e, x_e, y_e_point, npoints, b_half) * dy)

        summation = sum(sum_integral)
        torsion_int_distr.append(summation)

    # Make graphs
    ytab = []
    for n in range(npoints):
        ytab.append(n)

    return torsion_int_distr
