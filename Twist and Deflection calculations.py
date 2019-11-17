''' This script calculates the deflection and twist of a wing box cross section'''
# Keeping time and imports ------------------------------------------------------------------
import time
start_time = time.time()

import math
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import interpolate
from scipy import integrate

# Inputs ------------------------------------------------------------------------------------

# Given by the previous iterations
C_r = 3.695625413           # Root chord            [m]
C_t = 1.107118055           # Tip chord             [m]
b = 24.01371734             # Wing span             [m]
x_frontspar = 0.20          # Front spar position   [%]
x_rearspar = 0.60           # Rear spar position    [%]

# To be exactly determined
t_sheet = 0.001             # Sheet thickness       [m]
step = 0.1                  # Calculation accuracy  [-]
A_str = 0.001               # Area of a stringer    [m^2]

# Stringers
''' To insert stringers on your desired location, you first set the intervals at which
the stringers will vary, using distance_top. Then in the stringers_top array you input
the stringer location as % of the chord. With a spacing of [0.25, 0.5, 1] for example,
and locations [[0.2, 0.2, 0.2], [0.4, 0.4, 0], [0.6, 0, 0], [0.8, 0.8, 0.8]], there
will be 4 stringers located at [0.2, 0.4, 0.6, 0.8] percent of the chord up until 25%
of the span. Then until 50% it will look like [0.2, 0.4, 0.8] etc.'''

distance_top = np.array([0.25, 0.5, 1])             # Interval of stringer variance [%span]
stringers_top = np.array([[0.2, 0.2, 0.2], [0.4, 0.4, 0], [0.6, 0, 0], [0.8, 0.8, 0.8]])
topstr = sp.interpolate.interp1d(distance_top,stringers_top,kind='next',fill_value='extrapolate')

distance_bot = np.array([0.25, 0.5, 1])             # Interval of stringer variance [%span]
stringers_bot = np.array([[0.2, 0.2, 0.2], [0.4, 0.4, 0], [0.6, 0, 0], [0.8, 0.8, 0.8]])
botstr = sp.interpolate.interp1d(distance_bot,stringers_bot,kind='next',fill_value='extrapolate')

# Calculated
h_frontspar = 0.0908        # Front spar height     [%]     FILLER
h_rearspar = 0.0804         # Rear spar height      [%]     FILLER

M = np.array(0)             # Array of the moments  [Nm]    FILLER
E = 69 * 10**9              # E-modulus of material [Pa]    FILLER
G = 27 * 10**9              # G-modulus of material [Pa]    FILLER


# Calculations of length, area & centroid ---------------------------------------------------

# Creating the steps and interval. This will decide the level
#   of detail of the calculations
y = np.arange(0, b/2, step)
n_points = len(y)
C_y = C_r + C_t * (y/(b/2)) - C_r * (y/(b/2))   # Chord as function of y    [m]
ylst = y.tolist()

# Calculating the neutral axis, dividing the wing box plates into 4 parts
# Part length       [m]
L_I = (x_rearspar - x_frontspar) * C_y 
L_II = h_rearspar * C_y
L_IV = h_frontspar * C_y
L_III = np.sqrt( (L_IV-L_II)**2 + L_I**2)

# Part area         [m^2]
A_I = L_I * t_sheet
A_II = L_II * t_sheet
A_III = L_III * t_sheet
A_IV = L_IV * t_sheet

# Part centroid     [m]
z_I = L_IV
z_II = L_IV - 0.5*L_II
z_III = 0.5 * (L_IV - L_II)
z_IV = 0.5 * L_IV

# Calculating the moment of inertia for the stringers----------------------------------------

n_iterations = len(y)

# Top stringers
y_n = 0
n_strlist_top = []

# This loop computes the # of stringers at each y coord
for i in range(n_iterations):
    n_str_y_n = len(topstr(y_n/(b/2))[topstr(y_n/(b/2)) != 0])
    y_n += step
    n_strlist_top.append(n_str_y_n)

n_str_top = np.array(n_strlist_top)

# Bottom stringers
y_n = 0
n_strlist_bot = []

# This loop computes the # of stringers at each y coord
for i in range(n_iterations):
    n_str_y_n = len(botstr(y_n/(b/2))[botstr(y_n/(b/2)) != 0])
    y_n += step
    n_strlist_bot.append(n_str_y_n)

n_str_bot = np.array(n_strlist_bot)

# z-coordinate of neutral axis, as seen from the bottom of the front spar
z_n = (A_str*(n_str_top)*z_I + A_str*(n_str_bot)*z_III + A_I*z_I + A_II*z_II + A_III*z_III + A_IV*z_IV) / (A_I + A_II + A_III + A_IV + A_str*(n_str_top+n_str_bot))

# Stringer MOI
I_str_bot = n_str_bot * A_str * (z_I-z_n)**2
I_str_top = n_str_top * A_str * (z_I-z_n)**2

# Calculating the moment of inertia for the sheets-------------------------------------------


# Parts 1, 2 and 4
I_I = L_I * t_sheet**3 * (1/12) + A_I * (z_I-z_n)**2
I_II = t_sheet * L_II**3 * (1/12) + A_II * (z_II-z_n)**2
I_IV = t_sheet * L_IV**3 * (1/12) + A_IV * (z_IV-z_n)**2

# Part 3 is slanted
I_III_x = L_III * t_sheet**3 * (1/12)
I_III_y = t_sheet * L_III**3 * (1/12)
theta_III = np.tan((L_IV-L_II)/L_I)

I_III = I_III_x * np.cos(theta_III)**2 + I_III_y * np.sin(theta_III)**2 + A_III * (z_III-z_n)**2

# Final MOI         [m^4]
I = I_I + I_II + I_III + I_IV + I_str_top + I_str_bot
I_lst = np.array(I).tolist()
AMOI = sp.interpolate.interp1d(y,I_lst,kind='cubic',fill_value='extrapolate')
I_plot = np.array(AMOI(y))

def I_y(y):
    return AMOI(y)

# OUTPUT: Moment of inertia as a function of y: I_y(y) 
print('MOI done')

# Calculating the torsional constant --------------------------------------------------------


A_encl = L_I*L_II + (L_IV-L_II)*L_I*0.5         # Enclosed area [m^2]
ds_t = (L_I+L_II+L_III+L_IV) / t_sheet          # Integral of the reciprocal of t and length

J = 4 * A_encl**2 / ds_t                        # Torsional constant
T_C = sp.interpolate.interp1d(y,J,kind='cubic',fill_value='extrapolate')
J_plot = np.array(T_C(y))

def J_y(y):
    return T_C(y)

# OUTPUT: Torsional constant as a function of y: J_y(y)
print('J done')

# This section will be deleted once the moment and torque function is available -------------
''' The moment and torque functions can be inputted in this program as a list of values
Should it have another format, that doesn't have to be a problem but a little bit of
additional code will be required. '''

# Making a fake lift distribution
def L_prime(y):
    return 43.4

def L(y):
    return L_prime(0)*b/2


Mlst = []

for i in range(n_iterations):
    M = sp.integrate.quad(L, ylst[i], (b/2))
    Mlst.append(M[0])

# Fake moment distribution:
M = np.array(Mlst)
print('Fake calculations done')

# Fake torque distibution
T = 200
Tlst = []

for i in range(n_iterations):
    Tlst.append(T)

T = np.array(Tlst)

# Calculating deflection --------------------------------------------------------------------


# Retrieving -M/EI
MEI = -M/(E*I)

# First integration
MEIfunc = sp.interpolate.interp1d(y,MEI,kind='cubic',fill_value='extrapolate')
dvdylst = []

for i in range(n_iterations):
    dvdy = sp.integrate.quad(MEIfunc, ylst[i], (b/2))
    dvdylst.append(dvdy[0])

dvdy = np.array(dvdylst)

# Adding integration constant (because dvdy(0) = 0)
C_1 = -dvdy[0]
dvdy += C_1

print('First integration done')
# Second integration
dvdyfunc = sp.interpolate.interp1d(y,dvdy,kind='cubic',fill_value='extrapolate')
vlst = []

for i in range(n_iterations):
    v = sp.integrate.quad(dvdyfunc, ylst[i], (b/2))
    vlst.append(v[0])

v = np.array(vlst)

# Adding integration constant (because v(0) = 0)
C_2 = -v[0]
v += C_2
print('Second integration done')

# Calculating twist -------------------------------------------------------------------------


# Retrieving T/GJ
TGJ = T/(G*J)

# Integration
TGJfunc = sp.interpolate.interp1d(y,TGJ,kind='cubic',fill_value='extrapolate')
philst = []

for i in range(n_iterations):
    phi = sp.integrate.quad(TGJfunc, ylst[i], (b/2))
    philst.append(phi[0])

phi = np.array(philst)

# Adding integration constant (phi(0) = 0)

C_3 = -phi[0]
phi += C_3
print('Twist integration done')

# Statistics
v_max = -np.amin(v)
v_percentage = v_max/b * 100
phi_max = abs(-np.amin(phi) * 180 / math.pi)

print('The maximum deflection is', round(v_max, 4), '[m] or', round(v_percentage, 3), '[%] of the span.')
print('The maximum twist is', round(phi_max, 4),'[deg]')
print("Calculations took %s seconds." % round((time.time() - start_time), 2))
plt.plot(y, v)
plt.plot(y, dvdy)
plt.plot(y, I)
plt.plot(y, phi)
plt.show()
