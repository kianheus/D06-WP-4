from math import *
import numpy as np
import matplotlib.pyplot as plt

# General variables might be importable
S        =  57.7
CL_clean =  1.1
MTOW     = 30502 * 9.81
MLW      = 29142.1 * 9.81                                        # ASSUMED TO BE MTOW - 0.7*FUELWEIGHT!!
dCL      = 1.08


# at h = 0 
V_S0_takeoff   =  sqrt(2*MTOW/(1.225*(CL_clean+dCL)*S))        # stall speed full flaps at MTOW[m/s] "naming not very correct"
V_S0_landing   =  sqrt(2*MLW/(1.225*(CL_clean+dCL)*S))          #
V_S1_MTOW       =  sqrt(2*MTOW/(1.225*(CL_clean)*S))            # stall speed clean config  at MTOW[m/s]
V_S1_MLW        =  sqrt(2*MLW/(1.225*(CL_clean)*S))               # stall speed clean config at MLW [m/s]

Weight = MTOW/(0.454*9.81)
n_max = 2.1 + 2400/(Weight + 10000)
if n_max > 3.8:
    n_max = 3.8
elif n_max < 2.5:
    n_max = 2.5

V_A  = V_S1_MTOW * sqrt(n_max)

# at h = 35000 ft
M_C  = 0.77                                                     # Design cruise mach number [-]
c_sound = 296.508                                               # speed of sound [m/s]
V_C  = M_C * c_sound                                            # DEsign cruise speed [m/s]


V_D1 = V_C/0.8
V_D2 = (M_C/0.8)*c_sound
V_D3 = (M_C + 0.05)*c_sound
if max(V_D1,V_D2) == V_D1:
    V_D = V_D1                                                  # diving velocity if not bound by compress effects [m/s]
else:
    V_D = max(V_D2, V_D3)                                       # diving speed if bound by compress effefcts [m/s]
    

V_F_land        = 1.6 * V_S1_MTOW                                
V_F_approach    = 1.8 * V_S1_MLW
V_F_climb       = 1.8 * V_S0_landing
V_F = max(V_F_land, V_F_approach , V_F_climb)                   # Flapped speed [m/s] REQUIREMENTS VERRY VAGUE !!




n_maxflap = 2.0

V = np.arange(0,V_D,0.1)

n1 = (V/V_S0_takeoff)**2                                       # first parabolic part positive quadrantr
n1[n1 > 2] = 2                                                  # makes all values in array larger than 2 equal to 2
n2 = (V/V_S1_MTOW)**2                                           # second parabolic part to VA
n2[n2 > n_max] = n_max                                          # makes all values bigger than n_max equal to n_max

# bottom parabolic part assuming Clmin = -CLmax
# CS-25 says nmin = -1

n3 = -(V/V_S1_MTOW)**2
n3[n3 < -1] = -1
n_min = np.arange(-1,0,(1/(V_D - V_C)))                              # gives linear relation of n_min between V_D and V_C for integer speeds

print(V_S0_takeoff)
print(V_S1_MTOW)
print(V_S1_MLW)
print(V_S0_landing)
print(V_F)
print(V_A)

plt.plot(V,n1,'black')
plt.plot(V,n2,'black')
plt.plot(V[V<V_C],n3[V<V_C],'black')                                    # making the graph not continue after V_C
plt.plot(np.arange(V_C,V_D,1),n_min,'black')
plt.plot((V_D, V_D), (0, n_max), 'black')
plt.axhline(y=0,color='black')
plt.plot((V_A,V_A),(-1,n_max),color='black',linestyle='dashdot')
plt.plot((V_C,V_C),(-1,n_max),color='black',linestyle='--')
i = 0
while i < len(n1) and n1[i] < 2:
    V_X = V[i]
    i= i+1
plt.plot((V_X,V_X),(-1,2),color='black',linestyle=':')
plt.xlabel("V_EAS [m/s]")


plt.ylabel("Load factor [-]")
plt.show()







