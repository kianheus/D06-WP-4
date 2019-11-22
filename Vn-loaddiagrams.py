from math import *
import numpy as np
import matplotlib.pyplot as plt

# General variables might be importable
R        = 8.314510                                             # gasconstant [J/(mol*K)
M_air    = 0.0289645                                            # Molair mass air [kg/mol]
R_air    = R/M_air                                              # shorthand convention [J/kg*K]
g        = 9.80665                                              # gravitational acceleration [m/s^2]
rho_0    = 1.225                                                # density at sealevel [kg/m^3]
gamma    = 1.4                                                  # adiabatic index of air [-]
S        = 57.7                                                 # surface area wing [m^2]
CL_clean = 1.1                                                  #[-]
MTOW     = 30502 * g                                            #[N]
MLW      = 29142.1 * g                                          # ASSUMED TO BE MTOW - 0.7*FUELWEIGHT!![N]
OEW      = 21963 * g                                            # Operational empty weight [N]
dCL      = 1.08
M_C      = 0.77                                                 # Design cruise mach number [-]

W = np.linspace(OEW,MTOW,num = 10)
count = 0
fig, ax = plt.subplots(2,5, figsize=(15, 6))
ax = ax.ravel()
for i in W:
    for h in (1,2,3):
        if h == 1:                                                  # Sealevel conditions
            T      = 288.15                                         # Temperature in [K]
            rho    = 1.225                                          # density [kg/m^3]
            colour = 'black'
        elif h == 2:                                                # intermediate altitute 2000ft GENERICLY CHOSEN!!
            T      = 284.1876
            rho    = 1.1550945
            colour = 'red'
        elif h == 3:                                                # cruise altitude conditions 35000ft
            T      = 218.808
            rho    = 0.3795655
            colour = 'blue'
            
        V_S0            =  sqrt(2*i/(rho*(CL_clean+dCL)*S))*sqrt(rho/rho_0) 
        V_S0_takeoff    =  sqrt(2*MTOW/(rho*(CL_clean+dCL)*S))*sqrt(rho/rho_0)        # stall speed full flaps at MTOW[m/s] "naming not very correct"
        V_S0_landing    =  sqrt(2*MLW/(rho*(CL_clean+dCL)*S))*sqrt(rho/rho_0)         #
        V_S1            =  sqrt(2*i/(rho*(CL_clean)*S))*sqrt(rho/rho_0)
        V_S1_MTOW       =  sqrt(2*MTOW/(rho*(CL_clean)*S))*sqrt(rho/rho_0)            # stall speed clean config  at MTOW[m/s]
        V_S1_MLW        =  sqrt(2*MLW/(rho*(CL_clean)*S))*sqrt(rho/rho_0)             # stall speed clean config at MLW [m/s]

        Weight = MTOW/(0.454*g)                                                          # weight for the imperical suckers in [lb]
        n_max = 2.1 + 2400/(Weight + 10000)
        if n_max > 3.8:
            n_max = 3.8
        elif n_max < 2.5:
            n_max = 2.5

        V_A  = V_S1_MTOW * sqrt(n_max)
                                                            
        c_sound = sqrt(gamma * R_air * T)                                   # speed of sound [m/s]
        V_C  = M_C * sqrt(gamma * R_air * 218.808)* sqrt(0.37956557/rho_0)  # cruise speed constant for all plots [m/s]


        V_D1 = V_C/0.8
        V_D2 = (M_C/0.8)*c_sound * sqrt(rho/rho_0)
        V_D3 = (M_C + 0.05)*c_sound *sqrt(rho/rho_0)
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

        n1 = (V/V_S0)**2                                                # first parabolic part positive quadrantr
        n1[n1 > 2] = 2                                                  # makes all values in array larger than 2 equal to 2
        n2 = (V/V_S1)**2                                                # second parabolic part to VA
        n2[n2 > n_max] = n_max                                          # makes all values bigger than n_max equal to n_max

        # bottom parabolic part assuming Clmin = -CLmax
        # CS-25 says nmin = -1

        n3 = -(V/V_S1)**2
        n3[n3 < -1] = -1
        n_min = np.arange(-1,0,(1/(V_D - V_C)))                              # gives linear relation of n_min between V_D and V_C for integer speeds

        ax[count].plot(V,n1,colour)
        ax[count].plot(V,n2,colour)
        ax[count].plot(V[V<V_C],n3[V<V_C],colour)                                    # making the graph not continue after V_C
        ax[count].plot(np.arange(V_C,V_D,1),n_min, colour)
        ax[count].plot((V_D, V_D), (0, n_max),colour)
        ax[count].plot((V_A,V_A),(-1,n_max),color=colour,linestyle='dashdot')
        ax[count] .plot((V_C,V_C),(-1,n_max),color=colour,linestyle='--')

        i = 0
        while i < len(n1) and n1[i] < 2:
            V_X = V[i]
            i= i+1
        ax[count].plot((V_X,V_X),(-1,2),color= colour,linestyle=':')

    
    
    ax[count].axhline(y=0,color = 'black')
    count = count +1
    

fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10) = plt.subplots(10, sharex=True)    
plt.show()







