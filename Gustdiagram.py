from math import *

import numpy as np

import matplotlib.pyplot as plt



# General variables might be importable

R = 8.314510  # gasconstant [J/(mol*K)

M_air = 0.0289645  # Molair mass air [kg/mol]

R_air = R / M_air  # shorthand convention [J/kg*K]

g = 9.80665  # gravitational acceleration [m/s^2]

rho_0 = 1.225  # density at sealevel [kg/m^3]

gamma = 1.4  # adiabatic index of air [-]

S = 57.7  # surface area wing [m^2]

CL_clean = 1.1  # [-]

MTOW = 30502 * g  # [N]

MLW = 29142.1 * g  # ASSUMED TO BE MTOW - 0.7*FUELWEIGHT!![N]

MZFW = MTOW - 4533 * g

OEW = 21963 * g  # Operational empty weight [N]

dCL = 1.08

M_C = 0.77  # Design cruise mach number[-]

MAC = 2.63  # Mean aerodynamic chord [m]

altitude_cruise = 10668  # [m]

CL_a = 0.0762213  # CL_a at M=0 from excel calc with prandtl[-]



W = (OEW, MZFW, MTOW)



#fig, ax = plt.subplots(2, 5, figsize=(15, 6))

#ax = ax.ravel()

for i in W:

    for h in (1, 2, 3):

        if h == 1:  # Sealevel conditions

            T = 288.15  # Temperature in [K]

            rho = 1.225  # density [kg/m^3]

            colour = 'black'

            altitude = 0                # altitude in [m]

        elif h == 2:  # intermediate altitute 2000ft GENERICLY CHOSEN!!

            T = 284.1876

            rho = 1.1550945

            colour = 'red'

            altitude = 609.6                  # altitude in [m]

        elif h == 3:  # cruise altitude conditions 35000ft

            T = 218.808

            rho = 0.3795655

            colour = 'blue'

            altitude = 10668                # altitude in [m]



        #V_S0 = sqrt(2 * i / (rho * (CL_clean + dCL) * S)) * sqrt(rho / rho_0)

       #V_S0_takeoff = sqrt(2 * MTOW / (rho * (CL_clean + dCL) * S)) * sqrt(

            #rho / rho_0)  # stall speed full flaps at MTOW[m/s] "naming not very correct"

        #V_S0_landing = sqrt(2 * MLW / (rho * (CL_clean + dCL) * S)) * sqrt(rho / rho_0)  #

        V_S1 = sqrt(2 * i / (rho * (CL_clean) * S)) * sqrt(rho / rho_0)

        #V_S1_MTOW = sqrt(2 * MTOW / (rho * (CL_clean) * S)) * sqrt(

            #rho / rho_0)  # stall speed clean config  at MTOW[m/s]

        #V_S1_MLW = sqrt(2 * MLW / (rho * (CL_clean) * S)) * sqrt(rho / rho_0)  # stall speed clean config at MLW [m/s]



        #Weight = MTOW / (0.454 * g)  # weight for the imperical suckers in [lb]

        #n_max = 2.1 + 2400 / (Weight + 10000)

        #if n_max > 3.8:

            #n_max = 3.8

        #elif n_max < 2.5:

            #n_max = 2.5



        #V_A = V_S1_MTOW * sqrt(n_max)



        c_sound = sqrt(gamma * R_air * T)  # speed of sound [m/s]

        V_C = M_C * sqrt(gamma * R_air * 218.808) * sqrt(0.37956557 / rho_0)  # cruise speed constant for all plots [m/s]



        V_D1 = V_C / 0.8

        V_D2 = (M_C / 0.8) * c_sound * sqrt(rho / rho_0)

        V_D3 = (M_C + 0.05) * c_sound * sqrt(rho / rho_0)

        if max(V_D1, V_D2) == V_D1:

            V_D = V_D1  # diving velocity if not bound by compress effects [m/s]

        else:

            V_D = max(V_D2, V_D3)  # diving speed if bound by compress effefcts [m/s]



        #V_F_land = 1.6 * V_S1_MTOW

        #V_F_approach = 1.8 * V_S1_MLW

        #V_F_climb = 1.8 * V_S0_landing

        #V_F = max(V_F_land, V_F_approach, V_F_climb)  # Flapped speed [m/s] REQUIREMENTS VERRY VAGUE !!



        # New codecompared to V_n loaddiagram

        mu = (2 * i / S)/(rho * MAC * CL_a * g)

        K_G = (0.88 * mu) / (5.3 + mu)

        if altitude <= 4572:

            Uref = 17.07 - (3.66*altitude)/4572

        elif altitude >= 4572:

            Uref = 13.41 - (7.05*altitude)/18288

        
        V_B = V_S1 * sqrt(1 + (K_G * rho_0 * Uref * V_C * CL_a) / (2 * i / S))  # Design speed maximum gust intensity [m/s]

        V_TAS = np.arange(0.0001, V_B, 1) / sqrt(rho / rho_0)  # TAS that needed to be considred

        gustgrad = np.arange(9, max((107 + 1), (12.5 * MAC + 1)), 1)  # gust gradients needed to be considred [m]


        Fgz = 1 - (altitude_cruise / 76200)

        R1 = MLW / MTOW

        R2 = MZFW / MTOW

        Fgm = sqrt(R2 * tan(pi * R1 / 4))

        Fg = 0.5 * (Fgz + Fgm)

        loop = 0

        maxdns = []

        dns = []



        for henk in gustgrad:

            for joop in V_TAS:

                M = joop / c_sound

                CL_aM = CL_a / sqrt(1 - M ** 2)

                Uds = Uref * Fg * (henk / 107) ** (1 / 6)

                Ktime = 2 * i / (CL_aM * rho * joop * g * S)

                omega = pi * joop / henk

                t = np.arange(0, (2 * pi / omega), (2 * pi / omega) / 100)

                dns_p = Uds /(2*g) * (omega * np.sin(omega * t) + (np.exp(-t / Ktime)/Ktime - np.cos(omega * t)/Ktime - omega * np.sin(omega * t))/(1 + (omega * Ktime) ** -2))
                                      
                maxi = abs(max(dns_p,key=abs))

                maxdns.append(maxi)

                dns.append([dns_p,i,altitude,joop,V_D/sqrt(rho/rho_0)])









indexloop = maxdns.index(max(maxdns))

print(indexloop)

dns_p = dns[indexloop][0]

dns_n = -0.5 * dns_p

print(MTOW,MZFW,OEW)

print("max delta n occurs for altitude",dns[indexloop][2],"and the weight of",dns[indexloop][1])

print("at a speed of: ",dns[indexloop][3]," m/s")

print("V_D ",dns[indexloop][4])

plt.plot(t,dns_p,"b",label = "positive")

plt.plot(t,dns_n,"r",label = "negative")

plt.title("Delta Gustloading")

plt.xlabel("time in [s]")

plt.ylabel("delta n [-]")

plt.show()

print("Ready")
