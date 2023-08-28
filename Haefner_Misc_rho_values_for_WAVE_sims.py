"""
Matthew Haefner
SEA Lab
rho_sw value generator WAVE post-compression
3/7/22
"""

# Generates "post-compression" rho_sw values in column G of Haefner_WAVE_Simulations.xlsx
# Ultimately (in Haefner_WAVE_Simulations.xlsx) used to calculate "post-compression" S_f input values in column E of Haefner_WAVE_Simulations.xlsx

import math
import numpy as np

# General Parameters dictionary
p = {"g":9.81,             # gravitational acceleration [m/s^2]
     "rho_sw":1023.6,      # generic density of saltwater at 1 atm and 25 degrees celsius [kg/m^3]
     "rho_fw":996.9,       # generic density of freshwater at 1 atm and 25 degrees celsius [kg/m^3]
     "eta_hp":0.894,       # pump-side efficiency
     "eta_ht":0.894,       # pump-side efficiency [-]
     "T":25.0,             # temperature [degrees celsius]
     "S_sw":35.0,          # generic seawater salinity [g salt/kg seawater]
     "M_salt":58.44,       # molar mass of NaCl [g/mol]
     "P_p":14.696,         # permeate pressure, atmospheric for now [psi]
     "R":8.3145,           # universal gas constant [J/(mol*K)]
     "MW_salt_avg":30.548} # average molar mass of NaCl [g/mol]

# Converting salinity to density
def density_from_salinity(p,S,P):
    # Unpacking relevant parameters
    t = p["T"]             # celsius
    rho_fw = p["rho_fw"]   # kg/m^3
    P_o = p["P_p"]/145.038 # MPa (from psi)

    # Density at Atmospheric
    b1 = 8.020e2
    b2 = -2.001
    b3 = 1.677e-2
    b4 = -3.060e-5
    b5 = -1.613e-5

    S_kg_kg = S/1000 # kg/kg

    rho_sw_atm = rho_fw+((b1*S_kg_kg)+(b2*S_kg_kg*t)+(b3*S_kg_kg*(t**2))+(b4*S_kg_kg*(t**3))+(b5*(S_kg_kg**2)*(t**2))) # kg/m^3

    # Other Correction Term
    c1 = 5.0792e-4
    c2 = -3.4168e-6
    c3 = 5.6931e-8
    c4 = -3.7263e-10
    c5 = 1.4465e-12
    c6 = -1.7058e-15
    c7 = -1.3389e-6
    c8 = 4.8603e-9
    c9 = -6.8039e-13

    d1 = -1.1077e-6
    d2 = 5.5584e-9
    d3 = -4.2539e-11
    d4 = 8.3702e-9

    P = P/145.038 # MPa (from psi)

    Fp = math.exp(((P-P_o)*(c1+(c2*t)+(c3*(t**2))+(c4*(t**3))+(c5*(t**4))+(c6*(t**5))+(S*(d1+(d2*t)+(d3*(t**2))))))+((0.5*((P**2)-(P_o**2)))*(c7+(c8*t)+(c9*(t**3))+(d4*S))))

    rho_brine = rho_sw_atm*Fp # kg/m^3

    return rho_brine


# Setting input variable ranges
Q_f_min = 3.41 # m^3/hr
Q_f_max = 15.5 # m^3/hr

S_f_min = 35 # g/kg
S_f_max = 55 # g/kg

P_f_min = 300 # psi
P_f_max = 1200 # psi

# Input variable arrays
P_f_array = np.linspace(P_f_min,P_f_max,int(((P_f_max-P_f_min)/25)+1))
S_f_array = np.linspace(S_f_min,S_f_max,int(((S_f_max-S_f_min)/2)+1))
Q_f_array = np.linspace(Q_f_min,Q_f_max,int((Q_f_max-Q_f_min)+1))

full_factorial_array = np.zeros(3)

for i in P_f_array:
    for j in S_f_array:
        for k in Q_f_array:
            left_hand_side = (-18.222*i)+(282.748*j)+(200*k)+1334.934
            if left_hand_side >= 0:
                good_point = np.array([i,j,round(k,2)])
                full_factorial_array = np.vstack([full_factorial_array,good_point])
full_factorial_array = np.delete(full_factorial_array,0,0)

rho_sw_array = np.zeros(1)

for i in range(3363):
    rho_sw_array = np.vstack([rho_sw_array,round(density_from_salinity(p, full_factorial_array[i,1], full_factorial_array[i,0]),3)])
rho_sw_array = np.delete(rho_sw_array,0,0)
        