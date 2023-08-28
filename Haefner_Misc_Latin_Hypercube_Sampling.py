"""
Matthew Haefner
SEA Lab
Latin Hypercube Sampling for Sampling WAVE RO Design Space
3/3/23
"""

import lhsmdu
import numpy as np

# Latin Hypercube Sampling distribution
lhs_dist = lhsmdu.sample(3,50)

# Factor mins and maxes
Q_f_min = 3.41 # m^3/hr
Q_f_max = 15.5 # m^3/hr

S_f_min = 35 # g/kg
S_f_max = 55 # g/kg

P_f_min = 300 # psi
P_f_max = 1200 # psi

factor_lims = np.array([[Q_f_min,Q_f_max],[S_f_min,S_f_max],[P_f_min,P_f_max]])

lhs = np.zeros((3,50))

for i in range(lhs.shape[0]):
    for j in range(lhs.shape[1]):
        lhs[i,j] = (lhs_dist[i,j]*(factor_lims[i,1]-factor_lims[i,0]))+factor_lims[i,0]
        
