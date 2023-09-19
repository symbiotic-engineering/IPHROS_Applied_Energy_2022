IPHROS Applied Energy

The files in this repository are supplemental material for the paper submitted to Applied Energy entitled "Integrated Pumped Hydro Reverse Osmosis System Optimization Featuring Surrogate Model Development in Reverse Osmosis Modeling". They are used to generate the figures and results tables included in the paper, and any other values reported in the paper.


Haefner_Figure_6_Data.xlsx - Excel file that contains WAVE simulation results for the input configurations selected using the file Haefner_Misc_Latin_Hypercube_Sampling.py
Haefner_Figure_6_Generator.m - Matlab file that creates Figure 6 using the plane defined from Haefner_Figure_6_Hyperplane.py
Haefner_Figure_6_Hyperplane.py - Python file that calculates the normal vector of the hyperplane generated using an SVM classifier that separates WAVE simulations (from Haefner_Figure_6_Data.xlsx) that trigger RO design warnings from those that do not

Haefner_Figure_7_Table_4_Generator.m - MATLAB file that creates Figure 7 and generates the values in Table 4

Haefner_Figure_9_Data.xlsx - Excel file used to calculate the error between the predicted values of the permeate flowrate for an RO element and the actual values for a random 20% test set of Haefner_WAVE_Simulations_trimmed.xlsx
Haefner_Figure_9_Generator.m - Matlab file that creates Figure 9

Haefner_Figure_12_13_14_15_Table_5_Data.xlsx - Excel file containing the optimization results obtained from running Haefner_IPHROS_MOGA.m
Haefner_Figure_12_13_14_15_Table_5_Generator.m - Matlab file that generates the values and plots for Figures 12, 13, 14, and 15. It also generates the values for Table 5.

Haefner_Figure_16_Data.xlsx - Excel file used for the economic analysis found in Section 5.2
Haefner_Figure_16_Generator.m - Matlab file that creates Figure 16

Haefner_Figure_B17_Generator.jmpprj - JMP project that generates Figure B17, and contains the analysis used to determine that the fractional salt rejection rate is a function of the permeate flowrate

Haefner_IPHROS_MOGA - MATLAB file that runs the optimization of IPHROS using a multiobjective genetic algorithm
Haefner_NN_for_Qp.py - Python file that determines the architechure of the neural network illustrated in Figure 8, and generates the values for ro_weights.m. Both are done using AutoKeras
Haefner_WAVE_Simulations.xlsx - Excel file containing all WAVE simulation results
Haefner_WAVE_Simulations_Trimmed.xlsx - Excel file containing the results of WAVE simulation that do not trigger any design warnings

Haefner_Misc_Latin_Hypercube_Sampling.py - Python file that determines the Latin hypercube sampling of the WAVE simulations
Haefner_Misc_rho_values_for_WAVE_sims.py - Python file that generates the "post-compression" rho values in column G of Haefner_WAVE_Simulations.xlsx. These values are ultimately used (in Haefner_WAVE_Simulations.xlsx) to calculate the "post-compression" S_f input values in column E of Haefner_WAVE_Simulations.xlsx

ro_weights.mat - Matlab file that contains the neural network weights generated by AutoKeras via the file Haefner_NN_for_Qp.py
