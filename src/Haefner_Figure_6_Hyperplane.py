"""
Matthew Haefner
SEA Lab
Linear C-SVM for defining search region for WAVE full factorial
7/2/22
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from sklearn import svm

# Read in data, trim
df = pd.read_excel('..\data\Haefner_Figure_6_Data.xlsx')

data = df.to_numpy()
data_trimmed = data[19:68+1,0:31+1]

# Separate legal runs from illegal runs
good_data = np.zeros(data_trimmed.shape[1])
bad_data = np.zeros(data_trimmed.shape[1])

for i in range(data_trimmed.shape[0]):
    if data_trimmed[i,12] == 0:
        bad_data = np.vstack([bad_data,data_trimmed[i,:]])
    else:
        good_data = np.vstack([good_data,data_trimmed[i,:]])

good_data = np.delete(good_data,0,0)
bad_data= np.delete(bad_data,0,0)

# Extract plotting data        
Pf_good = good_data[:,1]
Pf_bad = bad_data[:,1]

Sf_good = good_data[:,2]
Sf_bad = bad_data[:,2]

Qf_good = good_data[:,5]
Qf_bad = bad_data[:,5]

# Specifying 3d plot
ax = plt.axes(projection = '3d')

# First plotting the good points
xdata_good = Pf_good
ydata_good = Sf_good
zdata_good = Qf_good
ax.scatter3D(xdata_good,ydata_good,zdata_good)

# Then plotting the bad points
xdata_bad = Pf_bad
ydata_bad = Sf_bad
zdata_bad = Qf_bad
ax.scatter3D(xdata_bad,ydata_bad,zdata_bad)

# Axis labels
ax.set_xlabel('P_f [psi]')
ax.set_ylabel('S_f [g/kg]')
ax.set_zlabel('Q_f [m^3/hr]')

# Change view angle
ax.view_init(40,-105)

#SVM stuff 
# Filling X and Y arrays
X = data_trimmed[:,1]
X = np.vstack([X,data_trimmed[:,2]])
X = np.vstack([X,data_trimmed[:,5]])
X = X.transpose()
Y = []
for i in range(data_trimmed.shape[0]): # 1 is bad, 0 is good
    if data_trimmed[i,12] == 0:
        Y = np.append(Y,[1])
    else:
        Y = np.append(Y,[0])

# Classification problem
model = svm.SVC(kernel='linear')
clf = model.fit(X, Y)

# The equation of the separating plane is given by all x so that np.dot(svc.coef_[0], x) + b = 0.
# Solve for w3 (z)
z = lambda x,y: (-clf.intercept_[0]-clf.coef_[0][0]*x -clf.coef_[0][1]*y) / clf.coef_[0][2]

# Generating mesh grid
x_range = np.linspace(300,1200,91)
y_range = np.linspace(35,55,21)
x,y = np.meshgrid(x_range,y_range)

# Creating figure and plotting
fig_svm = plt.figure()
ax_svm = fig_svm.add_subplot(111, projection='3d')
ax_svm.plot3D(X[Y==0,0], X[Y==0,1], X[Y==0,2],'ob')
ax_svm.plot3D(X[Y==1,0], X[Y==1,1], X[Y==1,2],'sr')
ax_svm.plot_surface(x, y, z(x,y))
ax_svm.view_init(30, 60)
#ax_svm.view_init(0, 70+180) # 30 60
#ax_svm.view_init(11, 90+180) # 30 60

# Axis labels
ax_svm.set_xlabel('P_f [psi]')
ax_svm.set_ylabel('S_f [g/kg]')
ax_svm.set_zlabel('Q_f [m^3/hr]')

#plt.show()

# Getting values of z(x,y) so I can define the plane
Qf_vals_mesh = z(x,y)

# Select 3 points
P1 = np.array([x_range[20],y_range[5],Qf_vals_mesh[5,20]])
P2 = np.array([x_range[40],y_range[9],Qf_vals_mesh[9,40]])
P3 = np.array([x_range[70],y_range[16],Qf_vals_mesh[16,70]])

# Make two vectors in the plane
V12 = P2-P1
V13 = P3-P1

# Normal vector (this will be the coefficients of the hyperplane)
NV = np.cross(V12,V13)
