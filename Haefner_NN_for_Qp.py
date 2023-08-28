"""
Matthew Haefner
SEA Lab
NN for Q_p
12/3/22
"""

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import math
import autokeras as ak
from string import ascii_uppercase
import itertools


# Note: BEFORE RUNNING, one must change the file path assigned to the variable "excelpath" in line 171


# Read in data, separate into training sets
# (for non-adjacent columns can use np.r_[1,3:5,etc.])
df = pd.read_excel('Haefner_WAVE_Simulations_Trimmed.xlsx')
x = df.iloc[:,1:4]
y = df.iloc[:,10]
x_train,x_test,y_train,y_test = train_test_split(x,y,test_size=0.2)


# AutoKeras Functions
reg = ak.StructuredDataRegressor(
    column_names=[
        "P_f [psi]",
        "S_f [g/kg]",
        "Q_f [m^3/hr]",
    ],
    column_types={"P_f [psi]": "numerical", "S_f [g/kg]": "numerical", "Q_f [m^3/hr]": "numerical"},
    loss = "mean_squared_error",
    max_trials=100,  # Default value; number of different models tried (I think)
    overwrite=True,
)

reg.fit(x_train, y_train, epochs=1000) # Default value; number of forward and back propagations done (I think)

predicted_y = reg.predict(x_test)

model = reg.export_model()
model.summary()

# Get the weights from the model, store them in a dictionary
d = {}
i = 1
for layer in model.layers:
    d["weights" + str(i)] = layer.get_weights()
    i = i+1

# Take the resulting dictionary and convert each entry into a data frame
def iter_all_strings():
    for size in itertools.count(1):
        for s in itertools.product(ascii_uppercase, repeat=size):
            yield "".join(s)

letter_list = []
for s in iter_all_strings():
    letter_list.append(s)
    if s == 'AMZ':
        break

df_list = []
#n_nodes_list = []
for i in range(1,len(d)+1):
    if d["weights" + str(i)] != []:
        for j in range(len(d["weights" + str(i)])):
            if i == 3 and j == 2:
                pass
            else:
                mat = np.matrix(d["weights" + str(i)][j])
                df_list.append(pd.DataFrame(mat,columns = letter_list[0:mat.shape[1]]))
                
# Next, need to verify the model (make sure that forward propagation matches model.predict())
# Function for normalizing the data
def normalization(A_df,df_list):
    # First, need to convert data frames to numpy arrays
    A = np.matrix(A_df)
    averages = np.matrix(df_list[0])
    variances = np.matrix(df_list[1])
    
    # Set up variables
    m = A.shape[0]
    n = A.shape[1]
    A_normalized = np.zeros([m,n])
    
    # Normalize each element of A
    for j in range(n):
        a = averages[0,j]
        v = variances[0,j]
        s = math.sqrt(v)
        for i in range(m):
            A_normalized[i,j] = (A[i,j]-a)/s
            
    return A_normalized

# Function for batch normalization
def batch_normalization(batch_in,gamma_df,beta_df,mean_df,std_df):
    # First must convert all data frames to matrics
    gamma = np.matrix(gamma_df)
    beta = np.matrix(beta_df)
    mean = np.matrix(mean_df)
    std = np.matrix(std_df)
    
    # Initialize variables
    n = batch_in.shape[1]
    batch_out = np.zeros([1,n])
    
    # Do batch normalization
    for i in range(n):
        batch_out[i] = (gamma[i] * (batch_in[i] - mean[i]) / math.sqrt((std[i]^2)+0.001)) + beta[i]
    return batch_out

# Forward propagation helper function
def single_forward_prop(a_i,Theta_df,bias_df):
    # Need to first convert theta and bias to numpy matrices
    Theta = np.matrix(Theta_df)
    bias = np.matrix(bias_df)
    
    # Set up variables
    num_from_nodes = Theta.shape[1]
    a_next = np.zeros([1,num_from_nodes])
    
    for j in range(num_from_nodes):
        T_j = Theta[:,j] # num_to_nodesx1
        z_next_j = np.dot(a_i,T_j)+bias[0,j] # 1x1
        a_next_j = max(z_next_j,0) # ReLU function, 1x1
        a_next[0,j] = a_next_j
        
    return a_next

# Forward propagation function
def forward_prop(x_df,df_list):
    # First need to normalize x_df (which also converts it to a matrix)
    x = normalization(x_df,df_list)
    
    # Set up variables
    m = x.shape[0] # number of samples
    y = np.zeros([m,1]) # one solution for each sample
    l = len(df_list)
    
    for i in range(m):
        a1_i = x[i,:] # one sample at a time, 1x3
        a = a1_i
        j = 2
        while j < l:
            if j == l-2:
                y[i,0] = np.dot(a,np.matrix(df_list[-2]))+np.matrix(df_list[-1])
                j = j+2
            elif all(df_list[j].columns.values==df_list[j+1].columns.values) and df_list[j].size == df_list[j+1].size and all(df_list[j+1].columns.values==df_list[j+2].columns.values) and df_list[j+1].size == df_list[j+2].size and all(df_list[j+2].columns.values==df_list[j+3].columns.values) and df_list[j+2].size == df_list[j+3].size:
                a = batch_normalization(a,df_list[j],df_list[j+1],df_list[j+2],df_list[j+3])
                j = j+4
            else:
                a = single_forward_prop(a, df_list[j], df_list[j+1])
                j = j+2
            
    return y         
    
# Actually veryifying things here
my_results = forward_prop(x_test,df_list)
if np.allclose(my_results.round(decimals=4),predicted_y.round(decimals=4)) == True:
    print('All is good!')
    # If all set, save each data frame to a sheet in an excel document
    # Excel path
    excelpath = r'C:\Users\mwh85\Documents\Applied Energy Paper\Analysis While Waiting for Feedback\RO_weights.xlsx'

    # Write my dataframes to difference sheets
    w = 1
    w_max = len(df_list)
    with pd.ExcelWriter(excelpath) as writer:
        while w <= w_max:
            df_list[w-1].to_excel(writer,sheet_name = 'Weights'+str(w))
            w = w+1
else:
    print('Something is off')