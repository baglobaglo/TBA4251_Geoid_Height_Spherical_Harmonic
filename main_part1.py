from part1_legendre_poly import create_legendre_poly_dict
import numpy as np
import math
import pandas as panda
import pygmt as ptgmt


a = 6378137.0000 #m
b = 6356752.3141 #m
f_power_one = 298.257222101
e_power_two = 0.006694380023
Radius = 6371000.7900 #m
GM_gravitational_parameter = 3986005*10**8 #m^3s^-2
gravity_small_gamma = 9.81 #m s^-2

#given values for J and J_hat
J_N_2 = 0.108263*10**-2
J_hat_N_2 = -(J_N_2)/np.sqrt(5)

def J_hat_N_n(n):
    #if(n == 0): return "Try with a different n"
    if(n == 2): return J_hat_N_2
    if(n % 2 == 0):
        return (-1)**(n/2)*((3*np.exp(n)*(1-(n/2)+(5/2)*(J_N_2 / np.exp(2))*n))/((n+1)*(n+3)*np.sqrt(2*n+1)))
    else: return 0
#print(J_hat_N_n(4))

def N_gravimetric(R, latitude, longitude, N_max, m):
    return #G*M/R*gamma

#test1 = create_legendre_poly_dict(math.sin(math.radians(30.0)), 180)
#print(test1[100, 89])

#Now we collect data from the different file(s)
#####column_coefficents = ["n", "m", "C", "S"]
#First model, GGM03S. Our variables is calles n, m, C and S
#but in the model txt file it is listed as L, M ,C S
GGM03S_model = panda.read_csv("./TBA4251_Geoid_Height_Spherical_Harmonic/data/GGM03S.txt", delim_whitespace=True, usecols=["L", "M", "C", "S"], index_col=False)
GGM03S_model_NMAX = 180
GGM03S_n = GGM03S_model["L"]
GGM03S_m = GGM03S_model["M"]
GGM03S_c = np.array(GGM03S_model["C"])
GGM03S_s = np.array(GGM03S_model["S"])

#Now we use multiindex, since MultiIndex in Pandas is a multi-level or hierarchical object that allows
#you to select more than one row and column in your index
#In that way it is more similar to the way we have it in our dictionary
GGM03S_multiindex = panda.MultiIndex.from_arrays([GGM03S_n, GGM03S_m], names=["n", "m"])

#Now we create the dataframe, which is similar to SQL spreadsheet or dictionary
GGM03S_dataframe= panda.DataFrame(np.transpose(np.array([GGM03S_c, GGM03S_s])), columns=["C", "S"], index = GGM03S_multiindex)

#test1 = GGM03S_dataframe.loc[40, 0]["C"]
#print(test1)

#Defining functions for C_nm and q_nm = S_nm
def C_nm(n,m):
    return GGM03S_dataframe.loc[n,m]["C"]
def S_q_nm(n,m):
    return GGM03S_dataframe.loc[n,m]["S"]

def R_nm(n,m):
    if(m == 0):
        C_nm(n,m) - J_hat_N_n(n)
    if(m != 0):
        C_nm(n,m)

#Now, since we have a large function for the N_gravimetric function, we split it into different pieces
#Inner sum

def N_gravimetric_inner_sum(n, longitude, P_legendre_poly):
    current_sum = 0
    for m in range(n+1):
        print(f"n: {n}, m: {m}, R_bar: {R_nm(n,m)}, q_nm: {S_q_nm(n,m)}")
        current_sum += (R_nm(n,m)*math.cos(math.radians(longitude)*m) + (S_q_nm(n,m)*math.sin(math.radians(longitude)*m)))*P_legendre_poly[(n, m)]
    return current_sum
print(N_gravimetric_inner_sum(40, 60, create_legendre_poly_dict(math.sin(math.radians(60)), 40)))