from task1_part1_legendre_poly import create_legendre_poly_dict
import numpy as np
import math
import pandas as panda
import pygmt as pygmt
from multiprocessing import Pool


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
        return (-1)**(n/2)*((3*(math.sqrt(e_power_two)**n)*(1-(n/2)+(5/2)*(J_N_2 / e_power_two)*n))/((n+1)*(n+3)*np.sqrt(2*n+1)))
    else: return 0
#print(J_hat_N_n(4))

#test1 = create_legendre_poly_dict(math.sin(math.radians(30.0)), 180)
#print(test1[100, 89])

#Now we collect data from the different file(s)

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

#Now we collect data from EGM2008 (SGG-UGM-2)
#Same as we did with GGM03S
EGM2008_model = panda.read_csv("./TBA4251_Geoid_Height_Spherical_Harmonic/data/SGG-UGM-2.txt", delim_whitespace=True, usecols=["L", "M", "C", "S"], index_col=False)
EGM2008_model_NMAX = 2190
EGM2008_n = EGM2008_model["L"]
EGM2008_m = EGM2008_model["M"]
EGM2008_c = np.array(EGM2008_model["C"])
EGM2008_s = np.array(EGM2008_model["S"])

EGM2008_multiindex = panda.MultiIndex.from_arrays([EGM2008_n, EGM2008_m], names=["n", "m"])
EGM2008_dataframe = panda.DataFrame(np.transpose(np.array([EGM2008_c, EGM2008_s])), columns=["C", "S"], index = EGM2008_multiindex)
#test2 = EGM2008_dataframe.loc[10, 1]["C"]
#print(test2)

#Defining functions for C_nm and q_nm = S_nm
def C_nm(n,m):
    return EGM2008_dataframe.loc[n, m]["C"]
    #return GGM03S_dataframe.loc[n,m]["C"]
def S_q_nm(n,m):
    return EGM2008_dataframe.loc[n, m]["S"]
    #return GGM03S_dataframe.loc[n,m]["S"]

def R_nm(n,m):
    if(m == 0):
        return C_nm(n,m) - J_hat_N_n(n)
    if(m != 0):
        return C_nm(n,m)
#Now, since we have a large function for the N_gravimetric function, we split it into different pieces
#Inner sum

def N_gravimetric_inner_sum(n, longitude, P_legendre_poly):
    current_sum = 0
    for m in range(n+1):
        current_sum += (R_nm(n,m)*math.cos(m * math.radians(longitude)) + S_q_nm(n,m)*math.sin(math.radians(longitude)*m)) * P_legendre_poly[(n, m)]
    return current_sum

#print(N_gravimetric_inner_sum(35, 20, create_legendre_poly_dict(0.5, 2190)))
#test = create_legendre_poly_dict(0.5, 180)
#print(test[(40, 3)])


def N_gravemetric_total_sum(latitude, longitude, n_max):
    build_legendre_poly = create_legendre_poly_dict(math.sin(math.radians(latitude)), n_max)
    outer_sum = 0
    for i in range(2, n_max+1):
        outer_sum += ((a/Radius)**i)*N_gravimetric_inner_sum(i, longitude, build_legendre_poly)
    total_sum = (GM_gravitational_parameter/(Radius*gravity_small_gamma))*outer_sum
    return total_sum

#Change to EGM2008_model_NMAX for EGM2008 model
print(N_gravemetric_total_sum(63, 10, EGM2008_model_NMAX))
#print(N_gravemetric_total_sum(61.6929259311394, 5.1957949286442, GGM03S_model_NMAX))


#Now we show the potensial results
def calc_geoid_for_GGM03():
    data = ['latitude,longitude,geoidheight']

    #This is for the world
    #latitudes = np.linspace(-90, 90, 361)
    #longitudes = np.linspace(-180, 180, 721)

    #This is for scandinavia ++
    latitudes = np.linspace(55, 73, 37)
    longitudes = np.linspace(-20, 40, 121)

    for i in latitudes:
        for j in longitudes:
            print(i, j)
            geoid_height_calculations = N_gravemetric_total_sum(i, j, GGM03S_model_NMAX)
            data.append(str(i) + ',' + str(j) + ',' + str(geoid_height_calculations))
    with open('geoid_calc_scandinacia_whole_Norway_GGM03S.csv', 'w') as new_file:
        new_file.write('\n'.join(data))

def calc_geoid_for_EGM2008():
    data = ['latitude,longitude,geoidheight']

    #This is for the world
    #latitudes = np.linspace(-90, 90, 361)
    #longitudes = np.linspace(-180, 180, 721)

    #This is for scandinavia ++
    latitudes = np.linspace(55, 73, 37)
    longitudes = np.linspace(-20, 40, 121)

    for i in latitudes:
        for j in longitudes:
            print(i, j)
            geoid_height_calculations = N_gravemetric_total_sum(i, j, EGM2008_model_NMAX)
            data.append(str(i) + ',' + str(j) + ',' + str(geoid_height_calculations))
    with open('geoid_calc_scandinavia_whole_Norway_EGM2008.csv', 'w') as new_file:
        new_file.write('\n'.join(data))

#calc_geoid_for_GGM03()
#calc_geoid_for_EGM2008()