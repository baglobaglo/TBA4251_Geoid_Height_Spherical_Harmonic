import numpy as np
from task1_part1_legendre_poly import create_legendre_poly_dict
import pandas as panda
#Second part

#first we define constants given

ro_ave = 5517
a = 6371
ro_w = 1000
love_number_MAX = 64

konstant_gldas = 123

ECCO_model = panda.read_csv("./TBA4251_Geoid_Height_Spherical_Harmonic/ECCO_data/2005.txt", delim_whitespace=True, usecols=["latitude", "longitude", "value"], index_col=False)
ECCO_model_latitude = ECCO_model["latitude"]
ECCO_model_longitude = ECCO_model["longitude"]
ECCO_model_value = ECCO_model["value"]

LOVE_number = panda.read_csv("./TBA4251_Geoid_Height_Spherical_Harmonic/data/load_love_number.txt", delim_whitespace=True, usecols=["n", "k"], index_col=False)
LOVE_number_value = LOVE_number["k"]

#Now we define the formula (4) given in the assignment

def ECCO_data(n):
    return ECCO_model_value[n-1]

def love_number(n):
    return LOVE_number_value[n-1]

def delta_R_nm_inner_sums(m, latitude, longitude, P_leg_poly):

    for n in range(1, m+1):
        inner_sum = ECCO_data(n) * np.cos(np.radians(n*longitude)) * P_leg_poly[(m, n)] * (((np.pi**2)*longitude*np.sin(np.radians(latitude)))/(32400))
    return inner_sum

def delta_R_nm(latitude, longitude, MAX_love):
    genereted_legendre_poly = create_legendre_poly_dict(np.cos(np.radians(latitude)), love_number_MAX)
    total_sum = 0
    for i in range(1, MAX_love+1):
            total_sum +=  (1/(4*(np.pi))) * ((1+love_number(i))/(2*i+1)) * ((3*ro_w)/(a*ro_ave)) * delta_R_nm_inner_sums(i, latitude, longitude, genereted_legendre_poly)
    return total_sum


def delta_q_nm_inner_sums(m, latitude, longitude, P_leg_poly):

    for n in range(1, m+1):
        inner_sum = ECCO_data(n) * np.sin(np.radians(n*longitude)) * P_leg_poly[(m, n)] * (((np.pi**2)*longitude*np.sin(np.radians(latitude)))/(32400))
    return inner_sum

def delta_q_nm(latitude, longitude, MAX_love):
    genereted_legendre_poly = create_legendre_poly_dict(np.cos(np.radians(latitude)), love_number_MAX)
    total_sum = 0
    for i in range(1, MAX_love+1):
            total_sum +=  (1/(4*(np.pi))) * ((1+love_number(i))/(2*i+1)) * ((3*ro_w)/(a*ro_ave)) * delta_R_nm_inner_sums(i, latitude, longitude, genereted_legendre_poly)
    return total_sum

print(delta_q_nm(170.5, -72.5, love_number_MAX))
print(delta_R_nm(170.5, -72.5, love_number_MAX))

def calc_coeff_for_ECCO_data():
    data = ['latitude,longitude,coeff_R_nm,coeff_q_nm']

    #This is for the world
    #latitudes = np.linspace(-90, 90, 361)
    #longitudes = np.linspace(-180, 180, 721)

    #This is for scandinavia ++
    latitudes = ECCO_model_latitude
    longitudes = ECCO_model_longitude

    for i in latitudes:
        for j in longitudes:
            print(i, j)
            r_nm_calculations = delta_R_nm(i, j, love_number_MAX)
            q_nm_calculations = delta_q_nm(i, j, love_number_MAX)

            data.append(str(i) + ',' + str(j) + ',' + str(r_nm_calculations) + ',' + str(q_nm_calculations))
    with open('ECCO_coeff_scandinavia.csv', 'w') as new_file:
        new_file.write('\n'.join(data))