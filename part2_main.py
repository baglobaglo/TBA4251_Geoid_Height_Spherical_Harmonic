import numpy as np
from part1_legendre_poly import create_legendre_poly_dict
import pandas as panda
#Second part

#first we define constants given

ro_ave = 5517
a = 6371
ro_w = 1000
love_number_MAX = 64

konstant_gldas = 123
#change 2006-20014 for different calculations
ECCO_model = panda.read_csv("./TBA4251_Geoid_Height_Spherical_Harmonic/ECCO_data/2006.txt", delim_whitespace=True, usecols=["latitude", "longitude", "value"], index_col=False)
ECCO_model_latitude = ECCO_model["latitude"]
ECCO_model_longitude = ECCO_model["longitude"]
ECCO_model_value = ECCO_model["value"]

LOVE_number = panda.read_csv("./TBA4251_Geoid_Height_Spherical_Harmonic/data/load_love_number.txt", delim_whitespace=True, usecols=["n", "k"], index_col=False)
LOVE_number_value = LOVE_number["k"]


calculate_coeffs_2005_ECCO = panda.read_csv("ECCO_coeff_2005.csv")
calculate_coeffs_2005_ECCO_R_nm = calculate_coeffs_2005_ECCO["coeff_R_nm"]
calculate_coeffs_2005_ECCO_q_nm = calculate_coeffs_2005_ECCO["coeff_q_nm"]
print(calculate_coeffs_2005_ECCO_q_nm)

#Now we define the formula (4) given in the assignment

def ECCO_data(n):
    return ECCO_model_value[n-1]

def love_number(n):
    return LOVE_number_value[n-1]

#divide the formula into two parts, one for R_nm and one for q_nm
#and for inner and out sum to get rid of complexity

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

#test
#print(delta_q_nm(170.5, -72.5, love_number_MAX))
#print(delta_R_nm(170.5, -72.5, love_number_MAX))

def calc_coeff_for_ECCO_data():
    data = ['latitude,longitude,coeff_R_nm,coeff_q_nm']
    latitudes = ECCO_model_latitude
    longitudes = ECCO_model_longitude

    for i in range(0, len(ECCO_model)):
        print(latitudes[i], longitudes[i])
        r_nm_calculations = delta_R_nm(latitudes[i], longitudes[i], love_number_MAX)
        q_nm_calculations = delta_q_nm(latitudes[i], longitudes[i], love_number_MAX)

        data.append(str(latitudes[i]) + ',' + str(longitudes[i]) + ',' + str(r_nm_calculations) + ',' + str(q_nm_calculations))
    with open('ECCO_coeff_2006.csv', 'w') as new_file:
        new_file.write('\n'.join(data))

#calc_coeff_for_ECCO_data()

def delta_sigma_calc_inner_sum(m, longitude, P_leg_poly):
    for n in range(1, m+1):
        innersum = (2*n +1)/(1+love_number(n)) * ( calculate_coeffs_2005_ECCO_R_nm[n-1]* np.cos(np.radians(n*longitude)) + calculate_coeffs_2005_ECCO_q_nm[n-1] * np.sin(np.radians(n*longitude))) * P_leg_poly[(m, n)]
    return innersum

def delta_sigma_calculations(latitude, longtitude, max_love):
    build_legendre_poly_dict = create_legendre_poly_dict(np.cos(np.radians(latitude)), max_love)
    sum = 0
    for n in range(1, love_number_MAX+1):
        sum += ((a*ro_ave)/3) * delta_sigma_calc_inner_sum(n, longtitude, build_legendre_poly_dict)
    return sum
#print(delta_sigma_calculations(170.5, -72.5, love_number_MAX))

#function which uses the coefficients calculated in the previous function to calculate the mass change
#will work for both GLDAS and ECCO data
def calc_massChange_for_ECCO_data():
    data = ['latitude,longitude,mass_change']
    latitudes = ECCO_model_latitude
    longitudes = ECCO_model_longitude

    for i in range(0, len(calculate_coeffs_2005_ECCO)):
        print(latitudes[i], longitudes[i])
        calculations_mass_change = delta_sigma_calculations(latitudes[i], longitudes[i], love_number_MAX)

        data.append(str(latitudes[i]) + ',' + str(longitudes[i]) + ',' + str(calculations_mass_change))
    with open('ECCO_mass_change_calc_2006.csv', 'w') as new_file:
        new_file.write('\n'.join(data))


calc_massChange_for_ECCO_data()