import numpy as np
from task1_part1_legendre_poly import create_legendre_poly_dict
#Second part

#first we define constants given

ro_ave = 5517
a = 6371
ro_w = 1000
love_number_MAX = 100

konstant_gldas = 123


#Now we define the formula (4) given in the assignment
def dA():
    return
def ECCO_data():
    return 3.21

def love_number(n):
    return n+1.0

def delta_R_nm_inner_sums(m, latitude, longitude, P_leg_poly):

    for n in range(1, m+1):
        inner_sum = np.cos(np.radians(n*longitude)) * P_leg_poly[(m, n)] * (((np.pi**2)*longitude*np.sin(np.radians(latitude)))/(32400))
    return inner_sum


def delta_R_nm(latitude, longitude, MAX_love):
    genereted_legendre_poly = create_legendre_poly_dict(np.cos(np.radians(latitude)), love_number_MAX)
    total_sum = 0
    for i in range(1, MAX_love+1):
            total_sum +=  (1/(4*(np.pi))) * ((1+love_number(i))/(2*i+1)) * ((3*ro_w)/(a*ro_ave)) * delta_R_nm_inner_sums(i, latitude, longitude, genereted_legendre_poly)
    return total_sum

print(delta_R_nm(63, 10, love_number_MAX))