import numpy as np
import task1_part1_legendre_poly
#Second part

#first we define constants given

ro_ave = 5517
a = 6371
ro_w = 1000

konstant_gldas = 123


#Now we define the formula (4) given in the assignment

def love_number(n):
    return

def delta_R_nm_inner_sums(n, m, latitude, longitude):
    sum = 0
    for i in range(n):
        for j in range(m):
            sum
    return sum


def delta_R_nm(n, m):
    total_sum = (1/(4*(np.pi))) * ((1+love_number(n))/(2*n+1)) * ((3*ro_w)/(a*ro_ave))
    return