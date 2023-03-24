import numpy as np
import math
import pandas as panda
import pygmt as ptgmt

#conda activate pygmt


#We use a dictionary to store all possible values for n and m
#in way that may remind you of how it is done i matlab with matrices

def first_initial_legendre_poly_values(sin_phi_t):
    legendre_poly_dict = {
        (0, 0): 1,
        (1, 0): sin_phi_t*math.sqrt(3),
        (1, 1): math.sqrt(3)*math.sqrt(1-sin_phi_t**2),
        (2, 0): math.sqrt(5)*((3/2)*sin_phi_t**2 - (1/2)),
        (2, 1): sin_phi_t * math.sqrt(15)*math.sqrt(1-sin_phi_t**2)
        }
    return legendre_poly_dict
#print(first_initial_legendre_poly_values(math.sin(math.radians(30.0))))

#Now we define the functionality and logic for the different legendre polynomials

def create_legendre_poly_dict(sin_phi_t, n):
    initial_legendre_poly = first_initial_legendre_poly_values(sin_phi_t)
    for i in range(2, n+1):
        for j in range(i+1):
            if (i>=3 and 1 <= j <= i-2):
                initial_legendre_poly[(i, j)] = -math.sqrt(((2*i+1) * (i+j-1) * (i-j-1))/((2*i-3) * (i+j) * (i-j))) * initial_legendre_poly[(i-2, j)] + sin_phi_t * math.sqrt(((2*i+1) * (2*i-1))/((i+j)*(i-j))) * initial_legendre_poly[(i-1, j)]
            if (i >= 2 and j == 0):
                initial_legendre_poly[(i, j)] = -((math.sqrt(2*i+1))/(i)) * (i-1)/(math.sqrt(2*i-3)) * initial_legendre_poly[(i-2, 0)] + sin_phi_t * ((math.sqrt(2*i+1))/(i)) * math.sqrt(2*i-1) * initial_legendre_poly[(i-1, 0)]
            if (i >= 2 and j == i):
                initial_legendre_poly[(i, j)] = math.sqrt((2*i+1)/(2*i)) * math.sqrt(1-sin_phi_t**2) * initial_legendre_poly[(i-1, i-1)]
            if (i >= 1 and j == i-1):
                initial_legendre_poly[(i,j)] = sin_phi_t * math.sqrt(2*i + 1) * initial_legendre_poly[(i-1, i-1)]
    return initial_legendre_poly

#build_it = create_legendre_poly_dict(math.sin(math.radians(30.0)), 180)

#print(build_it[100, 89])


"""
#Next we define the different cases of the Normalized Legendre's polynomials
# t is defined as sin(latitude)
def norm_legendre_poly_part1(n, t):
    #we try to do this iterative instead of recursice, since recursive functions
    #in python can be very slow. This function handles cases where m = 0 and n >=2
    p_previous_0 = 1
    p_previous_1 = t * math.sqrt(3)
    p = 0.0

    #using the recurison formula
    for i in range(2, n+1):
        p = (-(math.sqrt((2*i+1))/(i))*((i-1)/(math.sqrt(2*i-3))) *p_previous_0) + t*((math.sqrt(2*i+1))/(i))*(math.sqrt(2*i-1))*p_previous_1
        # update the previous two polynomials for the next iteration
        p_previous_0 = p_previous_1
        p_previous_1 = p
    return p
#print(norm_legendre_poly_part1(2, math.sin(math.radians(63.0))))

def norm_legendre_poly_part4(n, t):
    # for cases where n>=2; m = n
    if (n == 1):
        return math.sqrt(3)*math.sqrt(1-t**2)
    p_previous_0 = 1
    p_previous_1 = math.sqrt(3)*math.sqrt(1-t**2)
    p_n_n = 0

    for i in range(2, n+1):
        p_n_n = math.sqrt((2*i+1)/(2*i))*math.sqrt(1-t**2)*p_previous_1
        p_previous_0 = p_previous_1
        p_previous_1 = p_n_n
    return p_n_n

#print(norm_legendre_poly_part4(3, math.sin(math.radians(63.0))))

def norm_legendre_poly_part3(n, t):
    if (n == 1):
        return t*math.sqrt(3)
    if (n == 2):
        return t*math.sqrt(15)*math.sqrt(1-t**2)
    #for cases where n>=1; m = n - 1
    p_n_n_minus1 = 0
    for i in range(3, n+1):
        p_n_n_minus1 = t*math.sqrt(2*i+1)*norm_legendre_poly_part4(n-1, t)
    return p_n_n_minus1
#print(norm_legendre_poly_part3(3, math.sin(math.radians(63.0))))

#print(norm_legendre_poly_part4(2, math.sin(math.radians(63.0))))
#print(norm_legendre_poly_part4(3, math.sin(math.radians(63.0))))
#print(norm_legendre_poly_part4(4, math.sin(math.radians(63.0))))
#print(norm_legendre_poly_part4(5, math.sin(math.radians(63.0))))

#if (m == n-2):
        #p_n_m += -math.sqrt(((2*n+1)*(n+m-1)*(n-m-1))/((2*n-3)*(n+m)*(n-m)))*norm_legendre_poly_part4(n-2, t) + t*math.sqrt(((2*n+1)*(2*n-1))/((n+m)*(n-m)))*norm_legendre_poly_part3(n-1, t)
        #return p_n_m

def norm_legendre_poly_part2_3(n, m, t):
    #we try to do this next part iterative as well.
    #this function handles cases where n >= 3 and 1<= m <= n-2
    if(n == 1 & m == 1):
        return math.sqrt(3)*math.sqrt(1-t**2)
    if(n == 2 & m == 1):
        return t*math.sqrt(15)*math.sqrt(1-t**2)
    p_n_m = -(math.sqrt(((2*n+1)*(n+m-1)*(n-m-1))/((2*n-3)*(n+m)*(n-m))))*norm_legendre_poly_part2_3(n-2, m, t) + t*math.sqrt(((2*n+1)*(2*n-1))/((n+m)*(n-m)))*norm_legendre_poly_part2_3(n-1, m, t)
    return p_n_m

#print(norm_legendre_poly_part2_3(3, 1, math.sin(math.radians(30.0))))

def norm_legendre_poly_part2(n, m, t):
    if n == 0 and m == 0:
        return 1
    elif n == 1 and m == 0:
        return t*math.sqrt(3)
    elif n == 1 and m == 1:
        return math.sqrt(3)*math.sqrt(1 - t**2)
    elif n == 2 and m == 1:
        return t*math.sqrt(15)*math.sqrt(1-t**2)
    elif n >= 3 and 1 <= m <= n-2:
        return -(math.sqrt(((2*n+1)*(n+m-1)*(n-m-1))/((2*n-3)*(n+m)*(n-m))))*norm_legendre_poly_part2(n-2, m, t) \
                + t*math.sqrt(((2*n+1)*(2*n-1))/((n+m)*(n-m)))*norm_legendre_poly_part2(n-1, m, t)
    else:
        raise ValueError(f"Invalid values of n={n} and m={m}")
print(norm_legendre_poly_part2(17, 5, math.sin(math.radians(30))))  # Computes the value for n=3, m=1 and t=0.5

def noram_legendre_poly_part2_test(n, m, t):
    if n == 0 and m == 0:
        return 1
    elif n == 1 and m == 0:
        return t*math.sqrt(3)
    elif n == 1 and m == 1:
        return math.sqrt(3)*math.sqrt(1 - t**2)
    elif n == 2 and m == 1:
        return t*math.sqrt(15)*math.sqrt(1-t**2)
    p_n_m = 0
    for i in range(3, n+1):
        if i > 3 and 1<=m<=i-2:
            p_n_m += norm_legendre_poly_part2(i, m, t)
        if i > 3 and not (1<=m<=i-2):
            eps = 1e-15
            p_n_m += -(math.sqrt(((2*i+1)*(i+m-1)*(i-m-1))/((2*i-3)*(i+m)*(i-m)) + eps))*norm_legendre_poly_part3(i, t) \
                + t*math.sqrt(((2*i+1)*(2*i-1))/((i+m)*(i-m) + eps))*norm_legendre_poly_part4(i, t)
    return p_n_m
#print(noram_legendre_poly_part2_test(3, 1, math.sin(math.radians(30))))  # Computes the value for n=3, m=1 and t=0.5


def norm_legendre_poly_new_test(n, m, t):
    if n == 1 and m == 1:
        n_0 += math.sqrt(3)*math.sqrt(1-t**2)
        return n_0
    elif n == 2 and m == 1:
        n_0 += t*math.sqrt(15)*math.sqrt(1-t**2)
        return n_0
    elif m == n-1:
        n_0 += norm_legendre_poly_part3(n+1, t)
        return n_0
    elif m == n:
        n_0 += norm_legendre_poly_part4(n, t)
        return n_0
    else:
        n_2 = 1
        n_1 = math.sqrt(3)*math.sqrt(1-t**2)
        for k in range(2, n+1):
            n_0 = -(math.sqrt(((2*k+1)*(k+m-1)*(k-m-1))/((2*k-3)*(k+m)*(k-m))))*n_2 + t*math.sqrt(((2*k+1)*(2*k-1))/((k+m)*(k-m)))*n_1
            n_2 = n_1
            n_1 = n_0
        return n_0

#print(norm_legendre_poly_new_test(50, 34, math.sin(math.radians(30))))  # Computes the value for n=3, m=1 and t=0.5



# ----- COMBINE ALLE POLYNOMIALS ------- #
def combine_all_legendre_poly(n, m, t):
    if n >= 2 and m == 0:
        return norm_legendre_poly_part1(n, t)
    if n >= 3 and 1 <= m <= n-2:
        return norm_legendre_poly_part2(n, m, t)
    if n >= 1 and m == n - 1:
        return norm_legendre_poly_part3(n, t)
    if n >= 2 and m == n:
        return norm_legendre_poly_part4(n, t)


#print(combine_all_legendre_poly(100, 89, math.sin(math.radians(30))))
"""