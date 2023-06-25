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

#Now we define the functionality and logic for the different legendre polynomials given different n and m values

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

#build_it = create_legendre_poly_dict(0.5, 2190)
#print(build_it[100, 89])
