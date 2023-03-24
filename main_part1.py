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
