import numpy as np
import math


#Geoid Height Computations from the Global Gravity Models
#we start by defining some of the known and given constants

a = 6378137.0000 #m
b = 6356752.3141 #m
f_power_one = 298.257222101
e_power_two = 0.006694380023
Radius = 6371000.7900 #m
GM_gravitational_parameter = 3986005*10**8 #m^3s^-2
gravity_small_gamma = 9.81 #m s^-2


#next we try to define and structure the smaller function relying on eachother

#given values for J and J_hat
J_N_2 = 0.108263*10**-2
J_hat_N_2 = -(J_N_2)/np.sqrt(5)

def J_hat_N_n(n):
    if(n == 0): return "Try with a different n"
    if(n == 2): return J_hat_N_2
    elif(n % 2 == 0):
        return (-1)**(n/2)*((3*np.exp(n)*(1-(n/2)+(5/2)*(J_N_2 / np.exp(2))*n))/((n+1)*(n+3)*np.sqrt(2*n+1)))
    else: return 0
#print(J_hat_N_n(4))


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

def norm_legendre_poly_part2(n, m, t):
    #we try to do this next part iterative as well.
    #this function handles cases where n >= 3 and 1<= m <= n-2
    p_previous_0 = math.sqrt(3)*math.sqrt(1-t**2)
    p_previous_1 = t*math.sqrt(15)*math.sqrt(1-t**2)
    p_n_m = 0.0
    for i, j in range(n, n+2):
        n*m
    return

def norm_legendre_poly_part3(n, m, t):
    m = n - 1
    #for cases where n>=1; m = n - 1
    p_previous_0 = 1
    p_previous_1 = math.sqrt(3)*math.sqrt(1-t**2)
    p_n_n_minus1 = 0
    for i in range(1, n+1):
        p_n_n_minus1 = t*math.sqrt(2*i+1)*p_previous_0
        p_previous_0 = p_previous_1
        p_previous_1 = p_n_n_minus1

    return p_n_n_minus1

#print(norm_legendre_poly_part3(2, 0, math.sin(math.radians(63.0))))

def norm_legendre_poly_part4(n, t):
    # for cases where n>=2; m = n

    p_previous_0 = 1
    p_previous_1 = math.sqrt(3)*math.sqrt(1-t**2)
    p_n_n = 0

    for i in range(2, n+1):
        p_n_n = math.sqrt((2*i+1)/(2*i))*math.sqrt(1-t**2)*p_previous_1
        p_previous_0 = p_previous_1
        p_previous_1 = p_n_n

    return p_n_n

#print(norm_legendre_poly_part4(2, math.sin(math.radians(63.0))))
#print(norm_legendre_poly_part4(3, math.sin(math.radians(63.0))))
#print(norm_legendre_poly_part4(4, math.sin(math.radians(63.0))))
#print(norm_legendre_poly_part4(5, math.sin(math.radians(63.0))))


def combine_all_legendre_poly(n, m, t):
    return
