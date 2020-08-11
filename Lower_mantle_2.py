import pylab
from math import *
from scipy.constants import G
from scipy import optimize
import matplotlib.pyplot as plt
import numpy as np


def runkut(n, x, y, h):
#"Advances the solution of diff eqn defined by derivs from x to x+h"
    y0=y[:]
    k1=derivs(n, x, y)
    for i in range(1,n+1): y[i]=y0[i]+0.5*h*k1[i]
    k2=derivs(n, x+0.5*h, y)
    for i in range(1,n+1): y[i]=y0[i]+h*(0.2071067811*k1[i]+0.2928932188*k2[i])
    k3=derivs(n, x+0.5*h, y)
    for i in range(1,n+1): y[i]=y0[i]-h*(0.7071067811*k2[i]-1.7071067811*k3[i])
    k4=derivs(n, x+h, y)
    for i in range(1,n+1):1
    a=k1[i]+0.5857864376*k2[i]+3.4142135623*k3[i]+k4[i]
    y[i]=y0[i]+0.16666666667*h*a
	
    x+=h
    return (x,y)

#----------------------------------------------------------------------------

# B-M-Stixrude parameters for Mg-post-perovskite (1st-lower mantle) and Mg-perovskite (2nd-upper mantle)

K_0 = [252.58, 230.6]  
K_0_prime = [4.08, 4.0]
rho_0 = [4256.95, 4263.02]  
V_0 = [24.43, 24.49]    
theta_0 = [887.96, 838.5]  
gamma_0 = [1.584, 1.887]
beta = 1.1


p_ini = 130.0      # Initial pressure at the core-mantle boundary (Pcmb at r_cmb)
T_ini = 2528.90    # Initial temperature at the core-mantle boundary (Tcmb at r_cmb)
m_ini = 3.82965753582e+24
density = 5754.37
g_ini = 13.37
r_ini = 4370.90
q_ini = 0.00118

T_0 = 300.0 # Reference temperature in kelvins
n0 = 1.0
R = 8.314
epsilon = 7.38 * (10**(-11))

density = []


def integrate(f, a, b, n):

    h = float(b - a) / n
    result = (0.5 * f(a)) + (0.5 * f(b))

    for i in range(1, n):
        result += f(a + (i*h))
    result *= h

    return result


def differentiate(f, a, var='volume', h = 0.01):
    
    if var == 'volume':
        val = (a * V_0[1])
        t1 = f((val + h) / V_0[1])
        t2 = f((val - h) / V_0[1])
    elif var == 'temperature':
        t1 = f(a + h)
        t2 = f(a - h)

    return (t1 - t2) / (2.0 * h)


th = lambda t: (t**3.0) / (float(exp(t)+0.0000001) - 1.0)
c_v = lambda z: ((z**4.0) * exp(z)) / (((exp(z)+0.0000001) - 1.0)**2.0)



def f_BM_Stixrude(x, K_0, K_0_prime, theta_0, gamma_0, p, T):   

    f = (1.0 / 2.0) * ((x**(-2.0 / 3.0)) - 1.0)                # Birch-Murnaghan-Stixrude eqn. Calculated with x = V / V_0
    gamma = gamma_0 * (x**beta)
    theta_x = (theta_0**2.0) * (1 + (6.0 * gamma_0 * f) + (((-6.0 * gamma_0) + (18.0 * (gamma_0**2.0)) - (9.0 * beta * gamma_0))* (f**2.0)))
    if theta_x < 0:
        theta = ((-1.0 * theta_x))**0.5
    else:
        theta = (theta_x)**0.5

    ta = (3.0 * K_0 * f) * ((1.0 + (2.0 * f))**(5.0 / 2.0))
    tb = (1.0 + ((3.0 / 2.0) * (K_0_prime - 4.0) * f))

    V = x * V_0[1]

    "Calculates the thermal pressure for the BM-Stixrude regime."

    diff = ((T**4.0) * (integrate(th, 0.0, (theta/T), 600))) - ((T_0**4.0) * (integrate(th, 0.0, (theta/T_0), 600)))
    Pth = (((9.0 * gamma * n0 * R) / (V * (theta**3.0))) * diff)

    return (ta * tb) + (Pth / 1000.0)- p



def K_BM(x, K_0, K_0_prime, theta_0, gamma_0, T):

    f = (1.0 / 2.0) * ((x**(-2.0 / 3.0)) - 1.0)
    gamma = gamma_0 * (x**beta)
    theta_x = (theta_0**2.0) * (1 + (6.0 * gamma_0 * f) + (((-6.0 * gamma_0) + (18.0 * (gamma_0**2.0)) - (9.0 * beta * gamma_0))* (f**2.0)))
    theta = (theta_x)**0.5

    ext = ((1.0 + (2.0 * f))**(5.0 / 2.0))
    tc = K_0 + ((3.0 * K_0 * K_0_prime) - (5.0 * K_0)) * f
    td = (27.0 / 2.0) * ((K_0 * K_0_prime) - (4.0 * K_0)) * (f**2.0)
    K_i = ext * (tc + td)

    V = x * V_0[1]

    "Calculates the thermal pressure for the bulk modulus."

    diff_a = ((T**4.0) * (integrate(th, 0.0, (theta/T), 600))) - ((T_0**4.0) * (integrate(th, 0.0, (theta/T_0), 600)))
    E_har = ((9.0 * gamma * n0 * R) / (V * (theta**3.0))) * diff_a
    diff_b = ((T**4.0) * (integrate(c_v, 0.0, (theta/T), 600))) - ((T_0**4.0) * (integrate(c_v, 0.0, (theta/T_0), 600)))
    C_v = ((9.0 * (gamma**2.0) * n0 * R) / (V * (theta**3.0))) * diff_a

    K_th = ((gamma + 1.0 - beta) * E_har) - C_v

    return K_i + K_th


def p_th(x, K_0, K_0_prime, theta_0, gamma_0, T):   

    f = (1.0 / 2.0) * ((x**(-2.0 / 3.0)) - 1.0)                
    gamma = gamma_0 * (x**beta)
    theta_x = (theta_0**2.0) * (1 + (6.0 * gamma_0 * f) + (((-6.0 * gamma_0) + (18.0 * (gamma_0**2.0)) - (9.0 * beta * gamma_0))* (f**2.0)))
    theta = (theta_x)**0.5

    V = x * V_0[1]

    "Calculates the thermal pressure for the BM-Stixrude regime."

    diff = ((T**4.0) * (integrate(th, 0.0, (theta/T), 600))) - ((T_0**4.0) * (integrate(th, 0.0, (theta/T_0), 600)))
    Pth = (((9.0 * gamma * n0 * R) / (V * (theta**3.0))) * diff)
    #print "Therm", Pth

    return Pth / 1000.0

    
def derivs(n, x, y):
#"The function DERIVS calculates y' from x and y"
    dy=[0 for i in range(0,n+1)]

    if y[1] > 120.0:
        K_0 = [252.58, 230.6]  
        K_0_prime = [4.08, 4.0]
        rho_0 = [4256.95, 4263.02]  
        V_0 = [24.43, 24.49]    
        theta_0 = [887.96, 838.5]  
        gamma_0 = [1.584, 1.887]
        beta = 1.1   

    if y[1] < 120.0:     
        K_0 = [252.58, 246.252]  
        K_0_prime = [4.08, 4.078]
        rho_0 = [4256.95, 4258.67]  
        V_0 = [24.43, 24.51]    
        theta_0 = [887.96, 877.0]  
        gamma_0 = [1.584, 1.630]
        beta = 1.1 

    if y[1] < 90.00:     
        K_0 = [252.58, 254.668]  
        K_0_prime = [4.08, 4.099]
        rho_0 = [4256.95, 4279.37]  
        V_0 = [24.43, 24.52]    
        theta_0 = [887.96, 896.682]  
        gamma_0 = [1.584, 1.569]
        beta = 1.1  

    s2 = optimize.brentq(f_BM_Stixrude, 0.01, 10.0, args=(K_0[1], K_0_prime[1], theta_0[1], gamma_0[1], y[1], y[4]))
    V = s2 * V_0[1]
    rho = (1 / s2) * rho_0[1] 
    gamma_m = gamma_0[1] * (s2**beta)

    f = (1.0 / 2.0) * ((s2**(-2.0 / 3.0)) - 1.0)

    theta_x = (theta_0[1]**2.0) * (1 + (6.0 * gamma_0[1] * f) + (((-6.0 * gamma_0[1]) + (18.0 * (gamma_0[1]**2.0)) - (9.0 * beta * gamma_0[1]))* (f**2.0)))
    theta = (theta_x)**0.5

    tx = p_th(s2, K_0[1], K_0_prime[1], theta_0[1], gamma_0[1], y[4])
    ty = p_th(s2, K_0[1], K_0_prime[1], theta_0[1], gamma_0[1], (y[4] + 0.01))
    dydx_2 = (ty - tx) / 0.01
    Cv = dydx_2 * (V / gamma_m)

    fg = lambda xs: f_BM_Stixrude(xs, K_0[1], K_0_prime[1], theta_0[1], gamma_0[1], 0.0, y[4])
    dydx = differentiate(fg, s2, var='volume')
    K_T = dydx * (-1.0) * (s2 * V_0[1])
    alpha = (gamma_m * Cv) / (K_T * V)	
    K_s = K_T * (1.0 + (alpha * gamma_m * y[4]))

    print("density is", rho)
    density.append(rho)

    dy[1] = (-1.0) * (10**(-6.0))*(rho) * (y[3])                             #dy[1] is dP
    dy[2] = (10**9.0) * 4.0 * pi * (x**2) * (rho)                            #dy[2] is dm
    dy[3] = (10**3.0) * (4.0 * pi * G * (rho)) - ((2.0 * y[3])/ x )          #dy[3] is dg
    dy[4] = (-1.0) * (10**(-6.0)) * ((gamma_m * (rho) * y[3] * y[4]) / K_s)  #dy[4] is dT
    dy[5] = ((rho * epsilon) - ((2.0 * y[5]) / (x * 1.000)))                 #dy[5] is dq

    return dy

#----------------------------------------------------------------------------

N=10  # Step size

#x=r
#m=y[2] mass at centre is zero 
# Radius r is in Km (IMPORTANT)

mass = []
radius = []
pressure = []
gravity = []
Temperature = []
                                   
x=r_ini; y=[0.0, p_ini, m_ini, g_ini, T_ini, q_ini]   # Sets Boundary Conditions

P_term = ((y[4] - 800.0) * (-0.0017)) + 25

while y[1] > P_term:
        (x,y) = runkut(5, x, y, 1.0/N)
 
        print("mass is", y[2], "Kg") 
        print("radius is", x, "Km")
        print("pressure is", y[1], "GPa")
        print("gravity is", y[3])
        print("Temperature is", y[4], "Kelvin")
        print("Q is", (y[5] * 100))
                       
        mass.append(y[2])
        radius.append(x)
        pressure.append(y[1])
        gravity.append(y[3])
        Temperature.append(y[4])


import csv
with open('Proxima_lower_mantle_1B.txt', 'w+') as x:
    
    writer = csv.writer(x, delimiter='\t')
    writer.writerows(zip(pressure, Temperature, radius))

with open('Proxima_lower_mantle_2B.txt', 'w+') as x:
    
    writer = csv.writer(x, delimiter='\t')
    writer.writerows(zip(mass, gravity, radius))

with open('Proxima_lower_mantle_3B.txt', 'w+') as x:
    
    writer = csv.writer(x, delimiter='\t')
    writer.writerows(zip(density))
        

plt.plot(radius, pressure, color='blue')
plt.xlabel('Radius (Km)')
plt.ylabel('Pressure (GPa)')
plt.show()

#plt.plot(radius, density, color='green')
#plt.xlabel('Radius (Km)')
#plt.ylabel('Density (Kg/m$^3$)')
#plt.show()

plt.plot(radius, gravity, color='red')
plt.xlabel('Radius (Km)')
plt.ylabel('Gravity (m/s$^2$)')
plt.show()
    

print("Done")

