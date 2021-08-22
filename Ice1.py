import pylab
from math import *
from scipy.constants import G
from scipy import optimize
import matplotlib.pyplot as plt
import numpy as np

density = []


def runkut(n, x, y, h):
	"Advances the solution of diff eqn defined by derivs from x to x+h"
	y0=y[:]
	k1=derivs(n, x, y)
	for i in range(1,n+1): y[i]=y0[i]+0.5*h*k1[i]
	k2=derivs(n, x+0.5*h, y)
	for i in range(1,n+1): y[i]=y0[i]+h*(0.2071067811*k1[i]+0.2928932188*k2[i])
	k3=derivs(n, x+0.5*h, y)
	for i in range(1,n+1): y[i]=y0[i]-h*(0.7071067811*k2[i]-1.7071067811*k3[i])
	k4=derivs(n, x+h, y)
	for i in range(1,n+1):
		a=k1[i]+0.5857864376*k2[i]+3.4142135623*k3[i]+k4[i]
		y[i]=y0[i]+0.16666666667*h*a
	
	x+=h
	return (x,y)

#----------------------------------------------------------------------------

# Birch-Murnaghan EoS

K_0 = 10.9950                # GPa 
K_0_prime = 7.0            
rho_0 = 934.310           # Kg/m^3 
V_0 = [14.17, 12.30, 18.07]               # cm^3/mol 
q = 1.0
A = 3.348 * (10**(-8.0))
B = 1.559
C = -0.0040680
D = -2.051 * (10**(-5.0))
q_ini = 0.00083
theta_0 = 1470
n0 = 1.0
R = 8.314
epsilon = 7.38 * (10**(-11))
kappa = 3.480

p = 0.050 
T_i = 269.031
T_0 = 300.0
g_i = 12.29
m_ini = 1.98051879938e+25
r_ini = 10362.4

d_i = 1019.71


def integrate(f, a, b, n):

    h = float(b - a) / n
    result = (0.5 * f(a)) + (0.5 * f(b))

    for i in range(1, n):
        result += f(a + (i*h))
    result *= h

    return result


def differentiate(f, a, var='volume', h = 0.01):
    
    if var == 'volume':
        val = (a * V_0[2])
        t1 = f((val + h) / V_0[2])
        t2 = f((val - h) / V_0[2])
    elif var == 'temperature':
        t1 = f(a + h)
        t2 = f(a - h)

    return (t1 - t2) / (2.0 * h)

def f_M_ice(p, T):                 # Murnaghan eqn. 

    x = (-1.0) * (A / (B + 1.0)) * (T**(B + 1.0))
    rho_t = rho_0 * exp(x)
    K_T = K_0 + (C * T) + (D * (T**2.0))
    s = 1.0 + (p * (K_0_prime / K_T))
    y = 1 / K_0_prime
    rho_p = (rho_t / (s**y))

    return rho_p


def derivs(n, x, y):
	"The function DERIVS calculates y' from x and y"
	dy=[0 for i in range(0,n+1)]
 
        Cp = 1.2 #1000.00
        alpha = 4.0 * (10**(-5.0))
        rho = f_M_ice(y[1], y[4])	

        print "density is", rho
        density.append(rho)

	dy[1] = (-1.0) * (10**(-6.0))*(rho) * (y[3])                            #dy[1] is dP
	dy[2] = (10**9.0) * 4.0 * pi * (x**2) * (rho)                           #dy[2] is dm
        dy[3] = (10**3.0) * (4.0 * pi * G * (rho)) - ((2.0 * y[3])/ x )         #dy[3] is dg
        dy[4] = (-1.0) * (10**(5.0)) * (y[5] / kappa)                           #dy[4] is dT
        #dy[4] = (-1.0) * ((alpha * y[3] * y[4]) / Cp)                           #dy[4] is dT
        dy[5] = ((rho * epsilon) - ((2.0 * y[5]) / (x * 1.000)))                #dy[5] is dq

	return dy

#----------------------------------------------------------------------------

N=100  # Step size

"To get Isobaric heat capacity Cp, multiply the output of (alpha * K_T * V / gamma) by 1000 * (1000 / 18) where 18=molar mass of H20-ice"

#x=r
#m=y[2] mass at centre is zero 
# Radius r is in Km (IMPORTANT)

mass = []
radius = []
pressure = []
gravity = []
Temperature = []
                                   
x=r_ini; y=[0.0, p, m_ini, g_i, T_i, q_ini]   # Sets Boundary Conditions

while y[4] > 105.0:
        (x,y) = runkut(5, x, y, 1.0/N)
 
        print "mass is", y[2], "Kg" 
        print "radius is", x, "Km"
        print "pressure is", y[1], "GPa"
        print "gravity is", y[3]
        print "Temperature is", y[4]
        print "Q is", (y[5] * 100)
                       
        mass.append(y[2])
        radius.append(x)
        pressure.append(y[1])
        gravity.append(y[3])
        Temperature.append(y[4])


import csv
with open('Proxima_conductive_crust_1.txt', 'w+') as x:
    
    writer = csv.writer(x, delimiter='\t')
    writer.writerows(zip(pressure, Temperature, radius))

with open('Proxima_conductive_crust_2.txt', 'w+') as x:
    
    writer = csv.writer(x, delimiter='\t')
    writer.writerows(zip(mass, gravity, radius))

with open('Proxima_conductive_crust_3.txt', 'w+') as x:
    
    writer = csv.writer(x, delimiter='\t')
    writer.writerows(zip(density))
        

plt.plot(radius, pressure, color='blue')
plt.xlabel('Radius (Km)')
plt.ylabel('Pressure (GPa)')
plt.show()

plt.plot(radius, gravity, color='red')
plt.xlabel('Radius (Km)')
plt.ylabel('Gravity (m/s$^2$)')
plt.show()

plt.plot(pressure, Temperature, color='red')
plt.xlabel('Radius (Km)')
plt.ylabel('Temperature (K)')
plt.show()    

print "Done"

#dy[4] = (-1.0) * (10**(-6.0)) * ((gamma * (rho) * y[3] * y[4]) / K_s)  #dy[4] is dT
#rho = rho_p * exp((-1.0) * alpha * (y[4] - T_0))
