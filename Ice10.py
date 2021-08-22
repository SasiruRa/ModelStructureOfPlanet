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

K_0 = [14.05, 23.9, 162.8]                # GPa from 3BM for Ice-VI and Ice-VII, Ice-X
K_0_prime = [4.0, 4.2, 4.42]              # for from 3BM for Ice-VI and Ice-VII, Ice-X
rho_0 = [1271.4, 1460.0, 2318.6]          # Kg/m^3 from 3BM for Ice-VI and Ice-VII, Ice-X
V_0 = [14.17, 12.30, 7.77]                # cm^3/mol from 3BM for Ice-VI and Ice-VII, Ice-X
a_0 = (-1.0) * 4.2 * (10**-4.0)
a_1 = 1.58 * (10**-6.0)
b = [0.9, 1.1]                            # for 3BM, Vinet
q = 1.0

gamma_0 = 1.2
theta_0 = 1470
n0 = 1.0
R = 8.314

p = 120 # Initial pressure
T_i = 686.719
T_0 = 300.0
g_i = 15.38
m_ini = 1.23426219578e+25
r_ini = 6500


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


th = lambda t: (t**3.0) / (float(exp(t)+0.0000001) - 1.0)
c_v = lambda z: ((z**4.0) * exp(z)) / (((exp(z)+0.0000001) - 1.0)**2.0)


def f_BM(x, K_0, p):

    a = (3.0/2.0) * K_0
    b = (x**(7.0 / 3.0)) - (x**(5.0 / 3.0))

    return (a * b) - p


def f_3BM(x, K_0, K_0_prime, p):                               # 3rd order Birch-Murnaghan eqn. Calculated with x = rho / rho_0

    a = (3.0/2.0) * K_0
    b = (x**(7.0 / 3.0)) - (x**(5.0 / 3.0))
    c = 1.0 + (((3.0 / 4.0) * (K_0_prime - 4.0)) * (x**(2.0 / 3.0) - 1.0))

    return (a * b * c) - p


def f_3BM_Vol(x, K_0, K_0_prime, p, T):   

    f = (1.0 / 2.0) * ((x**(-2.0 / 3.0)) - 1.0)                # 3rd order Birch-Murnaghan eqn. Calculated with x = V / V_0
    gamma = gamma_0 * (x**q)
    theta = theta_0 * (x**(-1.0 * gamma))

    ta = (3.0 * K_0 * f) * ((1.0 + (2.0 * f))**(5.0 / 2.0))
    tb = (1.0 + ((3.0 / 2.0) * (K_0_prime - 4.0) * f))

    V = x * V_0[2]

    "Calculates the thermal pressure for the BM-Stixrude regime."

    diff = ((T**4.0) * (integrate(th, 0.0, (theta/T), 600))) - ((T_0**4.0) * (integrate(th, 0.0, (theta/T_0), 600)))
    Pth = (((9.0 * gamma * n0 * R) / (V * (theta**3.0))) * diff)

    return (ta * tb) + (Pth / 1000.0) - p



def derivs(n, x, y):
#"The function DERIVS calculates y' from x and y"
    dy=[0 for i in range(0,n+1)]
    alpha_0 = lambda T: (a_0 + (a_1 * T)) * (1.0 + ((K_0_prime[2] / K_0[2]) * y[1]))**(-1.0 * b[0])
    alpha = (a_0 + (a_1 * y[4])) * (1.0 + ((K_0_prime[2] / K_0[2]) * y[1]))**(-1.0 * b[0])
    s = optimize.brentq(f_3BM, 0.01, 100000.0, args=(K_0[2], K_0_prime[2], y[1]))
    rho_p = s * rho_0[2] 
    rho = rho_p * exp(integrate(alpha_0, T_0, y[4], 600))
    Vp = (1 / s) * V_0[2]
    V = Vp * exp(integrate(alpha_0, T_0, y[4], 600))

    gamma = gamma_0 * (s**(-1.0 * q))
    theta = theta_0 * (s**(1.0 * gamma))

    fg = lambda xs: f_3BM_Vol(xs, K_0[2], K_0_prime[2], 0.0, y[4])
    dydx = differentiate(fg, (1 / s), var='volume')
    K_T = dydx * (-1.0) * (V)	
    K_s = K_T * (1.0 + (alpha * gamma * y[4]))

    print("density is", rho)
    density.append(rho)

    dy[1] = (-1.0) * (10**(-6.0))*(rho) * (y[3])                            #dy[1] is dP
    dy[2] = (10**9.0) * 4.0 * pi * (x**2) * (rho)                           #dy[2] is dm
    dy[3] = (10**3.0) * (4.0 * pi * G * (rho)) - ((2.0 * y[3])/ x )         #dy[3] is dg
    dy[4] = (-1.0) * (10**(-6.0)) * ((gamma * (rho) * y[3] * y[4]) / K_s) #dy[4] is dT

    return dy

#----------------------------------------------------------------------------

N=10  # Step size

"To get Isobaric heat capacity Cp, multiply the output of (alpha * K_T * V / gamma) by 1000 * (1000 / 18) where 18=molar mass of H20-ice"

#x=r
#m=y[2] mass at centre is zero 
# Radius r is in Km (IMPORTANT)

mass = []
radius = []
pressure = []
gravity = []
Temperature = []
                                   
x=r_ini; y=[0.0, p, m_ini, g_i, T_i]   # Sets Boundary Conditions

while y[1] > 62.00:
    (x,y) = runkut(4, x, y, 1.0/N)
 
    print("mass is", y[2], "Kg")
    print("radius is", x, "Km")
    print("pressure is", y[1], "GPa")
    print("gravity is", y[3])
    print("Temperature is", y[4])
                       
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

#dy[4] = (-1.0) * (10**(-6.0)) * ((gamma * (rho) * y[3] * y[4]) / K_s)  #dy[4] is dT
#rho = rho_p * exp((-1.0) * alpha * (y[4] - T_0))
