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

K_0 = [14.05, 23.9, 2.22]                # GPa from 3BM for Ice-VI and Ice-VII, liquid H2O
K_0_prime = [4.0, 4.2, 4.0]             # for from 3BM for Ice-VI and Ice-VII, liquid H2O
rho_0 = [1271.4, 1460.0, 998.0]           # Kg/m^3 from 3BM for Ice-VI and Ice-VII, liquid H2O
V_0 = [14.17, 12.30, 18.07]               # cm^3/mol from 3BM for Ice-VI and Ice-VII, liquid H2O
a_0 = (-1.0) * 4.2 * (10**-4.0)
a_1 = 1.58 * (10**-6.0)
b = [0.9, 1.1]                            # for 3BM, Vinet
q = 1.0

theta_0 = 1470
n0 = 1.0
R = 8.314
gamma_0 = 0.14890
gamma_1 = -0.26150
gamma_2 = -6.5703 * (10**(-2.0))
gamma_3 = -1.590 * (10**(-2.0))
gamma_4 = -1.1839 * (10**(-3.0))

k = [0, 1, 2, 3, 4]
a1 = 0.119539337 * (10**7.0)
a2 = 0.808183159 * (10**5.0)
a3 = 0.333826860 * (10**4.0)
b1 = 0.3 * (10**1.0)
b2 = 0.257500 * (10**2.0)
b3 = 0.103750 * (10**3.0)

p = 12.13 
T_i = 488.708
T_0 = 300.0
g_i = 12.80
m_ini = 1.85855226669e+25
r_ini = 9839.50

d_i = 1959.47


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


def f_3BM(x, K_0, K_0_prime, p):                               # 3rd order Birch-Murnaghan eqn. Calculated with x = rho / rho_0

    a = (3.0/2.0) * K_0
    b = (x**(7.0 / 3.0)) - (x**(5.0 / 3.0))
    c = 1.0 + (((3.0 / 4.0) * (K_0_prime - 4.0)) * (x**(2.0 / 3.0) - 1.0))

    return (a * b * c) - p


def f_3BM_Vol(x, K_0, K_0_prime, p):                 # 3rd order Birch-Murnaghan eqn. Calculated with x = V / V_0

    V = x * V_0[2]
    f = (1.0 / 2.0) * ((x**(-2.0 / 3.0)) - 1.0)                

    ta = (3.0 * K_0 * f) * ((1.0 + (2.0 * f))**(5.0 / 2.0))
    tb = (1.0 + ((3.0 / 2.0) * (K_0_prime - 4.0) * f))

    return (ta * tb) - p



def derivs(n, x, y):
	"The function DERIVS calculates y' from x and y"
	dy=[0 for i in range(0,n+1)]
 
        Cv = 58.20 / 1000
        s = optimize.brentq(f_3BM_Vol, 0.001, 1.0, args=(K_0[2], K_0_prime[2], y[1]))
        rho = (1 / s) * rho_0[2] 
        V = s * V_0[2]

        k1 = (gamma_0 * (V - V_0[2])**(k[0])) + (gamma_1 * (V - V_0[2])**(k[1]))
        k2 = (gamma_2 * (V - V_0[2])**(k[2])) + (gamma_3 * (V - V_0[2])**(k[3]))
        k3 = gamma_4 * (V - V_0[2])**(k[4])
        gamma = k1 + k2 + k3
        
        fg = lambda xs: f_3BM_Vol(xs, K_0[2], K_0_prime[2], 0.0)
        dydx = differentiate(fg, s, var='volume')
        K_T = dydx * (-1.0) * (V)
        alpha = (gamma * Cv) / (K_T * V)
	
        K_s = K_T * (1.0 + (alpha * gamma * y[4]))

        print "density is", rho
        density.append(rho)

	dy[1] = (-1.0) * (10**(-6.0))*(rho) * (y[3])                            #dy[1] is dP
	dy[2] = (10**9.0) * 4.0 * pi * (x**2) * (rho)                           #dy[2] is dm
        dy[3] = (10**3.0) * (4.0 * pi * G * (rho)) - ((2.0 * y[3])/ x )         #dy[3] is dg
        dy[4] = (-1.0) * (10**(-6.0)) * ((gamma * (rho) * y[3] * y[4]) / K_s)    #dy[4] is dT  FLIP sign if T < 300 K

	return dy

#----------------------------------------------------------------------------

N=100  # Step size

"To get Isobaric heat capacity Cp, multiply the output of (alpha * K_T * V / gamma) by 1000 * (1000 / 18) where 18=molar mass of H20-ice"

mass = []
radius = []
pressure = []
gravity = []
Temperature = []
                                   
x=r_ini; y=[0.0, p, m_ini, g_i, T_i]   # Sets Boundary Conditions

theta_T = y[4] / 273.160
x1 = a1 * (1.0 - (theta_T**b1))
x2 = a2 * (1.0 - (theta_T**b2))
x3 = a3 * (1.0 - (theta_T**b3))
P_freeze = (1.0 + x1 + x2 + x3) * (611.657 * (10**(-7.0)))
print "ice", P_freeze

while y[1] > 0.050:
        (x,y) = runkut(4, x, y, 1.0/N)
 
        print "mass is", y[2], "Kg" 
        print "radius is", x, "Km"
        print "pressure is", y[1], "GPa"
        print "gravity is", y[3]
        print "Temperature is", y[4]
                       
        mass.append(y[2])
        radius.append(x)
        pressure.append(y[1])
        gravity.append(y[3])
        Temperature.append(y[4])


import csv
with open('Proxima_conductive_Upmantle_1.txt', 'w+') as x:
    
    writer = csv.writer(x, delimiter='\t')
    writer.writerows(zip(pressure, Temperature, radius))

with open('Proxima_conductive_Upmantle_2.txt', 'w+') as x:
    
    writer = csv.writer(x, delimiter='\t')
    writer.writerows(zip(mass, gravity, radius))

with open('Proxima_conductive_Upmantle_3.txt', 'w+') as x:
    
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
