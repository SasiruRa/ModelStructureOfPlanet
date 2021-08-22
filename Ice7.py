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

K_0 = [20.15, 23.9, 2.308]                # GPa for BM, 3BM, Vinet
K_0_prime = [4.0, 4.2, 4.532]             # for BM, 3BM, Vinet
rho_0 = [1442.4, 1460.0, 1463.0]          # Kg/m^3 for BM, 3BM, Vinet
V_0 = [12.49, 12.30]    # cm^3/mol
a_0 = (-1.0) * 4.2 * (10**-4.0)
a_1 = 1.58 * (10**-6.0)
b = [0.9, 1.1]                            # for 3BM, Vinet
q = 1.0

gamma_0 = 1.2
theta_0 = 1470
n0 = 1.0
R = 8.314
epsilon = 7.38 * (10**(-11))

p = 21.59604994419263  # Initial pressure
T_i = 641.1206673108625 # Initial temperature
T_0 = 300.0 #does not change ??
g_i = 9.935932325692756
m_ini = 3.4363650773737544e+24
r_ini = 4804.5000009999485
q_ini = 0.0008742777006077136

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
    V = x * V_0[1]

    "Calculates the thermal pressure for the BM-Stixrude regime."

    diff = ((T**4.0) * (integrate(th, 0.0, (theta/T), 600))) - ((T_0**4.0) * (integrate(th, 0.0, (theta/T_0), 600)))
    Pth = (((9.0 * gamma * n0 * R) / (V * (theta**3.0))) * diff)

    return (ta * tb) + (Pth / 1000.0) - p


def derivs(n, x, y):
#"The function DERIVS calculates y' from x and y"
    dy=[0 for i in range(0,n+1)]

    alpha_0 = lambda T: (a_0 + (a_1 * T)) * (1.0 + ((K_0_prime[1] / K_0[1]) * y[1]))**(-1.0 * b[0])
    alpha = (a_0 + (a_1 * y[4])) * (1.0 + ((K_0_prime[1] / K_0[1]) * y[1]))**(-1.0 * b[0])
    s = optimize.brentq(f_3BM, 0.01, 100000.0, args=(K_0[1], K_0_prime[1], y[1]))
    rho_p = s * rho_0[1] 
    rho = rho_p * exp(integrate(alpha_0, T_0, y[4], 600))
    Vp = (1 / s) * V_0[1]
    V = Vp * exp(integrate(alpha_0, T_0, y[4], 600))

    gamma = gamma_0 * (s**(-1.0 * q))
    theta = theta_0 * (s**(1.0 * gamma))

    fg = lambda xs: f_3BM_Vol(xs, K_0[1], K_0_prime[1], 0.0, y[4])
    dydx = differentiate(fg, (1 / s), var='volume')
    K_T = dydx * (-1.0) * (V)	
    K_s = K_T * (1.0 + (alpha * gamma * y[4]))

    print("density is", rho)
    density.append(rho)

    dy[1] = (-1.0) * (10**(-6.0))*(rho) * (y[3])                            #dy[1] is dP
    dy[2] = (10**9.0) * 4.0 * pi * (x**2) * (rho)                           #dy[2] is dm
    dy[3] = (10**3.0) * (4.0 * pi * G * (rho)) - ((2.0 * y[3])/ x )         #dy[3] is dg
    dy[4] = (-1.0) * (10**(-6.0)) * ((gamma * (rho) * y[3] * y[4]) / K_s)   #dy[4] is dT
    dy[5] = ((rho * epsilon) - ((2.0 * y[5]) / (x * 1.000)))                #dy[5] is dq

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
                                   
x=r_ini; y=[0.0, p, m_ini, g_i, T_i, q_ini]   # Sets Boundary Conditions

theta_t = y[4] / 355.0
num_1 = 0.173683 * 10.0 * (1.0 - (theta_t**(-1.0)))
num_2 = 0.544606 * (10.0**(-1.0)) * (1.0 - (theta_t**(5.0)))
num_3 = 0.806106 * (10.0**(-7.0)) * (1.0 - (theta_t**(22.0)))
p_melt = 2.216 * (exp(num_1 - num_2 + num_3))


while y[1] > p_melt:
    (x,y) = runkut(5, x, y, 1.0/N)
 
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
with open('Proxima_upper_mantle_1.txt', 'w+') as x:
    
    writer = csv.writer(x, delimiter='\t')
    writer.writerows(zip(pressure, Temperature, radius))

with open('Proxima_upper_mantle_2.txt', 'w+') as x:
    
    writer = csv.writer(x, delimiter='\t')
    writer.writerows(zip(mass, gravity, radius))

with open('Proxima_upper_mantle_3.txt', 'w+') as x:
    
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
