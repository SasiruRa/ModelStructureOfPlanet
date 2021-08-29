import pylab
from math import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import G
from scipy import optimize


# Bouchet EoS for Iron (Fe)

K_0 = 253.84  # GPa
K_0_prime = 4.719
rho_0 = 8878.4  # Kg/m^3
V_0 = 6.290    # cm^3/mol
z = 26.0  # Atomic number of iron
theta_0 = 44.7  # (K)
gamma_0 = 1.408
gamma_inf = 0.827
beta = 0.826

n0 = 6.0
a_0 = 2.121 * (10**-4.0)
m0 = 1.891
Kb = 1.38    #* 10**(-23) Boltzmann constant
Na = 6.022   #* 10**(23)  Avogadro constant
R = 8.314

epsilon = 7.38 * (10**(-11))

"----------------------Work horse for the 4th order Runge-Kutta algorithm------------------------------------------"


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

"-----------------------------------Differential function----------------------------------------------------------"

density = []

def differentiate(f, a, var='volume', h = 0.01):
    
    if var == 'volume':
        val = (a * V_0)
        t1 = f((val + h) / V_0)
        t2 = f((val - h) / V_0)
    elif var == 'temperature':
        t1 = f(a + h)
        t2 = f(a - h)

    return (t1 - t2) / (2.0 * h)


"--------------------------Calculate Bouchet Pressure w.r.t. density------------------------------------------------"
 
def f_Bouchet(x, n, m, K_0, K_0_prime, T, p):

    GPa = 10**9.0

    p_fg0 = 1003.6 * (z / V_0)**(5.0 / 3.0)
    c0 = (-1.0) * np.log((3.0 * K_0) / p_fg0)
    c2 = (3.0 / 2.0) * (K_0_prime - 3) - c0

    gamma = gamma_inf + ((gamma_0 - gamma_inf) * (x**(-beta)))
    theta = theta_0 * (x**gamma_inf) * exp(((1 - x**(-beta)) * (gamma_0 - gamma_inf)) / beta)

    try:
        lx = float(exp(theta / T))
        ax = theta / (lx - 1.0) 

        #print ax
    except OverflowError:
        ax = float('inf')
    
    p_har = (3.0 * gamma * n0 * Kb * Na) * (rho_0 * x) * ((theta / 2.0) + ax)
    p_ae = (3.0 / 2.0) * (a_0 * m0 * n0 * Kb * Na) * (rho_0 * x) * (x**(- m0)) * (T**2.0)
    px = p_har / GPa
    py = p_ae / GPa

    return (3.0 * K_0 * x**n) * (1 - x**(-m)) * (1 + c2*(x**(-m)) * (1 - x**(-m))) * exp(c0 * (1 - x**(-m))) + px + py - p


"--------------------------Calculate Bouchet Pressure w.r.t. volume-------------------------------------------------"

def f_Bouchet_vol(x, n, m, K_0, K_0_prime, p, T):

    p_fg0 = 1003.6 * (z / V_0)**(5.0 / 3.0)
    c0 = (-1.0) * np.log((3.0 * K_0) / p_fg0)
    c2 = (3.0 / 2.0) * (K_0_prime - 3) - c0

    gamma = gamma_inf + ((gamma_0 - gamma_inf) * (x**(beta)))
    theta = theta_0 * (x**(-gamma_inf)) * exp(((1 - x**(beta)) * (gamma_0 - gamma_inf)) / beta)

    an = float(exp(theta / T))
    if an == 1.0:
        ay = float(exp(theta / T)+0.00001)
    else:
        ay = float(exp(theta / T))

    try:
        ax = theta / (ay - 1.0)
        #print exp(theta / T)
    except OverflowError:
        ax = float('inf')

    p_har = ((3.0 * gamma * n0 * R) * ((theta / 2.0) + ax)) / (V_0 * x)
    p_ae = (3.0 / (2.0 * V_0 * x)) * (a_0 * m0) * (x**(m0)) * (T**2.0)
    px = p_har / 17925.178
    py = p_ae / 360.0
 
    return (3.0 * K_0 * x**(-n)) * (1 - x**(m)) * (1 + c2*(x**(m)) * (1 - x**(m))) * exp(c0 * (1 - x**(m))) + px + py - p


"--------------------------Calculate the Pth w.r.t. volume----------------------------------------------------------"

def p_th(x, T):
        

    gamma = gamma_inf + ((gamma_0 - gamma_inf) * (x**(beta)))
    theta = theta_0 * (x**(-gamma_inf)) * exp(((1 - x**(beta)) * (gamma_0 - gamma_inf)) / beta)

    an = float(exp(theta / T))
    if an == 1.0:
        ay = float(exp(theta / T)+0.00001)
    else:
        ay = float(exp(theta / T))

    try:
        ax = theta / (ay - 1.0)
    except OverflowError:
        ax = float('inf')

    p_har = ((3.0 * gamma * n0 * R) * ((theta / 2.0) + ax)) / (V_0 * x)
    p_ae = (3.0 / (2.0 * V_0 * x)) * (a_0 * m0) * (x**(m0)) * (T**2.0)
    px = p_har / 17925.178
    py = p_ae / 360.0

    return px + py


"--------------------------Solve the differential equations---------------------------------------------------------"

p = 450.0  # Initial pressure at the centre of the planet in GPa (Pc at r = 0)
T = 4500.0 # Initial temperature at the core in Kelvins (Tc at r = 0)
    
def derivs(n, x, y):
    "The function DERIVS calculates y from x and y"
    dy=[0 for i in range(0,n+1)]

    h_x=0.01
    s = optimize.brentq(f_Bouchet, 0.1, 100000000.0, args=(5.0/3.0, 1.0/3.0, K_0, K_0_prime, y[4], y[1]))
    rho = s * rho_0 
    gamma = gamma_inf + ((gamma_0 - gamma_inf) * (s**(-beta)))
    gamma_n = gamma_inf + ((gamma_0 - gamma_inf) * ((1.0 / s)**(beta)))

    f = lambda xs: f_Bouchet_vol(xs, 5.0/3.0, 1.0/3.0, K_0, K_0_prime, 0.0, y[4])
    dydx = differentiate(f, (1.0 / s), var='volume')
    K_T = dydx * (-1.0) * ((1.0 / s) * V_0)

    tx = p_th((1.0 / s), y[4])
    ty = p_th((1.0 / s), (y[4] + h_x))
    dydx_2 = (ty - tx) / h_x
    Cv = dydx_2 * (((1.0 / s) * V_0) / gamma_n)

    alpha = (gamma_n * Cv) / (K_T * ((1.0 / s) * V_0))
    K_s = K_T * (1.0 + (alpha * gamma_n * y[4]))
    print("density is", rho)
    density.append(rho)

    dy[1] = (-1.0) * (10**(-6.0))*(rho) * (y[3])                                  #dy[1] is dP
    dy[2] = (10**9.0) * 4.0 * pi * (x**2) * (rho)                                 #dy[2] is dm
    dy[3] = (10**3.0) * (4.0 * pi * G * (rho)) - ((2.0 * y[3])/ x )               #dy[3] is dg
    dy[4] = (-1.0) * (10**(-6.0)) * ((gamma_n * (rho) * y[3] * y[4]) / K_s)       #dy[4] is dT
    dy[5] = ((rho * epsilon) - ((2.0 * y[5]) / (x * 1.000)))                      #dy[5] is dq

    return dy

"----------------------------------End solving-----------------------------------------------------------------------"  

N=10  # Step size

#x=r
#m=y[2] mass at centre is zero 
# Radius r is in Km (IMPORTANT)

mass = []
radius = []
pressure = []
gravity = []
Temperature = []

q_i = 0.0000001                                   
x=0.000001; y=[0.0, p, 0.0, 0.0, T, q_i]   # Sets Boundary Conditions

ms=1.911187753123466e+24
#ms=(5.972e+24)*0.32
while y[2] < ms:
    (x,y) = runkut(5, x, y, 1.0/N)
         
    print("mass is", y[2], "Kg") 
    print("radius is", x, "Km")
    print("pressure is", y[1], "GPa")
    print("gravity is", y[3])
    print("Temperature is", y[4], "Kelvin")
    print("Q is", (y[5]))
                       
    mass.append(y[2])
    radius.append(x)
    pressure.append(y[1])
    gravity.append(y[3])
    Temperature.append(y[4])


import csv
with open('Proxima_Core_1.txt', 'w+') as x:
    
    writer = csv.writer(x, delimiter='\t')
    writer.writerows(zip(pressure, Temperature, radius))

with open('Proxima_Core_2.txt', 'w+') as x:
    
    writer = csv.writer(x, delimiter='\t')
    writer.writerows(zip(mass, gravity, radius))

with open('Proxima_Core_3.txt', 'w+') as x:
    
    writer = csv.writer(x, delimiter='\t')
    writer.writerows(zip(density))


        
plt.plot(radius, pressure, color='blue')
plt.xlabel('Radius (Km)')
plt.ylabel('Pressure (GPa)')
plt.show()

plt.plot(radius, gravity, color='green')
plt.xlabel('Radius (Km)')
plt.ylabel('Gravity (m/s$^2$)')
plt.show()

plt.plot(radius, Temperature, color='red')
plt.xlabel('Radius (Km)')
plt.ylabel('Temperature (K)')
plt.show()
    
print("Core is Done")
