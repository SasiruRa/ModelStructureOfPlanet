import pylab
from math import *
from scipy.constants import G
from scipy import optimize
import matplotlib.pyplot as plt
import numpy as np


def runkut(n, x, y, h):
	#Advances the solution of diff eqn defined by derivs from x to x+h
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

# B-M-Stixrude parameters for Mg-post-perovskite (1st-lower mantle) and Mg-perovskite (2nd-upper mantle)

K_0 = [230.6, 251]  # GPa
K_0_prime = [4.0, 4.1]
rho_0 = [4263.02, 4150.20]  # Kg/m^3
V_0 = [24.49, 24.45]    # cm^3/mol
theta_0 = [838.5, 905.0]  # (K)
gamma_0 = [1.887, 1.570]
#gamma_inf = 0.827
beta = 1.1
epsilon = 7.38 * (10**(-11))


p_ini = 180.50214631875014  
T_ini = 3116.3729304230701
m_ini = 1.194476183249915e+24
density = 12236.107440693842
g_ini = 10.22785774970102
r_ini = 2791.8000009985567
q_ini = 0.0009000266999994488

T_0 = 300.0 # Reference temperature in kelvins
n0 = 1.0
R = 8.314
kappa = 12.00

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
        val = (a * V_0[0])
        t1 = f((val + h) / V_0[0])
        t2 = f((val - h) / V_0[0])
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
    theta = (theta_x)**0.5

    ta = (3.0 * K_0 * f) * ((1.0 + (2.0 * f))**(5.0 / 2.0))
    tb = (1.0 + ((3.0 / 2.0) * (K_0_prime - 4.0) * f))

    V = x * V_0[0]

    #"Calculates the thermal pressure for the BM-Stixrude regime."

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

    V = x * V_0[0]

    #"Calculates the thermal pressure for the bulk modulus."

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

    V = x * V_0[0]

    #"Calculates the thermal pressure for the BM-Stixrude regime."

    diff = ((T**4.0) * (integrate(th, 0.0, (theta/T), 600))) - ((T_0**4.0) * (integrate(th, 0.0, (theta/T_0), 600)))
    Pth = (((9.0 * gamma * n0 * R) / (V * (theta**3.0))) * diff)
    #print "Therm", Pth

    return Pth / 1000.0

    
def derivs(n, x, y):
	#"The function DERIVS calculates y' from x and y"
    dy=[0 for i in range(0,n+1)]
    s2 = optimize.brentq(f_BM_Stixrude, 0.01, 10000.0, args=(K_0[0], K_0_prime[0], theta_0[0], gamma_0[0], y[1], y[4]))
    V = s2 * V_0[0]
    rho = (1 / s2) * rho_0[0] 
    gamma_m = gamma_0[0] * (s2**beta)

    f = (1.0 / 2.0) * ((s2**(-2.0 / 3.0)) - 1.0)

    theta_x = (theta_0[0]**2.0) * (1 + (6.0 * gamma_0[0] * f) + (((-6.0 * gamma_0[0]) + (18.0 * (gamma_0[0]**2.0)) - (9.0 * beta * gamma_0[0]))* (f**2.0)))
    theta = (theta_x)**0.5

    tx = p_th(s2, K_0[0], K_0_prime[0], theta_0[0], gamma_0[0], y[4])
    ty = p_th(s2, K_0[0], K_0_prime[0], theta_0[0], gamma_0[0], (y[4] + 0.01))
    dydx_2 = (ty - tx) / 0.01
    Cv = dydx_2 * (V / gamma_m)

    fg = lambda xs: f_BM_Stixrude(xs, K_0[0], K_0_prime[0], theta_0[0], gamma_0[0], 0.0, y[4])
    dydx = differentiate(fg, s2, var='volume')
    K_T = dydx * (-1.0) * (s2 * V_0[0])
    alpha = (gamma_m * Cv) / (K_T * V)	
    K_s = K_T * (1.0 + (alpha * gamma_m * y[4]))

    print("density is", rho)
    density.append(rho)

    dy[1] = (-1.0) * (10**(-6.0)) * (rho) * (y[3])                           #dy[1] is dP
    dy[2] = (10**9.0) * 4.0 * pi * (x**2) * (rho)                            #dy[2] is dm
    dy[3] = (10**3.0) * (4.0 * pi * G * (rho)) - ((2.0 * y[3])/ x )          #dy[3] is dg
    dy[4] = (-1.0) * (10**(5.0)) * (y[5] / kappa)                            #dy[4] is dT
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

while x < r_ini+150:
        (x,y) = runkut(5, x, y, 1.0/N)
 
        print("mass is", y[2], "Kg")
        print("radius is", x, "Km")
        print("pressure is", y[1], "GPa")
        print("gravity is", y[3])
        print("Temperature is", y[4], "Kelvin")
        print("Q is", y[5])
                       
        mass.append(y[2])
        radius.append(x)
        pressure.append(y[1])
        gravity.append(y[3])
        Temperature.append(y[4])


import csv
with open('Proxima_conductive_mantle_1.txt', 'w+') as x:
    
    writer = csv.writer(x, delimiter='\t')
    writer.writerows(zip(pressure, Temperature, radius))

with open('Proxima_conductive_mantle_2.txt', 'w+') as x:
    
    writer = csv.writer(x, delimiter='\t')
    writer.writerows(zip(mass, gravity, radius))

with open('Proxima_conductive_mantle_3.txt', 'w+') as x:
    
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

print("CMB is Done")