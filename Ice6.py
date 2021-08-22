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

K_0 = [14.05, 23.9]                # GPa from 3BM for Ice-VI and Ice-VII
K_0_prime = [4.0, 4.2]             # for from 3BM for Ice-VI and Ice-VII
rho_0 = [1271.4, 1460.0]           # Kg/m^3 from 3BM for Ice-VI and Ice-VII
V_0 = [14.17, 12.30]               # cm^3/mol from 3BM for Ice-VI and Ice-VII
a_0 = (-1.0) * 4.2 * (10**-4.0)
a_1 = 1.58 * (10**-6.0)
b = [0.9, 1.1]                            # for 3BM, Vinet for Ice-VII only
q = 1.0

gamma_0 = 1.2                             # These parameters are valid only for Ice-VII
theta_0 = 1470
n0 = 1.0
R = 8.314

p = 2.8819040374596   # Initial pressure at lowest level (r = r_ini)
T_i = 139.27603284258745 
T_0 = 300.0
g_i = 10.062185357162917
m_ini = 1.93231438063e+25
r_ini = 6427.700001005854
T_a = (-1.0) * 25  # Temperature in degrees centigrade


def integrate(f, a, b, n):

    h = float(b - a) / n
    result = (0.5 * f(a)) + (0.5 * f(b))

    for i in range(1, n):
        result += f(a + (i*h))
    result *= h

    return result


def differentiate(f, a, h = 0.01):
    
    t1 = f(a)
    t2 = f(a + h)

    return (t2 - t1) / h


th = lambda t: (t**3.0) / (float(exp(t)+0.0000001) - 1.0)


def f_3BM(x, K_0, K_0_prime, p):

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

    return (ta * tb) + (Pth / 1000.0)- p



def Volume_t(x, p, t):

    alpha_p = lambda T: (a_0 + (a_1 * T)) * (1.0 + ((K_0_prime[1] / K_0[1]) * p))**(-1.0 * b[0])
    Vp = x * V_0[1]
    V = Vp * exp(integrate(alpha_p, T_0, t, 600))
    return V


def V_IceVI(p, t):

    v_6 = (7.792 * (10**-4.0)) - ((4.884 * (10**-8.0)) * p) + ((2.488 * (10**-12.0)) * (p**2.0))
    f_6 = (t - T_a) * ((-3.9 * 10**(-4.0)) + ((7.85 * 10**(-7.0)) * (t + T_a)))
    g_6 = (1.0 + ((4.43 * 10**(-4.0)) * p))**(-0.55)

    return v_6 * exp(f_6 * g_6)


def derivs(n, x, y):
	#"The function DERIVS calculates y' from x and y"
	dy=[0 for i in range(0,n+1)]
	ref_p = 1.0
	a_vii = 72.49
	b_vii = 3.0240 / 1000.0
	c_vii = 1.442 * (10**6.0)

	alpha_0 = lambda T: (a_0 + (a_1 * T)) * (1.0 + ((K_0_prime[0] / K_0[0]) * y[1]))**(-1.0 * b[0])
	alpha = (a_0 + (a_1 * y[4])) * (1.0 + ((K_0_prime[0] / K_0[0]) * y[1]))**(-1.0 * b[0])
	s = optimize.brentq(f_3BM, 0.01, 100000.0, args=(K_0[0], K_0_prime[0], y[1]))
	rho_p = s * rho_0[0] 
	rho = rho_p * exp(integrate(alpha_0, T_0, y[4], 600))
	Vp = (1 / s) * V_0[0]
	V = Vp * exp(integrate(alpha_0, T_0, y[4], 600))

    #"Find second derivative of V: d^2V/dT^2 for Ice-VII" 
	def sec_deriv(g, prs, trs):

        h_x = 0.01
		tx = Volume_t(g, prs, (trs + h_x))
		ty = 2.0 * Volume_t(g, prs, trs)
		tz = Volume_t(g, prs, (trs - h_x))
		derivt = (tx - ty + tz) / (h_x**2.0)

	return derivt
         

	d2VdT2 = sec_deriv((1 / s), y[1], y[4])
        #"Integrate sec_deriv w.r.t P"
	sec_int = lambda pr: sec_deriv((1 / s), pr, y[4])
	integral = integrate(sec_int, ref_p, (y[1] * (10**4.0)), 600)

	Cp_0 = a_vii + (b_vii * y[4]) - (c_vii / (y[4]**2.0))             # Cp is in Joules per mole. Multiply by 55.56 for Joules per Kg.
	Cp_vii_a = Cp_0 - (y[4] * integral)
	Cp_vii = Cp_vii_a / 18
	print("Ice-VII Capacities", (Cp_0 / 18), Cp_vii)
		#"Finding derivatives of V6 and V7"
	deriv_VI_a = lambda temp_vi: V_IceVI((y[1] * 10), temp_vi)
	deriv_VI = differentiate(deriv_VI_a, (y[4] - 273))
	deriv_VII_a = lambda temp_vii: Volume_t((1 / s), (y[1] * (10**4.0)), temp_vii)
	deriv_VII = differentiate(deriv_VII_a, y[4])

	dPdT = 0.01250
	delta_V67 = V_IceVI((y[1] * 10), (y[4] - 273.0)) - Volume_t((1 / s), (y[1] * (10**4.0)), y[4])
	Cp0_VI = Cp_vii_a + (y[4] * (deriv_VI - deriv_VII) * dPdT) + (delta_V67 * dPdT)
	print("Kbar Capacity Ice-VI", (Cp0_VI / 18))
	"Find second derivative of V: d^2V/dT^2 for Ice-VI" 
	def sec_deriv_VI(prs_vi, trs_vi):

	    h_x = 0.01
	    tx = V_IceVI(prs_vi, (trs_vi + h_x))
	    ty = 2.0 * V_IceVI(prs_vi, trs_vi)
	    tz = V_IceVI(prs_vi, (trs_vi - h_x))
	    derivt_VI = (tx - ty + tz) / (h_x**2.0)

	    return derivt_VI
	P_kbar = 21.05 + (0.0125 * (y[4] - 273.0))
	#"Integrate sec_deriv_VI w.r.t P where P is in Kbars and T in centigrade"
	sec_int_VI = lambda prs: sec_deriv_VI(prs, (y[4] - 273.0))
	integral2 = integrate(sec_int_VI, P_kbar, (y[1] * (10**1.0)), 600)
	Cp_vi_a = Cp0_VI - (y[4] * integral2)
	Cp_vi = Cp_vi_a / 18
	           
	print("Total Capacity Ice-VI", Cp_vi)

	print("density is", rho)
	density.append(rho)

	dy[1] = (-1.0) * (10**(-6.0))*(rho) * (y[3])                           #dy[1] is dP
	dy[2] = (10**9.0) * 4.0 * pi * (x**2) * (rho)                          #dy[2] is dm
	dy[3] = (10**3.0) * (4.0 * pi * G * (rho)) - ((2.0 * y[3])/ x )        #dy[3] is dg
	dy[4] = (-1.0) * ((alpha * y[3] * y[4]) / Cp_vi)                       #dy[4] is dT

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
                                   
x=r_ini; y=[0.0, p, m_ini, g_i, T_i]   # Sets Boundary Conditions

theta_T = y[4] / 273.310
p_melt = 0.63240 * (1.0 - (1.07476 * (1.0 - (theta_T**4.6))))


while y[1] > p_melt:
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
    
print("Done")
