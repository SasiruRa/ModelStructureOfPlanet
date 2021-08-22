import matplotlib.pyplot as plt
import matplotlib
import numpy as np

file1 = open('Proxima_Core_1.txt')
file2 = open('Proxima_Core_2.txt')
file3 = open('Proxima_lower_mantle_1.txt')
file4 = open('Proxima_lower_mantle_2.txt')
file5 = open('Proxima_conductive_mantle_1.txt')
file6 = open('Proxima_conductive_mantle_2.txt')
file7 = open('Proxima_upper_mantle_1.txt')
file8 = open('Proxima_upper_mantle_2.txt')
file9 = open('Proxima_conductive_Upmantle_1.txt')
file10 = open('Proxima_conductive_Upmantle_2.txt')
file11 = open('Proxima_conductive_crust_1.txt')
file12 = open('Proxima_conductive_crust_2.txt')

file13 = open('Proxima_lower_mantle_1B.txt')
file14 = open('Proxima_lower_mantle_2B.txt')

file15 = open('Proxima_Core_3.txt')
file16 = open('Proxima_conductive_mantle_3.txt')
file17 = open('Proxima_lower_mantle_3.txt')
file18 = open('Proxima_lower_mantle_3B.txt')
file19 = open('Proxima_upper_mantle_3.txt')
file20 = open('Proxima_conductive_Upmantle_3.txt')
file21 = open('Proxima_conductive_crust_3.txt')

lst1 = []
lst2 = []
lst3 = []
lst4 = []
lst5 = []
lst6 = []
lst7 = []
lst8 = []
lst9 = []
lst10 = []
lst11 = []
lst12 = []

lst13 = []
lst14 = []

lst15 = []
lst16 = []
lst17 = []
lst18 = []
lst19 = []
lst20 = []
lst21 = []


for line in file1:
    lst1 += [line.split()]

for line2 in file2:
    lst2 += [line2.split()]

for line3 in file3:
    lst3 += [line3.split()]

for line4 in file4:
    lst4 += [line4.split()]

for line5 in file5:
    lst5 += [line5.split()]

for line6 in file6:
    lst6 += [line6.split()]

for line7 in file7:
    lst7 += [line7.split()]

for line8 in file8:
    lst8 += [line8.split()]

for line9 in file9:
    lst9 += [line9.split()]

for line10 in file10:
    lst10 += [line10.split()]

for line11 in file11:
    lst11 += [line11.split()]

for line12 in file12:
    lst12 += [line12.split()]

for line13 in file13:
    lst13 += [line13.split()]

for line14 in file14:
    lst14 += [line14.split()]

for line15 in file15:
    lst15 += [line15.split()]

for line16 in file16:
    lst16 += [line16.split()]

for line17 in file17:
    lst17 += [line17.split()]

for line18 in file18:
    lst18 += [line18.split()]

for line19 in file19:
    lst19 += [line19.split()]

for line20 in file20:
    lst20 += [line20.split()]

for line21 in file21:
    lst21 += [line21.split()]






R1 = []
T1 = []
P1 = []
g1 = []
R2 = []
T2 = []
P2 = []
g2 = []

R4 = []
T4 = []
P4 = []
g4 = []

R3 = []
T3 = []
P3 = []
g3 = []

R5 = []
T5 = []
P5 = []
g5 = []
R6 = []
T6 = []
P6 = []
g6 = []

R7 = []
T7 = []
P7 = []
g7 = []

d1 = []; d2 = []; d3 = []; d4 = []; d5 = []; d6 = []; d7 = []

for i in range(0, len(lst1), 1):
    P1.append(float(lst1[i][0]))
    T1.append(float(lst1[i][1]))
    R1.append(float(lst1[i][2]))

for i in range(0, len(lst2), 1):
    g1.append(float(lst2[i][1]))


for i in range(0, len(lst3), 1):
    P3.append(float(lst3[i][0]))
    T3.append(float(lst3[i][1]))
    R3.append(float(lst3[i][2]))

for i in range(0, len(lst4), 1):
    g3.append(float(lst4[i][1]))


for i in range(0, len(lst5), 1):
    P2.append(float(lst5[i][0]))
    T2.append(float(lst5[i][1]))
    R2.append(float(lst5[i][2]))

for i in range(0, len(lst6), 1):
    g2.append(float(lst6[i][1]))

for i in range(0, len(lst7), 1):
    P4.append(float(lst7[i][0]))
    T4.append(float(lst7[i][1]))
    R4.append(float(lst7[i][2]))

for i in range(0, len(lst8), 1):
    g4.append(float(lst8[i][1]))

for i in range(0, len(lst9), 1):
    P5.append(float(lst9[i][0]))
    T5.append(float(lst9[i][1]))
    R5.append(float(lst9[i][2]))

for i in range(0, len(lst10), 1):
    g5.append(float(lst10[i][1]))

for i in range(0, len(lst11), 1):
    P6.append(float(lst11[i][0]))
    T6.append(float(lst11[i][1]))
    R6.append(float(lst11[i][2]))

for i in range(0, len(lst12), 1):
    g6.append(float(lst12[i][1]))

for i in range(0, len(lst13), 1):
    P7.append(float(lst13[i][0]))
    T7.append(float(lst13[i][1]))
    R7.append(float(lst13[i][2]))

for i in range(0, len(lst14), 1):
    g7.append(float(lst14[i][1]))

for i in range(0, len(lst15), 4):
    d1.append(float(lst15[i][0]))

for i in range(0, len(lst16), 4):
    d2.append(float(lst16[i][0]))

for i in range(0, len(lst17), 4):
    d3.append(float(lst17[i][0]))

for i in range(0, len(lst18), 4):
    d4.append(float(lst18[i][0]))

for i in range(0, len(lst19), 4):
    d5.append(float(lst19[i][0]))

for i in range(0, len(lst20), 4):
    d6.append(float(lst20[i][0]))

for i in range(0, len(lst21), 4):
    d7.append(float(lst21[i][0]))



R1 = np.array(R1)
P1 = np.array(P1)
T1 = np.array(T1)
g1 = np.array(g1)
R2 = np.array(R2)
P2 = np.array(P2)
T2 = np.array(T2)
g2 = np.array(g2)

Rt = np.concatenate((R1, R2, R3, R7, R4, R5, R6), axis=0)
Pt = np.concatenate((P1, P2, P3, P7, P4, P5, P6), axis=0)
Tt = np.concatenate((T1, T2, T3, T7, T4, T5, T6), axis=0)
gt = np.concatenate((g1, g2, g3, g7, g4, g5, g6), axis=0)
dt = np.concatenate((d1, d2, d3, d4, d5, d6, d7), axis=0)


import csv
with open('Density.txt', 'w+') as x:
    
    writer = csv.writer(x, delimiter='\t')
    writer.writerows(zip(Rt, dt))

with open('Temp.txt', 'w+') as x:
    
    writer = csv.writer(x, delimiter='\t')
    writer.writerows(zip(Rt, Tt))


matplotlib.rcParams.update({'font.size': 20})

plt.plot(Rt, Pt, color='blue', ms=8.0)
ax = plt.gca()
plt.xlabel('Radius (Km)')
plt.ylabel('Pressure (GPa)')
plt.text(7000, 1200, 'Herath et al. 2019')
plt.legend()
plt.show()

plt.plot(Rt, gt, color='green', ms=8.0)
ax = plt.gca()
plt.xlabel('Radius (Km)')
plt.ylabel('Gravity (m/s$^2$)')
plt.text(7000, 12, 'Herath et al. 2019')
plt.legend()
plt.show()

plt.plot(Rt, Tt, color='red', ms=8.0)
ax = plt.gca()
plt.xlabel('Radius (Km)')
plt.ylabel('Temperature (Kelvins)')
plt.text(7000, 6000, 'Herath et al. 2019')
plt.legend()
plt.show()

plt.plot(Rt, dt, color='orange', ms=8.0)
ax = plt.gca()
plt.xlabel('Radius (Km)')
plt.ylabel('Density (Kg/m$^3$)')
plt.text(7000, 6000, 'Herath et al. 2019')
plt.legend()
plt.show()


print (len(lst15), len(lst16), len(lst17), len(lst18), len(lst19), len(lst20), len(lst21))
    
