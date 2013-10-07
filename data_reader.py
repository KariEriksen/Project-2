#!/usr/bin/env python
from math import exp, log10 
import os
import numpy as np
from matplotlib.pyplot import *

# Reading file data.dat, 

n = 200
rho_max = 5
omega = 0.01

os.system("g++ Plot_wavefunction.cpp -larmadillo")
os.system("./a.out %d %d %d" % (n, omega, rho_max)) 

infile = open('data.dat', 'r')
#n = int(infile.readline())
#omega = float(infile.readline())
#rho = float(infile.readline())

vec1 = []
vec2 = []
vec3 = []
abs_vec1 = []
abs_vec2 = []
abs_vec3 = []

for line in infile:
    
    if line == ' ':
        break
    else:

        value = line.split()
        vec1.append(float(value[0]))
        vec2.append(float(value[1]))
        vec3.append(float(value[2]))
        absvec1 = (float(value[0]))**2
        abs_vec1.append(absvec1)
        absvec2 = (float(value[1]))**2
        abs_vec2.append(absvec2)
        absvec3 = (float(value[2]))**2
        abs_vec3.append(absvec3)

infile.close()

rho = np.linspace(0, n, n)
plot(rho, vec1)
hold("on")
plot(rho, vec2)
plot(rho, vec3)

title("Wavefunction with n = %d, omega = %d, $rho_{max}$ = %d" % (n, omega, rho_max))

legend(("Lowest state", "Second lowest", "Third lowest"), loc='upper right')

xlabel("Rho")
ylabel("Wavefunction")
show()
"""
figure(2)
plot(rho, abs_vec1)
hold("on")
plot(rho, abs_vec2)
plot(rho, abs_vec3)
show()
"""
