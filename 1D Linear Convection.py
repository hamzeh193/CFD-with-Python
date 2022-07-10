import numpy  #Array lib
from matplotlib import pyplot  #plotting library
import time, sys     #some utilities

nx = 40
dx = 2 / (nx-1)
nt = 25
dt = .025
c = 1    #assume wavespeed of c = 1

#setting up the initial conditions

u = numpy.ones(nx)   # we define an array all equal to one with the nx size
u[int(0.5 / dx):int(1 / dx+1)] = 2  #setting u=2 between 0.5 and 1 for I.C

#pyplot.plot(numpy.linspace(0, 2, nx), u);

un = numpy.ones(nx) #initialize a temporary array

for n in range(nt):  #loop for values of n from 0 to nt, so it will run nt times
    un = u.copy() ##copy the existing values of u into un
    for i in range(1, nx): ## you can try commenting this line and...
    #for i in range(nx): ## ... uncommenting this line and see what happens!
        u[i] = un[i] - c * dt / dx * (un[i] - un[i-1])

pyplot.plot(numpy.linspace(0, 2, nx), u)
pyplot.show()
