import numpy
import sympy

from sympy import init_printing
init_printing(use_latex=True)

#Initial settup

x, nu, t = sympy.symbols('x nu t') # turn them to the symbolic variables
phi = (sympy.exp(-(x - 4*t)**2 / (4 * nu * (t+1))) + sympy.exp(-(x - 4 * t - 2 * sympy.pi)**2 / (4 * nu * (t+1))))

phiprime = phi.diff(x)

from sympy.utilities.lambdify import lambdify

u = -2 * nu * (phiprime / phi) + 4
ufunc = lambdify((t, x, nu), u)   # turn the symbolic variables to the scaler. this is a function

from matplotlib import pyplot
#%matplotlib inline this line is only for jupyter

###variable declaration
nx = 101
nt = 100
dx = 2 * numpy.pi / (nx - 1)
nu = 0.07
dt = dx * nu #why?

x = numpy.linspace(0, 2 * numpy.pi, nx)
un = numpy.empty(nx)
t = 0

u = numpy.asarray([ufunc(t, x0, nu) for x0 in x])

pyplot.figure(figsize=(10, 10), dpi=100) #resolution and height and width of the figure
pyplot.plot(x, u, marker='o', lw=2) #line width and the marker
pyplot.xlim([0, 2 * numpy.pi])
pyplot.ylim([0, 10]);


for n in range(nt):
    un = u.copy()
    for i in range(1, nx-1):
        u[i] = un[i] - un[i] * dt / dx *(un[i] - un[i-1]) + nu * dt / dx**2 *\
                (un[i+1] - 2 * un[i] + un[i-1])
    u[0] = un[0] - un[0] * dt / dx * (un[0] - un[-2]) + nu * dt / dx**2 *\
                (un[1] - 2 * un[0] + un[-2])
    u[-1] = u[0]

u_analytical = numpy.asarray([ufunc(nt * dt, xi, nu) for xi in x])
pyplot.figure(figsize=(11, 7), dpi=100)
pyplot.plot(x,u, marker='o', lw=2, label='Computational')
pyplot.plot(x, u_analytical, label='Analytical')
pyplot.xlim([0, 2 * numpy.pi])
pyplot.ylim([0, 10])
pyplot.legend();
pyplot.show()
