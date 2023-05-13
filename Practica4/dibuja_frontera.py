from numpy import *
from matplotlib.pyplot import *
    


def locfron(rho, sigma):
# Dibuja la frontera de la region de estabilidad absoluta 
# de un metodo multipaso.
# rho y sigma son los coeficientes de los polinomios caracteristicos
# ordenados de mayor a menor grado '''
    theta = arange(0, 2.*pi, 0.01)
    numer = polyval(rho, exp(theta*1j)) # rho(e^{theta*i})
    denom = polyval(sigma, exp(theta*1j)) # sigma(e^{theta*i})
    mu = numer/denom
    x = real(mu)
    y = imag(mu)
    plot(x, y)
    grid(True)
    axis('equal')

figure('AB')
# Ejemplo: AB3 y_{k+1} - y_k  = h/12*(23*f_k - 16*f_{k-1} + 5*f_{k-2})
rho = array([1., -1., 0.,0.]) # primero
sigma = array([0., 23., -16., 5.])/23. # segundo
locfron(rho,sigma)

#################  region estabilidad RK



def locfronRK(dR, N):
# Localizacion de la frontera de un metodo RK
#  Devidada de la funcion R
    Npoints = 5000
    T = 2*N*pi
    h = 2*N*pi/Npoints
    z = zeros(Npoints +1 , dtype = complex)
    z[0] = 0.
    t = 0
    for k in range(len(z)-1):
        z[k+1] = z[k]+ h*1j*exp(1j*t)/dR(z[k])
        t = t + h
    x = real(z)
    y = imag(z)
    plot(x,y)
    grid(True)
    axis('equal')

figure('RK explicitos')

# Euler: función de estabilidad R = 1 + z

def dREuler(z):
    return 1.
locfronRK(dREuler, 1.)

# RK2 explicitos: función de estabilidad R = 1 + z + z**2/2

def dRK2exp(z):
    return 1. + z
locfronRK(dRK2exp, 2.)



