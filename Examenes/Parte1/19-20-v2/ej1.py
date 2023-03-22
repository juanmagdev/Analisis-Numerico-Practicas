# from pylab import *
from time import perf_counter
from matplotlib.pyplot import *
from numpy import *

ep =1

def f(t, y):
    f1 = - (y[0]**3 - 3*y[0] + y[1])*(ep)**(-1)
    f2 = y[0]
    return(array([f1,f2]))


def rk4Sis(a, b, fun, N, y0):
    """Implementacion del metodo de rk4 para sistemas en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""

    h = (b-a)/N  # paso de malla
    t = zeros(N+1)  # inicializacion del vector de nodos
    y = zeros([len(y0), N+1])  # inicializacion del vector de resultados
    t[0] = a  # nodo inicial
    y[:, 0] = y0  # valor inicial
    k1 = zeros([len(y0), 1])
    k2 = zeros([len(y0), 1])
    k3 = zeros([len(y0), 1])
    k4 = zeros([len(y0), 1])

    # Metodo de rk4
    for k in range(N):
        k1 = fun(t[k], y[:, k])
        k2 = fun(t[k]+h/2, y[:, k]+h/2*k1)
        k3 = fun(t[k]+h/2, y[:, k]+h/2*k2)
        k4 = fun(t[k]+h, y[:, k]+h*k3)
        t[k+1] = t[k]+h
        y[:, k+1] = y[:, k] + h/6*(k1 + 2*k2 + 2*k3 + k4)

    return (t, y)


# Datos del problema
a = 0 # extremo inferior del intervalo
b = 30 # extremo superior del intervalo
y0 = array([1.7,0.3]) # condicion inicial
N = 150


tini = perf_counter()
(t1, y1) = rk4Sis(a, b, f, N, y0) 
tfin = perf_counter()

figure('Ejercicio 1a')
subplot(211)
plot(t1,y1[0],t1,y1[1]) 
legend(['(1)','(3)'])

subplot(212)
plot(y1[0],y1[1]) 
legend(['(1)','(3)'])

gcf().suptitle("Ejercicio 1a")
show()

# Apartado b
ep =0.4

def f(t, y):
    f1 = - (y[0]**3 - 3*y[0] + y[1])*(ep)**(-1)
    f2 = y[0]
    return(array([f1,f2]))

N = 150


tini = perf_counter()
(t1, y1) = rk4Sis(a, b, f, N, y0) 
tfin = perf_counter()

figure('Ejercicio 1b1')
subplot(211)
plot(t1,y1[0],t1,y1[1]) 
legend(['(1)','(3)'])

subplot(212)
plot(y1[0],y1[1]) 
legend(['(1)','(3)'])

gcf().suptitle("Ejercicio 1b1")
show()


ep =0.4

def f(t, y):
    f1 = - (y[0]**3 - 3*y[0] + y[1])*(ep)**(-1)
    f2 = y[0]
    return(array([f1,f2]))

N = 300


tini = perf_counter()
(t1, y1) = rk4Sis(a, b, f, N, y0) 
tfin = perf_counter()

figure('Ejercicio 1b2')
subplot(211)
plot(t1,y1[0],t1,y1[1]) 
legend(['(1)','(3)'])

subplot(212)
plot(y1[0],y1[1]) 
legend(['(1)','(3)'])

gcf().suptitle("Ejercicio 1b2")
show()