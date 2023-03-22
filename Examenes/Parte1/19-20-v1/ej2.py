from numpy import *  
from matplotlib.pyplot import *
# from pylab import *
from time import perf_counter

### Ejercicio 2
def f3(t,y,mu):
    return((1-mu*sin(t)*y)*y)


def euler(a, b, fun, N, y0, mu):
    """Implementacion del metodo de Euler en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""
    
    h = (b-a)/N # paso de malla
    t = zeros(N+1) # inicializacion del vector de nodos
    y = zeros(N+1) # inicializacion del vector de resultados
    t[0] = a # nodo inicial
    y[0] = y0 # valor inicial

    # Metodo de Euler
    for k in range(N):
        y[k+1] = y[k]+h*fun(t[k], y[k],mu)
        t[k+1] = t[k]+h
    
    return (t, y)

# Datos del problema
a = 0. # extremo inferior del intervalo
b = 2 # extremo superior del intervalo
N = 20 # numero de particiones
y0 = 3 # condicion inicial
mu = [1,10,20]  


figure('Ejercicio 2a')


for K in mu:
    tini = perf_counter()
    (t, y) = euler(a, b, f3, N, y0, K) # llamada al metodo de Euler
    tfin=perf_counter()

    plot(t,y)

legend(['mu = 1','mu = 10','mu = 20'])


gcf().suptitle("Ejercicio 2a")
show()


### Apartado b
def ptoFijo(y,k,h,fun):
    z = zeros(2)
    z[0] = y[k]
    z[1] = y[k]
    for n in range(100):
        z = array([y[k] + h*(1/2 *fun(t[k],z[0],K) - 1/2 *fun(t[k]+h, z[1],K)),
                   y[k] + h*(1/2 *fun(t[k],z[0],K) + 1/2 *fun(t[k]+h, z[1],K))])

    return(z)

def rk(a, b, fun, N, y0, K):
    h = (b-a)/N # paso de malla
    t = zeros(N+1) # inicializacion del vector de nodos
    y = zeros(N+1) # inicializacion del vector de resultados
    t[0] = a # nodo inicial
    y[0] = y0 # valor inicial

    for k in range(N):
        z = ptoFijo(y,k,h,fun)
        y[k+1] = y[k] +h*(1/2 * fun(t[k],z[0],K) + 1/2 * fun(t[k]+ h, z[1],K))
        t[k+1] = t[k]+h
        print(y[k+1])
    return(t,y)


figure('Ejercicio 2b')


for K in mu:
    tini = perf_counter()
    (t, y) = rk(a, b, f3, N, y0, K) 
    tfin=perf_counter()

    plot(t,y)

gcf().suptitle("Ejercicio 2b")
show()
