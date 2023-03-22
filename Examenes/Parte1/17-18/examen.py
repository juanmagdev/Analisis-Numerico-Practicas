from time import perf_counter
from matplotlib.pyplot import *
from numpy import *

def f(t, y):
    f1 = y[1]
    f2 = -20*y[1] - 101*y[0]
    return array([f1, f2])

def exacta(t):
    return exp(-10*t) * (cos(t))

def rk4Sis(a, b, fun, N, y0):
    """Implementacion del metodo de rk4 para sistemas en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""

    h = 7/N  # paso de malla
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
malla = [25, 100, 400]
y0 = array([1, -10])
a = 0
b = 7

for N in malla:
    tini = perf_counter()           # tiempo inicial

    (t, y) = rk4Sis(a, b, f, N, y0) # llamada al metodo de RK4

    tfin=perf_counter()             # tiempo final

    ye = exacta(t) # calculo de la solucion exacta

    # Dibujamos las soluciones
    plot(t, y[0], '-*') # dibuja la solucion aproximada
    xlabel('t')
    ylabel('y')
    legend(['RK', 'exacta'])
    grid(True)

    # Calculo del error cometido
    error = max(abs(y[0]-ye))

    # Resultados
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Error: ' + str(error))
    print('Paso de malla: ' + str((b-a)/N))
    print('-----')

plot(t, ye, 'k') # dibuja la solucion exacta

gcf().suptitle("Ejercicio 1")
show() # muestra la grafica


def eulerImplicito(a,b,fun,N,y0):
    h = (b-a)/N # paso de malla
    t = zeros(N+1) # inicializacion del vector de nodos
    y = zeros([len(y0), N+1]) # inicializacion del vector de resultados
    t[0] = a # nodo inicial
    y[:,0] = y0 # valor inicial
    
    for k in range(N):
        y[1,k+1] = (y[1,k] -101*h*y[0,k])/(1+20*h+101*h*h)
        y[0,k+1] = y[0,k] +h*y[1,k+1]
        t[k+1] = t[k] + h
    
    return(t,y)


figure('Apartado 3')
for N in malla:
    tini = perf_counter()
    (t, y) = eulerImplicito(a, b, f, N, y0) 
    tfin=perf_counter()
    ye = exacta(t)
    plot(t, y[0]) 
    error= max(abs(y[0,:]-ye))
    print('-----')
    print('N = '+str(N))
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Paso de malla: ' + str((b-a)/N))
    print('Error : '+str(error))
    if(N>malla[0]):
        print('Cociente de errores : '+str(errorOld/error))
    errorOld= error

plot(t,ye)
show()


