# -*- coding: utf-8 -*-
"""
@author: Juan Manuel García Delgado
"""

from numpy import *  
from matplotlib.pyplot import *
# from pylab import *
from time import perf_counter

# Ejemplo de resolucion de un problema de valor inicial

def f(t, y):
    """Funcion que define la ecuacion diferencial"""
    return 0.5*(t**2 - y)

def exacta(t):
    """Solucion exacta del problema de valor inicial"""
    return t**2 - 4*t + 8 - 7.*exp(-0.5*t)

def euler(a, b, fun, N, y0):
    """Implementacion del metodo de Euler en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""
    
    h = (b-a)/N # paso de malla
    t = zeros(N+1) # inicializacion del vector de nodos
    y = zeros(N+1) # inicializacion del vector de resultados
    t[0] = a # nodo inicial
    y[0] = y0 # valor inicial

    # Metodo de Euler
    for k in range(N):
        y[k+1] = y[k]+h*fun(t[k], y[k])
        t[k+1] = t[k]+h
    
    return (t, y)

# Datos del problema
a = 0. # extremo inferior del intervalo
b = 10. # extremo superior del intervalo
N = 20 # numero de particiones
y0 = 1. # condicion inicial

tini = perf_counter()           # tiempo inicial

(t, y) = euler(a, b, f, N, y0) # llamada al metodo de Euler

tfin=perf_counter()             # tiempo final

ye = exacta(t) # calculo de la solucion exacta

# Dibujamos las soluciones
plot(t, y, '-*') # dibuja la solucion aproximada
plot(t, ye, 'k') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['Euler', 'exacta'])
grid(True)

# Calculo del error cometido
error = max(abs(y-ye))

# Resultados
print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
print('Error: ' + str(error))
print('Paso de malla: ' + str((b-a)/N))
print('-----')

gcf().suptitle("Ejemplo")
show() # muestra la grafica


# ------------------------------
# EJERCICIO 1a
malla=  [10,20,40,80,160] # numero de particiones

def mallaF(f, metodoUnipaso, exacta, malla):

    for N in malla:
        tini = perf_counter()
        (t, y) = metodoUnipaso(a, b, f, N, y0) # llamada al metodo de Euler
        tfin=perf_counter()
        ye = exacta(t) # calculo de la solucion exacta
        plot(t, y, '-*') # dibuja la solucion aproximada
        error = max(abs(y-ye))
        print('-----')
        print('Tiempo CPU: ' + str(tfin-tini))
        print('Error: ' + str(error))
        if N>malla[0]:
            cociente=errorold/error
            print('Cociente de errores: ' + str(cociente))
        errorold=error
        print('Paso de malla: ' + str((b-a)/N))
    

    print('-----')
    plot(t, ye, 'k') # dibuja la solucion exacta
    xlabel('t')
    ylabel('y')
    legend(['Euler, N=10','Euler, N=20','Euler, N=40','Euler, N=80','Euler, N=160' ,'exacta'])
    grid(True)
    show() # muestra la grafica

gcf().suptitle("Ejercicio 1a")
mallaF(f, euler, exacta, malla)


# ------------------------------
# EJERCICIO 1b
def f(t, y):
    """Funcion que define la ecuacion diferencial"""
    return 6-y/10

def exacta(t):
    """Solucion exacta del problema de valor inicial"""
    return 60*(1-exp(-t/10))

# Datos del problema
a = 0. # extremo inferior del intervalo
b = 20. # extremo superior del intervalo
N = 20 # numero de particiones
y0 = 0. # condicion inicial

tini = perf_counter()           # tiempo inicial

(t, y) = euler(a, b, f, N, y0) # llamada al metodo de Euler

tfin=perf_counter()             # tiempo final

ye = exacta(t) # calculo de la solucion exacta

# Dibujamos las soluciones
plot(t, y, '-*') # dibuja la solucion aproximada
plot(t, ye, 'k') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['Euler', 'exacta'])
grid(True)

# Calculo del error cometido
error = max(abs(y-ye))

# Resultados
print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
print('Error: ' + str(error))
print('Paso de malla: ' + str((b-a)/N))
print('-----')

gcf().suptitle("Ejercicio 1b")
show() # muestra la grafica

gcf().suptitle("Ejercicio 1b")
mallaF(f, euler, exacta, malla)

# ------------------------------
# EJERCICIO 1c
# Como cambia la velocidad, la nueva EDO es:
def f(t, y):
    """Funcion que define la ecuacion diferencial"""
    return 9-y/10

def exacta(t):
    """Solucion exacta del problema de valor inicial"""
    return 60*(1-exp(-t/10))

# Datos del problema
a = 0. # extremo inferior del intervalo
b = 20. # extremo superior del intervalo
N = 2000 # numero de particiones
y0 = 0. # condicion inicial

tini = perf_counter()           # tiempo inicial

(t, y) = euler(a, b, f, N, y0) # llamada al metodo de Euler

tfin=perf_counter()             # tiempo final

ye = exacta(t) # calculo de la solucion exacta

# Dibujamos las soluciones
plot(t, y, '-*') # dibuja la solucion aproximada
plot(t, ye, 'k') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['Euler', 'exacta'])
grid(True)

# Calculo del error cometido
error = max(abs(y-ye))

# Resultados
print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
print('Error: ' + str(error))
print('Paso de malla: ' + str((b-a)/N))
print('-----')

gcf().suptitle("Ejercicio 1c")
show() # muestra la grafica

# ------------------------------
# EJERCICIO 2
def taylor2(a, b, fun, d1fun, N, y0):
    """Implementacion del metodo de Taylor de orden 2 en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""
    
    h = (b-a)/N # paso de malla
    t = zeros(N+1) # inicializacion del vector de nodos
    y = zeros(N+1) # inicializacion del vector de resultados
    t[0] = a # nodo inicial
    y[0] = y0 # valor inicial

    # Metodo de Euler
    for k in range(N):
        y[k+1] = y[k]+h*fun(t[k], y[k])+h**2/2*d1fun(t[k], y[k])
        t[k+1] = t[k]+h
    
    return (t, y)

def taylor3(a, b, fun, d1fun, d2fun, N, y0):
    """Implementacion del metodo de Taylor de orden 2 en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""
    
    h = (b-a)/N # paso de malla
    t = zeros(N+1) # inicializacion del vector de nodos
    y = zeros(N+1) # inicializacion del vector de resultados
    t[0] = a # nodo inicial
    y[0] = y0 # valor inicial

    # Metodo de Euler
    for k in range(N):
        y[k+1] = y[k]+h*fun(t[k], y[k])+h**2/2*d1fun(t[k], y[k])+h**3/6*d2fun(t[k], y[k])
        t[k+1] = t[k]+h
    
    return (t, y)



def f(t, y):
    """Funcion que define la ecuacion diferencial"""
    return 6-y/10
# primera derivada de la funcion f
def f1(t, y):
    """Funcion que define la ecuacion diferencial"""
    return -1/10
# segunada derivada de la funcion f
def f2(t, y):
    """Funcion que define la ecuacion diferencial"""
    return 0

def exacta(t):
    """Solucion exacta del problema de valor inicial"""
    return 60*(1-exp(-t/10))

def mallaFTaylor2(f, d1,exacta, malla):
    for N in malla:
        tini = perf_counter()
        (t, y) = taylor2(a, b, f, d1, N, y0) # llamada al metodo de Euler
        tfin=perf_counter()
        ye = exacta(t) # calculo de la solucion exacta
        plot(t, y, '-*') # dibuja la solucion aproximada
        error = max(abs(y-ye))
        print('-----')
        print('Tiempo CPU: ' + str(tfin-tini))
        print('Error: ' + str(error))
        if N>malla[0]:
            cociente=errorold/error
            print('Cociente de errores: ' + str(cociente))
        errorold=error
        print('Paso de malla: ' + str((b-a)/N))
    
    print('-----')
    plot(t, ye, 'k') # dibuja la solucion exacta
    xlabel('t')
    ylabel('y')
    legend(['Euler, N=10','Euler, N=20','Euler, N=40','Euler, N=80','Euler, N=160' ,'exacta'])
    grid(True)
    show() # muestra la grafica

def mallaFTaylor3(f, d1, d2, exacta, malla):
    for N in malla:
        tini = perf_counter()
        (t, y) = taylor3(a, b, f, d1, d2, N, y0) # llamada al metodo de Euler
        tfin=perf_counter()
        ye = exacta(t) # calculo de la solucion exacta
        plot(t, y, '-*') # dibuja la solucion aproximada
        error = max(abs(y-ye))
        print('-----')
        print('Tiempo CPU: ' + str(tfin-tini))
        print('Error: ' + str(error))
        if N>malla[0]:
            cociente=errorold/error
            print('Cociente de errores: ' + str(cociente))
        errorold=error
        print('Paso de malla: ' + str((b-a)/N))
    
    print('-----')
    plot(t, ye, 'k') # dibuja la solucion exacta
    xlabel('t')
    ylabel('y')
    legend(['Euler, N=10','Euler, N=20','Euler, N=40','Euler, N=80','Euler, N=160' ,'exacta'])
    grid(True)
    show() # muestra la grafica
   

gcf().suptitle("Ejercicio 2 Taylor 2")
mallaFTaylor2(f,f1, exacta, malla)
gcf().suptitle("Ejercicio 2 Taylor 3")
mallaFTaylor3(f,f1,f2, exacta, malla)

# ------------------------------
# EJERCICIO 3
def Heun(a, b, fun, N, y0):
    """Implementacion del metodo de Heun en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""
    
    h = (b-a)/N # paso de malla
    t = zeros(N+1) # inicializacion del vector de nodos
    y = zeros(N+1) # inicializacion del vector de resultados
    t[0] = a # nodo inicial
    y[0] = y0 # valor inicial

    # Metodo de Euler
    for k in range(N):
        t[k+1] = t[k]+h
        ystar = y[k]+h*fun(t[k], y[k])
        y[k+1] = y[k]+0.5*h*(fun(t[k],y[k])+fun(t[k+1],ystar))
        
    return (t, y)


def PuntoMedio(a, b, fun, N, y0):
    """Implementacion del metodo de punto medio en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""
    
    h = (b-a)/N # paso de malla
    t = zeros(N+1) # inicializacion del vector de nodos
    y = zeros(N+1) # inicializacion del vector de resultados
    t[0] = a # nodo inicial
    y[0] = y0 # valor inicial

    # Metodo de punto medio
    for k in range(N):
        tstar = t[k]+h*0.5
        ystar = y[k]+0.5*h*fun(t[k], y[k])
        y[k+1] = y[k]+h*fun(tstar,ystar)
        t[k+1] = t[k]+h
    return (t, y)


def RK4(a, b, fun, N, y0):
    """Implementacion del metodo de punto medio en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""
    
    h = (b-a)/N # paso de malla
    t = zeros(N+1) # inicializacion del vector de nodos
    y = zeros(N+1) # inicializacion del vector de resultados
    t[0] = a # nodo inicial
    y[0] = y0 # valor inicial

    # Metodo de punto medio
    for k in range(N):
        k1=fun(t[k],y[k])
        tstar= t[k]+h/2
        k2=fun(tstar,y[k]+0.5*h*k1)
        k3=fun(tstar,y[k]+0.5*h*k2)
        t[k+1]=t[k] +h
        k4=fun(t[k+1],y[k]+h*k3)
        y[k+1]=y[k]+(h/6) *(k1+2*k2+2*k3+k4)
    return (t, y)

gcf().suptitle("Ejercicio 3 Heun")
mallaF(f, Heun, exacta, malla)
gcf().suptitle("Ejercicio 3 Punto Medio")
mallaF(f, PuntoMedio, exacta, malla)
gcf().suptitle("Ejercicio 3 RK4")
mallaF(f, RK4, exacta, malla)

# ------------------------------
#EJERCICIO 4 Modelo depredador/presa de Lotka-Volterra.


# ______________________________
# EJERCICIO 5
def f(t,y):
    f1= y[1] 
    f2= -20*y[1] -101*y[0]
    return array([f1,f2])


def eulerSis(a,b,fun,N,y0):
    h = (b-a)/N
    t = zeros(N+1)
    y = zeros([len(y0),N+1])
    t[0] = a
    y[:,0] = y0
    
    for k in range (N):
        y[:,k+1] = y[:,k] + h*fun(t[k],y[:,k])
        t[k+1] = t[k] + h
    return (t,y)


def exacta(t):
    return(exp(-10*t)*cos(t))


# Datos del problema
a = 0. # extremo inferior del intervalo
b = 7. # extremo superior del intervalo
y0 = array([1,-10]) # condicion inicial

malla=  [40,80,160,320,640] # numero de particiones

print('Método de Euler para Sistemas')
figure('Método de Euler para Sistemas')
for N in malla:
    tini = perf_counter()
    (t, y) = eulerSis(a, b, f, N, y0) # llamada al metodo de Euler
    tfin=perf_counter()
    ye=exacta(t)
    plot(t, y[0,:]) 
    error= max(abs(y[0]-ye))
    
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Paso de malla: ' + str((b-a)/N))
    print('Error : '+str(error))
    if(N>malla[0]):
        cociente = errorOld/error
        print('Cociente de errores: '+str(cociente))
    errorOld=error

gcf().suptitle("Ejercicio 5 Euler para Sistemas")
plot(t, exacta(t))
show()



### Usando el método de Heun

def heunSis(a, b, fun, N, y0):
    """Implementacion del metodo de eun para sistemas en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""
    
    h = (b-a)/N # paso de malla
    t = zeros(N+1) # inicializacion del vector de nodos
    y = zeros([len(y0), N+1]) # inicializacion del vector de resultados
    t[0] = a # nodo inicial
    y[:,0] = y0 # valor inicial
    ystar = zeros([len(y0),1])

    # Metodo de Heun
    for k in range(N):
        ystar[:,0] = y[:,k] + h*fun(t[k],y[:,k])
        t[k+1] = t[k]+h
        y[:, k+1] = y[:,k] + h*1/2*(fun(t[k], y[:,k]) + fun(t[k+1], ystar[:,0]))
        
    
    return (t, y)
print('-----')
print('Método de Heun')
figure('Método de Heun')
for N in malla:
    tini = perf_counter()
    (t, y) = heunSis(a, b, f, N, y0) # llamada al metodo de Euler
    tfin=perf_counter()
    ye=exacta(t)
    plot(t, y[0,:]) 
    error= max(abs(y[0]-ye))
    
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Paso de malla: ' + str((b-a)/N))
    print('Error : '+str(error))
    if(N>malla[0]):
        cociente = errorOld/error
        print('Cociente de errores: '+str(cociente))
    errorOld=error

gcf().suptitle("Ejercicio 5 Heun para Sistemas")
plot(t, exacta(t))
show()