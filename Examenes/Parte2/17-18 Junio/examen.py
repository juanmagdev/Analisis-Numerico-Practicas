from time import perf_counter
from matplotlib.pyplot import *
from numpy import *
from numpy.linalg import eig

def AB4Sis(a, b, fun, N, y0):
    y = zeros([len(y0), N+1])
    t = zeros(N+1)
    f = zeros([len(y0), N+1])
    t[0] = a
    h = (b-a) / N
    y[:, 0] = y0
    f[:, 0] = fun(a, y[:, 0])

    ## Metodo de Runge-Kutta de orden 4 Sistemas
    for k in range(3):
        k1 = fun(t[k], y[:, k])
        k2 = fun(t[k]+h/2, y[:, k]+h/2*k1)
        k3 = fun(t[k]+h/2, y[:, k]+h/2*k2)
        k4 = fun(t[k]+h, y[:, k]+h*k3)
        t[k+1] = t[k]+h
        y[:, k+1] = y[:, k] + h/6*(k1 + 2*k2 + 2*k3 + k4)
    
    for k in range(3, N):
        y[:, k+1] = y[:, k] + h/24 * (55*f[:, k] - 59*f[:, k-1] + 37*f[:, k-2] - 9*f[:, k-3])
        t[k+1] = t[k] + h
        f[:, k+1] = fun(t[k+1], y[:, k+1])
        
    return (t, y)

def rk4(a, b, fun, N, y0):
    h = (b-a)/N
    t = zeros(N+1)
    y = zeros(N+1)
    t[0] = a 
    y[0] = y0 

    for k in range(N):
        t[k+1] = t[k] + h
        k1 = fun(t[k], y[k])
        k2 = fun(t[k] + h/2, y[k] + h/2 * k1)
        k3 = fun(t[k] + h/2, y[k] + h/2 * k2)
        k4 = fun(t[k+1], y[k] + h*k3)
        y[k+1] = y[k] + h/6 *(k1 + 2*k2 + 2*k3 + k4)
    
    return (t, y)

def localizar_frontera(rho, sigma):
    theta = arange(0, 2.*pi, 0.01)
    numer = polyval(rho, exp(theta*1j))
    denom = polyval(sigma, exp(theta*1j))
    mu = numer/denom
    x, y = real(mu), imag(mu)
    plot(x, y)
    grid(True)
    axis('equal')

# Ejercicio 4

def f(t, y):
    """Definicion del sistema de ecuaciones diferenciales"""
    f1 = y[1]
    f2 = -20*y[1] - 101*y[0]
    return array([f1, f2])

def exacta(t):
    """Definicion de la solucion exacta"""
    return exp(-10*t) * (cos(t))

# Datos del problema
malla = [25, 100, 400]
y0 = array([1, -10])
a = 0
b = 7

for N in malla:
    tini = perf_counter()           # tiempo inicial

    (t, y) = AB4Sis(a, b, f, N, y0) # llamada al metodo de AB4 Sistemas

    tfin=perf_counter()             # tiempo final

    ye = exacta(t) # calculo de la solucion exacta

    # Dibujamos las soluciones
    plot(t, y[0], '-*') # dibuja la solucion aproximada

    # Calculo del error cometido
    error = max(abs(y[0]-ye))

    # Resultados
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Error: ' + str(error))
    print('Paso de malla: ' + str((b-a)/N))
    print('-----')

plot(t, ye, 'k') # dibuja la solucion exacta

xlabel('t')
ylabel('y')
grid(True)

gcf().suptitle("Ejercicio 3")
leyenda = ['RK, N=' + str(N) for N in malla]
leyenda.append('Exacta')
legend(leyenda)
show() # muestra la grafica


# Ejercicio 5
localizar_frontera(array([1, -1, 0, 0, 0]), array([0, 55, -59, 37, - 9])/24) # AB4Sis
title("Frontera AB4")


A = array([[0, 1], [-101, -20]])
print("Matriz A:\n", A)
print("Autovalores de A:", eig(A)[0])

rel, iml = -10, 1
# Obtenemos h para llevar las raíces a la frontera
h = 0.0299 ## estimar h
plot(h*rel, h*iml, ".", h*rel, -h*iml, ".")
show()

