from time import perf_counter
from matplotlib.pyplot import *
from numpy import *
from numpy.linalg import eig, norm

# Métodos
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

def localizar_frontera(rho, sigma):
    theta = arange(0, 2.*pi, 0.01)
    numer = polyval(rho, exp(theta*1j))
    denom = polyval(sigma, exp(theta*1j))
    mu = numer/denom
    x, y = real(mu), imag(mu)
    plot(x, y)
    grid(True)
    axis('equal')

# Ejericio 1
e = 15
d = 3
def f(t, y):
    """Definicion del sistema de ecuaciones diferenciales"""
    f1 = y[1]
    f2 = -d*y[1] - e*y[0]
    return array([f1, f2])

# Datos
a = 0
b = 20
N = 200
y0 = array([1, 0])

# Solucion

tini = perf_counter()           # tiempo inicial

(t, y) = AB4Sis(a, b, f, N, y0) # llamada al metodo de AB4 Sistemas

tfin=perf_counter()             # tiempo final


# Dibujamos las soluciones
figure("Ejercicio 1a")
subplot(2, 1, 1)
plot(t, y[0], '-*') # dibuja la solucion aproximada

# Resultados
print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
print('Paso de malla: ' + str((b-a)/N))
print('-----')

xlabel('t')
ylabel('y')
grid(True)



gcf().suptitle("Ejercicio 1a")

legend(['RK, N=' + str(N)])
subplot(2, 1, 2)
plot(y[0],y[1])

xlabel('t')
ylabel('y')
grid(True)


show() # muestra la grafica

localizar_frontera(array([1, -1, 0, 0, 0]), array([0, 55, -59, 37, - 9])/24) # AB4

title("Frontera AB4")
# show()

A = array([[0, 1], [-e, -d]])
print("Matriz A:\n", A)
print("Autovalores de A:", eig(A)[0])

rel, iml = -0.5, 3.1225
# Obtenemos h para llevar las raíces a la frontera
h = 0.1
plot(h*rel, h*iml, ".", h*rel, -h*iml, ".")
show()

# Como los puntos de los autovalores se encuentran dentro de la región que delimita la 
# frontera de estabilidad, eso indica que los autovalores están dentro de la región 
# de estabilidad absoluta del método. Esto es una buena señal, ya que significa que 
# el método numérico utilizado es estable para los autovalores correspondientes a la
# matriz del sistema.

# Ejercicio 2
# Método AM4
# def AM4Sis(a, b, fun, N, y0):
#     y = zeros([len(y0), N+1])
#     t = zeros(N+1)
#     f = zeros([len(y0), N+1])
#     t[0] = a
#     h = (b-a) / N
#     y[:, 0] = y0
#     f[:, 0] = fun(a, y[:, 0])
#     iteraciones = 0


#     ## Metodo de Runge-Kutta de orden 4 Sistemas
#     for k in range(3):
#         k1 = fun(t[k], y[:, k])
#         k2 = fun(t[k]+h/2, y[:, k]+h/2*k1)
#         k3 = fun(t[k]+h/2, y[:, k]+h/2*k2)
#         k4 = fun(t[k]+h, y[:, k]+h*k3)
#         t[k+1] = t[k]+h
#         y[:, k+1] = y[:, k] + h/6*(k1 + 2*k2 + 2*k3 + k4)
    
#     for k in range(3, N):
#         z = y[k] + h/24 * (55*f[:, k] - 59*f[:, k-1] + 37*f[:, k-2] - 9*f[:, k-3])
#         contador = 0
#         distancia = 1 + 1e-12
#         C = y[k] + h/720 * (251*f[:, k] + 646*f[:, k-1] - 264*f[:, k-2] + 106*f[:, k-3] - 19*f[:, k-4])
#         while distancia > 1e-12 and contador < 100:
#             z_nuevo = y[:, k] + h/720 * (251*fun(t[k+1], z) + 646*f[:, k] - 264*f[:, k-1] + 106*f[:, k-2] - 19*f[:, k-3])
#             distancia = norm(z_nuevo - z)
#             z = z_nuevo
#             contador += 1
#         if contador == 100:
#             print("No se ha alcanzado la convergencia en el paso", k)
#         iteraciones = max(iteraciones, contador)
#         y[:, k+1] = z
#         t[k+1] = t[k] + h
#         f[:, k+1] = fun(t[k+1], y[:, k+1])
    
#     return (t, y, iteraciones)


import numpy as np

def AM4Sis(a, b, fun, N, y0):
    y = np.zeros([len(y0), N+1])
    t = np.zeros(N+1)
    f = np.zeros([len(y0), N+1])
    t[0] = a
    h = (b-a) / N
    y[:, 0] = y0
    f[:, 0] = fun(a, y[:, 0])
    iteraciones = 0

    ## Método de Runge-Kutta de orden 4 para sistemas
    for k in range(3):
        k1 = fun(t[k], y[:, k])
        k2 = fun(t[k]+h/2, y[:, k]+h/2*k1)
        k3 = fun(t[k]+h/2, y[:, k]+h/2*k2)
        k4 = fun(t[k]+h, y[:, k]+h*k3)
        t[k+1] = t[k] + h
        y[:, k+1] = y[:, k] + h/6*(k1 + 2*k2 + 2*k3 + k4)
    
    for k in range(3, N):
        z = y[:, k] + h/24 * (55*f[:, k] - 59*f[:, k-1] + 37*f[:, k-2] - 9*f[:, k-3])
        contador = 0
        distancia = 1 + 1e-12
        C = y[:, k] + h/720 * (251*fun(t[k+1], z) + 646*f[:, k] - 264*f[:, k-1] + 106*f[:, k-2] - 19*f[:, k-3])
        while distancia > 1e-12 and contador < 200:
            z_nuevo = y[:, k] + h/720 * (251*fun(t[k+1], z) + 646*fun(t[k], y[:, k]) - 264*fun(t[k-1], y[:, k-1]) + 106*fun(t[k-2], y[:, k-2]) - 19*fun(t[k-3], y[:, k-3]))
            distancia = np.linalg.norm(z_nuevo - z)
            z = z_nuevo
            contador += 1
        if contador == 200:
            print("No se ha alcanzado la convergencia en el paso", k)
        iteraciones = max(iteraciones, contador)
        y[:, k+1] = z
        t[k+1] = t[k] + h
        f[:, k+1] = fun(t[k+1], y[:, k+1])
    
    return (t, y, iteraciones)


e = 15
d = 3
def f(t, y):
    """Definicion del sistema de ecuaciones diferenciales"""
    f1 = y[1]
    f2 = -d*y[1] - e*y[0]
    return array([f1, f2])

# Datos
a = 0
b = 20
N = 200
y0 = array([1, 0])

# Solucion

tini = perf_counter()           # tiempo inicial

(t, y, i) = AM4Sis(a, b, f, N, y0) # llamada al metodo de AB4 Sistemas

tfin=perf_counter()             # tiempo final


# Dibujamos las soluciones
figure("Ejercicio 1a")
subplot(2, 1, 1)
plot(t, y[0], '-*') # dibuja la solucion aproximada

# Resultados
print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
print('Paso de malla: ' + str((b-a)/N))
print('-----')

xlabel('t')
ylabel('y')
grid(True)



gcf().suptitle("Ejercicio 1a")

legend(['RK, N=' + str(N)])
subplot(2, 1, 2)
plot(y[0],y[1])

xlabel('t')
ylabel('y')
grid(True)


show() # muestra la grafica

localizar_frontera(array([1, -1, 0, 0, 0]), array([251, 646, -264, 106, -19])/720) # AM4

title("Frontera AB4")
# show()

A = array([[0, 1], [-e, -d]])
print("Matriz A:\n", A)
print("Autovalores de A:", eig(A)[0])

rel, iml = -0.5, 3.1225
# Obtenemos h para llevar las raíces a la frontera
h = 0.1
plot(h*rel, h*iml, ".", h*rel, -h*iml, ".")
show()
