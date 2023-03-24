# Autor: Juan Manuel Garcia Delgado
# Doble Grado Ingieneria Informática y Matemáticas

from time import perf_counter
from matplotlib.pyplot import *
from numpy import *

# Ejercicio 1: Un cohete...
M = 10
F = 100
g = 9.81
C = 0.02

# y'' = -g + F/M - C/M * y'*absy'
def f(t, y):
    """Funcion que devuelve el lado derecho de la ecuacion diferencial"""
    f1 = y[1]
    f2 = -g + F/M - C/M * abs(y[1])*y[1]
    return array([f1, f2])


def rk4Sis(a, b, fun, N, y0):
    """Implementacion del metodo de rk4 para sistemas en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""

    h = 0.01  # paso de malla
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
        if (y[0, k] > 1000):
            break
        if (t[k] > 1000):
            break
    
        k1 = fun(t[k], y[:, k])
        k2 = fun(t[k]+h/2, y[:, k]+h/2*k1)
        k3 = fun(t[k]+h/2, y[:, k]+h/2*k2)
        k4 = fun(t[k]+h, y[:, k]+h*k3)
        t[k+1] = t[k]+h
        y[:, k+1] = y[:, k] + h/6*(k1 + 2*k2 + 2*k3 + k4)

    return (t, y)


# Datos del problema
a = 0
b = 1000
y0 = array([0, 0])


# Creamos la malla de puntos, para ver la evolucion de las aproximaciones
malla = [10, 100, 1000, 10000]

figure("Ejercicio 1")
for N in malla:
    tini = perf_counter()           # tiempo inicial

    (t, y) = rk4Sis(a, b, f, N, y0) # llamada al metodo de RK4

    tfin=perf_counter()             # tiempo final
    

    subplot(2, 1, 1)
    title("Altura frente al tiempo")
    # Dibujamos la altura frente al tiempo
    plot(t, y[0], '-*') # dibuja la solucion aproximada
    xlabel('t')
    ylabel('y')
    grid(True)

    subplot(2, 1, 2)
    title("Velocidad frente al tiempo")
    # Dibujamos la velocidad frente al tiempo
    plot(t, y[1], '-*') # dibuja la solucion aproximada
    xlabel('t')
    ylabel('dy/dt')
    grid(True)

    # Calculo del error cometido
    # Resultados
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Paso de malla: ' + str((b-a)/N))
    print('-----')

leyenda = ['RK, N=' + str(N) for N in malla]
legend(leyenda)
subplots_adjust(hspace=0.8)

print("El tiempo en alcanzar los 1000 metros es de " + str(t[len(t)-1]) + " segundos.")
print("La velocidad en alcanzar los 1000 metros es de " + str(y[1, len(t)-2]) + " m/s.")

show()

# Ejercicio 2: 
def rk45Sis(a, b, fun, y0, h0, tol):

    hmin = 1.e-5  # paso de malla minimo
    hmax = 0.5  # paso de malla maximo

    # coeficientes RK
    q = 6  # orden del metodo mas uno
    A = zeros([q, q])
    A[1, 0] = 1/4
    A[2, 0] = 3/32
    A[2, 1] = 9/32
    A[3, 0] = 1932/2197
    A[3, 1] = -7200/2197
    A[3, 2] = 7296/2197
    A[4, 0] = 439/216
    A[4, 1] = -8
    A[4, 2] = 3680/513
    A[4, 3] = -845/4104
    A[5, 0] = -8/27
    A[5, 1] = 2
    A[5, 2] = -3544/2565
    A[5, 3] = 1859/4104
    A[5, 4] = -11/40

    B = zeros(q)
    B[0] = 25/216
    B[2] = 1408/2565
    B[3] = 2197/4104
    B[4] = -1/5

    BB = zeros(q)
    BB[0] = 16/135
    BB[2] = 6656/12825
    BB[3] = 28561/56430
    BB[4] = -9/50
    BB[5] = 2/55

    C = zeros(q)
    for i in range(q):
        C[i] = sum(A[i, :])

    # inicializacion de variables
    t = array([a])  # nodos
    y = y0.reshape(len(y0), 1)   # soluciones
    h = array([h0])  # pasos de malla
    K = zeros([len(y0), q])
    k = 0  # contador de iteraciones

    while (t[k] < b):
        if (y[0, k] > 1000):
            break
        if (t[k] > 1000):
            break

        h[k] = min(h[k], b-t[k])  # ajuste del ultimo paso de malla
        for i in range(q):
            K[:, i] = fun(t[k]+C[i]*h[k], y[:, k]+h[k]
                          * dot(A[i, :], transpose(K)))

        incrlow = dot(B, transpose(K))  # metodo de orden 4
        incrhigh = dot(BB, transpose(K))  # metodo de orden 5

        error = linalg.norm(h[k]*(incrhigh-incrlow), inf)  # estimacion del error
        y = column_stack((y, y[:, k]+h[k]*incrlow))
        t = append(t, t[k]+h[k])  # t_(k+1)
        hnew = 0.9*h[k]*abs(tol/error)**(1./5)  # h_(k+1)
        hnew = min(max(hnew, hmin), hmax)  # hmin <= h_(k+1) <= hmax
        h = append(h, hnew)
        k += 1

    return (t, y, h)

# Datos del problema

tol = 1.e-8
h0 = 0.01
a = 0
b = 1000

tini = perf_counter()
(t, y, h) = rk45Sis(a, b, f, y0, h0, tol) 
tfin = perf_counter()


# Resultados
print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
print("El tiempo en alcanzar los 1000 metros es de " + str(t[len(t)-1]) + " segundos.")
print("La velocidad en alcanzar los 1000 metros es de " + str(y[1, len(t)-2]) + " m/s.")
print('-----')

figure("Ejercicio 2")

subplot(3, 1, 1)
title("Altura frente al tiempo")
plot(t, y[0], '-*') # dibuja la solucion aproximada
xlabel('t')
ylabel('y')
legend(['RK45'])

subplot(3, 1, 2)
title("Velocidad frente al tiempo")
plot(t, y[1], '-*') # dibuja la solucion aproximada
xlabel('t')
ylabel('dy/dt')
legend(['RK45'])

subplot(313)
title("Pasos hk")
plot(t[:-1], h[:-1]) 

subplots_adjust(hspace=0.8)
show()

# Podemos observar que las graficas de altura de ambos metodos difieren un poco,
# Seria util tener la ecuacion exacta para poder compararlas.
# El primer metodo detecta que la altura de 1000m se alcanza en 100 segundos, mientras que el segundo metodo lo hace en 130 segundos.
# Sin embargo, la velocidad a los 1000m es muy similar en ambos metodos, 9.35 m/s y 9.65 m/s respectivamente.
# En ambas graficas podemos ver que la altura es creciente, pero la velocidad va disminuyendo, lo cual significa que en algun momento 
# la velocidad sera negativa, y el cohete comenzara a caer. 
# No se llega a ver el momento de caida por la condicion de calcular solo a los 1000m.