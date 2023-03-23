from time import perf_counter
from matplotlib.pyplot import *
from numpy import *

# y′ = 1 − 10y,
# y(0) = 1.
def f(t, y):
    return 1 - 10*y

# EDO de variables separables, resuelta a papel:
def exacta(t):
    return (9*exp(-10*t)  +1)/10

def rk(a, b, fun, N, y0):
    """Implementacion del metodo de RK4 en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""

    h = (b-a)/N  # paso de malla
    t = zeros(N+1)  # inicializacion del vector de nodos
    y = zeros(N+1)  # inicializacion del vector de resultados
    t[0] = a  # nodo inicial
    y[0] = y0  # valor inicial

    # Metodo de RK
    for k in range(N):
        # t1= t[k] +h/3
        # y1 = y[k] + 1/3*h*fun(y1, y1)
        # y1 = y[k] + 1/3*h*(1-10*y1)
        # y1 = y[k] + 1/3*h - 1/3*h*10*y1
        # ...

        k1 = fun(t[k] +h/3, (y[k] + h/3)/1  + 10*h/3)
        t[k+1] = t[k] + h
        y[k+1] = y[k]+(h) * (k1)
    return (t, y)


# apartado b)
# Datos del problema
a = 0   
b  =2
malla = [20, 40, 80, 160, 320, 640, 1280]
y0 = 1

for N in malla:
    tini = perf_counter()           # tiempo inicial

    (t, y) = rk(a, b, f, N, y0) # llamada al metodo de RK

    tfin=perf_counter()             # tiempo final

    ye = exacta(t) # calculo de la solucion exacta

    # Dibujamos las soluciones
    plot(t, y, '-*') # dibuja la solucion aproximada
    xlabel('t')
    ylabel('y')
    # legend(['RK', 'exacta'])
    grid(True)

    # Calculo del error cometido
    error = max(abs(y-ye))

    # Resultados
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Error: ' + str(error))
    print('Paso de malla: ' + str((b-a)/N))

leyenda = ['RK, N=' + str(N) for N in malla]
leyenda.append('Exacta')
legend(leyenda)
plot(t, ye, 'k') # dibuja la solucion exacta
show()


# Ejercicio 2

# y′′ + 4y′ + 29y = 0
def f(t, y):
    f1 = y[1]
    f2 = -4*y[1] - 29*y[0]
    return array([f1, f2])

# y(t) = e−2t cos(5t)
def exacta(t):
    return exp(-2*t) * cos(5*t)

def rk12Sis(a, b, fun, y0, h0, tol):

    hmin = 1.e-5  # paso de malla minimo
    hmax = 1.  # paso de malla maximo

    # coeficientes RK
    q = 2  # numero de etapas
    p = 1  # orden del método menos preciso
    A = zeros([q, q])
    A[1, 0] = 1/2.

    B = zeros(q)
    B[0] = 1.

    BB = zeros(q)
    BB[1] = 1

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
        h[k] = min(h[k], b-t[k])  # ajuste del ultimo paso de malla
        for i in range(q):
            K[:, i] = fun(t[k]+C[i]*h[k], y[:, k]+h[k]
                          * dot(A[i, :], transpose(K)))

        incrlow = dot(B, transpose(K))  # metodo de orden 4
        incrhigh = dot(BB, transpose(K))  # metodo de orden 5

        error = linalg.norm(h[k]*(incrhigh-incrlow), inf)  # estimacion del error
        y = column_stack((y, y[:, k]+h[k]*incrlow))
        t = append(t, t[k]+h[k])  # t_(k+1)
        hnew = 0.9*h[k]*abs(tol/error)**(1./(p+1))  # h_(k+1)
        hnew = min(max(hnew, hmin), hmax)  # hmin <= h_(k+1) <= hmax
        h = append(h, hnew)
        k += 1

    return (t, y, h)


# Datos del problema
a = 0
b = 5
y0 = array([1, -2])
h0 = 0.003
tol = 1.e-5


tini = perf_counter()
(t, y, h) = rk12Sis(a, b, f, y0, h0, tol) 
tfin = perf_counter()

ye = exacta(t) # calculo de la solucion exacta

# Calculo del error cometido
error = max(abs(y[0]-ye))

# Resultados
print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
print('Error: ' + str(error))
print('Paso de malla: ' + str((b-a)/N))
print('-----')

subplot(2, 1, 1)
title("exacta y RK12")
plot(t, y[0], '-*') # dibuja la solucion aproximada
plot(t, ye, 'k') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['RK12', 'exacta'])

subplot(212)
title("Pasos hk")
plot(t[:-1], h[:-1]) 

subplots_adjust(hspace=0.8)
show()