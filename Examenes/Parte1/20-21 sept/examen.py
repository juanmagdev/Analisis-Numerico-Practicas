from time import perf_counter
from matplotlib.pyplot import *
from numpy import *


# y[0] = x
# y[1] = y
def f(t, y):
    f1 = 1/y[0] - (exp(t**2)/t**2) * y[1] - t
    f2 = 1/y[1] - exp(t**2) - 2*t*exp(-t**2)
    return array([f1, f2])

def exacta(t):
    f1 = 1/t
    f2 = exp(-t**2)
    return array([f1, f2])

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
a = 1
b = 2
N = 1000
y0 = array([1, exp(-1)])

tini = perf_counter()           # tiempo inicial

(t, y) = rk4Sis(a, b, f, N, y0) # llamada al metodo de RK4

tfin=perf_counter()             # tiempo final

ye = exacta(t) # calculo de la solucion exacta

# Dibujamos las soluciones
figure("1a")
subplot(2, 1, 1)
title("X frente a T")
plot(t, y[0], '-*') # dibuja la solucion aproximada
plot(t, ye[0], 'k') # dibuja la solucion exacta
xlabel('t')
ylabel('x')
legend(['RK', 'exacta'])
grid(True)

subplot(2, 1, 2)
title("Y frente a T")
plot(t, y[1], '-*') # dibuja la solucion aproximada
plot(t, ye[1], 'k') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['RK', 'exacta'])
grid(True)

subplots_adjust(hspace=0.8)

# Calculo del error cometido
error = max(abs(y[0]-ye[0]))

# Resultados
print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
print('Error: ' + str(error))
print('Paso de malla: ' + str((b-a)/N))
print('-----')


# gcf().tight_layout()
show()


# Apartado b
malla = [20, 40, 80, 160, 320, 640]

for N in malla:
    tini = perf_counter()           # tiempo inicial

    (t, y) = rk4Sis(a, b, f, N, y0) # llamada al metodo de RK4

    tfin=perf_counter()             # tiempo final

    ye = exacta(t) # calculo de la solucion exacta

    # Dibujamos las soluciones
    figure("1a")
    subplot(2, 1, 1)
    title("X frente a T")
    plot(t, y[0], '-*') # dibuja la solucion aproximada
    plot(t, ye[0], 'k') # dibuja la solucion exacta
    xlabel('t')
    ylabel('y')
    legend(['RK', 'exacta'])
    grid(True)

    subplot(2, 1, 2)
    title("Y frente a T")
    plot(t, y[1], '-*') # dibuja la solucion aproximada
    plot(t, ye[1], 'k') # dibuja la solucion exacta
    xlabel('t')
    ylabel('y')
    legend(['RK', 'exacta'])
    grid(True)

    subplots_adjust(hspace=0.8)

    # Calculo del error cometido
    error = max(abs(y[0]-ye[0]))

    # Resultados
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Error: ' + str(error))
    print('Paso de malla: ' + str((b-a)/N))
    print('-----')

show()

def rk23Sis(a, b, fun, y0, h0, tol):

    hmin = (b-a)*1.e-5  # paso de malla minimo
    hmax = (b-a)/10.  # paso de malla maximo

    # coeficientes RK
    q = 3  # orden del metodo mas uno
    A = zeros([q, q])
    A[1, 0] = 1./2.
    A[2, 0] = -1.
    A[2, 1] = 2.

    B = zeros(q)
    B[1] = 1.

    BB = zeros(q)
    BB[0] = 1./6.
    BB[1] = 2./3.
    BB[2] = 1./6.

    C = zeros(q)
    for i in range(q):
        C[i] = sum(A[i, :])

    # inicializacion de variables
    t = array([a])  # nodos
    y = y0.reshape(len(y0), 1)  # soluciones
    h = array([h0])  # pasos de malla
    K = zeros([len(y0), 3])
    k = 0  # contador de iteraciones

    while (t[k] < b):
        h[k] = min(h[k], b-t[k])  # ajuste del ultimo paso de malla
        for i in range(3):
            K[:, i] = fun(t[k]+C[i]*h[k], y[:, k]+h[k]
                          * dot(A[i, :], transpose(K)))

        incrlow = dot(B, transpose(K))  # metodo de orden 2
        incrhigh = dot(BB, transpose(K))  # metodo de orden 3

        error = linalg.norm(h[k]*(incrhigh-incrlow), inf)  # estimacion del error
        y = column_stack((y, y[:, k]+h[k]*incrlow))
        t = append(t, t[k]+h[k])  # t_(k+1)
        hnew = 0.9*h[k]*abs(tol/error)**(1./q)  # h_(k+1)
        hnew = min(max(hnew, hmin), hmax)  # hmin <= h_(k+1) <= hmax
        h = append(h, hnew)
        k += 1

    return (t, y, h)

tol = 1.e-6
h0 = 0.01


tini = perf_counter()           # tiempo inicial

(t, y, h) = rk23Sis(a, b, f, y0, h0, tol) # llamada al metodo de RK4

tfin=perf_counter()             # tiempo final

ye = exacta(t) # calculo de la solucion exacta

# Dibujamos las soluciones
figure("2")
subplot(2, 2, 1)
title("X frente a T")
plot(t, y[0], '-*') # dibuja la solucion aproximada
plot(t, ye[0], 'k') # dibuja la solucion exacta
xlabel('t')
ylabel('x')
legend(['RK', 'exacta'])
grid(True)

subplot(2, 2, 2)
title("Y frente a T")
plot(t, y[1], '-*') # dibuja la solucion aproximada
plot(t, ye[1], 'k') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['RK', 'exacta'])
grid(True)

subplot(2, 2, 3)
title("Orbitas")
plot(y[0], y[1], '-*') # dibuja la solucion aproximada
plot(ye[0], ye[1], 'k') # dibuja la solucion exacta
xlabel('x')
ylabel('y')
legend(['RK', 'exacta'])
grid(True)

subplot(2, 2, 4)
title("Pasos hk")
plot(t, h, '-*') # dibuja la solucion aproximada
xlabel('t')
ylabel('h')
legend(['RK', 'exacta'])
grid(True)

subplots_adjust(hspace=0.8)
subplots_adjust(wspace=0.3)

# Calculo del error cometido
error = max(abs(y[0]-ye[0]))

# Resultados
print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
print('Error: ' + str(error))
print('Paso de malla: ' + str((b-a)/N))
print('-----')


# gcf().tight_layout()
show()