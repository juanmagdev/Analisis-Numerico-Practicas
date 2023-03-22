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