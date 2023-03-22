from time import perf_counter
from matplotlib.pyplot import *
from numpy import *

def f(t, y):
    return -9*y

def exacta(t):
    return exp(-9*t)

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
        k1 = fun(t[k] +h/3, y[k]/(1+3*h))
        k2 = fun(t[k] +(2*h)/3, y[k]/(1+3*h)+ (h*k1)/(3/1+3*h))
        t[k+1] = t[k] + h
        y[k+1] = y[k]+(h/2) * (k1+ k2)
    return (t, y)


# Datos del problema
a = 0
b  =2
malla = [10, 20, 40, 80, 160]
y0 = 1

for N in malla:
    tini = perf_counter()           # tiempo inicial

    (t, y) = rk(a, b, f, N, y0) # llamada al metodo de RK

    tfin=perf_counter()             # tiempo final

    ye = exacta(t) # calculo de la solucion exacta

    # Dibujamos las soluciones
    plot(t, y, '-*') # dibuja la solucion aproximada
    plot(t, ye, 'k') # dibuja la solucion exacta
    xlabel('t')
    ylabel('y')
    legend(['RK', 'exacta'])
    grid(True)

    # Calculo del error cometido
    error = max(abs(y-ye))

    # Resultados
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Error: ' + str(error))
    print('Paso de malla: ' + str((b-a)/N))
    print('-----')

gcf().suptitle("Ejercicio 1")
show() # muestra la grafica

# Ejercicio 2
def f(t, y):
    f1 = y[1]
    f2 = 2*(y[0]-t)*(y[1]-1)
    return(array([f1,f2]))

def exacta(t):
    return tan(t)+t

def rk45Sis(a, b, fun, y0, h0, tol):

    hmin = 1.e-5  # paso de malla minimo
    hmax = 0.1  # paso de malla maximo

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
a = 0 # extremo inferior del intervalo
b = 1.3 # extremo superior del intervalo
y0 = array([0,2]) # condicion inicial
h0 = 2.e-4 # paso de malla inicial
tol = 1.e-5 # tolerancia


tini = perf_counter()
(t, y, h) = rk45Sis(a, b, f, y0, h0, tol) 
tfin = perf_counter()

ye = exacta(t) # calculo de la solucion exacta

# Calculo del error cometido
# error = max(abs(y-ye))

# Resultados
print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
# print('Error: ' + str(error))
print('Paso de malla: ' + str((b-a)/N))
print('-----')

# Dibujamos las soluciones
figure('Ejercicio 1a')
subplot(211)
grid(True)
# Dibujamos las soluciones
plot(t, y[1], '-*') # dibuja la solucion aproximada
plot(t, ye, 'k') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['RK45', 'exacta'])


subplot(212)
plot(t, h) 

gcf().suptitle("Ejercicio 2")
show()

# La evolución de los pasos de tiempo se debe...