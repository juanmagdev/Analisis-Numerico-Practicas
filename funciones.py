from time import perf_counter
from matplotlib.pyplot import *
from numpy import *


def euler(a, b, fun, N, y0):
    """Implementacion del metodo de Euler en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""

    h = (b-a)/N  # paso de malla
    t = zeros(N+1)  # inicializacion del vector de nodos
    y = zeros(N+1)  # inicializacion del vector de resultados
    t[0] = a  # nodo inicial
    y[0] = y0  # valor inicial

    # Metodo de Eulerf
    for k in range(N):
        y[k+1] = y[k]+h*fun(t[k], y[k])
        t[k+1] = t[k]+h

    return (t, y)

# IMPLEMENTACION PARA EL EXAMEN DE 2018


def eulerImplicito2018(a, b, fun, N, y0):
    h = (b-a)/N  # paso de malla
    t = zeros(N+1)  # inicializacion del vector de nodos
    y = zeros([len(y0), N+1])  # inicializacion del vector de resultados
    t[0] = a  # nodo inicial
    y[:, 0] = y0  # valor inicial

    for k in range(N):
        y[1, k+1] = (y[1, k] - 101*h*y[0, k])/(1+20*h+101*h*h)
        y[0, k+1] = y[0, k] + h*y[1, k+1]
        t[k+1] = t[k] + h

    return (t, y)


def taylor2(a, b, N, y0):
    """Implementacion del metodo de Taylor2 en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""

    h = (b-a)/N  # paso de malla
    t = zeros(N+1)  # inicializacion del vector de nodos
    y = zeros(N+1)  # inicializacion del vector de resultados
    t[0] = a  # nodo inicial
    y[0] = y0  # valor inicial

    # Metodo de Taylor de 2º orden
    for k in range(N):
        dy = 0.5*(t[k]**2 - y[k])
        d2y = t[k]-0.5*dy
        y[k+1] = y[k]+h*dy + 0.5*h*h*d2y
        t[k+1] = t[k]+h

    return (t, y)


def taylor3(a, b, N, y0):
    """Implementacion del metodo de Taylor3 en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""

    h = (b-a)/N  # paso de malla
    t = zeros(N+1)  # inicializacion del vector de nodos
    y = zeros(N+1)  # inicializacion del vector de resultados
    t[0] = a  # nodo inicial
    y[0] = y0  # valor inicial

    # Metodo de Taylor de 3er orden
    for k in range(N):
        dy = 0.5*(t[k]**2 - y[k])
        d2y = t[k]-0.5*dy
        d3y = 1-0.5*d2y
        y[k+1] = y[k]+h*dy + 0.5*h*h*d2y + h*h*h/6 * d3y
        t[k+1] = t[k]+h

    return (t, y)


def puntoMedio(a, b, fun, N, y0):
    """Implementacion del metodo de punto medio en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""

    h = (b-a)/N  # paso de malla
    t = zeros(N+1)  # inicializacion del vector de nodos
    y = zeros(N+1)  # inicializacion del vector de resultados
    t[0] = a  # nodo inicial
    y[0] = y0  # valor inicial

    # Metodo de punto medio
    for k in range(N):
        tstar = t[k]+h*0.5
        ystar = y[k]+0.5*h*fun(t[k], y[k])
        y[k+1] = y[k]+h*fun(tstar, ystar)
        t[k+1] = t[k]+h
    return (t, y)


def heun(a, b, fun, N, y0):
    """Implementacion del metodo de Heun en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""

    h = (b-a)/N  # paso de malla
    t = zeros(N+1)  # inicializacion del vector de nodos
    y = zeros(N+1)  # inicializacion del vector de resultados
    t[0] = a  # nodo inicial
    y[0] = y0  # valor inicial

    # Metodo de Heun
    for k in range(N):
        t[k+1] = t[k]+h
        ystar = y[k]+h*fun(t[k], y[k])
        y[k+1] = y[k]+0.5*h*(fun(t[k], y[k])+fun(t[k+1], ystar))
    return (t, y)


def rk4(a, b, fun, N, y0):
    """Implementacion del metodo de RK4 en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""

    h = (b-a)/N  # paso de malla
    t = zeros(N+1)  # inicializacion del vector de nodos
    y = zeros(N+1)  # inicializacion del vector de resultados
    t[0] = a  # nodo inicial
    y[0] = y0  # valor inicial

    # Metodo de RK4
    for k in range(N):
        k1 = fun(t[k], y[k])
        tstar = t[k]+h/2
        k2 = fun(tstar, y[k]+0.5*h*k1)
        k3 = fun(tstar, y[k]+0.5*h*k2)
        t[k+1] = t[k] + h
        k4 = fun(t[k+1], y[k]+h*k3)
        y[k+1] = y[k]+(h/6) * (k1+2*k2+2*k3+k4)
    return (t, y)


def eulerSis(a, b, fun, N, y0):
    h = (b-a)/N
    t = zeros(N+1)
    y = zeros([len(y0), N+1])
    t[0] = a
    y[:, 0] = y0

    for k in range(N):
        y[:, k+1] = y[:, k] + h*fun(t[k], y[:, k])
        t[k+1] = t[k] + h

    return (t, y)


def heunSis(a, b, fun, N, y0):
    """Implementacion del metodo de eun para sistemas en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""

    h = (b-a)/N  # paso de malla
    t = zeros(N+1)  # inicializacion del vector de nodos
    y = zeros([len(y0), N+1])  # inicializacion del vector de resultados
    t[0] = a  # nodo inicial
    y[:, 0] = y0  # valor inicial
    ystar = zeros([len(y0), 1])

    # Metodo de Heun
    for k in range(N):
        ystar[:, 0] = y[:, k] + h*fun(t[k], y[:, k])
        t[k+1] = t[k]+h
        y[:, k+1] = y[:, k] + h*1/2 * \
            (fun(t[k], y[:, k]) + fun(t[k+1], ystar[:, 0]))

    return (t, y)


def puntoMedioSis(a, b, fun, N, y0):
    """Implementacion del metodo de pto medio para sistemas en el intervalo [a, b]
    usando N particiones y condicion inicial y0"""

    h = (b-a)/N  # paso de malla
    t = zeros(N+1)  # inicializacion del vector de nodos
    y = zeros([len(y0), N+1])  # inicializacion del vector de resultados
    t[0] = a  # nodo inicial
    y[:, 0] = y0  # valor inicial
    ystar = zeros([len(y0), 1])

    # Metodo de pto medio
    for k in range(N):
        ystar[:, 0] = y[:, k] + h/2*fun(t[k], y[:, k])
        t[k+1] = t[k]+h
        y[:, k+1] = y[:, k] + h*fun(t[k]+h/2, ystar[:, 0])

    return (t, y)


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


def rk23(a, b, fun, y0, h0, tol):
    """Implementacion del metodo encajado RK2(3)
    en el intervalo [a, b] con condicion inicial y0,
    paso inicial h0 y tolerancia tol"""

    hmin = (b-a)*1.e-5  # paso de malla minimo
    hmax = (b-a)/10.  # paso de malla maximo

    # coeficientes RK
    q = 3  # numero de etapas
    p = 2  # orden del método menos preciso
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
    y = array([y0])  # soluciones
    h = array([h0])  # pasos de malla
    K = zeros(3)
    k = 0  # contador de iteraciones

    while (t[k] < b):
        h[k] = min(h[k], b-t[k])  # ajuste del ultimo paso de malla
        for i in range(q):
            K[i] = fun(t[k]+C[i]*h[k], y[k]+h[k]*sum(A[i, :]*K))
        incrlow = sum(B*K)  # metodo de orden 2
        incrhigh = sum(BB*K)  # metodo de orden 3
        error = h[k]*(incrhigh-incrlow)  # estimacion del error
        y = append(y, y[k]+h[k]*incrlow)  # y_(k+1)
        t = append(t, t[k]+h[k])  # t_(k+1)
        hnew = 0.9*h[k]*abs(tol/error)**(1./(p+1))  # h_(k+1)
        hnew = min(max(hnew, hmin), hmax)  # hmin <= h_(k+1) <= hmax
        h = append(h, hnew)
        k += 1

    return (t, y, h)


def rk45(a, b, fun, y0, h0, tol):
    """Implementacion del metodo encajado RK2(3)
    en el intervalo [a, b] con condicion inicial y0,
    paso inicial h0 y tolerancia tol"""

    hmin = (b-a)*1.e-5  # paso de malla minimo
    hmax = (b-a)/10.  # paso de malla maximo

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
    y = array([y0])  # soluciones
    h = array([h0])  # pasos de malla
    K = zeros(6)
    k = 0  # contador de iteraciones

    while (t[k] < b):
        h[k] = min(h[k], b-t[k])  # ajuste del ultimo paso de malla
        for i in range(q):
            K[i] = fun(t[k]+C[i]*h[k], y[k]+h[k]*sum(A[i, :]*K))
        incrlow = sum(B*K)  # metodo de orden 2
        incrhigh = sum(BB*K)  # metodo de orden 3
        error = h[k]*(incrhigh-incrlow)  # estimacion del error
        y = append(y, y[k]+h[k]*incrlow)  # y_(k+1)
        t = append(t, t[k]+h[k])  # t_(k+1)
        hnew = 0.9*h[k]*abs(tol/error)**(1/5)  # h_(k+1)
        hnew = min(max(hnew, hmin), hmax)  # hmin <= h_(k+1) <= hmax
        h = append(h, hnew)
        k += 1

    print("Numero de iteraciones: "+str(k))
    return (t, y, h)


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


def rk23SisPIZARRA(a, b, fun, y0, h0, tol):

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
    y = zeros([len(y0), 1])  # soluciones
    y[:, 0] = y0
    # y = y0.reshape(len(y0),1)
    h = array([h0])  # pasos de malla
    K = zeros([len(y0), q])
    k = 0  # contador de iteraciones

    while (t[k] < b):
        h[k] = min(h[k], b-t[k])  # ajuste del ultimo paso de malla
        for i in range(3):
            K[:, i] = fun(t[k]+C[i]*h[k], y[:, k]+h[k]
                          * dot(A[i, :], transpose(K)))

        incrlow = dot(B, transpose(K))  # metodo de orden 2
        incrhigh = dot(BB, transpose(K))  # metodo de orden 3

        error = norm(h[k]*(incrhigh-incrlow), inf)  # estimacion del error
        y = column_stack((y, y[:, k]+h[k]*incrlow))
        t = append(t, t[k]+h[k])  # t_(k+1)
        hnew = 0.9*h[k]*abs(tol/error)**(1./q)  # h_(k+1)
        hnew = min(max(hnew, hmin), hmax)  # hmin <= h_(k+1) <= hmax
        h = append(h, hnew)
        k += 1

    return (t, y, h)


def rk45Sis(a, b, fun, y0, h0, tol):

    hmin = (b-a)*1.e-5  # paso de malla minimo
    hmax = (b-a)/10.  # paso de malla maximo

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
