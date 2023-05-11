from pylab import *
from time import perf_counter
from matplotlib.pyplot import *
from numpy import *

#######################################################
def AB2(a, b, fun, N, y0):
    y = zeros(N+1)
    t = zeros(N+1)
    f = zeros(N+1)
    t[0] = a
    h = (b-a) / N 
    y[0] = y0
    f[0] = fun(a, y[0])
    y[1] = y[0] + h*f[0]
    t[1] = a + h
    f[1] = fun(t[1], y[1])

    for k in range(1, N):
        y[k+1] = y[k] + h/2 * (3*f[k] - f[k-1])
        t[k+1] = t[k] + h
        f[k+1] = fun(t[k+1], y[k+1])
        
    return (t, y)

def AB3(a, b, fun, N, y0):
    y = zeros(N+1)
    t = zeros(N+1)
    f = zeros(N+1)
    t[0] = a
    h = (b-a) / N
    y[0] = y0
    f[0] = fun(a, y[0])

    for k in range(2):
        y[k+1] = y[k] + h*fun(t[k] + h/2, y[k] + h/2 * f[k])
        t[k+1] = t[k] + h
        f[k+1] = fun(t[k+1], y[k+1])
    
    for k in range(2, N):
        y[k+1] = y[k] + h/12 * (23*f[k] - 16*f[k-1] + 5*f[k-2])
        t[k+1] = t[k] + h
        f[k+1] = fun(t[k+1], y[k+1])
        
    return (t, y)

def AM3_punto_fijo(a, b, fun, N, y0):
    y = zeros(N+1)
    t = zeros(N+1)
    f = zeros(N+1)
    iteraciones = 0
    t[0] = a
    h = (b-a) / N
    y[0] = y0
    f[0] = fun(a, y[0])
    
    for k in range(2):
        t[k+1] = t[k] + h
        k1 = f[k]
        k2 = fun(t[k] + h/2, y[k] + h/2 * k1)
        k3 = fun(t[k] + h/2, y[k] + h/2 * k2)
        k4 = fun(t[k+1], y[k] + h*k3)
        y[k+1] = y[k] + h/6 *(k1 + 2*k2 + 2*k3 + k4)
        f[k+1] = fun(t[k+1], y[k+1])
    
    for k in range(2, N):
        z = y[k]
        contador = 0
        distancia = 1 + 1e-12
        C = y[k] + h/24 * (19*f[k] - 5*f[k-1] + f[k-2])
        t[k+1] = t[k] + h

        while distancia >= 1e-12 and contador < 200:
            z_nuevo = 9*h*fun(t[k+1], z)/24 + C
            distancia = abs(z - z_nuevo)
            z = z_nuevo
            contador += 1
        if contador == 200:
            print("El método no converge.")
        iteraciones = max(iteraciones, contador)
        y[k+1] = z
        f[k+1] = fun(t[k+1], y[k+1])

    return (t, y, iteraciones)

def AM3_newton(a, b, fun, dfun, N, y0):
    y = zeros(N+1)
    t = zeros(N+1)
    f = zeros(N+1)
    iteraciones = 0
    t[0] = a
    h = (b-a) / N
    y[0] = y0
    f[0] = fun(a, y[0])
    
    for k in range(2):
        t[k+1] = t[k] + h
        k1 = f[k]
        k2 = fun(t[k] + h/2, y[k] + h/2 * k1)
        k3 = fun(t[k] + h/2, y[k] + h/2 * k2)
        k4 = fun(t[k+1], y[k] + h*k3)
        y[k+1] = y[k] + h/6 *(k1 + 2*k2 + 2*k3 + k4)
        f[k+1] = fun(t[k+1], y[k+1])
    
    for k in range(2, N):
        z = y[k]
        contador = 0
        distancia = 1 + 1e-12
        C = y[k] + h/24 * (19*f[k] - 5*f[k-1] + f[k-2])
        t[k+1] = t[k] + h

        while distancia >= 1e-12 and contador < 200:
            z_nuevo = z - (z - 9*h*fun(t[k+1], z)/24 - C)/(1 - 9*h*dfun(t[k+1], z)/24)
            distancia = abs(z - z_nuevo)
            z = z_nuevo
            contador += 1
        if contador == 200:
            print("El método no converge.")
        iteraciones = max(iteraciones, contador)
        y[k+1] = z
        f[k+1] = fun(t[k+1], y[k+1])

    return (t, y, iteraciones)

def predictor_corrector(a, b, fun, N, y0):
    y = zeros(N+1)
    t = zeros(N+1)
    f = zeros(N+1)
    t[0] = a
    h = (b-a) / N
    y[0] = y0
    f[0] = fun(a, y[0])
    
    for k in range(2):
        t[k+1] = t[k] + h
        k1 = f[k]
        k2 = fun(t[k] + h/2, y[k] + h/2 * k1)
        k3 = fun(t[k] + h/2, y[k] + h/2 * k2)
        k4 = fun(t[k+1], y[k] + h*k3)
        y[k+1] = y[k] + h/6 *(k1 + 2*k2 + 2*k3 + k4)
        f[k+1] = fun(t[k+1], y[k+1])
    
    for k in range(2, N):
        t[k+1] = t[k] + h
        y_nuevo = y[k] + h/12 * (23*f[k] - 16*f[k-1] + 5*f[k-2])
        f_nuevo = fun(t[k+1], y_nuevo)
        y[k+1] = y[k] + h/24 * (9*f_nuevo + 19*f[k] - 5*f[k-1] + f[k-2])
        f[k+1] = fun(t[k+1], y[k+1])

    return (t, y)

def AB3_sistema(a, b, fun, N, y0):
    y = zeros([len(y0), N+1])
    t = zeros(N+1)
    f = zeros([len(y0), N+1])
    t[0] = a
    h = (b-a) / N
    y[:, 0] = y0
    f[:, 0] = fun(a, y[:, 0])

    for k in range(2):
        y[:, k+1] = y[:, k] + h*fun(t[k] + h/2, y[:, k] + h/2 * f[:, k])
        t[k+1] = t[k] + h
        f[:, k+1] = fun(t[k+1], y[:, k+1])
    
    for k in range(2, N):
        y[:, k+1] = y[:, k] + h/12 * (23*f[:, k] - 16*f[:, k-1] + 5*f[:, k-2])
        t[k+1] = t[k] + h
        f[:, k+1] = fun(t[k+1], y[:, k+1])
        
    return (t, y)

def AM3_punto_fijo_sistema(a, b, fun, N, y0):
    y = zeros([len(y0), N+1])
    t = zeros(N+1)
    f = zeros([len(y0), N+1])
    iteraciones = 0
    t[0] = a
    h = (b-a) / N
    y[:, 0] = y0
    f[:, 0] = fun(a, y[:, 0])
    
    for k in range(2):
        t[k+1] = t[k] + h
        k1 = f[:, k]
        k2 = fun(t[k] + h/2, y[:, k] + h/2 * k1)
        k3 = fun(t[k] + h/2, y[:, k] + h/2 * k2)
        k4 = fun(t[k+1], y[:, k] + h*k3)
        y[:, k+1] = y[:, k] + h/6 *(k1 + 2*k2 + 2*k3 + k4)
        f[:, k+1] = fun(t[k+1], y[:, k+1])
    
    for k in range(2, N):
        z = y[:, k]
        contador = 0
        distancia = 1 + 1e-12
        C = y[:, k] + h/24 * (19*f[:, k] - 5*f[:, k-1] + f[:, k-2])
        t[k+1] = t[k] + h

        while distancia >= 1e-12 and contador < 200:
            z_nuevo = 9*h*fun(t[k+1], z)/24 + C
            distancia = max(abs(z - z_nuevo))
            z = z_nuevo
            contador += 1
        if contador == 200:
            print("El método no converge.")
        iteraciones = max(iteraciones, contador)
        y[:, k+1] = z
        f[:, k+1] = fun(t[k+1], y[:, k+1])

    return (t, y, iteraciones)

def predictor_corrector_sistema(a, b, fun, N, y0):
    y = zeros([len(y0), N+1])
    t = zeros(N+1)
    f = zeros([len(y0), N+1])
    t[0] = a
    h = (b-a) / N
    y[:, 0] = y0
    f[:, 0] = fun(a, y[:, 0])
    
    for k in range(2):
        t[k+1] = t[k] + h
        k1 = f[:, k]
        k2 = fun(t[k] + h/2, y[:, k] + h/2 * k1)
        k3 = fun(t[k] + h/2, y[:, k] + h/2 * k2)
        k4 = fun(t[k+1], y[:, k] + h*k3)
        y[:, k+1] = y[:, k] + h/6 *(k1 + 2*k2 + 2*k3 + k4)
        f[:, k+1] = fun(t[k+1], y[:, k+1])
    
    for k in range(2, N):
        t[k+1] = t[k] + h
        y_nuevo = y[:, k] + h/12 * (23*f[:, k] - 16*f[:, k-1] + 5*f[:, k-2])
        f_nuevo = fun(t[k+1], y_nuevo)
        y[:, k+1] = y[:, k] + h/24 * (9*f_nuevo + 19*f[:, k] - 5*f[:, k-1] + f[:, k-2])
        f[:, k+1] = fun(t[k+1], y[:, k+1])

    return (t, y)

#######################################################
# Ejercicio 1
# Apartado a)
def f(t,y):
    return -y + exp(-t)*cos(t)

def exacta(t):
    return exp(-t)*sin(t)

a, b = 0, 5
N = 30

(t, y) = AB2(a, b, f, N, 0)

# Dibujamos las soluciones
ye = exacta(t)

figure("Ejercicio 1a")
title("Ejercicio 1a")

plot(t, y, '.')
plot(t, ye)
legend(["AB2", "Exacta"])
show()

h = (b-a) / N
error = max(abs(y - ye))
print("h =", h)
print("Error =", error)

# Apartado b)
def f(t,y):
    return -y + exp(-t)*cos(t)

def exacta(t):
    return exp(-t)*sin(t)

a, b = 0, 5
N = 30

(t, y) = AB3(a, b, f, N, 0)

# Dibujamos las soluciones
figure("Ejercicio 1b")
title("Ejercicio 1b")

ye = exacta(t)
plot(t, y, '.')
plot(t, ye)
legend(["AB3", "Exacta"])
show()

h = (b-a) / N
error = max(abs(y - ye))
print("h =", h)
print("Error =", error)

# Apartado c)
def punto_medio(a, b, fun, N, y0):    
    h = (b-a)/N
    t = zeros(N+1)
    y = zeros(N+1)
    t[0] = a
    y[0] = y0

    for k in range(N):
        t[k+1] = t[k] + h
        y[k+1] = y[k] + h * fun(t[k] + h/2, y[k] + h/2 * fun(t[k], y[k]))
    
    return (t, y)

def f(t,y):
    return -y + exp(-t)*cos(t)

def exacta(t):
    return exp(-t)*sin(t)

a, b = 0, 5
N = 3000
y0 = 0

(t, y) = punto_medio(a, b, f, N, y0)

figure("Ejercicio 1c")
title("Ejercicio 1c")

plot(t, y)
ye = exacta(t)
error_unipaso = max(abs(y - ye))
print("Error unipaso =", error_unipaso)

(t, y) = AB3(a, b, f, N, y0)
plot(t, y)
ye = exacta(t)
error_multipaso = max(abs(y - ye))
print("Error multipaso =", error_multipaso)

legend(["Unipaso", "Multipaso", "Exacta"])
show()



# Ejericio 2
# Apartado a)
def AM3_especifico(a, b, fun, N, y0):
    y = zeros(N+1)
    t = zeros(N+1)
    f = zeros(N+1)
    t[0] = a
    h = (b-a) / N
    y[0] = y0
    f[0] = fun(a, y[0])

    for k in range(2):
        t[k+1] = t[k] + h
        k1 = f[k]
        k2 = fun(t[k] + h/2, y[k] + h/2 * k1)
        k3 = fun(t[k] + h/2, y[k] + h/2 * k2)
        k4 = fun(t[k+1], y[k] + h*k3)
        y[k+1] = y[k] + h/6 *(k1 + 2*k2 + 2*k3 + k4)
        f[k+1] = fun(t[k+1], y[k+1])

    for k in range(2, N):
        C = y[k] + h/24 * (19*f[k] - 5*f[k-1] + f[k-2])
        t[k+1] = t[k] + h
        y[k+1] = (9*h*exp(-t[k+1])*cos(t[k+1])/24 + C) / (1 + 9*h/24)
        f[k+1] = fun(t[k+1], y[k+1])
    
    return (t, y)

def f(t,y):
    return -y + exp(-t)*cos(t)

def exacta(t):
    return exp(-t)*sin(t)

a, b = 0, 5
N = 30

(t, y) = AM3_especifico(a, b, f, N, 0)

figure("Ejercicio 2a")
title("Ejercicio 2a")
# Dibujamos las soluciones
ye = exacta(t)
plot(t, y, '.')
plot(t, ye)
legend(["AM3", "Exacta"])
show()

h = (b-a) / N
error = max(abs(y - ye))
print("h =", h)
print("Error =", error)

# Apartado b)
def f(t,y):
    return -y + exp(-t)*cos(t)

def exacta(t):
    return exp(-t)*sin(t)

a, b = 0, 5
N = 30

(t, y, iteraciones) = AM3_punto_fijo(a, b, f, N, 0)

# Dibujamos las soluciones
ye = exacta(t)

figure("Ejercicio 2b")
title("Ejercicio 2b")

plot(t, y, '.')
plot(t, ye)
legend(["AM3", "Exacta"])
show()

h = (b-a) / N
error = max(abs(y - ye))
print("h =", h)
print("Error =", error)
print("Iteraciones:", iteraciones)

# Apartado c)
def f(t, y):
    return -y + exp(-t)*cos(t)

def df(t, y):
    return -1

def exacta(t):
    return exp(-t)*sin(t)

a, b = 0, 5
N = 30

(t, y, iteraciones) = AM3_newton(a, b, f, df, N, 0)

# Dibujamos las soluciones
ye = exacta(t)

figure("Ejercicio 2c")
title("Ejercicio 2c")

plot(t, y, '.')
plot(t, ye)
legend(["AM3", "Exacta"])
show()

h = (b-a) / N
error = max(abs(y - ye))
print("h =", h)
print("Error =", error)
print("Iteraciones:", iteraciones)

# Apartado d)
def f(t, y):
    return 1 + y**2

def df(t, y):
    return 2*y

def exacta(t):
    return tan(t)

a, b = 0, 1
N = 320
y0 = 0

h = (b-a) / N
print("h =", h)

(t, y, iteraciones) = AM3_punto_fijo(a, b, f, N, y0)
ye = exacta(t)

figure("Ejercicio 2d")
title("Ejercicio 2d")

error = max(abs(y - ye))
print("Error punto fijo =", error)
print("Iteraciones punto fijo:", iteraciones)
plot(t, y, '.')

(t, y, iteraciones) = AM3_newton(a, b, f, df, N, y0)
ye = exacta(t)
error = max(abs(y - ye))
print("Error Newton =", error)
print("Iteraciones Newton:", iteraciones)
plot(t, y, '.')

plot(t, ye)
legend(["AM3 punto fijo", "AM3 Newton", "Exacta"])
show()




# Ejercicio 3
def f(t,y):
    return -y + exp(-t)*cos(t)

def exacta(t):
    return exp(-t)*sin(t)

a, b = 0, 5
N = 30

(t, y) = predictor_corrector(a, b, f, N, 0)

# Dibujamos las soluciones
ye = exacta(t)

figure("Ejercicio 3")
title("Ejercicio 3")

plot(t, y, '.')
plot(t, ye)
legend(["Predictor-corrector", "Exacta"])
show()

h = (b-a) / N
error = max(abs(y - ye))
print("h =", h)
print("Error =", error)




# Ejercicio 4
def f(t, y):
    (y1, y2) = (y[0], y[1])
    return array([0.25*y1 - 0.01*y1*y2, -y2 + 0.01*y1*y2])

a, b = 0, 20
N = 200
y0 = array([80, 30])

(t, y) = AB3_sistema(a, b, f, N, y0)
plot(t, y[0, :], t, y[1, :])
show()

(t, y, iteraciones) = AM3_punto_fijo_sistema(a, b, f, N, y0)
plot(t, y[0, :], t, y[1, :])
show()

(t, y) = predictor_corrector_sistema(a, b, f, N, y0)
plot(t, y[0, :], t, y[1, :])
show()