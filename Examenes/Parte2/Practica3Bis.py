from funciones import *
from time import perf_counter


# Ejercicio 1


def fun(t, y):
    return (-y + 2*sin(t))


def exacta(t):
    return ((pi + 1)*exp(-t) + sin(t) - cos(t))


y0 = pi
a = 0.0
b = 10

N = 50


tini = perf_counter()
(t, y) = AB3(a, b, fun, N, y0)
tfin = perf_counter()
ye = exacta(t)
# clf()
plot(t, y)

error = max(abs(y-ye))
tcpu = tfin-tini

print('---------------')
print('AB3')
print('---------------')
print('N = '+str(N))
print('Error= '+str(error))
print('Tiempo CPU= '+str(tcpu))
print('---------------')


tini = perf_counter()
(t, y, maxiter) = AM3_punto_fijo(a, b, fun, N, y0)
tfin = perf_counter()
ye = exacta(t)
# clf()
plot(t, y)

error = max(abs(y-ye))
tcpu = tfin-tini

print('---------------')
print('AM3')
print('---------------')
print('N = '+str(N))
print('Error= '+str(error))
print('maxiter : '+str(maxiter))
print('Tiempo CPU= '+str(tcpu))
print('---------------')


tini = perf_counter()
(t, y) = AB3(a, b, fun, N, y0)
tfin = perf_counter()
ye = exacta(t)
# clf()
plot(t, y)

error = max(abs(y-ye))
tcpu = tfin-tini

print('---------------')
print('ABM3')
print('---------------')
print('N = '+str(N))
print('Error= '+str(error))
print('Tiempo CPU= '+str(tcpu))
print('---------------')

legend(('AB3', 'AM3', 'ABM3'))
show()


# Ejercicio 2


def fun2(t, y):
    f1 = 3*y[0] - 2*y[1]
    f2 = -y[0] + 3*y[1] - 2*y[2]
    f3 = -y[1] + 3*y[2]
    return array([f1, f2, f3])


def exacta2(t):
    f1 = -1/4 * exp(5*t) + 3/2 * exp(3*t) - 1/4 * exp(t)
    f2 = 1/4 * exp(5*t) - 1/4 * exp(t)
    f3 = -1/8 * exp(5*t) - 3/4 * exp(3*t) - 1/8 * exp(t)
    return (array([f1, f2, f3]))


y0 = array([1, 0, -1])

a = 0.0
b = 1

N = 100

tini = perf_counter()
(t, y) = AB3_sistema(a, b, fun2, N, y0)
tfin = perf_counter()
ye = exacta2(t)


errorX = max(abs(y[0]-ye[0]))
errorY = max(abs(y[1]-ye[1]))
errorZ = max(abs(y[2]-ye[2]))
tcpu = tfin-tini

print('---------------')
print('AB3 Sistemas')
print('---------------')
print('N = '+str(N))
print('ErrorX= '+str(errorX))
print('ErrorY= '+str(errorY))
print('ErrorZ= '+str(errorZ))
print('Tiempo CPU= '+str(tcpu))
print('---------------')

tini = perf_counter()
(t, y, maxiter) = AM3_punto_fijo_sistema(a, b, fun2, N, y0)
tfin = perf_counter()
ye = exacta2(t)


errorX = max(abs(y[0]-ye[0]))
errorY = max(abs(y[1]-ye[1]))
errorZ = max(abs(y[2]-ye[2]))
tcpu = tfin-tini

print('---------------')
print('AM3 Sistemas')
print('---------------')
print('N = '+str(N))
print('ErrorX= '+str(errorX))
print('ErrorY= '+str(errorY))
print('ErrorZ= '+str(errorZ))
print('Maxiter = '+str(maxiter))
print('Tiempo CPU= '+str(tcpu))
print('---------------')


tini = perf_counter()
(t, y) = ABM3Sis(a, b, fun2, N, y0)  
tfin = perf_counter()
ye = exacta2(t)


errorX = max(abs(y[0]-ye[0]))
errorY = max(abs(y[1]-ye[1]))
errorZ = max(abs(y[2]-ye[2]))
tcpu = tfin-tini

print('---------------')
print('ABM3 Sistemas')
print('---------------')
print('N = '+str(N))
print('ErrorX= '+str(errorX))
print('ErrorY= '+str(errorY))
print('ErrorZ= '+str(errorZ))
print('Tiempo CPU= '+str(tcpu))
print('---------------')
