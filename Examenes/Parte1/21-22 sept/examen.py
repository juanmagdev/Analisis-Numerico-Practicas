from time import perf_counter
from matplotlib.pyplot import *
from numpy import * 


# y′ = 4 − 1000y,
# y(0) = 10. 

def f(t, y):
    return 4 - 1000*y

def exacta(t):
    return (1/250)*(2499*exp(-1000*t)  +1)



# q = 1 etapa
# t1 = tk + h/2
# y1 = yk + h/2 * f(t1, y1)
#    = yk + h/2*(4-1000*y1)
#    = yk + 4h/2 - 1000*h/2*y1
#    = yk + 2h - 500*h*y1

# y1 = (yk + 2h)/(1 + 500*h)
# y[k+1] = yk + h*f(t1, y1)
def implicito(a, b, fun, N, y0):
    h = (b-a)/N  # paso de malla
    t = zeros(N+1)  # inicializacion del vector de nodos
    y = zeros(N+1)  # inicializacion del vector de resultados
    t[0] = a  # nodo inicial
    y[0] = y0  # valor inicial

    for k in range(N):
        y[k+1] = y[k]+h*fun(t[k], (y[k]+2*h)/(1 + 500*h))
        t[k+1] = t[k]+h

    return (t, y)



malla = [40, 80, 160, 320, 640, 1280]
y0 = 10.
a= 0.
b = 10.

for N in malla:
    tini = perf_counter()           # tiempo inicial

    (t, y) = implicito(a, b, f, N, y0) # llamada al metodo de RK

    tfin=perf_counter()             # tiempo final

    ye = exacta(t) # calculo de la solucion exacta

    # Dibujamos las soluciones
    plot(t, y, '-*') # dibuja la solucion aproximada
    xlabel('t')
    ylabel('y')

    grid(True)

    # Calculo del error cometido
    error = max(abs(y-ye))

    # Resultados
    print('-----')
    print('Tiempo CPU: ' + str(tfin-tini))
    print('Error: ' + str(error))
    print('Paso de malla: ' + str((b-a)/N))
    print('-----')


plot(t, ye, 'k') # dibuja la solucion exacta
leyenda = ['RK, N=' + str(N) for N in malla]
leyenda.append('Exacta')
legend(leyenda)
show()