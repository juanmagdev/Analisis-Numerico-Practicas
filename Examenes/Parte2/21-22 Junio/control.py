from time import perf_counter
from matplotlib.pyplot import *
from numpy import *
from numpy.linalg import eig



# y'' + 4y' + 29y = 0
def f(t, y):
    return array([y[1], -4*y[1] - 29*y[0]])

# y(0) = 1, y'(0) = -2  
def exacta(t):
    return exp(-2*t)*cos(5*t)

def BDF2_sistema(a, b, fun, N, y0):
    y = zeros([len(y0), N+1])
    t = zeros(N+1)
    f = zeros([len(y0), N+1])
    iteraciones = 0
    t[0] = a
    h = (b-a) / N
    y[:, 0] = y0
    f[:, 0] = fun(a, y[:, 0])
    ystar = zeros([len(y0), 1])
    
    # Metodo de Heun
    for k in range(2):
        ystar[:, 0] = y[:, k] + h*fun(t[k], y[:, k])
        t[k+1] = t[k]+h
        y[:, k+1] = y[:, k] + h*1/2 * \
            (fun(t[k], y[:, k]) + fun(t[k+1], ystar[:, 0]))
    
    for k in range(2, N):
        z = y[:, k]
        contador = 0
        distancia = 1 + 1e-12
        t[k+1] = t[k] + h

        while distancia >= 1e-12 and contador < 200:
            z_nuevo = (4*y[:, k] - y[:, k-1] + 2*h*f[:, k])/3
            distancia = max(abs(z - z_nuevo))
            z = z_nuevo
            contador += 1
        if contador == 200:
            print("El método no converge.")
        iteraciones = max(iteraciones, contador)
        y[:, k+1] = z
        f[:, k+1] = fun(t[k+1], y[:, k+1])

    return (t, y, iteraciones)

# Datos del problema
a = 0
b = 5

N = 100
y0 = array([1, -2])

# Llamada al metodo de BDF2
tini = perf_counter()           # tiempo inicial

(t, y, i) = BDF2_sistema(a, b, f, N, y0) # llamada al metodo de RK4

tfin=perf_counter()             # tiempo final

ye = exacta(t) # calculo de la solucion exacta

# Dibujamos las soluciones
plot(t, y[0], '-*') # dibuja la solucion aproximada
xlabel('t')
ylabel('y')
legend(['BDF2', 'exacta'])
grid(True)

# Calculo del error cometido
error = max(abs(y[0]-ye))

# Resultados
print('-----')
print('Tiempo CPU: ' + str(tfin-tini))
print('Error: ' + str(error))
print('Paso de malla: ' + str((b-a)/N))
print('Iteraciones maximas: ' + str(i))
print('-----')

plot(t, ye, 'k') # dibuja la solucion exacta

gcf().suptitle("Ejercicio 3")
show() # muestra la grafica


# Ejercicio 4
# Lo hago solo para el metodo multipaso de BDF2


def localizar_frontera(rho, sigma):
    theta = arange(0, 2.*pi, 0.01)
    numer = polyval(rho, exp(theta*1j))
    denom = polyval(sigma, exp(theta*1j))
    mu = numer/denom
    x, y = real(mu), imag(mu)
    plot(x, y)
    grid(True)
    axis('equal')

A = array([[0, 1], [-29, -4]])
print("Matriz A:\n", A)
print("Autovalores de A:", eig(A)[0])

# El metodo es convergente si los autovalores de A estan en el semiplano izquierdo

localizar_frontera(array([1, -1, 0, 0, 0]), array([0, 0, 0, 0, 2/3]))

rel, iml = -2, 5
# Obtenemos h para llevar las raíces a la frontera
h = 0.1
plot(h*rel, h*iml, ".", h*rel, -h*iml, ".")
show()

