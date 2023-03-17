from numpy import *  
from matplotlib.pyplot import *
# from pylab import *
from time import perf_counter

# Ejercicio 1
# a)
g = 9.81
L = 0.25

def f1(t,y):
    f1 = y[1]
    f2 = -g/L * sin(y[0])
    return(array([f1,f2]))

def f2(t,y):
    f1 = y[1]
    f2 = -g/L * y[0]
    return(array([f1,f2]))

def rkSistemas45(a, b, fun, y0, h0, tol):
    
    hmin = 1.e-5 # paso de malla minimo
    hmax = 0.1 # paso de malla maximo

    
    # coeficientes RK
    q = 6 # orden del metodo mas uno
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
        C[i] = sum(A[i,:])
    
    # inicializacion de variables
    t = array([a]) # nodos
    y = y0 # soluciones
    h = array([h0]) # pasos de malla
    K = zeros([len(y0),q])
    k = 0 # contador de iteraciones
    
    
    
    while (t[k] < b):
        h[k] = min(h[k], b-t[k]) # ajuste del ultimo paso de malla
        for i in range(q):
            K[:,i] = fun(t[k]+C[i]*h[k], y[:,k]+h[k]*dot(A[i,:],transpose(K)))
        
        incrlow = dot(B,transpose(K)) # metodo de orden 4
        incrhigh = dot(BB,transpose(K)) # metodo de orden 5
            
        error = linalg.norm(h[k]*(incrhigh-incrlow),inf) # estimacion del error
        y = column_stack((y, y[:,k]+h[k]*incrlow))
        t = append(t, t[k]+h[k]); # t_(k+1)
        hnew = 0.9*h[k]*abs(tol/error)**(1./5) # h_(k+1)
        hnew = min(max(hnew,hmin),hmax) # hmin <= h_(k+1) <= hmax
        h = append(h, hnew)
        k += 1
        
    return (t, y, h)

# Datos del problema
a = 0 # extremo inferior del intervalo
b = 4 # extremo superior del intervalo
y0 = array([pi/6,0]) # condicion inicial
y0 = y0.reshape(2,1)
h0 = 0.03 #paso inicial
tol = 1.e-6 #tolerancia



tini = perf_counter()
(t1, y1, h1) = rkSistemas45(a, b, f1, y0, h0, tol) # llamada al metodo RK4(5)
(t2, y2, h2) = rkSistemas45(a, b, f2, y0, h0, tol) # llamada al metodo RK4(5)
tfin = perf_counter()

figure('Apartado a')
subplot(311)
plot(t1,y1[0],t2,y2[0]) 
legend(['(1)','(3)'])
subplot(312)
plot(y1[0],y1[1], y2[0],y2[1])
legend(['(1)','(3)'])
subplot(313)
plot(t1[:-1],h1[:-1],t2[:-1],h2[:-1])# se excluye el ultimo valor de h porque no se usa para avanzar
legend(['(1)','(3)'])

gcf().suptitle("Ejercicio 1a")
show()