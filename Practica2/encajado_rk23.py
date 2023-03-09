# --------- Ejercicio 1 ----------

# Resolucion del problema de valor inicial
# y'=f(t,y), y(t0)=y0,
# mediante el metodo RK2(3).

from doctest import NORMALIZE_WHITESPACE
from pylab import *
from time import perf_counter
from matplotlib.pyplot import *
from numpy import *

def f(t, y):
    """Funcion que define la ecuacion diferencial"""
    return -y+exp(-t)*cos(t)

def exacta(t):
    """Solucion exacta del problema de valor inicial"""
    return exp(-t)*sin(t)

def rk23(a, b, fun, y0, h0, tol):
    """Implementacion del metodo encajado RK2(3)
    en el intervalo [a, b] con condicion inicial y0,
    paso inicial h0 y tolerancia tol"""
    

    hmin = (b-a)*1.e-5 # paso de malla minimo
    hmax = (b-a)/10. # paso de malla maximo

    
    # coeficientes RK
    q = 3 # numero de etapas
    p = 2 # orden del método menos preciso
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
        C[i] = sum(A[i,:])
    
    # inicializacion de variables
    t = array([a]) # nodos
    y = array([y0]) # soluciones
    h = array([h0]) # pasos de malla
    K = zeros(3)
    k = 0 # contador de iteraciones
    
    while (t[k] < b):
        h[k] = min(h[k], b-t[k]) # ajuste del ultimo paso de malla
        for i in range(q):
            K[i] = fun(t[k]+C[i]*h[k], y[k]+h[k]*sum(A[i,:]*K))
        incrlow = sum(B*K) # metodo de orden 2
        incrhigh = sum(BB*K) # metodo de orden 3
        error = h[k]*(incrhigh-incrlow) # estimacion del error
        y = append(y, y[k]+h[k]*incrlow) # y_(k+1)
        t = append(t, t[k]+h[k]); # t_(k+1)
        hnew = 0.9*h[k]*abs(tol/error)**(1./(p+1)) # h_(k+1)
        hnew = min(max(hnew,hmin),hmax) # hmin <= h_(k+1) <= hmax
        h = append(h, hnew)
        k += 1
    

    return (t, y, h)


    
# Datos del problema
a = 0. # extremo inferior del intervalo
b = 30. # extremo superior del intervalo
y0 = 0. # condicion inicial
h0 = 0.1 #paso inicial
tol = 1.e-6 #tolerancia


tini = perf_counter()
(t, y, h) = rk23(a, b, f, y0, h0, tol) # llamada al metodo RK2(3)
tfin = perf_counter()


# calculo de la solucion exacta
te = linspace(a,b,200)
ye = exacta(te) 

# Dibujamos las soluciones
subplot(211)
plot(t, y, 'bo-') # dibuja la solucion aproximada
plot(te, ye, 'r') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['RK2(3)', 'exacta'])
grid(True)
subplot(212)
plot(t[:-1],h[:-1],'*-')# se excluye el ultimo valor de h porque no se usa para avanzar
xlabel('t')
ylabel('h')
legend(['pasos usados'])
gcf().suptitle("Ejemplo")
show()
# Calculo del error cometido
error = max(abs(y-exacta(t)))
hn = min(h[:-2]) # minimo valor utilizado del paso de malla
hm = max(h[:-2]) # maximo valor utilizado del paso de malla
# se elimina el ultimo h calculado porque no se usa y el penúltimo porque se ha ajustado para terminar en b

# Resultados
print('-----')
print("Error: " + str(error))
print("No. de iteraciones: " + str(len(y)))
print('Tiempo CPU: ' + str(tfin-tini))
print("Paso de malla minimo: " + str(hn))
print("Paso de malla maximo: " + str(hm))
print(sum(h))
print('-----')






# --------- Ejercicio 2 ----------
def f(t, y):
    """Funcion que define la ecuacion diferencial"""
    return -y+2*sin(t)

def exacta(t):
    """Solucion exacta del problema de valor inicial"""
    return (pi+1)*exp(-t) + sin(t)-cos(t)



def rk45(a, b, fun, y0, h0, tol):
    """Implementacion del metodo encajado RK2(3)
    en el intervalo [a, b] con condicion inicial y0,
    paso inicial h0 y tolerancia tol"""
    

    hmin = (b-a)*1.e-5 # paso de malla minimo
    hmax = (b-a)/10. # paso de malla maximo

    # definicion de los coeficientes del metodo (tablero)
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
    y = array([y0]) # soluciones
    h = array([h0]) # pasos de malla
    K = zeros(6)
    k = 0 # contador de iteraciones
    
    while (t[k] < b):
        h[k] = min(h[k], b-t[k]) # ajuste del ultimo paso de malla
        for i in range(q):
            K[i] = fun(t[k]+C[i]*h[k], y[k]+h[k]*sum(A[i,:]*K))
        incrlow = sum(B*K) # metodo de orden 4
        incrhigh = sum(BB*K) # metodo de orden 5
        error = h[k]*(incrhigh-incrlow) # estimacion del error
        y = append(y, y[k]+h[k]*incrlow) # y_(k+1)
        t = append(t, t[k]+h[k]); # t_(k+1)
        hnew = 0.9*h[k]*abs(tol/error)**(1/5) # h_(k+1)
        hnew = min(max(hnew,hmin),hmax) # hmin <= h_(k+1) <= hmax
        h = append(h, hnew)
        k += 1
    
    # print("Numero de iteraciones: "+str(k))
    return (t, y, h)




# Datos del problema
a = 0. # extremo inferior del intervalo
b = 10. # extremo superior del intervalo
y0 = pi # condicion inicial
h0 = (b-a)/2000 #paso inicial
tol = 1.e-6 #tolerancia


tini = perf_counter()
(t, y, h) = rk45(a, b, f, y0, h0, tol) # llamada al metodo RK2(3)
tfin = perf_counter()


# calculo de la solucion exacta
te = linspace(a,b,200)
ye = exacta(te) 

# Dibujamos las soluciones
subplot(211)
plot(t, y, 'bo-') # dibuja la solucion aproximada
plot(te, ye, 'r') # dibuja la solucion exacta
xlabel('t')
ylabel('y')
legend(['RK4(5)', 'exacta'])
grid(True)
subplot(212)
plot(t[:-1],h[:-1],'*-')# se excluye el ultimo valor de h porque no se usa para avanzar
xlabel('t')
ylabel('h')
legend(['pasos usados'])
gcf().suptitle("Ejercicio 2")
show()
# Calculo del error cometido
error = max(abs(y-exacta(t)))
hn = min(h[:-2]) # minimo valor utilizado del paso de malla
hm = max(h[:-2]) # maximo valor utilizado del paso de malla
# se elimina el ultimo h calculado porque no se usa y el penúltimo porque se ha ajustado para terminar en b

# Resultados
print('-----')
print("Error: " + str(error))
print("No. de iteraciones: " + str(len(y)))
print('Tiempo CPU: ' + str(tfin-tini))
print("Paso de malla minimo: " + str(hn))
print("Paso de malla maximo: " + str(hm))
print(sum(h))
print('-----')





# --------- Ejercicio 3 ----------
def f(t,y):
    f1= 3*y[0] - 2*y[1]
    f2= -y[0] +3*y[1] -2*y[2]
    f3= -y[1] +3*y[2]
    return array([f1,f2,f3])

def exacta(t):
    f1=-1/4 *exp(5*t) +3/2 *exp(3*t) -1/4 * exp(t)
    f2=1/4 *exp(5*t) -1/4 * exp(t)
    f3=-1/8 *exp(5*t) -3/4 *exp(3*t) -1/8 * exp(t)
    return(array([f1,f2,f3]))

def rkSistemas23(a, b, fun, y0, h0, tol):
    
    hmin = (b-a)*1.e-5 # paso de malla minimo
    hmax = (b-a)/10. # paso de malla maximo

    
    # coeficientes RK
    q = 3 # orden del metodo mas uno
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
        C[i] = sum(A[i,:])
    
    # inicializacion de variables
    t = array([a]) # nodos
    y = y0 # soluciones
    h = array([h0]) # pasos de malla
    K = zeros([len(y0),3])
    k = 0 # contador de iteraciones
    
    while (t[k] < b):
        h[k] = min(h[k], b-t[k]) # ajuste del ultimo paso de malla
        for i in range(3):
            K[:,i] = fun(t[k]+C[i]*h[k], y[:,k]+h[k]*dot(A[i,:],transpose(K)))
        
        incrlow = dot(B,transpose(K)) # metodo de orden 2
        incrhigh = dot(BB,transpose(K)) # metodo de orden 3
            
        error = linalg.norm(h[k]*(incrhigh-incrlow),inf) # estimacion del error
        y = column_stack((y, y[:,k]+h[k]*incrlow))
        t = append(t, t[k]+h[k]); # t_(k+1)
        hnew = 0.9*h[k]*abs(tol/error)**(1./q) # h_(k+1)
        hnew = min(max(hnew,hmin),hmax) # hmin <= h_(k+1) <= hmax
        h = append(h, hnew)
        k += 1

    return (t, y, h)




# Datos del problema
a = 0. # extremo inferior del intervalo
b = 1. # extremo superior del intervalo
y0 = array([1,0,-1]) # condicion inicial
y0 = y0.reshape(3,1)
h0 = (b-a)/2000 #paso inicial
tol = 1.e-6 #tolerancia


tini = perf_counter()
(t, y, h) = rkSistemas23(a, b, f, y0, h0, tol) # llamada al metodo RK2(3)
tfin = perf_counter()


# calculo de la solucion exacta
te = linspace(a,b,200)
ye = exacta(te) 


subplot(211)
plot(t,y[0],t,y[1],t,y[2]) 
subplot(212)
plot(y[0],y[1])
print('-----')
# Calculo del error cometido
error1= max(abs(y[0]-exacta(t)[0]))
error2= max(abs(y[1]-exacta(t)[1]))
error3= max(abs(y[2]-exacta(t)[2]))
hn = min(h[:-2]) # minimo valor utilizado del paso de malla
hm = max(h[:-2]) # maximo valor utilizado del paso de malla
# se elimina el ultimo h calculado porque no se usa y el penúltimo porque se ha ajustado para terminar en b

gcf().suptitle("Ejercicio 3")
show()

# Resultados
print('-----')

print("No. de iteraciones: " + str(len(y[0])-1))
print(error1)
print(error2)
print(error3)
print('Tiempo CPU: ' + str(tfin-tini))
print("Paso de malla minimo: " + str(hn))
print("Paso de malla maximo: " + str(hm))
print(sum(h))
print('-----')