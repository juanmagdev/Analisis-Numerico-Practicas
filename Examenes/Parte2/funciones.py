from time import perf_counter
from matplotlib.pyplot import *
from numpy import *
from numpy.linalg import eig

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

#####################################################

# dR función que toma un número complejo como entrada y devuelve otro número complejo
def localizar_frontera_RK(dR, N):
    Npoints = 5000
    h = 2*N*pi/Npoints
    z = zeros(Npoints +1 , dtype = complex)
    z[0] = 0
    t = 0
    for k in range(len(z)-1):
        z[k+1] = z[k]+ h*1j*exp(1j*t)/dR(z[k])
        t = t + h
    x, y = real(z), imag(z)
    plot(x, y)
    grid(True)
    axis('equal')

def localizar_frontera(rho, sigma):
    theta = arange(0, 2.*pi, 0.01)
    numer = polyval(rho, exp(theta*1j))
    denom = polyval(sigma, exp(theta*1j))
    mu = numer/denom
    x, y = real(mu), imag(mu)
    plot(x, y)
    grid(True)
    axis('equal')

def euler(a, b, fun, N, y0):
    h = (b-a)/N
    t = zeros(N+1)
    y = zeros(N+1)
    t[0] = a
    y[0] = y0

    for k in range(N):
        t[k+1] = t[k] + h
        y[k+1] = y[k] + h*fun(t[k], y[k])
    
    return (t, y)

def rk4(a, b, fun, N, y0):
    h = (b-a)/N
    t = zeros(N+1)
    y = zeros(N+1)
    t[0] = a 
    y[0] = y0 

    for k in range(N):
        t[k+1] = t[k] + h
        k1 = fun(t[k], y[k])
        k2 = fun(t[k] + h/2, y[k] + h/2 * k1)
        k3 = fun(t[k] + h/2, y[k] + h/2 * k2)
        k4 = fun(t[k+1], y[k] + h*k3)
        y[k+1] = y[k] + h/6 *(k1 + 2*k2 + 2*k3 + k4)
    
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

def AB4(a, b, fun, N, y0):
    y = zeros(N+1)
    t = zeros(N+1)
    f = zeros(N+1)
    t[0] = a
    h = (b-a) / N
    y[0] = y0
    f[0] = fun(a, y[0])

    for k in range(3):
        y[k+1] = y[k] + h*fun(t[k] + h/2, y[k] + h/2 * f[k])
        t[k+1] = t[k] + h
        f[k+1] = fun(t[k+1], y[k+1])
    
    for k in range(3, N):
        y[k+1] = y[k] + h/24 * (55*f[k] - 59*f[k-1] + 37*f[k-2] - 9*f[k-3])
        t[k+1] = t[k] + h
        f[k+1] = fun(t[k+1], y[k+1])
        
    return (t, y)

def AM3(a, b, fun, N, y0):
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

def AB4Sis(a, b, fun, N, y0):
    y = zeros([len(y0), N+1])
    t = zeros(N+1)
    f = zeros([len(y0), N+1])
    t[0] = a
    h = (b-a) / N
    y[:, 0] = y0
    f[:, 0] = fun(a, y[:, 0])

    ## Metodo de Runge-Kutta de orden 4 Sistemas
    for k in range(3):
        k1 = fun(t[k], y[:, k])
        k2 = fun(t[k]+h/2, y[:, k]+h/2*k1)
        k3 = fun(t[k]+h/2, y[:, k]+h/2*k2)
        k4 = fun(t[k]+h, y[:, k]+h*k3)
        t[k+1] = t[k]+h
        y[:, k+1] = y[:, k] + h/6*(k1 + 2*k2 + 2*k3 + k4)
    
    for k in range(3, N):
        y[:, k+1] = y[:, k] + h/24 * (55*f[:, k] - 59*f[:, k-1] + 37*f[:, k-2] - 9*f[:, k-3])
        t[k+1] = t[k] + h
        f[:, k+1] = fun(t[k+1], y[:, k+1])
        
    return (t, y)

def ABM3(a,b,N,y0,fun):
    
    y = zeros(N+1)
    t = zeros(N+1)
    f = zeros(N+1)
    t[0] = a
    h = (b-a)/float(N) 
    y[0] = y0
    f[0] = fun(a,y[0])

    for k in range(2):  ### Arranco con RK4 estandar
        k1 = f[k]
        k2 = fun(t[k]+0.5*h, y[k]+0.5*h*k1)
        k3 = fun(t[k]+0.5*h, y[k]+0.5*h*k2)
        k4 = fun(t[k]+h, y[k]+h*k3)
        y[k+1] = y[k] + h*(k1+2*k2+2*k3+k4)/6
        t[k+1] = t[k]+h
        f[k+1] = fun(t[k+1],y[k+1])
    for k in range(2,N):
        Ck = y[k] +h/24*(19*f[k] - 5*f[k-1] + f[k-2])
        t[k+1] = t[k]+h
       
        
        z0 = y[k]+h/12*(23*f[k] - 16*f[k-1] + 5*f[k-2])
        

        y[k+1] = h*9/24*fun(t[k+1],z0) + Ck
        f[k+1] = fun(t[k+1],y[k+1])
    

    return(t,y)

def ABM3Sis(a,b,fun,N,y0):
    
    tol = 1.e-12
    stop = 200
    
    y = zeros([len(y0),N+1])
    t = zeros(N+1)
    f = zeros([len(y0),N+1])
    t[0] = a
    h = (b-a)/float(N) 
    y[:,0] = y0
    f[:,0] = fun(a,y0)

    for k in range(2):  ### Arranco con RK4 estandar
        k1 = fun(t[k],y[:,k])
        k2 = fun(t[k]+0.5*h, y[:,k]+0.5*h*k1)
        k3 = fun(t[k]+0.5*h, y[:,k]+0.5*h*k2)
        k4 = fun(t[k]+h, y[:,k]+h*k3)
        y[:,k+1] = y[:,k] + h*(k1+2*k2+2*k3+k4)/6
        t[k+1] = t[k]+h
        f[:,k+1] = fun(t[k+1],y[:,k+1])

    for k in range(2,N):
        Ck = y[:,k] +h/24*(19*f[:,k] - 5*f[:,k-1] + f[:,k-2])
        t[k+1] = t[k]+h

               
        z0 = y[:,k]+h/12*(23*f[:,k] - 16*f[:,k-1] + 5*f[:,k-2])
        

        y[:,k+1] = h*9/24*fun(t[k+1],z0) + Ck
        f[:,k+1] = fun(t[k+1],y[:,k+1])
    
    return(t,y)

# Otros
def BDF5sis(a,b,N,y0,fun):
    tol = 1.e-8
    stop = 200
    
    y = zeros([len(y0),N+1])
    t = zeros(N+1)
    f = zeros([len(y0),N+1])
    t[0] = a
    h = (b-a)/float(N) 
    y[:,0] = y0
    f[:,0] = fun(a,y0)
    l = zeros(N+1)
    for k in range(4):  ### Arranco con el metodo RK5 del RK45
        y0 = y0.reshape(len(y0),1)
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
        B[0] = 16/135
        B[2] = 6656/12825
        B[3] = 28561/56430
        B[4] = -9/50
        B[5] = 2/55
        
   
    
        C = zeros(q)
        for i in range(q):
            C[i] = sum(A[i,:]) 
            
        K = zeros([len(y0),q])
        
        for i in range(q):
            K[:,i] = fun(t[k]+C[i]*h, y[:,k]+h*dot(A[i,:],transpose(K)))
    
        incrlow = dot(B,transpose(K)) # metodo de orden 5

        y[:,k+1] =  y[:,k]+h*incrlow

        t[k+1] = t[k]+h
        f[:,k+1] = fun(t[k+1],y[:,k+1]) # Aqui acaba el RK5
            
            
    for k in range(4,N):
        Ck = 5*y[:,k] - 5*y[:,k-1] + 10/3*y[:,k-2] - 5/4*y[:,k-3] +1/5*y[:,k-4]
        t[k+1] = t[k]+h
        
        zold = y[:,k]
        znew = (h*fun(t[k+1],zold) + Ck)*60/137
        
        while(l[k]<stop and max(abs(znew -zold))>=tol):
            zold = znew
            znew = (h*fun(t[k+1],zold) + Ck)*60/137 
            l[k] +=1
        if (l[k] >= stop):
            print('El algoritmo de punto fijo no ha convergido')
        y[:,k+1] = znew
        f[:,k+1] = fun(t[k+1],y[:,k+1])
    
    maxiter = max(l)
    return(t,y,maxiter)

