from time import perf_counter
from matplotlib.pyplot import *
from numpy import *
from numpy.linalg import eig

# x'' + u(x^2 - 1)x' + x = 0
u = 2
def fun(t, y):
    return array([y[1], -y[0] - u*(y[0]**2 - 1)*y[1]])

# condición inicial x(0) = 2, x'(0) = 0
y0 = array([2, 0])
