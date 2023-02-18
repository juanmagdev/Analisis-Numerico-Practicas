# Práctica 1: Métodos unipaso para problemas de valor inicial

Esta carpeta contiene el código y la documentación correspondiente a la Práctica 1 del curso de Análisis Numérico. En esta práctica, se implementan diversos métodos numéricos unipaso para la resolución de problemas de valor inicial (PVI) de ecuaciones diferenciales ordinarias (EDOs).

Los métodos implementados en esta práctica incluyen el método de Euler, el método de Euler mejorado (también conocido como método de Heun), el método de Runge-Kutta de segundo orden y el método de Runge-Kutta de cuarto orden.

Además del código fuente de los métodos implementados, esta carpeta también contiene un archivo practica1.md (este archivo) que proporciona una breve introducción a la práctica, así como un resumen de los métodos utilizados.

El objetivo de esta práctica es que el estudiante adquiera conocimientos prácticos sobre la implementación y el uso de métodos numéricos para la resolución de EDOs, así como sobre la interpretación y análisis de los resultados obtenidos.



## Métodos implementados

Los métodos implementados en esta práctica son los siguientes:

### Método de Euler

El método de Euler es uno de los métodos numéricos más simples para la resolución de EDOs. En este método, se aproxima la solución en cada punto como la solución anterior más el producto del tamaño del paso y la derivada en ese punto. La fórmula matemática del método de Euler es la siguiente:

$$y_{n+1} = y_n + hf(t_n,y_n)$$

Este método es de orden 1 y tiene un error de truncamiento local proporcional al cuadrado del tamaño del paso. Es un método útil para entender el comportamiento cualitativo de las soluciones, pero no es muy preciso y puede requerir tamaños de paso muy pequeños para obtener resultados precisos.

### Método de Taylor

El método de Taylor es un método numérico para la resolución de EDOs que se basa en la expansión de la solución en series de Taylor alrededor de un punto dado. En este método, se utiliza la serie de Taylor truncada para aproximar la solución en cada punto. La fórmula matemática del método de Taylor de orden $k$ es la siguiente:

$$y_{n+1} = y_n + h\sum_{i=1}^{k}\frac{f^{(i-1)}(t_n,y_n)}{i!}$$

Este método es de orden $k$ y tiene un error de truncamiento local proporcional al tamaño del paso elevado a la potencia $k+1$. Es un método más preciso que el método de Euler, pero requiere más cálculos por paso.

### Método de Runge-Kutta de orden 2

El método de Runge-Kutta de orden 2 (también conocido como método del punto medio) es un método numérico para la resolución de EDOs que utiliza dos evaluaciones de la función $f$. En este método, se aproxima la solución en cada punto como la solución anterior más el producto del tamaño del paso y una combinación lineal de las dos evaluaciones de la función $f$ en ese punto y en el punto medio del intervalo. La fórmula matemática del método de Runge-Kutta de orden 2 es la siguiente:

$$k_1 = f(t_n,y_n)$$

$$k_2 = f(t_n + h/2,y_n + hk_1/2)$$

$$y_{n+1} = y_n + hk_2$$

Este método es de orden 2 y tiene un error de truncamiento local proporcional al tamaño del paso elevado al cubo. Es un método más preciso que el método de Euler y requiere menos cálculos por paso que el método de Taylor.

### Método de Heun
El método de Heun, también conocido como método del punto medio mejorado, es una extensión del método del punto medio que utiliza la información sobre la derivada en ambos extremos del intervalo de tiempo para mejorar la precisión de la aproximación. La idea principal del método de Heun es aproximar el valor de la función en el punto medio del intervalo de tiempo utilizando una aproximación lineal de la función. Luego se utiliza este valor para estimar la pendiente de la función en ese punto y se utiliza esta pendiente para obtener una aproximación del valor de la función al final del intervalo de tiempo. La fórmula matemática del método de Heun es la siguiente:

$$\begin{aligned}
k_1 &= f(t_n, y_n), \\
k_2 &= f(t_n + h, y_n + hk_1), \\
y_{n+1} &= y_n + \frac{h}{2}(k_1 + k_2).
\end{aligned}$$

En esta fórmula, $k_1$ es la pendiente de la función en el inicio del intervalo de tiempo, $k_2$ es la pendiente de la función en el punto medio del intervalo de tiempo, $h$ es el tamaño del paso de tiempo y $y_n$ y $y_{n+1}$ son los valores aproximados de la función en los puntos inicial y final del intervalo de tiempo, respectivamente. El método de Heun es un método de segundo orden, lo que significa que su error local es proporcional a $h^3$, y su error global es proporcional a $h^2$.

### Método de Runge-Kutta de orden 4 (RK4)
El método de Runge-Kutta de orden 4 (RK4) es uno de los métodos numéricos más populares para la solución de ecuaciones diferenciales ordinarias (EDO) de primer orden. Este método utiliza cuatro evaluaciones de la función en cada paso de tiempo para obtener una aproximación más precisa del valor de la función en el siguiente punto de tiempo. La fórmula matemática del método de RK4 es la siguiente:

$$\begin{aligned}
k_1 &= f(t_n, y_n), \\
k_2 &= f(t_n + \frac{h}{2}, y_n + \frac{h}{2}k_1), \\
k_3 &= f(t_n + \frac{h}{2}, y_n + \frac{h}{2}k_2), \\
k_4 &= f(t_n + h, y_n + hk_3), \\
y_{n+1} &= y_n + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4).
\end{aligned}$$

En esta fórmula, $k_1$, $k_2$, $k_3$ y $k_4$ son las constantes de pendiente que se utilizan para aproximar la solución en el intervalo de tiempo. El valor de $y_{n+1}$ se calcula a partir de la combinación lineal ponderada de estas pendientes. RK4 es un método de cuarto orden, lo que significa que su error local es proporcional a $h^5$, y su error global es proporcional a $h^4$.

En resumen, el método de RK4 es un método numérico de alto orden y es más preciso que los métodos de Euler, Heun y Runge-Kutta de segundo orden. Sin embargo, este método también es más costoso computacionalmente debido a las cuatro evaluaciones de la función en cada paso de tiempo.

