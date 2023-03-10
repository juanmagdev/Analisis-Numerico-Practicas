# Práctica 2: Métodos Runge-Kutta encajados
Esta carpeta contiene el código y la documentación correspondiente a la Práctica 1 del curso de Análisis Numérico. En esta práctica, se implementan dos métodos encajados:  RK2(3) y RK4(5) en sus versiones para EDOs y para Sistemas de Ecuaciones Diferenciales. 

Además del código fuente de los métodos implementados, esta carpeta también contiene un archivo README.md (este archivo) que proporciona una breve introducción a la práctica, así como un resumen del método Runge Kutta encajado.

## Método Runge-Kutta encajado (RK2(3))
El método de Runge-Kutta encajado RK2(3) es un método numérico para resolver ecuaciones diferenciales ordinarias (EDO) de primer orden. Como su nombre indica, este método combina dos métodos de Runge-Kutta de segundo orden (RK2) para obtener una solución de tercer orden (RK3).

El método RK2(3) se llama "encajado" porque uno de los métodos de RK2 se utiliza para estimar la solución en el siguiente paso del tiempo, mientras que el otro método de RK2 se utiliza para estimar el error en la estimación. Si el error es mayor que una cierta tolerancia, se repite el paso con un tamaño de paso de tiempo más pequeño hasta que se alcanza la tolerancia deseada.

El método RK2(3) tiene una precisión de tercer orden y es más preciso que los métodos de RK2 convencionales. Sin embargo, el método también es más costoso computacionalmente debido a la necesidad de realizar cálculos adicionales para estimar el error. En general, el método RK2(3) se utiliza cuando se requiere una mayor precisión que los métodos de RK2 convencionales, pero cuando el costo computacional adicional no es un problema.


## Generalización del método RKp(p+1)
Idea: Usamos 2 métodos Rk de orden $p$ y $p + 1$:

Una vez calculado $y_k$, para estimar el error:
 - Se calcula las etapas (comunes a los dos métodos):
    $$y_k^{(i)} = y_k + h_k\sum_{j=1}^s a_{ij}f(t_k^{(j)}, y_k^{j}), j= 1, ..., s$$

    $$t_k^{(i)} = t_k + c_ih $$

 - La primera función es: 
  
  $$ \Phi(t_k, y_k, h_k) = \sum_{i = 1}^s b_if(t_k^{(i)}, y_k^{(i)}) $$

 - La segunda función es:

   $$ \Phi^*(t_k, y_k, h_k) = \sum_{i = 1}^s b_i^*f(t_k^{(i)}, y_k^{(i)}) $$

 - Por tanto, el error es:
  
  $$ \~{\epsilon}_k = h_k\sum_{i = 1}^s(b_i^*-b_i)f(t_k^{(i)}, y_k^{(i)}) = h_k(\Phi^*(t_k, y_k, h_k) - \Phi(t_k, y_k, h_k))$$


El nuevo paso a considerar, usando el error estimado seria:

$$ h_{k+1} = \frac{\epsilon}{|\~{\epsilon}_k|}^{\frac{1}{p+1}} $$

De igual forma que el método Runge-Kutta habitual, la solución aproximada viene dada por:

$$ y_{k+1} = y_k + h_k\Phi(t_k, y_k;h_k)$$