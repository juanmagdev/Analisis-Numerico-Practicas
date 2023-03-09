# Práctica 2: Métodos RK encajados
Esta carpeta contiene el código y la documentación correspondiente a la Práctica 1 del curso de Análisis Numérico. En esta práctica, se implementan dos métodos encajados:  RK2(3) y RK4(5) en sus versiones para EDOs y para Sistemas de Ecuaciones Diferenciales. 

Además del código fuente de los métodos implementados, esta carpeta también contiene un archivo README.md (este archivo) que proporciona una breve introducción a la práctica, así como un resumen del método Runge Kutta encajado.

## Método RungeKutta encajado (RK2(3))
El método de Runge-Kutta encajado RK2(3) es un método numérico para resolver ecuaciones diferenciales ordinarias (EDO) de primer orden. Como su nombre indica, este método combina dos métodos de Runge-Kutta de segundo orden (RK2) para obtener una solución de tercer orden (RK3).

El método RK2(3) se llama "encajado" porque uno de los métodos de RK2 se utiliza para estimar la solución en el siguiente paso del tiempo, mientras que el otro método de RK2 se utiliza para estimar el error en la estimación. Si el error es mayor que una cierta tolerancia, se repite el paso con un tamaño de paso de tiempo más pequeño hasta que se alcanza la tolerancia deseada.

El método RK2(3) tiene una precisión de tercer orden y es más preciso que los métodos de RK2 convencionales. Sin embargo, el método también es más costoso computacionalmente debido a la necesidad de realizar cálculos adicionales para estimar el error. En general, el método RK2(3) se utiliza cuando se requiere una mayor precisión que los métodos de RK2 convencionales, pero cuando el costo computacional adicional no es un problema.