# Práctica 3: Métodos de Aproximación de Ecuaciones Diferenciales

En esta práctica, exploraremos varios métodos numéricos utilizados para aproximar ecuaciones diferenciales. En particular, nos centraremos en los métodos de Adams-Bashforth de dos pasos, métodos de Adams-Moulton, métodos Predictor-Corrector y una generalización para tres pasos.

## Métodos de Adams-Bashforth de dos pasos

Los métodos de Adams-Bashforth de dos pasos son fórmulas explícitas para la aproximación de ecuaciones diferenciales. Estos métodos utilizan información de los puntos anteriores para estimar el siguiente valor. La fórmula general para los métodos de Adams-Bashforth de dos pasos es la siguiente:

$$y_{n+1} = y_n + h * (3/2 * f(t_n, y_n) - 1/2 * f(t_{n-1}, y_{n-1})) $$

donde $y_n$ y $y_{n-1}$ son los valores aproximados de la solución en los pasos anteriores, $h$ es el tamaño del paso y $f(t, y)$ es la función que define la ecuación diferencial.

**Ventajas:**
- Son métodos explícitos y fáciles de implementar.
- No requieren la solución de sistemas de ecuaciones no lineales.

**Desventajas:**
- Son menos precisos en comparación con los métodos implícitos.
- Pueden volverse inestables para ecuaciones con coeficientes grandes.

## Métodos de Adams-Moulton

Los métodos de Adams-Moulton son fórmulas implícitas utilizadas para aproximar ecuaciones diferenciales. Estos métodos utilizan información de los puntos anteriores, incluyendo el punto actual, para estimar el siguiente valor. La fórmula general para los métodos de Adams-Moulton es la siguiente:

$$y_{n+1} = y_n + h * (α * f(t_{n+1}, y_{n+1}) + β * f(t_n, y_n))$$

donde $y_{n+1}$ es el valor aproximado de la solución en el siguiente paso, $α$ y $β$ son coeficientes determinados por el método específico y las condiciones iniciales, $h$ es el tamaño del paso y $f(t, y)$ es la función que define la ecuación diferencial.

**Ventajas:**
- Proporcionan una mayor precisión que los métodos explícitos.
- Son adecuados para ecuaciones con coeficientes grandes.

**Desventajas:**
- Requieren la solución de sistemas de ecuaciones no lineales en cada paso.
- Son más complejos de implementar que los métodos explícitos.

## Métodos Predictor-Corrector
Los métodos Predictor-Corrector combinan un método predictor y un método corrector para aproximar ecuaciones diferenciales. El método predictor proporciona una estimación inicial del siguiente valor, y luego se utiliza un método corrector para mejorar esta estimación. La fórmula general para los métodos Predictor-Corrector es la siguiente:

$$y_{n+1}^* = y_n + h * α * f(t_n, y_n)$$
$$y_{n+1} = y_n + h * [(α/2) * (f(t_n, y_n) + f(t_{n+1}, y_{n+1}^*)]$$

donde $y_{n+1}$ es el valor aproximado de la solución en el siguiente paso, α es un coeficiente determinado por el método específico, h es el tamaño del paso, $f(t, y)$ es la función que define la ecuación diferencial y $y_{n+1}^*$ es la estimación inicial obtenida por el método predictor.

**Ventajas**:
 - Proporcionan una mayor precisión que los métodos explícitos.
 - Son menos propensos a la inestabilidad que los métodos explícitos.
 - No requieren la solución de sistemas de ecuaciones no lineales en cada paso.

**Desventajas**:

 - Pueden requerir un mayor costo computacional debido al cálculo adicional del método corrector.
 - La precisión puede verse afectada por la precisión del método predictor utilizado.

## Generalización para tres pasos
Una generalización de los métodos de Adams-Bashforth y Adams-Moulton para tres pasos puede lograrse utilizando la información de los dos pasos anteriores para estimar el siguiente valor. La fórmula general para esta generalización es la siguiente:

### Método de Adams-Bashforth de tres pasos:
$$y_{n+1} = y_n + h * (23/12 * f(t_n, y_n) - 4/3 * f(t_{n-1}, y_{n-1}) + 5/12 * f(t_{n-2}, y_{n-2}))$$

### Método de Adams-Moulton de tres pasos:
$$y_{n+1} = y_n + h * (5/12 * f(t_{n+1}, y_{n+1}) + 2/3 * f(t_n, y_n) - 1/12 * f(t_{n-1}, y_{n-1}))$$

Estas generalizaciones proporcionan una mayor precisión al utilizar más puntos anteriores en la aproximación de la solución.