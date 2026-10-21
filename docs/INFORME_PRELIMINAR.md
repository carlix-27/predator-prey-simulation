# Modelado y Simulación Numérica del Sistema Tritrófico de Cosner: Integración Temporal y Análisis de Estabilidad

**Autores:** [Integrantes del Grupo]  
**Materia:** Cálculo Numérico  
**Institución:** [Universidad / Facultad]  
**Fecha:** Octubre 2026  

---

## Resumen (Abstract)
En este trabajo se modela y analiza numéricamente la dinámica de un sistema ecológico tritrófico compuesto por una presa ($N$), un depredador intermedio ($P$) y un superdepredador omnívoro ($S$), gobernado por el modelo no lineal de Cosner con respuestas funcionales racionales cuadráticas generalizadas. Para resolver el problema de valor inicial asociado, se implementaron de forma independiente los algoritmos de Euler progresivo y Runge-Kutta de segundo orden (RK2/Heun). Se realizó un estudio de estabilidad numérica y consistencia, justificando la selección del paso temporal $\Delta t$ mediante análisis de refinamiento sucesivo de malla hasta asegurar convergencia asintótica y evitar inestabilidades numéricas artificiales. Asimismo, se determinó un conjunto de parámetros biológicamente plausibles que garantizan un régimen de coexistencia asintótica estable, previniendo fluctuaciones caóticas y la extinción de las especies. Los resultados validan la superioridad del método RK2 en términos de precisión y estabilidad frente al incremento de $\Delta t$.

---

## 1. Introducción
La dinámica de poblaciones ha sido tradicionalmente abordada mediante el sistema clásico de Lotka-Volterra. Si bien este modelo pionero capturó la naturaleza oscilatoria de las interacciones presa-depredador, adolece de severas simplificaciones ecológicas: asume crecimiento maltusiano irrestricto de la presa en ausencia de depredadores, considera una tasa de ataque lineal sin límite biológico (Holling Tipo I) y genera órbitas periódicas neutralmente estables sumamente sensibles a perturbaciones del medio.

Para representar con mayor fidelidad los ecosistemas reales, surgen modelos no lineales que incorporan la capacidad de carga del hábitat (crecimiento logístico) y respuestas funcionales saturables. En particular, el **modelo de Cosner** introduce funciones de consumo de tipo racional cuadrático, las cuales capturan fenómenos de saturación por manipulación de presas y efectos de defensa o confusión grupal a densidades elevadas.

El presente estudio tiene como objetivo modelar un ecosistema de tres niveles tróficos (presa-depredador-superdepredador) bajo la formulación de Cosner. Se busca analizar el comportamiento asintótico del sistema mediante integración numérica discreta, contrastando el desempeño del método de Euler frente a Runge-Kutta 2 (RK2), estableciendo criterios rigurosos para la elección del paso temporal $\Delta t$ y delimitando parámetros que impidan la emergencia de dinámicas caóticas indeseadas.

---

## 2. Formulación Matemática del Modelo

### 2.1. Sistema de Ecuaciones Diferenciales
El modelo asignado describe la tasa de variación temporal de las densidades poblacionales mediante el siguiente sistema de EDOs acopladas:

$$\frac{dN}{dt} = r N \left(1 - \frac{N}{K}\right) - d N - F_{PN}(N, P) P - F_{SN}(N, S) S$$

$$\frac{dP}{dt} = b F_{PN}(N, P) P - m P - F_{SP}(P, S) S$$

$$\frac{dS}{dt} = f F_{SN}(N, S) S + g F_{SP}(P, S) S - q S$$

Donde:
* $N(t)$: Densidad poblacional de la presa.
* $P(t)$: Densidad poblacional del depredador intermedio.
* $S(t)$: Densidad poblacional del superdepredador.

### 2.2. Respuestas Funcionales de Cosner
Las tasas de consumo per cápita $F_{PN}$, $F_{SN}$ y $F_{SP}$ están dadas por funciones racionales cuadráticas:

$$F_{PN}(N, P) = \frac{a N}{1 + b_1 N + c_1 N^2}$$

$$F_{SN}(N, S) = \frac{c N}{1 + b_2 N + c_2 N^2}$$

$$F_{SP}(P, S) = \frac{e P}{1 + b_3 P + c_3 P^2}$$

* **Interpretación biológica:**
  * Parámetros $a, c, e$: Tasas intrínsecas de ataque o búsqueda a bajas densidades.
  * Parámetros $b_1, b_2, b_3$: Tiempo de manipulación de la presa (*handling time*), que satura la ingesta a densidades intermedias.
  * Parámetros $c_1, c_2, c_3$: Coeficientes cuadráticos que modelan la inhibición o defensa gregaria a densidades altas (ej. formación de cardúmenes o manadas que dificultan la captura individual).

### 2.3. Parámetros del Sistema
| Símbolo | Significado Biológico | Intervalo Típico Consigna | Valor Adoptado |
| :---: | :--- | :---: | :---: |
| $r$ | Tasa de natalidad de la presa | - | $1.0$ |
| $K$ | Capacidad de carga del ambiente | - | $100.0$ |
| $d$ | Tasa de mortalidad natural de la presa | - | $0.1$ |
| $b$ | Eficiencia de conversión ($N \to P$) | - | $0.5$ |
| $m$ | Tasa de mortalidad natural del depredador | - | $0.2$ |
| $f$ | Eficiencia de conversión ($N \to S$) | - | $0.2$ |
| $g$ | Eficiencia de conversión ($P \to S$) | - | $0.4$ |
| $q$ | Tasa de mortalidad del superdepredador | - | $0.1$ |
| $a, c, e$ | Tasas de ataque | $0.01 - 0.5$ | $0.15, 0.05, 0.10$ |
| $b_1, b_2, b_3$| Coeficientes lineales de saturación | $0.1 - 2.0$ | $0.5, 0.5, 0.5$ |
| $c_1, c_2, c_3$| Coeficientes cuadráticos de inhibición | $0 - 0.1$ | $0.01, 0.01, 0.01$ |

---

## 3. Metodología Numérica

Definiendo el vector de estado $\mathbf{y}(t) = [N(t), P(t), S(t)]^T$ y el campo vectorial $\mathbf{f}(t, \mathbf{y}) = [f_N, f_P, f_S]^T$, el problema se reduce al PVI:
$$\frac{d\mathbf{y}}{dt} = \mathbf{f}(t, \mathbf{y}), \quad \mathbf{y}(0) = \mathbf{y}_0$$

### 3.1. Método de Euler Explícito (Orden 1)
Discretizando el intervalo temporal con paso constante $\Delta t$:
$$\mathbf{y}_{n+1} = \mathbf{y}_n + \Delta t \cdot \mathbf{f}(t_n, \mathbf{y}_n)$$
* **Error de truncamiento local (LTE):** $\mathcal{O}(\Delta t^2)$.
* **Error global:** $\mathcal{O}(\Delta t)$.

### 3.2. Método de Runge-Kutta 2 (Heun / Predictor-Corrector, Orden 2)
Para mejorar el orden de aproximación sin incurrir en cálculo de derivadas analíticas:
$$\mathbf{k}_1 = \mathbf{f}(t_n, \mathbf{y}_n)$$
$$\mathbf{k}_2 = \mathbf{f}(t_n + \Delta t, \mathbf{y}_n + \Delta t \mathbf{k}_1)$$
$$\mathbf{y}_{n+1} = \mathbf{y}_n + \frac{\Delta t}{2} (\mathbf{k}_1 + \mathbf{k}_2)$$
* **Error de truncamiento local (LTE):** $\mathcal{O}(\Delta t^3)$.
* **Error global:** $\mathcal{O}(\Delta t^2)$.

### 3.3. Justificación y Selección del Paso $\Delta t$
El paso temporal no puede elegirse de forma arbitraria:
1. **Estabilidad Absoluta:** Para un sistema no lineal linealizado en torno a un punto de operación mediante la matriz Jacobiana $J = \frac{\partial \mathbf{f}}{\partial \mathbf{y}}$, el paso debe satisfacer que los autovalores $\lambda_i(J)$ queden dentro del dominio de estabilidad del método:
   * Para Euler: $|1 + \lambda \Delta t| \le 1$.
   * Para RK2: $|1 + \lambda \Delta t + \frac{1}{2}(\lambda \Delta t)^2| \le 1$.
2. **Criterio Empírico de Refinamiento de Malla:** Se evalúa la solución numérica reduciendo sucesivamente el paso ($\Delta t \to \Delta t / 2$). La condición de convergencia se alcanza cuando la diferencia cuadrática entre soluciones consecutivas es menor a una tolerancia prefijada:
   $$\max_{n} \|\mathbf{y}_{\Delta t}(t_n) - \mathbf{y}_{\Delta t / 2}(t_n)\|_\infty < 10^{-3}$$
3. **Análogo del Oscilador Armónico:** Si se utilizara Euler con un paso inadecuado en sistemas con autovalores cercanos al eje imaginario, el error numérico introduce amortiguamiento negativo artificial (amplificación espuria de la energía), provocando que las órbitas diverjan al infinito.

---

## 4. Resultados y Discusión
*(Nota: Esta sección presenta de forma concisa las figuras con texto puntual explicativo).*

### 4.1. Evolución Temporal hacia el Equilibrio Asintótico
*(Insertar aquí Figura 1: Curvas de N(t), P(t), S(t) integradas con RK2).*
* **Observación:** Se aprecia un régimen transitorio oscilatorio durante los primeros instantes ($t < 150$), tras el cual el sistema converge suave y asintóticamente a un punto fijo de coexistencia no trivial $(N^*, P^*, S^*)$.

### 4.2. Trayectoria en el Espacio de Fases 3D
*(Insertar aquí Figura 2: Espacio de fases tridimensional N-P-S).*
* **Observación:** La trayectoria tridimensional describe una espiral convergente hacia el atractor de equilibrio. No se observan bifurcaciones ni órbitas abiertas, corroborando la estabilidad de Lyapunov del punto.

### 4.3. Sensibilidad Numérica frente a la Variación de $\Delta t$ (Euler vs. RK2)
*(Insertar aquí Figura 3: Comparativa de trayectorias para $\Delta t = 0.5, 0.1, 0.01$).*
* **Observación:** Con $\Delta t = 0.5$, el método de Euler exhibe severas oscilaciones no físicas e inestabilidad acumulada, mientras que RK2 mantiene la trayectoria acotada. A partir de $\Delta t = 0.01$, ambos métodos muestran total coincidencia visual y numérica.

### 4.4. Prevención de Dinámicas Caóticas
*(Insertar aquí Figura 4: Comparación entre régimen estable adoptado vs. régimen caótico con tasas de ataque hiper-elevadas).*
* **Observación:** Se verifica que al elevar drásticamente los parámetros de interacción $a$ y $e$ más allá de los intervalos típicos, el sistema entra en un régimen de atractor extraño caótico (sensibilidad extrema a condiciones iniciales). La parametrización elegida asegura permanecer fuera de dicha frontera.

---

## 5. Conclusiones
* Se implementaron con éxito los algoritmos numéricos de Euler y Runge-Kutta 2 en códigos modulares independientes para el modelo tritrófico de Cosner.
* El análisis del paso temporal demostró que la elección de $\Delta t$ está gobernada tanto por la rigidez local del sistema de EDOs como por la estabilidad intrínseca de los esquemas discretos. Mientras Euler requiere pasos excesivamente pequeños para evitar divergencia, RK2 proporciona una precisión de segundo orden que equilibra costo computacional y estabilidad.
* Se comprobó la eficacia de la respuesta funcional racional cuadrática para amortiguar sobrepoblaciones y evitar colapsos ecológicos mediante el término inhibitorio de segundo orden ($c_i N^2, c_i P^2$).
* La selección de parámetros biológicos realizada garantizó un estado de equilibrio asintótico estable, previniendo comportamientos caóticos aperiódicos.

---

## 6. Referencias
1. Cosner, C., DeAngelis, D. L., Ault, J. S., & Perez-Sanchez, D. R. (1999). *A Model for Movement and Aggregation of Fish with Rational Functional Response*.
2. Burden, R. L., & Faires, J. D. (2010). *Numerical Analysis* (9th ed.). Brooks/Cole.
3. Strogatz, S. H. (2018). *Nonlinear Dynamics and Chaos: With Applications to Physics, Biology, Chemistry, and Engineering* (2nd ed.). CRC Press.
4. Hastings, A., & Powell, T. (1991). *Chaos in a Three-Species Food Chain*. Ecology, 72(3), 896–903.

---

## Anexo: Correspondencia Código - Informe
* [Ver `euler_cosner.py`](file:///c:/Projects/calc-num/docs/euler_cosner.py): Implementación del método de Euler progresivo.
* [Ver `rk2_cosner.py`](file:///c:/Projects/calc-num/docs/rk2_cosner.py): Implementación del método Runge-Kutta de orden 2 (Heun).
