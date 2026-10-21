# Fuentes Bibliográficas y Marco Teórico de Referencia

Este compendio de fuentes proporciona el sustento teórico y matemático necesario para redactar el informe, justificar las decisiones de diseño numérico y preparar la defensa oral (10 minutos de charla + 3 minutos de preguntas).

---

## 1. Modelado Matemático y Sistemas Predador-Presa Tritróficos

1. **Cosner, C., DeAngelis, D. L., Ault, J. S., & Perez-Sanchez, D. R. (1999).**
   * *A Model for Movement and Aggregation of Fish with Rational Functional Response.*
   * **Relevancia:** Origen y justificación biológica de las respuestas funcionales racionales del tipo $\frac{a N}{1 + b N + c N^2}$. Explica la no-monotonía (inhibición por alta densidad, formación de cardúmenes o defensa grupal) frente a las funciones de Holling Tipo II y III estándar.
   * **Concepto clave para la defensa:** ¿Por qué el término cuadrático $c N^2$ en el denominador? Porque a densidades extremadamente altas de presa, se produce interferencia o defensa colectiva (efecto *group defense* o confusión del depredador), reduciendo la tasa de captura per cápita.

2. **Hastings, A., & Powell, T. (1991).**
   * *Chaos in a Three-Species Food Chain.* Ecology, 72(3), 896–903.
   * **Relevancia:** Paper seminal que demuestra que los sistemas depredador-presa de 3 niveles tróficos (presa, depredador, superdepredador) exhiben bifurcaciones hacia el caos mediante cascadas de duplicación de período (atractor de Hastings-Powell).
   * **Concepto clave para la defensa:** Explica por qué en 2 especies (teorema de Poincaré-Bendixson) el caos en tiempo continuo es imposible, pero en 3 especies (dimensión $\ge 3$) la dinámica caótica emerge de forma natural si los parámetros de ataque no están acotados.

3. **Murray, J. D. (2002).**
   * *Mathematical Biology I: An Introduction* (3rd ed.). Springer-Verlag.
   * **Relevancia:** Capítulos 3 y 4: Dinámica de poblaciones continuas, modelo logístico de Verhulst, capacidad de carga $K$, depredación con respuesta funcional y análisis cualitativo en el plano de fases.

4. **Strogatz, S. H. (2018).**
   * *Nonlinear Dynamics and Chaos: With Applications to Physics, Biology, Chemistry, and Engineering* (2nd ed.). CRC Press.
   * **Relevancia:** Excelente referencia pedagógica para linealización alrededor de puntos fijos, cálculo de autovalores del Jacobiano, estabilidad asintótica vs. órbitas periódicas, y teoría del atractor de Lorenz / efecto mariposa.

---

## 2. Análisis Numérico y Métodos de Ecuaciones Diferenciales Ordinarias (EDOs)

1. **Burden, R. L., & Faires, J. D. (2010).**
   * *Numerical Analysis* (9th ed.). Brooks/Cole, Cengage Learning.
   * **Capítulo 5: Problemas de Valor Inicial para EDOs.**
   * **Relevancia:**
     * Deducción de la fórmula de Taylor para el método de Euler.
     * Demostración de que el error de truncamiento local (LTE) de Euler es $\mathcal{O}(\Delta t^2)$ y el error global es $\mathcal{O}(\Delta t)$.
     * Familia de métodos de Runge-Kutta de orden 2 (método de Heun y método del Punto Medio) con LTE $\mathcal{O}(\Delta t^3)$ y error global $\mathcal{O}(\Delta t^2)$.
     * Análisis de estabilidad absoluta para la ecuación test de Dahlquist $y' = \lambda y$.

2. **Hairer, E., Nørsett, S. P., & Wanner, G. (1993).**
   * *Solving Ordinary Differential Equations I: Nonstiff Problems.* Springer-Verlag.
   * **Relevancia:** Referencia clásica mundial para comprender la propagación de errores acumulativos en sistemas de EDOs no lineales y el control adaptativo del paso temporal.

3. **Chapra, S. C., & Canale, R. P. (2015).**
   * *Numerical Methods for Engineers* (7th ed.). McGraw-Hill Education.
   * **Relevancia:** Análisis comparativo del oscilador armónico simple resuelto por Euler vs RK2. Demuestra analíticamente que Euler produce ganancia artificial de energía (las trayectorias se abren en espiral divergente debido a que la región de estabilidad no contiene el eje imaginario puro), mientras que RK2 reduce drásticamente esta dispersión numérica.

---

## 3. Justificación de Criterios Numéricos Clave (Preguntas de Examen)

### ¿Por qué justificar el $\Delta t$?
* **Cuestión:** Si $\Delta t$ es demasiado grande, los métodos explícitos se vuelven inestables numéricamente (los autovalores del Jacobiano multiplicados por $\Delta t$ caen fuera de la región de estabilidad absoluta $|1 + \lambda \Delta t| \le 1$ para Euler). Las poblaciones computadas divergen a $\pm \infty$ o dan valores biológicamente absurdos (negativos).
* **Solución de doble paso:** Se toma un $\Delta t_0$, se computa la solución hasta $t_{max}$, luego se repite con $\Delta t_1 = \Delta t_0 / 2$. Si $\|\mathbf{y}_{\Delta t_1}(t) - \mathbf{y}_{\Delta t_0}(t)\| < \text{Tolerancia}$, el paso queda validado tanto por estabilidad como por convergencia de malla.

### ¿Cómo evitar el caos?
* El informe de la consigna señala: *"Aparecen oscilaciones caóticas al combinarse con interacciones no locales y un retraso de mayor duración. Un alcance elevado de las interacciones no locales conduce al caos."*
* Para garantizar un estado asintótico de equilibrio:
  1. Mantener las tasas de ataque $a, c, e$ dentro de rangos moderados (ej. $0.05 - 0.20$).
  2. Ajustar los coeficientes de saturación $b_i$ y $c_i$ para que no ocurran saltos abruptos en el consumo.
  3. Asegurar que las tasas de mortalidad $m, q$ estén equilibradas con la eficiencia de conversión $b, f, g$, logrando que la traza del Jacobiano en el punto de coexistencia tenga parte real negativa ($\text{Re}(\lambda_i) < 0$).
