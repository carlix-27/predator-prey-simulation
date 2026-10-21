# Guía Estratégica para la Defensa Oral (10 Minutos + Preguntas)

Esta guía organiza los 10 minutos de exposición grupal y prepara las respuestas clave para la ronda de 2-3 minutos de preguntas del docente.

---

## 1. Cronograma de Exposición (Total: 10 Minutos)

```
[00:00 - 02:00]  1. Introducción y Contexto del Modelo (Cosner vs Lotka-Volterra)
[02:00 - 05:00]  2. Metodología Numérica: Euler vs. RK2 y Justificación de dt
[05:00 - 08:00]  3. Resultados Gráficos y Régimen Asintótico (Coexistencia y Caos)
[08:00 - 10:00]  4. Trazabilidad del Código y Conclusiones Finales
```

---

### Minutos 0 a 2: El Modelo Ecológico (Cosner Tritrófico)
* **Mensaje clave:** No nos quedamos en Lotka-Volterra clásico porque es biológicamente irreal (crecimiento infinito de presas, depredación sin límite, órbitas cerradas no estructuralmente estables).
* **Qué destacar de Cosner:**
  * Tres niveles: Presa ($N$), Depredador ($P$) y Superdepredador ($S$).
  * Respuesta racional cuadrática $F(N) = \frac{aN}{1 + bN + cN^2}$.
  * El término cuadrático $cN^2$ modela la *defensa grupal* / *confusión del depredador* a densidades altas. Evita explosiones no realistas de biomasa.

---

### Minutos 2 a 5: Métodos Numéricos y Justificación de $dt$
* **Comparativa de Métodos:**
  * **Euler Explícito:** Orden 1. Error de truncamiento local $\mathcal{O}(\Delta t^2)$, global $\mathcal{O}(\Delta t)$. Muy sensible al paso.
  * **RK2 (Heun):** Orden 2. Error local $\mathcal{O}(\Delta t^3)$, global $\mathcal{O}(\Delta t^2)$. Predictor-corrector que promedia pendientes.
* **Justificación del paso $dt$ (¡Pregunta fija de examen!):**
  * *¿Por qué no usar un $dt$ grande?* Porque el Jacobiano del sistema tiene autovalores propios del acoplamiento. Si $\Delta t$ supera el límite de estabilidad de Dahlquist, el método genera amplificación numérica artificial (análogo al oscilador armónico que gana energía espuria con Euler).
  * *Criterio adoptado:* Reducción sucesiva de paso ($\Delta t \to \Delta t/2$) hasta obtener coincidencia numérica asintótica. Se adoptó $\Delta t = 0.05$ por balance óptimo entre precisión y costo.

---

### Minutos 5 a 8: Resultados Gráficos y Evitación del Caos
* **Mostrar Figura 1 (Evolución temporal):** Destacar que tras el transitorio inicial las tres poblaciones convergen asintóticamente a valores constantes (equilibrio no trivial).
* **Mostrar Figura 2 (Espacio de fases 3D):** Explicar la trayectoria en espiral hacia el atractor de punto fijo.
* **El tema del Caos:** Explicar que en sistemas de dimensión 3 con acoplamiento no local o retardos pueden surgir atractores extraños caóticos (efecto mariposa / sensibilidad a condiciones iniciales). Explicar que se calibraron las tasas de ataque dentro de los rangos típicos de la consigna ($a, c, e \le 0.15$) para asegurar estabilidad asintótica de Lyapunov y evitar el caos.

---

### Minutos 8 a 10: Trazabilidad del Código y Cierre
* Mostrar que el código `euler_cosner.py` y `rk2_cosner.py` tienen nombres de variables exactamente idénticos a los del informe (`r`, `K`, `F_PN`, `F_SN`, etc.) y que cada bloque numérico está referenciado con su sección correspondiente.

---

## 2. Banco de Preguntas Típicas del Jurado ("¿Por qué?")

### P1: "¿Por qué el método de Euler no es adecuado para pasos de tiempo moderados o sistemas oscilatorios?"
> **Respuesta:** "Porque la región de estabilidad absoluta de Euler explícito es el círculo $|1 + \lambda \Delta t| \le 1$. En sistemas presa-depredador con componentes oscilatorias (autovalores con parte imaginaria pura o dominante), Euler siempre queda fuera del disco de estabilidad si $\Delta t$ no es infinitamente pequeño, sumando error en cada ciclo y creando amplificación artificial de las poblaciones."

### P2: "¿Cómo justificaron matemáticamente la elección del $\Delta t$?"
> **Respuesta:** "Realizamos un test de refinamiento de malla o doble paso: computamos la solución con pasos decrecientes ($\Delta t = 0.5, 0.1, 0.05, 0.01$). Observamos que con $\Delta t = 0.5$ Euler oscila e diverge, mientras que a partir de $\Delta t = 0.05$ la diferencia relativa con $\Delta t = 0.01$ es inferior al $0.1\%$, logrando convergencia numérica tanto en Euler como en RK2."

### P3: "¿Por qué en este modelo tritrófico puede haber caos mientras que en Lotka-Volterra clásico de 2 especies no?"
> **Respuesta:** "Por el Teorema de Poincaré-Bendixson: en el plano bidimensional continuo ($\mathbb{R}^2$), las trayectorias no pueden cruzarse y los únicos conjuntos límite posibles son puntos fijos o ciclos límites periódicos. Al agregar un tercer nivel trófico (superdepredador), el espacio de estados pasa a ser tridimensional ($\mathbb{R}^3$), lo cual permite el estiramiento y plegado del espacio de fases característico del caos determinista."
