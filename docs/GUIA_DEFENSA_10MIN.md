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
  * *Criterio adoptado:* Doble paso, $D(\Delta t) = \max_n \|\mathbf{y}_{\Delta t} - \mathbf{y}_{\Delta t/2}\|_\infty < 10^{-3}$. Con RK2 y $\Delta t = 0.05$ se obtiene $D = 1.5\times10^{-4}$ con 5000 pasos. Euler con $h=0.001$ da $D=5.4532\times10^{-4}$ y requiere 250 000 evaluaciones para resolver; validar añade 500 000. Para RK2 las cifras son 10 000 y 20 000.
  * *Órdenes verificados:* pendientes log-log de 1.04 (Euler) y 2.00 (RK2), es decir, al dividir $\Delta t$ por 2 el error baja a la mitad o a la cuarta parte.

---

### Minutos 5 a 8: Resultados Gráficos y Evitación del Caos
* **Mostrar Figura 1 (Evolución temporal):** Destacar que tras el transitorio oscilatorio (período ≈ 61) las tres poblaciones convergen a la coexistencia $(N^*, P^*, S^*) = (5.105,\ 4.308,\ 1.878)$, un foco estable con autovalores $-0.1433$ y $-0.0386 \pm 0.1033i$.
* **Mostrar Figura 2 (Espacio de fases 3D):** Explicar la trayectoria en espiral hacia el atractor de punto fijo.
* **Mostrar Figura 3 (Euler vs. RK2):** Con $\Delta t = 2$, Euler amplifica las oscilaciones y reduce mucho la densidad del superdepredador ($S(250) \approx 0.0091$); eso no demuestra extinción. RK2 sigue la referencia.
* **Análisis complementario:** el cruce de autovalores cerca de $K=11.18$ sugiere Hopf, pero no prueba la estabilidad de un ciclo. Una nube de máximos indica irregularidad; requiere comprobar transitorios y refinamientos antes de afirmar caos.

> **Ojo:** la parametrización original ($K = 100$, $m = 0.2$, $q = 0.1$, $b_1 = b_3 = 0.5$) extinguía a $P$ y $S$, porque $b \cdot \max F_{PN} = 0.107 < m = 0.2$. Ver la Sección 2.4 del informe.

---

### Minutos 8 a 10: Trazabilidad del Código y Cierre
* Mostrar que el código `euler_cosner.py` y `rk2_cosner.py` tienen nombres de variables exactamente idénticos a los del informe (`r`, `K`, `F_PN`, `F_SN`, etc.) y que cada bloque numérico está referenciado con su sección correspondiente.

---

## 2. Banco de Preguntas Típicas del Jurado ("¿Por qué?")

### P1: "¿Por qué el método de Euler no es adecuado para pasos de tiempo moderados o sistemas oscilatorios?"
> **Respuesta:** "Porque la región de estabilidad absoluta de Euler explícito es el círculo $|1 + \lambda \Delta t| \le 1$. Para un autovalor $\lambda = \alpha + i\beta$ el paso máximo es $-2\alpha/|\lambda|^2$, que tiende a cero cuando el modo es casi conservativo ($\alpha \approx 0$). En nuestro equilibrio interior eso da $\Delta t_{\max} = 6.35$. En cambio, en el borde $N$–$P$ (autovalores $-0.0137 \pm 0.247i$) Euler con $\Delta t = 2$ tiene $|1 + \lambda\Delta t| = 1.091$, amplifica las oscilaciones un 9% por paso y termina extinguiendo numéricamente al superdepredador (Figura 3)."

### P2: "¿Cómo justificaron matemáticamente la elección del $\Delta t$?"
> **Respuesta:** "Con el criterio de doble paso: integramos con $\Delta t$ y $\Delta t/2$ y medimos $D(\Delta t) = \max_n \|\mathbf{y}_{\Delta t} - \mathbf{y}_{\Delta t/2}\|_\infty$ para $\Delta t$ entre 0.5 y 0.01. RK2 cumple $D < 10^{-3}$ desde $\Delta t = 0.1$, y elegimos $\Delta t = 0.05$ ($D = 1.5\times10^{-4}$) para tener margen. La pendiente log-log confirma los órdenes teóricos: 1.04 para Euler y 2.00 para RK2 (Figura 4)."

### P2b: "¿Por qué cambiaron los parámetros respecto del planteo inicial?"
> **Respuesta:** "Porque con $m = 0.2$ el depredador no puede compensar su mortalidad: el máximo de $b\,F_{PN}$ es $0.5 \cdot 0.15/(0.5 + 2\sqrt{0.01}) = 0.107 < 0.2$, y el sistema converge a $(90, 0, 0)$. Ajustamos $b_1$, $b_3$ (dentro de rango), $m$, $q$ y $K$ (sin rango en la consigna). Con $K = 100$ el equilibrio de coexistencia sería inestable por la paradoja del enriquecimiento; el umbral es $K_H \approx 11.18$."

### P2c: "¿Por qué no iteran el corrector de Heun hasta que converja?"
> **Respuesta:** "Porque no hace falta. Hacemos un paso con una sola corrección: predictor de Euler, $\mathbf{k}_2$ en el punto predicho y $\mathbf{y}_{n+1} = \mathbf{y}_n + \frac{\Delta t}{2}(\mathbf{k}_1 + \mathbf{k}_2)$. Si iteráramos el corrector, la iteración de punto fijo convergería a la regla del trapecio implícita, que también es de orden 2. Iterar no sube el orden, pero modifica el costo y puede modificar la estabilidad. Con una corrección ya obtenemos orden 2 (medimos 2.00) con dos evaluaciones por paso y el método sigue siendo explícito."

### P3: "¿Por qué en este modelo tritrófico puede haber caos mientras que en Lotka-Volterra clásico de 2 especies no?"
> **Respuesta:** "Por el Teorema de Poincaré-Bendixson: para EDO autónomas planas, bajo las hipótesis del teorema, un conjunto límite compacto sin equilibrios es una órbita periódica. Al agregar un tercer nivel trófico (superdepredador), el espacio de estados pasa a ser tridimensional ($\mathbb{R}^3$), lo cual permite el estiramiento y plegado del espacio de fases característico del caos determinista."
