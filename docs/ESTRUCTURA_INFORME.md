# Estructura y Plan Maestro del Informe: Modelo Tritrófico de Cosner

## 1. Portada y Metadatos
* **Título sugerido:** *Modelado y Simulación Numérica del Sistema Tritrófico Presa-Depredador-Superdepredador con Respuesta Funcional Racional Generalizada (Cosner): Análisis de Estabilidad y Métodos de Integración Temporal (Euler vs. RK2).*
* **Autores:** [Nombres y Apellidos de los integrantes del grupo]
* **Cátedra:** Cálculo Numérico / Métodos Numéricos
* **Fecha de entrega:** Octubre 2026

---

## 2. Estructura de Secciones (Formato Paper)

```
├── 1. Resumen (Abstract)
├── 2. Introducción y Contexto Biológico
├── 3. Modelo Matemático: Cosner Tritrófico
│   ├── 3.1. Ecuaciones Diferenciales Acopladas
│   ├── 3.2. Respuesta Funcional de Cosner (Racional Cuadrática)
│   └── 3.3. Puntos de Equilibrio y Estabilidad Cualitativa
├── 4. Metodología Numérica y Análisis de Convergencia
│   ├── 4.1. Método de Euler Progresivo (Primer Orden)
│   ├── 4.2. Método de Runge-Kutta de Orden 2 (RK2 / Heun)
│   ├── 4.3. Error Local y Global de Truncamiento
│   ├── 4.4. Criterio de Selección y Justificación del Paso dt
│   └── 4.5. Fenómenos de Caos, Rigidez Numérica e Inestabilidad
├── 5. Resultados Numéricos y Discusión
│   ├── 5.1. Selección de Parámetros Biológicamente Viables (Equilibrio Estable)
│   ├── 5.2. Evolución Temporal de Poblaciones (N(t), P(t), S(t))
│   ├── 5.3. Espacio de Fases 3D y Planos de Fase 2D
│   ├── 5.4. Estudio de Sensibilidad al dt y Comparación Euler vs. RK2
│   └── 5.5. Exploración de Zonas de Bifurcación / Dinámica Caótica
├── 6. Conclusiones
├── 7. Referencias Bibliográficas
└── Anexo: Código Fuente (Trazabilidad línea a línea con el informe)
```

---

## 3. Desglose Contenido por Sección

### Sección 1: Resumen (Abstract)
* **Objetivo:** Un único párrafo (150 - 250 palabras) auto-contenido.
* **Puntos clave a incluir:**
  1. Planteo del problema tritrófico (Presa $N$, Depredador $P$, Superdepredador $S$) con respuesta de Cosner.
  2. Implementación de los métodos de integración temporal: Euler explícito y Runge-Kutta 2 (RK2).
  3. Criterio de selección del paso $\Delta t$ basado en estabilidad numérica y precisión asintótica.
  4. Resultados principales: coexistencia estable vs. oscilaciones, evitando regiones caóticas numéricas y biológicas.

---

### Sección 2: Introducción
* **Antecedentes:** Evolución desde el modelo maltusiano exponencial $\frac{dN}{dt} = rN$, pasando por Verhulst (logístico con capacidad de carga $K$) hasta Lotka-Volterra clásico (interacciones bilineales $NP$).
* **Limitaciones de Lotka-Volterra clásico:** Supone que un depredador puede consumir infinitas presas sin saturación biológica (respuesta funcional lineal tipo Holling I) y es neutralmente estable (órbitas conservativas no robustas ante perturbaciones).
* **Motivación del modelo asignado:** Cadena trófica de tres niveles (presa, depredador intermedio y superdepredador omnívoro) con saturación e interferencia (respuesta racional generalizada de Cosner).
* **Objetivo del trabajo:** Resolver el sistema mediante algoritmos discretos, comparar la eficiencia de Euler vs RK2, justificar rigurosamente la elección del paso $\Delta t$ y la parametrización ecológica para alcanzar un equilibrio asintótico.

---

### Sección 3: Formulación Matemática del Modelo de Cosner
* **Sistema de EDOs acopladas:**
  $$\frac{dN}{dt} = r N \left(1 - \frac{N}{K}\right) - d N - F_{PN}(N, P) P - F_{SN}(N, S) S$$
  $$\frac{dP}{dt} = b F_{PN}(N, P) P - m P - F_{SP}(P, S) S$$
  $$\frac{dS}{dt} = f F_{SN}(N, S) S + g F_{SP}(P, S) S - q S$$

* **Respuestas Funcionales no lineales (Cosner):**
  $$F_{PN}(N, P) = \frac{a N}{1 + b_1 N + c_1 N^2}$$
  $$F_{SN}(N, S) = \frac{c N}{1 + b_2 N + c_2 N^2}$$
  $$F_{SP}(P, S) = \frac{e P}{1 + b_3 P + c_3 P^2}$$
  * *Interpretación biológica:* A densidades bajas, el consumo crece con la densidad; a densidades medias se satura (tiempo de manipulación $b_i$); a densidades altas el término cuadrático $c_i N^2$ modela inhibición grupal, defensa de grupo o evasión cooperativa.
* **Tabla de Parámetros:** Definición conceptual, unidades físicas/temporales y rangos admisibles según consigna.

---

### Sección 4: Metodología Numérica y Análisis de Estabilidad
* **Formulación discreta del problema de Cauchy / PVI:**
  $$\mathbf{y}'(t) = \mathbf{f}(t, \mathbf{y}), \quad \mathbf{y}(t_0) = \mathbf{y}_0, \quad \mathbf{y} = [N, P, S]^T$$
* **Algoritmo 1: Euler Explícito**
  $$\mathbf{y}_{n+1} = \mathbf{y}_n + \Delta t \cdot \mathbf{f}(t_n, \mathbf{y}_n)$$
  * Error local: $\mathcal{O}(\Delta t^2)$, error global acumulado: $\mathcal{O}(\Delta t)$.
  * Región de estabilidad lineal (condición CFL sobre el Jacobiano $J = \frac{\partial \mathbf{f}}{\partial \mathbf{y}}$).
* **Algoritmo 2: Runge-Kutta 2 (Heun / Punto Medio)**
  $$\mathbf{k}_1 = \mathbf{f}(t_n, \mathbf{y}_n)$$
  $$\mathbf{k}_2 = \mathbf{f}(t_n + \Delta t, \mathbf{y}_n + \Delta t \mathbf{k}_1)$$
  $$\mathbf{y}_{n+1} = \mathbf{y}_n + \frac{\Delta t}{2} (\mathbf{k}_1 + \mathbf{k}_2)$$
  * Error local: $\mathcal{O}(\Delta t^3)$, error global: $\mathcal{O}(\Delta t^2)$.
* **Criterio de Selección de $\Delta t$:**
  * Justificación de convergencia: Criterio de doblado de paso (Runge / Richardson). Reducción sucesiva de $\Delta t \to \Delta t/2$ hasta que la diferencia relativa $\|\mathbf{y}_{\Delta t}(t_{max}) - \mathbf{y}_{\Delta t/2}(t_{max})\| < \epsilon_{tol}$.
  * Comparación conceptual con el oscilador armónico simple: demostrar por qué Euler amplifica la energía artificialmente si $\Delta t$ no es suficientemente pequeño y cómo RK2 amortigua el residuo.

---

### Sección 5: Resultados y Discusión
* **Estructura visual recomendada (gráficos con poco texto explicativo):**
  1. **Figura 1: Dinámica Temporal hacia el Equilibrio Asintótico** ($N(t), P(t), S(t)$ vs. $t$). Muestra cómo las 3 poblaciones se estabilizan en valores no nulos (coexistencia).
  2. **Figura 2: Diagrama de Fases Tridimensional (3D)** en el espacio $(N, P, S)$ mostrando la convergencia espiralada al atractor de punto fijo.
  3. **Figura 3: Análisis de Convergencia vs $\Delta t$ (Euler vs RK2)**. Gráfico de discrepancia o trayectorias con $\Delta t \in \{0.5, 0.1, 0.01, 0.001\}$. Evidencia empírica de inestabilidad de Euler para $\Delta t$ grande frente a la robustez de RK2.
  4. **Figura 4: Régimen de Caos / Bifurcación (Opcional demostrativo):** Demostración de qué ocurre cuando se eleva la tasa de ataque o interacciones no locales (órbitas no periódicas aperiódicas tipo Atractor Extraño) y por qué se calibraron los parámetros para operar fuera de esa zona.

---

### Sección 6: Conclusiones
* Síntesis de hallazgos: consistencia biológica del modelo de Cosner, verificación de órdenes de convergencia numérica, superioridad de RK2 en relación precisión/costo computacional y delimitación de condiciones estables.

---

### Sección 7: Referencias y Anexos
* Formato estándar IEEE/APA.
* Anexo con scripts limpios, modulares y profusamente comentados (`euler_cosner.py` y `rk2_cosner.py`).
