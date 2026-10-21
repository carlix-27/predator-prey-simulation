# Modelo Predador-Presa Tritrófico de Cosner (Cálculo Numérico)

Simulación numérica y análisis de estabilidad del modelo ecológico de tres niveles tróficos (Presa $N$, Depredador $P$, Superdepredador $S$) con respuesta funcional racional cuadrática generalizada (Cosner), desarrollado para la cátedra de **Cálculo Numérico**.

---

## 📚 Estructura de Documentación (`docs/`)

* **[`docs/ESTRUCTURA_INFORME.md`](docs/ESTRUCTURA_INFORME.md)**: Plan maestro y estructura detallada del informe en formato paper.
* **[`docs/INFORME_PRELIMINAR.md`](docs/INFORME_PRELIMINAR.md)**: Borrador completo del informe (Abstract, Introducción, Modelo, Metodología Numérica, Resultados, Conclusiones y Referencias).
* **[`docs/FUENTES_Y_BIBLIOGRAFIA.md`](docs/FUENTES_Y_BIBLIOGRAFIA.md)**: Compendio de papers y libros de referencia (Cosner, Hastings-Powell, Strogatz, Burden-Faires).
* **[`docs/GUIA_DEFENSA_10MIN.md`](docs/GUIA_DEFENSA_10MIN.md)**: Planificación de la exposición oral de 10 minutos y banco de preguntas típicas de defensa.

---

## 💻 Implementación Numérica

* **[`docs/euler_cosner.py`](docs/euler_cosner.py)**: Algoritmo de Euler explícito (primer orden, $\mathcal{O}(\Delta t)$) con evaluación de estabilidad y sensibilidad al paso $\Delta t$.
* **[`docs/rk2_cosner.py`](docs/rk2_cosner.py)**: Algoritmo Runge-Kutta 2 / Heun (segundo orden, $\mathcal{O}(\Delta t^2)$) con generación de curvas asintóticas y diagrama de fases tridimensional $(N, P, S)$.

### Ejecución
```bash
# Integración con Euler Explícito
python docs/euler_cosner.py

# Integración con Runge-Kutta 2 (Heun)
python docs/rk2_cosner.py
```
