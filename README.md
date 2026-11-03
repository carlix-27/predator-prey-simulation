# Modelo Predador-Presa Tritrófico de Cosner (Cálculo Numérico)

Simulación numérica y análisis de estabilidad del modelo ecológico de tres niveles tróficos (Presa $N$, Depredador $P$, Superdepredador $S$) con respuesta funcional racional cuadrática generalizada (Cosner), desarrollado para la cátedra de **Cálculo Numérico**.

---

## 📄 Informe

* **[`docs/INFORME_PRELIMINAR.md`](docs/INFORME_PRELIMINAR.md)**: informe completo (fuente). GitHub muestra las ecuaciones y las figuras; el Anexo A incluye el código fuente.
* **[`docs/INFORME_EDITABLE.docx`](docs/INFORME_EDITABLE.docx)**: versión Word con ecuaciones nativas editables. También se abre en Google Docs (*Archivo → Abrir → Subir*).
* **[`docs/INFORME.pdf`](docs/INFORME.pdf)**: versión para imprimir o entregar.
* **[`docs/INFORME_EDITABLE.html`](docs/INFORME_EDITABLE.html)**: versión web autocontenida, que funciona sin conexión.
* **[`docs/figuras/`](docs/figuras/)**: figuras 1 a 6 del informe.

## 📚 Documentación de apoyo (`docs/`)

* **[`docs/ESTRUCTURA_INFORME.md`](docs/ESTRUCTURA_INFORME.md)**: plan maestro y estructura del informe en formato paper.
* **[`docs/FUENTES_Y_BIBLIOGRAFIA.md`](docs/FUENTES_Y_BIBLIOGRAFIA.md)**: compendio de papers y libros de referencia (Cosner, Hastings-Powell, Rosenzweig, Strogatz, Burden-Faires).
* **[`docs/GUIA_DEFENSA_10MIN.md`](docs/GUIA_DEFENSA_10MIN.md)**: planificación de la exposición oral de 10 minutos y banco de preguntas típicas de defensa.

---

## 💻 Implementación Numérica

* **[`docs/euler_cosner.py`](docs/euler_cosner.py)**: Euler explícito (primer orden, $\mathcal{O}(\Delta t)$).
* **[`docs/rk2_cosner.py`](docs/rk2_cosner.py)**: Runge-Kutta 2 / Heun (segundo orden, $\mathcal{O}(\Delta t^2)$).
* **[`docs/generar_figuras.py`](docs/generar_figuras.py)**: reutiliza ambos integradores y genera todas las figuras y valores del informe. Incluye equilibrio por Newton, autovalores, paso máximo estable, tabla de doble paso, bifurcación de Hopf y exponentes de Lyapunov.
* **[`docs/construir_informe.py`](docs/construir_informe.py)**: inserta el código como anexo y exporta el informe a `.docx`, `.html` y `.pdf`.

### Ejecución
```bash
pip install -r requirements.txt

# Integración individual con cada método
python docs/euler_cosner.py
python docs/rk2_cosner.py

# Figuras y resultados del informe (~5 min por el diagrama de bifurcación)
python docs/generar_figuras.py

# Regenerar INFORME_EDITABLE.docx / .html e INFORME.pdf (el PDF usa Edge o Chrome)
python docs/construir_informe.py
```
