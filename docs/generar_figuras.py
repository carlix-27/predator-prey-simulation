"""
===============================================================================
Trabajo Práctico: Cálculo Numérico - Modelo Tritrófico de Cosner
Archivo: generar_figuras.py
Propósito: Generar todas las figuras y tablas numéricas del informe
===============================================================================
Descripción:
    Reutiliza el campo vectorial y los integradores implementados en
    euler_cosner.py y rk2_cosner.py (no reimplementa ningún método) y produce:

        Figura 1: Evolución temporal N(t), P(t), S(t) con RK2 (dt elegido).
        Figura 2: Espacio de fases 3D (N, P, S).
        Figura 3: Sensibilidad al paso dt: Euler vs. RK2.
        Figura 4: Orden de convergencia observado (criterio de doble paso).
        Figura 5: Diagrama de bifurcación respecto de la capacidad de carga K.
        Figura 6: Régimen caótico: atractor extraño y sensibilidad a y0.

    Además imprime en consola el punto de equilibrio, los autovalores del
    Jacobiano, los pasos máximos estables, la bifurcación de Hopf, los
    exponentes de Lyapunov y la tabla del criterio de doble paso del informe.

Ejecución:
    python docs/generar_figuras.py      (las imágenes se guardan en docs/figuras/)
===============================================================================
"""

import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from euler_cosner import resolver_euler
from rk2_cosner import sistema_cosner, resolver_rk2


# =============================================================================
# 1. CONFIGURACIÓN (Tabla de la Sección 2.3 del Informe)
# =============================================================================

PARAMETROS = {
    'r': 1.0,           # Natalidad de la presa
    'K': 10.0,          # Capacidad de carga
    'd': 0.1,           # Mortalidad de la presa
    'b': 0.5,           # Eficiencia conversión presa -> depredador
    'm': 0.1,           # Mortalidad depredador
    'f': 0.2,           # Eficiencia conversión presa -> superdepredador
    'g': 0.4,           # Eficiencia conversión depredador -> superdepredador
    'q': 0.12,          # Mortalidad superdepredador
    'a': 0.15,          # Tasa de ataque F_PN (0.01 - 0.5)
    'c_ataque': 0.05,   # Tasa de ataque F_SN (0.01 - 0.5)
    'e': 0.10,          # Tasa de ataque F_SP (0.01 - 0.5)
    'b1': 0.1,          # Coeficiente lineal N en F_PN (0.1 - 2.0)
    'b2': 0.5,          # Coeficiente lineal N en F_SN (0.1 - 2.0)
    'b3': 0.1,          # Coeficiente lineal P en F_SP (0.1 - 2.0)
    'c1': 0.01,         # Coeficiente cuadrático N^2 en F_PN (0.0 - 0.1)
    'c2': 0.01,         # Coeficiente cuadrático N^2 en F_SN (0.0 - 0.1)
    'c3': 0.01          # Coeficiente cuadrático P^2 en F_SP (0.0 - 0.1)
}

ESTADO_INICIAL = np.array([8.0, 2.0, 1.0])   # [N0, P0, S0]
T_INICIO = 0.0
T_FINAL = 250.0
DT_ELEGIDO = 0.05

CARPETA_FIGURAS = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figuras')
COLORES = {'N': 'tab:blue', 'P': 'tab:orange', 'S': 'tab:green'}


# =============================================================================
# 2. HERRAMIENTAS DE ANÁLISIS (Secciones 3.3 y 4 del Informe)
# =============================================================================

def jacobiano(y, params, h=1e-6):
    """
    Aproxima la matriz Jacobiana J = df/dy por diferencias centradas:
        J[:, j] = (f(y + h e_j) - f(y - h e_j)) / (2h)
    """
    J = np.zeros((3, 3))
    for j in range(3):
        e_j = np.zeros(3)
        e_j[j] = h
        J[:, j] = (sistema_cosner(0.0, y + e_j, params) - sistema_cosner(0.0, y - e_j, params)) / (2.0 * h)
    return J


def equilibrio_newton(y_semilla, params, tol=1e-12, max_iter=50):
    """
    Refina un punto de equilibrio f(y*) = 0 con el método de Newton-Raphson:
        J(y_k) * delta = -f(y_k),   y_{k+1} = y_k + delta
    """
    y = np.array(y_semilla, dtype=float)
    for _ in range(max_iter):
        delta = np.linalg.solve(jacobiano(y, params), -sistema_cosner(0.0, y, params))
        y = y + delta
        if np.max(np.abs(delta)) < tol:
            break
    return y


def paso_maximo_estable(autovalores, factor_amplificacion, dt_max=50.0, n=500001):
    """
    Mayor dt tal que |R(lambda * dt)| <= 1 para todos los autovalores (ecuación test
    de Dahlquist). R(z) = 1 + z para Euler y R(z) = 1 + z + z^2/2 para RK2.
    """
    pasos = np.linspace(0.0, dt_max, n)[1:]
    inestable = np.zeros_like(pasos, dtype=bool)
    for lam in autovalores:
        inestable |= np.abs(factor_amplificacion(lam * pasos)) > 1.0
    return pasos[np.argmax(inestable)] if inestable.any() else np.inf


def max_parte_real(params, y_semilla):
    y = equilibrio_newton(y_semilla, params)
    return np.max(np.linalg.eigvals(jacobiano(y, params)).real)


def K_critico_hopf(params, K_estable, K_inestable, y_semilla, tol=1e-4):
    """
    Bisección sobre K del máximo Re(lambda) en el equilibrio interior:
    el valor donde cruza el cero es la bifurcación de Hopf.
    """
    while K_inestable - K_estable > tol:
        K_medio = 0.5 * (K_estable + K_inestable)
        if max_parte_real(dict(params, K=K_medio), y_semilla) < 0.0:
            K_estable = K_medio
        else:
            K_inestable = K_medio
    return 0.5 * (K_estable + K_inestable)


def diferencia_doble_paso(resolver, dt, params, y0=ESTADO_INICIAL, t_max=T_FINAL):
    """
    Criterio de doble paso (Runge / Richardson):
        D(dt) = max_n || y_dt(t_n) - y_{dt/2}(t_n) ||_inf
    Se compara únicamente en los instantes comunes a ambas mallas.
    """
    t_grueso, y_grueso = resolver(sistema_cosner, y0, T_INICIO, t_max, dt, params)
    t_fino, y_fino = resolver(sistema_cosner, y0, T_INICIO, t_max, dt / 2.0, params)
    indices = np.searchsorted(t_fino, t_grueso)
    indices = np.minimum(indices, len(t_fino) - 1)
    previos = np.maximum(indices - 1, 0)
    indices = np.where(abs(t_fino[previos] - t_grueso) < abs(t_fino[indices] - t_grueso), previos, indices)
    validos = np.isclose(t_grueso, t_fino[indices], rtol=1e-12, atol=1e-12)
    return np.max(np.abs(y_grueso[validos] - y_fino[indices[validos]]))


def guardar(fig, nombre):
    os.makedirs(CARPETA_FIGURAS, exist_ok=True)
    ruta = os.path.join(CARPETA_FIGURAS, nombre)
    fig.savefig(ruta, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  -> {ruta}")


# =============================================================================
# 3. FIGURAS DEL INFORME (Sección 4)
# =============================================================================

def figura_1_evolucion(t, sol):
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(t, sol[:, 0], label='Presa $N(t)$', color=COLORES['N'], lw=2)
    ax.plot(t, sol[:, 1], label='Depredador $P(t)$', color=COLORES['P'], lw=2)
    ax.plot(t, sol[:, 2], label='Superdepredador $S(t)$', color=COLORES['S'], lw=2)
    ax.set_xlabel('Tiempo $t$')
    ax.set_ylabel('Densidad poblacional')
    ax.set_title(f'Evolución temporal del sistema de Cosner (RK2, $\\Delta t$ = {DT_ELEGIDO})')
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.legend()
    guardar(fig, 'fig1_evolucion_temporal.png')


def figura_2_fases_3d(sol):
    fig = plt.figure(figsize=(7.5, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(sol[:, 0], sol[:, 1], sol[:, 2], color='purple', lw=1.2)
    ax.scatter(*sol[0], color='red', s=40, label='Condición inicial $y_0$')
    ax.scatter(*sol[-1], color='green', s=60, label='Equilibrio $(N^*, P^*, S^*)$')
    ax.set_xlabel('Presa $N$')
    ax.set_ylabel('Depredador $P$')
    ax.set_zlabel('Superdepredador $S$')
    ax.set_title('Espacio de fases 3D (RK2)')
    ax.legend(loc='upper left')
    guardar(fig, 'fig2_fases_3d.png')


def figura_3_euler_vs_rk2(params):
    pasos = [2.0, 0.5, 0.1]
    t_ref, sol_ref = resolver_rk2(sistema_cosner, ESTADO_INICIAL, T_INICIO, T_FINAL, 0.01, params)
    fig, ejes = plt.subplots(2, 2, figsize=(11, 7), sharex=True, sharey='row')
    metodos = [('Euler explícito', resolver_euler), ('RK2 (Heun)', resolver_rk2)]
    for col, (nombre, resolver) in enumerate(metodos):
        soluciones = {dt: resolver(sistema_cosner, ESTADO_INICIAL, T_INICIO, T_FINAL, dt, params) for dt in pasos}
        for fila, (indice, especie) in enumerate([(1, 'depredador $P(t)$'), (2, 'superdepredador $S(t)$')]):
            ax = ejes[fila, col]
            ax.plot(t_ref, sol_ref[:, indice], color='black', lw=2.5, alpha=0.35, label='Referencia (RK2, $\\Delta t$ = 0.01)')
            for dt in pasos:
                t, sol = soluciones[dt]
                ax.plot(t, sol[:, indice], lw=1.3, label=f'$\\Delta t$ = {dt}')
            ax.set_title(f'{nombre}: {especie}')
            ax.grid(True, linestyle='--', alpha=0.6)
    for ax in ejes[1]:
        ax.set_xlabel('Tiempo $t$')
    for ax in ejes[:, 0]:
        ax.set_ylabel('Densidad poblacional')
    ejes[0, 1].legend(loc='upper right', fontsize=9)
    guardar(fig, 'fig3_euler_vs_rk2_dt.png')


def figura_4_orden_convergencia(params):
    pasos = np.array([0.5, 0.25, 0.1, 0.05, 0.025, 0.01])
    D_euler = np.array([diferencia_doble_paso(resolver_euler, dt, params) for dt in pasos])
    D_rk2 = np.array([diferencia_doble_paso(resolver_rk2, dt, params) for dt in pasos])

    print("\n  Tabla de doble paso  D(dt) = max_n ||y_dt - y_dt/2||_inf")
    print(f"  {'dt':>7} | {'Euler':>11} | {'RK2':>11} | {'orden local Euler':>15} | {'RK2':>6}")
    for i, (dt, de, dr) in enumerate(zip(pasos, D_euler, D_rk2)):
        if i + 1 < len(pasos):
            r = pasos[i] / pasos[i + 1]
            ordenes = f"{np.log(de / D_euler[i + 1]) / np.log(r):>15.2f} | {np.log(dr / D_rk2[i + 1]) / np.log(r):>6.2f}"
        else:
            ordenes = f"{'-':>15} | {'-':>6}"
        print(f"  {dt:>7} | {de:>11.3e} | {dr:>11.3e} | {ordenes}")
    # Pendiente ajustada en el rango asintótico (dt <= 0.1)
    asintotico = pasos <= 0.1
    orden_euler = np.polyfit(np.log(pasos[asintotico]), np.log(D_euler[asintotico]), 1)[0]
    orden_rk2 = np.polyfit(np.log(pasos[asintotico]), np.log(D_rk2[asintotico]), 1)[0]
    print(f"  Orden observado (pendiente log-log, dt <= 0.1): Euler = {orden_euler:.2f}, RK2 = {orden_rk2:.2f}")

    fig, ax = plt.subplots(figsize=(7, 4.8))
    ax.loglog(pasos, D_euler, 'o-', label=f'Euler (pendiente {orden_euler:.2f})')
    ax.loglog(pasos, D_rk2, 's-', label=f'RK2 (pendiente {orden_rk2:.2f})')
    ax.loglog(pasos, D_euler[-1] * (pasos / pasos[-1]), 'k--', lw=0.8, label='Referencia $\\mathcal{O}(\\Delta t)$')
    ax.loglog(pasos, D_rk2[-1] * (pasos / pasos[-1]) ** 2, 'k:', lw=0.8, label='Referencia $\\mathcal{O}(\\Delta t^2)$')
    ax.axhline(1e-3, color='gray', lw=0.8)
    ax.text(pasos[0], 1.3e-3, 'tolerancia $10^{-3}$', ha='right', color='gray', fontsize=9)
    ax.set_xlabel('Paso temporal $\\Delta t$')
    ax.set_ylabel('$D(\\Delta t)=\\max_n\\|y_{\\Delta t}-y_{\\Delta t/2}\\|_\\infty$')
    ax.set_title('Criterio de doble paso y orden de convergencia')
    ax.grid(True, which='both', linestyle='--', alpha=0.5)
    ax.legend()
    guardar(fig, 'fig4_orden_convergencia.png')


def maximos_locales(x):
    """Máximos locales de una serie temporal muestreada."""
    centro = x[1:-1]
    return centro[(centro > x[:-2]) & (centro >= x[2:])]


def exponente_lyapunov(params, y0=ESTADO_INICIAL, dt=DT_ELEGIDO, t_transitorio=2000.0,
                       t_promedio=4000.0, tau=2.0, d0=1e-8):
    """
    Máximo exponente de Lyapunov por el método de dos trayectorias (Benettin):
    se integra con RK2 una trayectoria de referencia y otra separada d0, se acumula
    log(d / d0) cada tau unidades de tiempo y se renormaliza la separación a d0.
        lambda_max > 0  ->  divergencia exponencial de trayectorias vecinas (caos)
        lambda_max ~ 0  ->  ciclo límite (órbita periódica)
        lambda_max < 0  ->  equilibrio estable
    """
    _, sol = resolver_rk2(sistema_cosner, y0, 0.0, t_transitorio, dt, params)
    y = sol[-1]
    v = np.ones(3) / np.sqrt(3.0) * d0
    suma = 0.0
    n = int(t_promedio / tau)
    for _ in range(n):
        _, ref = resolver_rk2(sistema_cosner, y, 0.0, tau, dt, params)
        _, pert = resolver_rk2(sistema_cosner, y + v, 0.0, tau, dt, params)
        y, separacion = ref[-1], pert[-1] - ref[-1]
        d = np.linalg.norm(separacion)
        suma += np.log(d / d0)
        v = separacion * (d0 / d)
    return suma / (n * tau)


def figura_5_bifurcacion_K(params, K_hopf):
    """
    Diagrama de bifurcación respecto de la capacidad de carga K: para cada K se
    descarta el transitorio y se grafican los máximos locales de N(t) y S(t).
        1 punto   -> equilibrio estable (se grafica el valor de equilibrio)
        p puntos  -> ciclo límite de período p
        nube      -> régimen aperiódico (caótico)
    """
    valores_K = np.arange(9.0, 15.001, 0.025)
    t_transitorio, t_total = 2000.0, 3000.0
    fig, ejes = plt.subplots(1, 2, figsize=(11, 4.3))
    for K in valores_K:
        t, sol = resolver_rk2(sistema_cosner, ESTADO_INICIAL, T_INICIO, t_total, DT_ELEGIDO, dict(params, K=K))
        cola = sol[t > t_transitorio]
        for ax, indice, color in [(ejes[0], 0, COLORES['N']), (ejes[1], 2, COLORES['S'])]:
            maximos = maximos_locales(cola[:, indice])
            if maximos.size == 0 or np.ptp(cola[:, indice]) < 1e-4:
                maximos = cola[-1:, indice]
            ax.plot(np.full(maximos.size, K), maximos, ',', color=color, alpha=0.8)
            ax.plot(np.full(maximos.size, K), maximos, '.', color=color, ms=1.2, alpha=0.6)
    for ax, especie in zip(ejes, ['presa $N$', 'superdepredador $S$']):
        ax.axvline(K_hopf, color='red', ls='--', lw=1, label=f'Hopf $K_H \\approx$ {K_hopf:.2f}')
        ax.axvline(params['K'], color='gray', ls=':', lw=1.2, label=f'Valor adoptado $K$ = {params["K"]:g}')
        ax.set_xlabel('Capacidad de carga $K$')
        ax.set_ylabel(f'Máximos locales de {especie}')
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(loc='upper left', fontsize=9)
    ejes[0].set_title('Diagrama de bifurcación: presa')
    ejes[1].set_title('Diagrama de bifurcación: superdepredador')
    guardar(fig, 'fig5_bifurcacion_K.png')


def figura_6_caos(params, K_caos, lyap):
    """
    Régimen caótico (K = K_caos): atractor extraño y sensibilidad a condiciones iniciales.
    """
    p = dict(params, K=K_caos)
    t, sol = resolver_rk2(sistema_cosner, ESTADO_INICIAL, T_INICIO, 3000.0, DT_ELEGIDO, p)
    perturbado = ESTADO_INICIAL + np.array([1e-6, 0.0, 0.0])
    _, sol_p = resolver_rk2(sistema_cosner, perturbado, T_INICIO, 3000.0, DT_ELEGIDO, p)
    separacion = np.max(np.abs(sol - sol_p), axis=1)

    fig = plt.figure(figsize=(13, 4.2))
    ax3d = fig.add_subplot(1, 3, 1, projection='3d')
    cola = t > 1500.0
    ax3d.plot(sol[cola, 0], sol[cola, 1], sol[cola, 2], color='purple', lw=0.4)
    ax3d.set_xlabel('$N$')
    ax3d.set_ylabel('$P$')
    ax3d.set_zlabel('$S$')
    ax3d.set_title(f'Atractor extraño ($K$ = {K_caos:g})')

    ax = fig.add_subplot(1, 3, 2)
    ventana = (t >= 1500.0) & (t <= 2500.0)
    ax.plot(t[ventana], sol[ventana, 0], lw=1.2, label='$y_0$')
    ax.plot(t[ventana], sol_p[ventana, 0], '--', lw=1.2, label='$y_0 + (10^{-6}, 0, 0)$')
    ax.set_xlabel('Tiempo $t$')
    ax.set_ylabel('Presa $N(t)$')
    ax.set_title('Sensibilidad a condiciones iniciales')
    ax.set_ylim(0.0, 13.0)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='upper center', ncol=2, fontsize=9)

    ax = fig.add_subplot(1, 3, 3)
    ax.semilogy(t, separacion, lw=0.8, color='tab:red')
    # Recta de referencia con pendiente lambda, anclada en la media geométrica de la separación en t ~ 400
    ancla = np.exp(np.mean(np.log(separacion[(t > 300.0) & (t < 500.0)])))
    t_ref = np.array([400.0, 1800.0])
    ax.semilogy(t_ref, ancla * np.exp(lyap * (t_ref - 400.0)), 'k--', lw=1.2,
                label=f'$e^{{\\lambda t}}$, $\\lambda$ = {lyap:.4f}')
    ax.set_xlabel('Tiempo $t$')
    ax.set_ylabel('$\\|y(t) - \\tilde{y}(t)\\|_\\infty$')
    ax.set_title('Divergencia de trayectorias vecinas')
    ax.grid(True, which='both', linestyle='--', alpha=0.4)
    ax.legend(loc='lower right', fontsize=9)
    fig.tight_layout()
    guardar(fig, 'fig6_regimen_caotico.png')


# =============================================================================
# 4. EJECUCIÓN
# =============================================================================

if __name__ == '__main__':
    t, sol = resolver_rk2(sistema_cosner, ESTADO_INICIAL, T_INICIO, T_FINAL, DT_ELEGIDO, PARAMETROS)
    equilibrio = equilibrio_newton(sol[-1], PARAMETROS)
    autovalores = np.linalg.eigvals(jacobiano(equilibrio, PARAMETROS))
    print(f"Equilibrio interior (Newton) (N*, P*, S*) = {np.round(equilibrio, 4)}")
    print(f"Residuo ||f(y*)||_inf = {np.max(np.abs(sistema_cosner(0.0, equilibrio, PARAMETROS))):.2e}")
    print(f"Distancia ||y(T_FINAL) - y*||_inf = {np.max(np.abs(sol[-1] - equilibrio)):.2e}")
    print(f"Autovalores del Jacobiano en y*: {np.round(autovalores, 4)}")

    dt_euler = paso_maximo_estable(autovalores, lambda z: 1 + z)
    dt_rk2 = paso_maximo_estable(autovalores, lambda z: 1 + z + z ** 2 / 2)
    print(f"Paso máximo linealmente estable en y*: Euler = {dt_euler:.2f}, RK2 = {dt_rk2:.2f}")

    # Equilibrio de borde sin superdepredador (S = 0 es invariante)
    borde = equilibrio_newton([1.6, 5.9, 0.0], PARAMETROS)
    autovalores_borde = np.linalg.eigvals(jacobiano(borde, PARAMETROS))
    print(f"Equilibrio de borde (N, P, 0) = {np.round(borde, 4)}, autovalores: {np.round(autovalores_borde, 4)}")
    print(f"  |1 + lambda*dt| de Euler con dt = 2 en el borde: {np.round(np.abs(1 + 2.0 * autovalores_borde), 3)}")

    K_hopf = K_critico_hopf(PARAMETROS, 10.0, 12.0, equilibrio)
    print(f"Bifurcación de Hopf en K_H = {K_hopf:.3f}")

    print("\nMáximo exponente de Lyapunov (RK2, dt = 0.05):")
    exponentes = {}
    for K in [10.0, 11.5, 12.0, 12.25, 13.5]:
        exponentes[K] = exponente_lyapunov(dict(PARAMETROS, K=K))
        print(f"  K = {K:>5}: lambda_max = {exponentes[K]:+.4f}")

    print("\nGenerando figuras:")
    figura_1_evolucion(t, sol)
    figura_2_fases_3d(sol)
    figura_3_euler_vs_rk2(PARAMETROS)
    figura_4_orden_convergencia(PARAMETROS)
    figura_5_bifurcacion_K(PARAMETROS, K_hopf)
    figura_6_caos(PARAMETROS, 13.5, exponentes[13.5])
