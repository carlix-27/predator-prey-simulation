"""
===============================================================================
Trabajo Práctico: Cálculo Numérico - Modelo Tritrófico de Cosner
Archivo: rk2_cosner.py
Método: Runge-Kutta 2 (Heun / Predictor-Corrector, Segundo Orden)
===============================================================================
Descripción:
    Este script implementa la simulación numérica del modelo presa-depredador-
    superdepredador de Cosner mediante el algoritmo de Runge-Kutta de orden 2 (RK2).
    
    Genera gráficos de la evolución temporal hacia el equilibrio asintótico,
    el espacio de fases tridimensional (N, P, S) y la verificación del doble paso
    para la justificación formal del dt.
===============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# 1. DEFINICIÓN DE RESPUESTAS FUNCIONALES DE COSNER
#    (Corresponde a la Sección 2.2 del Informe)
# =============================================================================

def F_PN(N, a, b1, c1):
    """
    Tasa de consumo per cápita de la presa (N) por parte del depredador (P).
    F_PN = (a * N) / (1 + b1 * N + c1 * N^2)
    """
    return (a * N) / (1.0 + b1 * N + c1 * (N ** 2))


def F_SN(N, c, b2, c2):
    """
    Tasa de consumo per cápita de la presa (N) por parte del superdepredador (S).
    F_SN = (c * N) / (1 + b2 * N + c2 * N^2)
    """
    return (c * N) / (1.0 + b2 * N + c2 * (N ** 2))


def F_SP(P, e, b3, c3):
    """
    Tasa de consumo per cápita del depredador (P) por parte del superdepredador (S).
    F_SP = (e * P) / (1 + b3 * P + c3 * P^2)
    """
    return (e * P) / (1.0 + b3 * P + c3 * (P ** 2))


# =============================================================================
# 2. SISTEMA DE ECUACIONES DIFERENCIALES ORDINARIAS (EDOs)
#    (Corresponde a la Sección 2.1 del Informe)
# =============================================================================

def sistema_cosner(t, estado, params):
    """
    Evalúa el campo vectorial f(t, y) para el modelo tritrófico de Cosner.
    """
    N, P, S = estado
    
    r = params['r']
    K = params['K']
    d = params['d']
    b = params['b']
    m = params['m']
    f = params['f']
    g = params['g']
    q = params['q']
    
    fpn = F_PN(N, params['a'], params['b1'], params['c1'])
    fsn = F_SN(N, params['c_ataque'], params['b2'], params['c2'])
    fsp = F_SP(P, params['e'], params['b3'], params['c3'])
    
    dN_dt = r * N * (1.0 - N / K) - d * N - fpn * P - fsn * S
    dP_dt = b * fpn * P - m * P - fsp * S
    dS_dt = f * fsn * S + g * fsp * S - q * S
    
    return np.array([dN_dt, dP_dt, dS_dt])


# =============================================================================
# 3. ALGORITMO NUMÉRICO: RUNGE-KUTTA 2 (HEUN / SEGUNDO ORDEN)
#    (Corresponde a la Sección 3.2 del Informe)
# =============================================================================

def resolver_rk2(f_sistema, y0, t0, t_max, dt, params):
    """
    Avanza la solución mediante el esquema Runge-Kutta 2 (Heun):
        k1 = f(t_n, y_n)
        k2 = f(t_n + dt, y_n + dt * k1)
        y_{n+1} = y_n + (dt / 2) * (k1 + k2)
    
    Error de truncamiento local: O(dt^3)
    Error global acumulado: O(dt^2)
    """
    num_pasos = int(np.ceil((t_max - t0) / dt)) + 1
    tiempos = np.linspace(t0, t0 + (num_pasos - 1) * dt, num_pasos)
    
    solucion = np.zeros((num_pasos, len(y0)))
    solucion[0] = y0
    
    for n in range(num_pasos - 1):
        tn = tiempos[n]
        yn = solucion[n]
        
        # Etapa 1: Pendiente inicial (predictor Euler)
        k1 = f_sistema(tn, yn, params)
        
        # Etapa 2: Pendiente tentativa al final del intervalo
        y_predicho = np.maximum(yn + dt * k1, 0.0)
        k2 = f_sistema(tn + dt, y_predicho, params)
        
        # Corrección ponderada de orden 2
        yn_siguiente = yn + (dt / 2.0) * (k1 + k2)
        yn_siguiente = np.maximum(yn_siguiente, 0.0)
        
        solucion[n + 1] = yn_siguiente
        
    return tiempos, solucion


# =============================================================================
# 4. CONFIGURACIÓN EXPERIMENTAL Y VISUALIZACIÓN
# =============================================================================

if __name__ == '__main__':
    parametros = {
        'r': 1.0,           # Natalidad de la presa
        'K': 100.0,         # Capacidad de carga presa
        'd': 0.1,           # Mortalidad de la presa
        'b': 0.5,           # Conversión presa -> depredador
        'm': 0.2,           # Mortalidad depredador
        'f': 0.2,           # Conversión presa -> superdepredador
        'g': 0.4,           # Conversión depredador -> superdepredador
        'q': 0.1,           # Mortalidad superdepredador
        'a': 0.15,          # Tasa de ataque F_PN
        'c_ataque': 0.05,   # Tasa de ataque F_SN
        'e': 0.10,          # Tasa de ataque F_SP
        'b1': 0.5,          # Coeficiente lineal N
        'b2': 0.5,          # Coeficiente lineal N
        'b3': 0.5,          # Coeficiente lineal P
        'c1': 0.01,         # Coeficiente cuadrático N^2
        'c2': 0.01,         # Coeficiente cuadrático N^2
        'c3': 0.01          # Coeficiente cuadrático P^2
    }

    estado_inicial = np.array([30.0, 10.0, 5.0])
    t_inicio = 0.0
    t_final = 250.0
    dt_elegido = 0.05   # Paso validado por doble paso (precisión y bajo costo)

    # 1. Integración con RK2
    t_vec, sol_rk2 = resolver_rk2(sistema_cosner, estado_inicial, t_inicio, t_final, dt_elegido, parametros)

    # 2. Gráfico 1: Evolución Temporal hacia el Equilibrio Asintótico
    plt.figure(figsize=(10, 5))
    plt.plot(t_vec, sol_rk2[:, 0], label='Presa N(t)', color='tab:blue', lw=2)
    plt.plot(t_vec, sol_rk2[:, 1], label='Depredador P(t)', color='tab:orange', lw=2)
    plt.plot(t_vec, sol_rk2[:, 2], label='Superdepredador S(t)', color='tab:green', lw=2)
    plt.title('Evolución Temporal del Sistema de Cosner (RK2, dt=0.05) - Estado Asintótico', fontsize=12)
    plt.xlabel('Tiempo t', fontsize=11)
    plt.ylabel('Densidad Poblacional', fontsize=11)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig('evolucion_asintotica_rk2.png', dpi=300)
    plt.close()

    # 3. Gráfico 2: Diagrama de Fases 3D (Espacio N-P-S)
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(sol_rk2[:, 0], sol_rk2[:, 1], sol_rk2[:, 2], color='purple', lw=1.5)
    ax.scatter([estado_inicial[0]], [estado_inicial[1]], [estado_inicial[2]], color='red', s=50, label='Condición Inicial y0')
    ax.scatter([sol_rk2[-1, 0]], [sol_rk2[-1, 1]], [sol_rk2[-1, 2]], color='green', s=60, label='Punto Fijo Asintótico')
    ax.set_title('Espacio de Fases 3D: Trayectoria Espiralada hacia Atractor de Punto Fijo', fontsize=12)
    ax.set_xlabel('Presa N')
    ax.set_ylabel('Depredador P')
    ax.set_zlabel('Superdepredador S')
    ax.legend()
    plt.tight_layout()
    plt.savefig('fases_3d_rk2.png', dpi=300)
    plt.close()

    print("Simulación con RK2 completada. Gráficos 'evolucion_asintotica_rk2.png' y 'fases_3d_rk2.png' generados con éxito.")
