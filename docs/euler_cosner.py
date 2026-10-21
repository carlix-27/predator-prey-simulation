"""
===============================================================================
Trabajo Práctico: Cálculo Numérico - Modelo Tritrófico de Cosner
Archivo: euler_cosner.py
Método: Euler Progresivo (Primer Orden)
===============================================================================
Descripción:
    Este script implementa la simulación numérica del modelo presa-depredador-
    superdepredador propuesto por Cosner con respuesta funcional racional cuadrática,
    utilizando el método de integración temporal de Euler Explícito.
    
    El código evalúa la convergencia y estabilidad frente a la variación del paso dt
    y grafica la evolución temporal de las densidades poblacionales N(t), P(t), S(t).
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
    Respuesta funcional racional cuadrática generalizada.
    """
    return (a * N) / (1.0 + b1 * N + c1 * (N ** 2))


def F_SN(N, c, b2, c2):
    """
    Tasa de consumo per cápita de la presa (N) por parte del superdepredador (S).
    """
    return (c * N) / (1.0 + b2 * N + c2 * (N ** 2))


def F_SP(P, e, b3, c3):
    """
    Tasa de consumo per cápita del depredador (P) por parte del superdepredador (S).
    """
    return (e * P) / (1.0 + b3 * P + c3 * (P ** 2))


# =============================================================================
# 2. SISTEMA DE ECUACIONES DIFERENCIALES ORDINARIAS (EDOs)
#    (Corresponde a la Sección 2.1 del Informe)
# =============================================================================

def sistema_cosner(t, estado, params):
    """
    Calcula las derivadas temporales dN/dt, dP/dt, dS/dt según el modelo de Cosner.
    
    Parámetros:
        t: Tiempo actual (float).
        estado: Vector de estado [N, P, S] con densidades poblacionales.
        params: Diccionario con todos los coeficientes biológicos del modelo.
    
    Retorna:
        np.array con [dN_dt, dP_dt, dS_dt].
    """
    N, P, S = estado
    
    # Extraer parámetros ecológicos
    r = params['r']      # Tasa de natalidad intrínseca de la presa
    K = params['K']      # Capacidad de carga del medio para la presa
    d = params['d']      # Tasa de mortalidad natural de la presa
    b = params['b']      # Eficiencia de conversión de presa en depredador
    m = params['m']      # Tasa de mortalidad natural del depredador
    f = params['f']      # Eficiencia de conversión de presa en superdepredador
    g = params['g']      # Eficiencia de conversión de depredador en superdepredador
    q = params['q']      # Tasa de mortalidad natural del superdepredador
    
    # Evaluar tasas de ataque no lineales
    fpn = F_PN(N, params['a'], params['b1'], params['c1'])
    fsn = F_SN(N, params['c_ataque'], params['b2'], params['c2'])
    fsp = F_SP(P, params['e'], params['b3'], params['c3'])
    
    # Ecuaciones de variación poblacional
    dN_dt = r * N * (1.0 - N / K) - d * N - fpn * P - fsn * S
    dP_dt = b * fpn * P - m * P - fsp * S
    dS_dt = f * fsn * S + g * fsp * S - q * S
    
    return np.array([dN_dt, dP_dt, dS_dt])


# =============================================================================
# 3. ALGORITMO NUMÉRICO: MÉTODO DE EULER EXPLÍCITO
#    (Corresponde a la Sección 3.1 del Informe)
# =============================================================================

def resolver_euler(f_sistema, y0, t0, t_max, dt, params):
    """
    Avanza la solución paso a paso mediante el esquema de Euler explícito:
        y_{n+1} = y_n + dt * f(t_n, y_n)
    
    Retorna:
        tiempos: Vector con la discretización temporal.
        solucion: Matriz de dimensiones (pasos, 3) con [N(t), P(t), S(t)].
    """
    num_pasos = int(np.ceil((t_max - t0) / dt)) + 1
    tiempos = np.linspace(t0, t0 + (num_pasos - 1) * dt, num_pasos)
    
    # Reserva de memoria para almacenar la trayectoria
    solucion = np.zeros((num_pasos, len(y0)))
    solucion[0] = y0
    
    for n in range(num_pasos - 1):
        tn = tiempos[n]
        yn = solucion[n]
        
        # Paso de avance de Euler
        derivadas = f_sistema(tn, yn, params)
        yn_siguiente = yn + dt * derivadas
        
        # Control físico: las poblaciones no pueden ser negativas
        yn_siguiente = np.maximum(yn_siguiente, 0.0)
        solucion[n + 1] = yn_siguiente
        
    return tiempos, solucion


# =============================================================================
# 4. CONFIGURACIÓN EXPERIMENTAL Y PARÁMETROS BIOLÓGICOS
#    (Garantizan estabilidad asintótica y evitan zonas de caos)
# =============================================================================

if __name__ == '__main__':
    parametros = {
        'r': 1.0,           # Natalidad de la presa
        'K': 100.0,         # Capacidad de carga
        'd': 0.1,           # Mortalidad de la presa
        'b': 0.5,           # Eficiencia conversión presa -> depredador
        'm': 0.2,           # Mortalidad depredador
        'f': 0.2,           # Eficiencia conversión presa -> superdepredador
        'g': 0.4,           # Eficiencia conversión depredador -> superdepredador
        'q': 0.1,           # Mortalidad superdepredador
        'a': 0.15,          # Tasa de ataque depredador -> presa (0.01 - 0.5)
        'c_ataque': 0.05,   # Tasa de ataque superdepredador -> presa (0.01 - 0.5)
        'e': 0.10,          # Tasa de ataque superdepredador -> depredador (0.01 - 0.5)
        'b1': 0.5,          # Coeficiente lineal N en F_PN (0.1 - 2.0)
        'b2': 0.5,          # Coeficiente lineal N en F_SN (0.1 - 2.0)
        'b3': 0.5,          # Coeficiente lineal P en F_SP (0.1 - 2.0)
        'c1': 0.01,         # Coeficiente cuadrático N^2 en F_PN (0.0 - 0.1)
        'c2': 0.01,         # Coeficiente cuadrático N^2 en F_SN (0.0 - 0.1)
        'c3': 0.01          # Coeficiente cuadrático P^2 en F_SP (0.0 - 0.1)
    }

    # Condiciones iniciales: Poblaciones iniciales no nulas
    estado_inicial = np.array([30.0, 10.0, 5.0])  # [N0, P0, S0]
    t_inicio = 0.0
    t_final = 250.0   # Tiempo suficiente para verificar estado asintótico

    # Justificación de dt: Comparación de estabilidad
    pasos_dt = [0.5, 0.1, 0.01]
    
    plt.figure(figsize=(12, 7))
    
    for dt in pasos_dt:
        t_vec, sol = resolver_euler(sistema_cosner, estado_inicial, t_inicio, t_final, dt, parametros)
        plt.plot(t_vec, sol[:, 0], label=f'Presa N(t) [dt={dt}]', alpha=0.8)
    
    plt.title('Sensibilidad del Método de Euler frente al Paso dt en el Modelo de Cosner', fontsize=13)
    plt.xlabel('Tiempo t', fontsize=11)
    plt.ylabel('Densidad de Población', fontsize=11)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig('comparacion_euler_dt.png', dpi=300)
    plt.close()
    print("Simulación con Euler completada. Gráfico guardado en 'comparacion_euler_dt.png'.")
