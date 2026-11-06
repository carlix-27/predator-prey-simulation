"""Refinamiento y sensibilidad temporal del exponente máximo para K=13.5."""
from generar_figuras import PARAMETROS, exponente_lyapunov

if __name__ == '__main__':
    for h in (0.05, 0.025, 0.0125):
        for transitorio, promedio in ((2000, 4000), (4000, 8000)):
            valor = exponente_lyapunov(dict(PARAMETROS, K=13.5), dt=h,
                                       t_transitorio=transitorio, t_promedio=promedio)
            print(f'h={h:g} transitorio={transitorio} promedio={promedio}: {valor:.6f}', flush=True)
