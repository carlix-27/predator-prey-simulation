"""
===============================================================================
Trabajo Práctico: Cálculo Numérico - Modelo Tritrófico de Cosner
Archivo: construir_presentacion.py
Propósito: Generar PRESENTACION.html a partir de presentacion_fuente.html
===============================================================================
Descripción:
    1. Convierte las fórmulas LaTeX \\( ... \\) y \\[ ... \\] a MathML con pandoc,
       de modo que se vean sin conexión a internet.
    2. Embebe las figuras de docs/figuras/ como data URI.
    3. Inserta el bucle temporal real de resolver_rk2 (rk2_cosner.py) con
       resaltado de sintaxis, para que la diapositiva de código no quede
       desactualizada.

Requisitos:
    pip install pypandoc_binary

Ejecución (después de python docs/generar_figuras.py):
    python docs/construir_presentacion.py
===============================================================================
"""

import base64
import html
import re
from pathlib import Path

import pypandoc

DOCS = Path(__file__).resolve().parent
FUENTE = DOCS / 'presentacion_fuente.html'
SALIDA = DOCS / 'PRESENTACION.html'

BLOQUES_INTACTOS = re.compile(r'(<script\b.*?</script>|<style\b.*?</style>)', re.S)
FORMULA = re.compile(r'\\\[(.+?)\\\]|\\\((.+?)\\\)', re.S)
FIGURA = re.compile(r'src="(figuras/[^"]+\.png)"')
COMENTARIO_FUENTE = re.compile(r'<!-- Fuente de la presentación\..*?-->', re.S)


def formulas_a_mathml(texto):
    """Convierte todas las fórmulas fuera de <script> y <style> en una sola llamada a pandoc."""
    partes = BLOQUES_INTACTOS.split(texto)
    formulas = []
    for i in range(0, len(partes), 2):
        for m in FORMULA.finditer(partes[i]):
            formulas.append(('$$' + m.group(1).strip() + '$$') if m.group(1) is not None else ('$' + m.group(2).strip() + '$'))
    markdown = '\n\n'.join(f'FORMULA{i} {f}' for i, f in enumerate(formulas))
    salida = pypandoc.convert_text(markdown, 'html5', format='markdown', extra_args=['--mathml'])
    mathml = {int(n): cuerpo.strip() for n, cuerpo in re.findall(r'<p>FORMULA(\d+)\s*(.*?)</p>', salida, re.S)}
    faltan = [formulas[i] for i in range(len(formulas)) if not mathml.get(i, '').startswith('<math')]
    if faltan:
        raise ValueError(f"pandoc no pudo convertir {len(faltan)} fórmulas, por ejemplo: {faltan[:3]}")

    contador = iter(range(len(formulas)))
    for i in range(0, len(partes), 2):
        partes[i] = FORMULA.sub(lambda m: mathml[next(contador)], partes[i])
    return ''.join(partes), len(formulas)


def embeber_figuras(texto):
    def a_data_uri(m):
        datos = base64.b64encode((DOCS / m.group(1)).read_bytes()).decode('ascii')
        return f'src="data:image/png;base64,{datos}"'
    return FIGURA.sub(a_data_uri, texto)


def resaltar_python(codigo):
    """Resaltado mínimo: comentarios, cadenas, palabras clave y números."""
    patron = re.compile(r"('(?:[^'\\]|\\.)*')|\b(for|in|range|def|return|if|else|while|import)\b|\b(\d+\.?\d*)\b")
    clases = {1: 'st', 2: 'kw', 3: 'nu'}
    lineas = []
    for linea in codigo.split('\n'):
        cuerpo, _, comentario = linea.partition('#')
        cuerpo = patron.sub(lambda m: next(f'<span class="{clases[g]}">{m.group(g)}</span>'
                                           for g in clases if m.group(g) is not None), html.escape(cuerpo, quote=False))
        if comentario:
            cuerpo += f'<span class="cm">#{html.escape(comentario, quote=False)}</span>'
        lineas.append(cuerpo)
    return '\n'.join(lineas)


def extracto_rk2():
    """Bucle temporal de resolver_rk2, tal como está en rk2_cosner.py."""
    lineas = (DOCS / 'rk2_cosner.py').read_text(encoding='utf-8').split('\n')
    inicio = next(i for i, l in enumerate(lineas) if l.strip().startswith('for n in range(pasos):'))
    fin = next(i for i in range(inicio, len(lineas)) if lineas[i].strip().startswith('solucion[n + 1] = siguiente'))
    bloque = [l[4:] if l.startswith('    ') else l for l in lineas[inicio:fin + 1]]
    return resaltar_python('\n'.join(l.rstrip() for l in bloque))


if __name__ == '__main__':
    texto = FUENTE.read_text(encoding='utf-8')
    texto = COMENTARIO_FUENTE.sub('<!-- Generado por construir_presentacion.py a partir de presentacion_fuente.html. No editar a mano. -->', texto)
    texto, n = formulas_a_mathml(texto)
    texto = embeber_figuras(texto)
    texto = texto.replace('<!--CODIGO:rk2-->', extracto_rk2())
    SALIDA.write_text(texto, encoding='utf-8')
    print(f"{n} fórmulas convertidas a MathML -> {SALIDA} ({SALIDA.stat().st_size / 1024:.0f} KB)")
