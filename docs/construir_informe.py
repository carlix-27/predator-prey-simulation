"""
===============================================================================
Trabajo Práctico: Cálculo Numérico - Modelo Tritrófico de Cosner
Archivo: construir_informe.py
Propósito: Generar el informe final a partir de INFORME_PRELIMINAR.md
===============================================================================
Descripción:
    1. Reescribe el Anexo A del Markdown con el contenido actual de los scripts
       (euler_cosner.py, rk2_cosner.py, generar_figuras.py).
    2. Genera INFORME_EDITABLE.docx: ecuaciones nativas de Word (editables y
       compatibles con Google Docs), figuras embebidas y código como anexo.
    3. Genera INFORME_EDITABLE.html autocontenido (MathML, imágenes embebidas,
       sin dependencias de internet).
    4. Genera INFORME.pdf imprimiendo el HTML con Microsoft Edge o Google Chrome
       en modo headless (se omite si no hay ninguno instalado).

Requisitos:
    pip install pypandoc_binary      (incluye el ejecutable de pandoc)

Ejecución (después de python docs/generar_figuras.py):
    python docs/construir_informe.py
===============================================================================
"""

import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import pypandoc

DOCS = Path(__file__).resolve().parent
FUENTE = DOCS / 'INFORME_PRELIMINAR.md'
PLANTILLA_WORD = DOCS / 'plantilla_word.docx'
SALIDA_DOCX = DOCS / 'INFORME_EDITABLE.docx'
SALIDA_HTML = DOCS / 'INFORME_EDITABLE.html'
SALIDA_PDF = DOCS / 'INFORME.pdf'

ANEXO = [
    ('A.1', 'euler_cosner.py', 'Método de Euler explícito'),
    ('A.2', 'rk2_cosner.py', 'Método de Runge-Kutta 2 (Heun)'),
    ('A.3', 'generar_figuras.py', 'Análisis, figuras y tablas del informe'),
]
MARCA_INICIO = '<!-- ANEXO:INICIO'
MARCA_FIN = '<!-- ANEXO:FIN -->'
TITULO_ANEXO = '## Anexo A: Código Fuente'

SALTO_DE_PAGINA_WORD = '```{=openxml}\n<w:p><w:r><w:br w:type="page"/></w:r></w:p>\n```\n\n'

ESTILO_HTML = """
body { font-family: Georgia, 'Times New Roman', serif; line-height: 1.55; color: #222;
       max-width: 820px; margin: 2em auto; padding: 0 1.2em; background: #fff; }
header, header p { text-align: center; margin-bottom: 2em; }
h1.title { font-family: Arial, sans-serif; font-size: 1.6em; margin-bottom: .3em; }
p.subtitle { font-family: Arial, sans-serif; font-size: 1.1em; color: #444; margin: 0; }
p.author, p.date { font-family: Arial, sans-serif; color: #555; margin: .25em 0; }
h2 { font-family: Arial, sans-serif; color: #1a365d; border-bottom: 1px solid #cbd5e0;
     padding-bottom: .2em; margin-top: 1.8em; }
h3 { font-family: Arial, sans-serif; color: #2b6cb0; margin-top: 1.3em; }
h2, h3 { break-after: avoid; page-break-after: avoid; }
p, li { text-align: justify; }
figure { margin: 1.2em 0; text-align: center; break-inside: avoid; }
figure img { max-width: 100%; max-height: 9.5cm; height: auto; }
figcaption { font-family: Arial, sans-serif; font-size: .88em; color: #333; font-style: italic;
             margin-top: .4em; }
table { border-collapse: collapse; margin: 1em auto; font-family: Arial, sans-serif; font-size: .9em; }
th, td { border: 1px solid #cbd5e1; padding: .35em .7em; }
th { background: #f1f5f9; }
math[display="block"] { margin: .8em 0; }
pre, code { font-family: Consolas, 'Courier New', monospace; }
pre { font-size: .72em; background: #f5f7fa; border: 1px solid #e2e8f0; padding: .8em;
      white-space: pre-wrap; word-break: break-word; }
h2[id^="anexo"] { break-before: page; }
@page { size: letter; margin: 2.2cm; }
@media print { body { margin: 0; max-width: none; } }
"""


def actualizar_anexo(texto):
    """Reemplaza el contenido entre las marcas del anexo por el código actual."""
    bloques = []
    for numero, archivo, descripcion in ANEXO:
        codigo = (DOCS / archivo).read_text(encoding='utf-8').rstrip()
        bloques.append(f"### {numero}. `{archivo}`: {descripcion}\n\n```python\n{codigo}\n```\n")
    inicio = texto.index(MARCA_INICIO)
    fin_linea_inicio = texto.index('\n', inicio)
    fin = texto.index(MARCA_FIN)
    return texto[:fin_linea_inicio + 1] + '\n' + '\n'.join(bloques) + '\n' + texto[fin:]


def preparar_para_pandoc(texto):
    """
    En el Markdown (pensado para GitHub) cada figura va seguida de un párrafo en
    cursiva con su epígrafe. Para Word/HTML ese epígrafe pasa a ser el caption
    real de la figura, de modo que quede centrado y unido a la imagen.
    """
    return re.sub(r'!\[[^\]]*\]\(([^)]+)\)\n\n\*(.+?)\*\n', r'![\2](\1)\n', texto)


def ejecutar_pandoc(texto, formato, salida, extra):
    pypandoc.convert_text(
        texto, formato, format='markdown', outputfile=str(salida),
        extra_args=['--standalone', f'--resource-path={DOCS}', '--syntax-highlighting=tango'] + extra,
    )
    print(f"  -> {salida}")


def buscar_navegador():
    candidatos = [
        r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe',
        r'C:\Program Files\Microsoft\Edge\Application\msedge.exe',
        r'C:\Program Files\Google\Chrome\Application\chrome.exe',
        '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
        '/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge',
        '/Applications/Chromium.app/Contents/MacOS/Chromium',
        shutil.which('msedge'), shutil.which('google-chrome'), shutil.which('chromium'),
    ]
    return next((c for c in candidatos if c and os.path.exists(c)), None)


def imprimir_pdf(html, pdf):
    navegador = buscar_navegador()
    if navegador is None:
        print("  (PDF omitido: no se encontró Edge ni Chrome)")
        return
    with tempfile.TemporaryDirectory() as perfil:
        subprocess.run([navegador, '--headless', '--disable-gpu', '--no-pdf-header-footer',
                        f'--user-data-dir={perfil}', f'--print-to-pdf={pdf}', html.as_uri()],
                       check=True, capture_output=True, timeout=180)
    print(f"  -> {pdf}")


if __name__ == '__main__':
    texto = actualizar_anexo(FUENTE.read_text(encoding='utf-8'))
    FUENTE.write_text(texto, encoding='utf-8')
    print(f"Anexo actualizado en {FUENTE}")

    texto_pandoc = preparar_para_pandoc(texto)
    print("Generando documentos:")

    texto_word = texto_pandoc.replace(TITULO_ANEXO, SALTO_DE_PAGINA_WORD + TITULO_ANEXO)
    ejecutar_pandoc(texto_word, 'docx', SALIDA_DOCX, [f'--reference-doc={PLANTILLA_WORD}'])

    with tempfile.NamedTemporaryFile('w', suffix='.css', delete=False, encoding='utf-8') as css:
        css.write(ESTILO_HTML)
    try:
        ejecutar_pandoc(texto_pandoc, 'html5', SALIDA_HTML, ['--embed-resources', '--mathml', f'--css={css.name}'])
    finally:
        os.unlink(css.name)

    imprimir_pdf(SALIDA_HTML, SALIDA_PDF)
