# -*- coding: utf-8 -*-
"""
diseno_integrado.py
-------------------
Etapa 2 completa del flujo de diseño de una edificación:

  1. Secciones por rigidez: diseno_automatizado_v2.disenar_edificio_automatico
     (deriva E.030 <= 0.007).
  2. Geometría mínima de vigas E.060 21.5.1.3: bw >= 250 mm y bw >= 0.3 h.
  3. Diseño del refuerzo E.060 (fuerzas_diseno.disenar_refuerzo_edificio).
  4. Si el refuerzo no es posible, se agranda la sección responsable y se
     repite desde 3 (el modelo se vuelve a analizar, por lo que la deriva y las
     fuerzas se recalculan con las nuevas secciones):
       - acero de viga excede el máximo  -> h_viga + 5 cm (y bw si hace falta)
       - columna fuerte–viga débil / rho > 6% -> b_col y h_col + 5 cm
  Las secciones solo crecen, así que la deriva sigue cumpliendo (se verifica).

Devuelve las secciones finales, el resultado del modelo elástico, el refuerzo
y el historial de ajustes.
"""
import io
import math
import contextlib

from diseno_automatizado_v2 import disenar_edificio_automatico
from fuerzas_diseno import disenar_refuerzo_edificio

B_MIN_VIGA = 0.25           # m, E.060 21.5.1.3


def _r5(x):
    return round(math.ceil(round(x / 0.05, 6)) * 0.05, 2)


def _ajustar_ancho_viga(b, h):
    return max(b, B_MIN_VIGA, _r5(0.3 * h))


def disenar_edificio_completo(geo, sitio, max_ajustes=12, deriva_limite=0.007,
                              verbose=False):
    """geo: n_pisos, h_piso, n_vanos_x, lx, n_vanos_y, ly
       sitio: fc, CV, CM_losa, CM_acabados, CM_tabiques, Z, U, S, Tp, TL, R"""
    with contextlib.redirect_stdout(io.StringIO()):
        sec = disenar_edificio_automatico(**geo, **sitio, verbose=False)
    b_v, h_v, b_c, h_c = sec['b_viga'], sec['h_viga'], sec['b_col'], sec['h_col']
    historial = [('deriva', b_v, h_v, b_c, h_c)]
    b_v = _ajustar_ancho_viga(b_v, h_v)
    if b_v != sec['b_viga']:
        historial.append(('ancho_min_viga', b_v, h_v, b_c, h_c))

    d = None
    for _ in range(max_ajustes):
        p = dict(**geo, **sitio, b_viga=b_v, h_viga=h_v, b_col=b_c, h_col=h_c)
        with contextlib.redirect_stdout(io.StringIO()):
            d = disenar_refuerzo_edificio(p, verbose=False)
        falla_v = any(not r['ok'] for r in d['vigas'].values())
        falla_c = any(not r['ok'] for r in d['columnas'].values())
        deriva_ok = d['modelo']['deriva_max'] <= deriva_limite + 1e-9
        if not (falla_v or falla_c) and deriva_ok:
            break
        if falla_v:
            h_v = round(h_v + 0.05, 2)
            b_v = _ajustar_ancho_viga(b_v, h_v)
        if falla_c or not deriva_ok:
            b_c, h_c = round(b_c + 0.05, 2), round(h_c + 0.05, 2)
        historial.append(('resistencia', b_v, h_v, b_c, h_c))
        if verbose:
            print(f"  ajuste -> viga {b_v:.2f}x{h_v:.2f}, col {b_c:.2f}x{h_c:.2f}")

    return {
        'convergio_deriva': sec['convergio'],
        'iteraciones_deriva': sec['iteraciones'],
        'b_viga': b_v, 'h_viga': h_v, 'b_col': b_c, 'h_col': h_c,
        'modelo': d['modelo'], 'vigas': d['vigas'], 'columnas': d['columnas'],
        'refuerzo_ok': d['ok'], 'deriva_ok': d['modelo']['deriva_max'] <= deriva_limite + 1e-9,
        'n_ajustes': len(historial) - 1, 'historial': historial,
    }
