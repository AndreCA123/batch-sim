# -*- coding: utf-8 -*-
"""
rutas.py
--------
Centralización de rutas del proyecto con pathlib para compatibilidad
multiplataforma (Linux, Windows y Google Colab).
"""
from pathlib import Path

# Directorio raíz del proyecto
PROYECTO = Path(__file__).resolve().parent.parent

# Subdirectorios principales
CODIGO = PROYECTO / "codigo"
DATOS = PROYECTO / "datos"
DATOS_IDA_PARTES = DATOS / "ida_partes"
RESULTADOS = PROYECTO / "resultados"
REGISTROS = PROYECTO / "registros_sismicos"
DATOS_SISMICOS = PROYECTO / "datos_sismicos"
MODELOS_ETABS = PROYECTO / "modelos_etabs_150"
FIGURAS = PROYECTO / "figuras_reporte"

# Crear carpetas de salida si no existen
for _dir in [DATOS, DATOS_IDA_PARTES, RESULTADOS, MODELOS_ETABS, FIGURAS]:
    _dir.mkdir(parents=True, exist_ok=True)
