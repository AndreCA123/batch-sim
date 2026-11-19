"""
=============================================================================
 POBLACION PARAMETRICA - MUESTREO LHS ESTRATIFICADO POR ZONA x SUELO
=============================================================================
Paso 6 de la hoja de ruta. Genera las N edificaciones de tu poblacion de
entrenamiento y las diseña automaticamente (llamando al modulo del paso 5),
guardando resultados en un CSV que se puede RETOMAR si la sesion de Colab
se corta a la mitad (no se pierde el trabajo ya hecho).

COMO FUNCIONA EL MUESTREO (recordando lo que conversamos):
  1. Las variables de SITIO (zona, suelo) no se muestrean: se tratan como
     6 estratos fijos (2 zonas x 3 suelos), y la poblacion se reparte en
     partes iguales entre ellos.
  2. Dentro de cada estrato, las variables de DISENO (pisos, vanos, luz,
     f'c, uso) se generan con hipercubo latino (LHS): cada variable se
     divide en tantos intervalos como edificaciones le toquen a ese
     estrato, se saca un valor de cada intervalo, y se combinan sin
     repetir intervalo -- as[i] se cubre bien el rango de cada variable
     sin necesitar miles de combinaciones.
  3. Simplificacion pactada: se usa el MISMO numero de vanos y la MISMA
     luz en direccion X e Y (edificios de planta cuadrada/rectangular
     uniforme), igual que en Kazemi et al. Si mas adelante quieres
     vanos/luces distintos por direccion, se puede ampliar el muestreo
     a 7 dimensiones en vez de 5.

QUE HACE CADA ARCHIVO CSV QUE SE GENERA:
  - poblacion_150.csv    : la DEFINICION de las 150 edificaciones (que
    parametros le toco a cada una). Se genera UNA sola vez.
  - resultados_diseno.csv: el RESULTADO de diseñar cada una (secciones
    finales, periodos, derivas). Se va llenando fila por fila, y el
    script SALTA las filas que ya existen si lo vuelves a correr -- por
    eso puedes cortar la sesion de Colab a media corrida sin perder nada.
=============================================================================
"""

import os
import csv
import numpy as np
from scipy.stats import qmc

from diseno_automatizado_v2 import disenar_edificio_automatico

# -----------------------------------------------------------------------
# ESPACIO DE DISENO (Tabla 2.1 de tu plan de tesis)
# -----------------------------------------------------------------------
RANGO_PISOS = (2, 10)
RANGO_VANOS = (2, 5)
RANGO_LUZ = (3.5, 7.0)          # m
NIVELES_FC = [210.0, 280.0]     # kg/cm2

NIVELES_USO = [
    {'nombre': 'vivienda', 'CV': 0.200},
    {'nombre': 'oficina', 'CV': 0.250},
    {'nombre': 'comercio', 'CV': 0.500},
]

CM_LOSA, CM_ACABADOS, CM_TABIQUES = 0.300, 0.100, 0.150   # tonf/m2, Tabla 3
R_SISMICO = 8.0
U_USO = 1.0
H_PISO = 3.0   # m (simplificacion: altura de piso uniforme)

# -----------------------------------------------------------------------
# ESTRATOS DE SITIO (Tabla 2.2 de tu plan de tesis)
# -----------------------------------------------------------------------
ESTRATOS_SITIO = [
    {'zona': 4, 'suelo': 'S1', 'Z': 0.45, 'S': 1.00, 'Tp': 0.40, 'TL': 2.50},
    {'zona': 4, 'suelo': 'S2', 'Z': 0.45, 'S': 1.05, 'Tp': 0.60, 'TL': 2.00},
    {'zona': 4, 'suelo': 'S3', 'Z': 0.45, 'S': 1.10, 'Tp': 1.00, 'TL': 1.60},
    {'zona': 3, 'suelo': 'S1', 'Z': 0.35, 'S': 1.00, 'Tp': 0.40, 'TL': 2.50},
    {'zona': 3, 'suelo': 'S2', 'Z': 0.35, 'S': 1.15, 'Tp': 0.60, 'TL': 2.00},
    {'zona': 3, 'suelo': 'S3', 'Z': 0.35, 'S': 1.20, 'Tp': 1.00, 'TL': 1.60},
]


def generar_poblacion(n_total=150, semilla=42):
    """Genera la definicion de las N edificaciones (sin diseñarlas todavia).
    Devuelve una lista de diccionarios, uno por edificacion."""
    n_estratos = len(ESTRATOS_SITIO)
    n_por_estrato = n_total // n_estratos
    resto = n_total - n_por_estrato * n_estratos   # reparte lo que sobra

    poblacion = []
    id_edificio = 1

    for idx_estrato, estrato in enumerate(ESTRATOS_SITIO):
        n_este_estrato = n_por_estrato + (1 if idx_estrato < resto else 0)

        # Una semilla distinta por estrato para que no se repita el mismo
        # patron de muestreo en cada combinacion zona-suelo
        sampler = qmc.LatinHypercube(d=5, seed=semilla + idx_estrato)
        muestra = sampler.random(n=n_este_estrato)   # valores en [0,1), forma (n,5)

        for fila in muestra:
            u_pisos, u_vanos, u_luz, u_fc, u_uso = fila

            n_pisos = min(RANGO_PISOS[1],
                          int(RANGO_PISOS[0] + u_pisos * (RANGO_PISOS[1] - RANGO_PISOS[0] + 1)))
            n_vanos = min(RANGO_VANOS[1],
                          int(RANGO_VANOS[0] + u_vanos * (RANGO_VANOS[1] - RANGO_VANOS[0] + 1)))
            luz = round(RANGO_LUZ[0] + u_luz * (RANGO_LUZ[1] - RANGO_LUZ[0]), 2)
            fc = NIVELES_FC[0] if u_fc < 0.5 else NIVELES_FC[1]
            idx_uso = min(int(u_uso * len(NIVELES_USO)), len(NIVELES_USO) - 1)
            uso = NIVELES_USO[idx_uso]

            poblacion.append({
                'id': id_edificio,
                'n_pisos': n_pisos,
                'n_vanos_x': n_vanos, 'lx': luz,
                'n_vanos_y': n_vanos, 'ly': luz,
                'fc': fc,
                'uso': uso['nombre'], 'CV': uso['CV'],
                'zona': estrato['zona'], 'suelo': estrato['suelo'],
                'Z': estrato['Z'], 'S': estrato['S'],
                'Tp': estrato['Tp'], 'TL': estrato['TL'],
            })
            id_edificio += 1

    return poblacion


def guardar_poblacion(poblacion, ruta_csv):
    campos = list(poblacion[0].keys())
    with open(ruta_csv, 'w', newline='') as f:
        escritor = csv.DictWriter(f, fieldnames=campos)
        escritor.writeheader()
        escritor.writerows(poblacion)


def cargar_poblacion(ruta_csv):
    with open(ruta_csv, newline='') as f:
        lector = csv.DictReader(f)
        poblacion = []
        for fila in lector:
            fila['id'] = int(fila['id'])
            fila['n_pisos'] = int(fila['n_pisos'])
            fila['n_vanos_x'] = fila['n_vanos_y'] = int(fila['n_vanos_x'])
            fila['lx'] = fila['ly'] = float(fila['lx'])
            fila['fc'] = float(fila['fc'])
            fila['CV'] = float(fila['CV'])
            fila['zona'] = int(fila['zona'])
            fila['Z'] = float(fila['Z']); fila['S'] = float(fila['S'])
            fila['Tp'] = float(fila['Tp']); fila['TL'] = float(fila['TL'])
            poblacion.append(fila)
    return poblacion


def ejecutar_diseno_poblacion(poblacion, ruta_resultados, verbose_disenos=False):
    """Diseña cada edificacion de la poblacion y guarda los resultados.
    RETOMABLE: si el archivo ya tiene filas de una corrida anterior, las
    edificaciones ya procesadas se saltan."""
    ids_procesados = set()
    if os.path.exists(ruta_resultados):
        with open(ruta_resultados, newline='') as f:
            for fila in csv.DictReader(f):
                ids_procesados.add(int(fila['id']))

    escribir_encabezado = not os.path.exists(ruta_resultados)
    campos = ['id', 'n_pisos', 'n_vanos', 'luz', 'fc', 'uso', 'CV', 'zona', 'suelo',
              'b_viga', 'h_viga', 'b_col', 'h_col',
              'T1', 'T2', 'T3', 'Vx', 'Vy', 'deriva_max', 'convergio', 'iteraciones']

    n_saltados = len(ids_procesados)
    n_pendientes = len(poblacion) - n_saltados
    print(f"Edificaciones ya procesadas (se saltan): {n_saltados}")
    print(f"Edificaciones pendientes en esta corrida: {n_pendientes}")

    with open(ruta_resultados, 'a', newline='') as f:
        escritor = csv.DictWriter(f, fieldnames=campos)
        if escribir_encabezado:
            escritor.writeheader()

        for edif in poblacion:
            if edif['id'] in ids_procesados:
                continue

            print(f"[{edif['id']}/{len(poblacion)}] Diseñando: {edif['n_pisos']} pisos, "
                  f"{edif['n_vanos_x']} vanos x {edif['lx']} m, f'c={edif['fc']}, "
                  f"{edif['uso']}, Zona {edif['zona']}/{edif['suelo']}...", end=' ')

            resultado = disenar_edificio_automatico(
                n_pisos=edif['n_pisos'], h_piso=H_PISO,
                n_vanos_x=edif['n_vanos_x'], lx=edif['lx'],
                n_vanos_y=edif['n_vanos_y'], ly=edif['ly'],
                fc=edif['fc'], CV=edif['CV'],
                CM_losa=CM_LOSA, CM_acabados=CM_ACABADOS, CM_tabiques=CM_TABIQUES,
                Z=edif['Z'], U=U_USO, S=edif['S'], Tp=edif['Tp'], TL=edif['TL'],
                R=R_SISMICO, verbose=verbose_disenos,
            )

            rf = resultado['resultado_final']
            fila = {
                'id': edif['id'], 'n_pisos': edif['n_pisos'],
                'n_vanos': edif['n_vanos_x'], 'luz': edif['lx'],
                'fc': edif['fc'], 'uso': edif['uso'], 'CV': edif['CV'],
                'zona': edif['zona'], 'suelo': edif['suelo'],
                'b_viga': resultado['b_viga'], 'h_viga': resultado['h_viga'],
                'b_col': resultado['b_col'], 'h_col': resultado['h_col'],
                'T1': round(rf['periodos'][0], 4), 'T2': round(rf['periodos'][1], 4),
                'T3': round(rf['periodos'][2], 4),
                'Vx': round(rf['Vx'], 2), 'Vy': round(rf['Vy'], 2),
                'deriva_max': round(rf['deriva_max'], 5),
                'convergio': resultado['convergio'], 'iteraciones': resultado['iteraciones'],
            }
            escritor.writerow(fila)
            f.flush()   # escribe en disco DE INMEDIATO -- esto es lo que hace
                        # posible retomar si la sesion se corta a media corrida

            estado = "OK" if resultado['convergio'] else "*** NO CONVERGIO ***"
            print(f"{estado} ({resultado['iteraciones']} iter., "
                  f"col {resultado['b_col']:.2f}x{resultado['h_col']:.2f} m)")


# =============================================================================
# PRUEBA: poblacion pequeña (12 edificaciones, 2 por estrato) para verificar
# que todo el flujo corre de punta a punta antes de lanzar las 150 reales.
# =============================================================================
if __name__ == '__main__':
    poblacion = generar_poblacion(n_total=12, semilla=42)
    guardar_poblacion(poblacion, 'poblacion_prueba.csv')
    print(f"Poblacion de prueba generada: {len(poblacion)} edificaciones")
    print("-" * 74)
    for e in poblacion:
        print(f"  id={e['id']:2d}  pisos={e['n_pisos']:2d}  vanos={e['n_vanos_x']}  "
              f"luz={e['lx']:.2f} m  f'c={e['fc']:.0f}  uso={e['uso']:9s}  "
              f"Zona {e['zona']}/{e['suelo']}")
    print("-" * 74)

    ejecutar_diseno_poblacion(poblacion, 'resultados_prueba.csv')
    print("-" * 74)
    print("Listo. Revisa resultados_prueba.csv")
