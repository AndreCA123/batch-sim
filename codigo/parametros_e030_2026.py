# -*- coding: utf-8 -*-
"""
parametros_e030_2026.py
-----------------------
Parámetros de la Norma Técnica E.030 "Diseño Sismorresistente" modificada por
la R.M. N° 183-2026-VIVIENDA (El Peruano, 3 de mayo de 2026).

  * Tabla N° 1  : factores de zona Z.
  * Tablas N° 4 y 5: S, TP y TL. Para S2 y S3 la norma da intervalos que se
    interpolan con Vs30; sin ese dato se toma el mayor S del intervalo,
    TP = 0,6 s / TL = 2,0 s (S2) y TP = 0,9 s / TL = 1,6 s (S3) (nota (*)).
    Es la opción adoptada en la población paramétrica (no se dispone de Vs30
    por edificio).
  * Tabla N° 6  : factor C con la nueva rama de periodos cortos
    (T < 0,2 TP: C = 1 + 7,5 T/TP). Para el cortante estático se usa C = 2,5
    en 0 <= T <= TP (art. 18.3 y 34.1).

Zonificación de la provincia de Arequipa (Anexo II): 21 distritos en Zona 3
(Arequipa, Cayma, Cerro Colorado, J. L. Bustamante y Rivero, Paucarpata,
Socabaya, Yanahuara, etc.) y 8 en Zona 4 (La Joya, Uchumayo, Vítor,
Yarabamba, Polobaya y los tres Siguas).
"""

Z_ZONA = {4: 0.45, 3: 0.35, 2: 0.25, 1: 0.10}

# S por zona y perfil (valor superior del intervalo cuando no hay Vs30)
S_SUELO = {
    4: {'S0': 0.80, 'S1': 1.00, 'S2': 1.10, 'S3': 1.20},
    3: {'S0': 0.80, 'S1': 1.00, 'S2': 1.15, 'S3': 1.20, 'S4': 1.30},
}
TP_TL = {'S0': (0.3, 3.0), 'S1': (0.4, 2.5), 'S2': (0.6, 2.0),
         'S3': (0.9, 1.6), 'S4': (1.2, 1.6)}

# Estratos de sitio de la tesis (mismo orden que poblacion_lhs.ESTRATOS_SITIO,
# para conservar las semillas del muestreo LHS)
ESTRATOS_SITIO_2026 = []
for zona in (4, 3):
    for suelo in ('S1', 'S2', 'S3'):
        tp, tl = TP_TL[suelo]
        ESTRATOS_SITIO_2026.append({'zona': zona, 'suelo': suelo, 'Z': Z_ZONA[zona],
                                    'S': S_SUELO[zona][suelo], 'Tp': tp, 'TL': tl})


def C_estatico(T, Tp, TL):
    """Factor C para el cortante estático (art. 18.3 / 34.1)."""
    if T <= Tp:
        return 2.5
    if T < TL:
        return 2.5 * Tp / T
    return 2.5 * Tp * TL / T ** 2


def C_espectral(T, Tp, TL):
    """Factor C completo de la Tabla N° 6 (incluye rama de periodos cortos)."""
    if T < 0.2 * Tp:
        return 1.0 + 7.5 * T / Tp
    return C_estatico(T, Tp, TL)


def espectro_e030(T, Z, U, S, Tp, TL, R=1.0, g=1.0):
    """Sa = Z·U·C·S/R · g (art. 41). Con R = 1 es el espectro elástico que se usa
    para escalar registros (art. 47.5)."""
    return Z * U * C_espectral(T, Tp, TL) * S / R * g


if __name__ == '__main__':
    print("Estratos E.030-2026:")
    for e in ESTRATOS_SITIO_2026:
        print(" ", e)
