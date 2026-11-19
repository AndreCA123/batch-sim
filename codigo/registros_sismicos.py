"""
=============================================================================
 REGISTROS SISMICOS - LECTURA, ESPECTRO DE RESPUESTA Y ESCALADO A LA E.030
=============================================================================
Paso 7: antes del IDA necesitas un conjunto de acelerogramas ya escalados
al espectro de diseno de cada combinacion zona-suelo de tu poblacion.
Este modulo hace tres cosas:

  1. LEE un acelerograma desde archivo (dos formatos: PEER NGA .AT2, que es
     el formato estandar mas usado en investigacion sismica, y un CSV
     generico de dos columnas tiempo-aceleracion, por si el formato de
     REDACIS/CISMID no calza exacto con el .AT2 -- ver nota abajo).

  2. CALCULA el espectro de respuesta elastico (pseudo-aceleracion) del
     registro, mediante integracion de Newmark (metodo de la aceleracion
     lineal), que es el metodo estandar de la dinamica estructural para
     esto -- no es una funcion "magica" de una libreria, se calcula
     resolviendo la ecuacion de movimiento de un oscilador de 1 grado de
     libertad para cada periodo de interes.

  3. ESCALA el registro para que la media de su espectro, en el rango de
     periodos 0.2*T1 a 1.5*T1, no sea menor que la media del espectro de
     diseno E.030 en ese mismo rango -- el criterio de escalado habitual
     de la norma peruana.

IMPORTANTE SOBRE DE DONDE SACAR LOS REGISTROS REALES:
  La fuente correcta para tu tesis es REDACIS-CISMID (Universidad Nacional
  de Ingenieria): https://www.cismid.uni.edu.pe/ceois/redacis/red/
  Tienen estaciones DENTRO de Arequipa (UNSA, SENCICO-Arequipa) y el
  acelerograma del sismo del 23 de junio de 2001 (Mw 8.4), el evento de
  subduccion mas fuerte registrado en el sur del Peru en decadas -- el
  candidato obligado para tu conjunto de registros. No pude verificar el
  formato EXACTO de los archivos que entrega su portal (requiere
  navegacion interactiva), asi que este modulo trae DOS lectores: prueba
  primero con leer_peer_at2() y si el archivo no calza, usa
  leer_registro_csv() ajustando las columnas. Avisame si tienes un
  archivo de muestra y ajusto el lector exacto.

EL EJEMPLO AL FINAL usa un acelerograma SINTETICO (generado con ruido
filtrado), solo para comprobar que las funciones calculan correctamente.
NO es un registro real -- no lo uses para tu IDA final.
=============================================================================
"""

import numpy as np


# =============================================================================
# 1. LECTURA DE REGISTROS
# =============================================================================
def leer_peer_at2(ruta):
    """Lee un acelerograma en formato PEER NGA (.AT2), el formato estandar
    en investigacion sismica. Estructura: 4 lineas de encabezado (la 4ta
    trae NPTS y DT), luego los valores de aceleracion en varias columnas
    por linea. Devuelve (tiempo, aceleracion) en (s, g)."""
    # latin-1: algunos encabezados NGA-Sub traen caracteres no UTF-8 (p. ej. Ñ)
    with open(ruta, 'r', encoding='latin-1') as f:
        lineas = f.readlines()

    # La linea 4 (indice 3) trae algo como:
    # "NPTS= 4000, DT= .0200 SEC"
    linea_info = lineas[3]
    npts = int(linea_info.split('NPTS=')[1].split(',')[0].strip())
    dt = float(linea_info.split('DT=')[1].strip().split()[0])

    valores = []
    for linea in lineas[4:]:
        valores.extend(float(v) for v in linea.split())
    valores = np.array(valores[:npts])

    tiempo = np.arange(npts) * dt
    return tiempo, valores, dt


def leer_cismid_redacis(ruta, componente='EW'):
    """Lee un acelerograma en el formato REAL que entrega REDACIS-CISMID
    (verificado contra un archivo descargado del sismo del 23/06/2001,
    estacion MOQ001). Estructura real del archivo:

      - Encabezado de texto libre con metadatos de estacion y sismo
        (codificado en latin-1, no UTF-8 -- trae tildes y enies).
      - Una linea "NUMBER OF SAMPLES : <n>" con el numero de puntos.
      - Una linea "5. ACCELERATION DATA" que marca el inicio de la tabla.
      - Le sigue una linea de titulos de columna: "T   EW   NS   UD".
      - Luego <n> filas con 4 columnas: tiempo(s), EW(cm/s2), NS(cm/s2),
        UD/vertical(cm/s2).

    Los archivos UPTAC/TAC2 de 2014 traen ademas "STANDARD FILE OF
    ACCELERATION DATA" al inicio de la cabecera, por eso la tabla se ubica
    por la fila de titulos "T EW NS UD" y no por ese texto.

    componente: 'EW', 'NS' o 'UD' -- que componente horizontal (o
    vertical) devolver. Devuelve (tiempo, aceleracion_en_g, dt).
    """
    with open(ruta, encoding='latin-1') as f:
        lineas = f.readlines()

    # Ubicar la fila de titulos de columna (T EW NS UD); los datos empiezan debajo
    idx_titulos = next((i for i, l in enumerate(lineas)
                        if l.split() and l.split()[0] == 'T'
                        and set(l.split()[1:]) <= {'EW', 'NS', 'UD', 'NE', 'NW', 'Z'}), None)
    if idx_titulos is None:
        raise ValueError(f'{ruta}: no se encontro la fila de titulos "T EW NS UD"')
    titulos = lineas[idx_titulos].split()
    if componente not in titulos:
        raise ValueError(f'{ruta}: no hay columna {componente} (columnas: {titulos})')
    idx_columna = titulos.index(componente)   # columna por nombre, no por posicion
    idx_datos = idx_titulos + 1

    tiempo, acel = [], []
    for linea in lineas[idx_datos:]:
        partes = linea.split()
        if len(partes) < 4:
            continue   # fin de archivo o linea en blanco
        tiempo.append(float(partes[0]))
        acel.append(float(partes[idx_columna]))

    tiempo = np.array(tiempo)
    acel_g = np.array(acel) / 981.0   # cm/s2 -> g
    dt = tiempo[1] - tiempo[0]
    return tiempo, acel_g, dt


def leer_registro_csv(ruta, columna_tiempo=0, columna_acel=1,
                       delimitador=',', tiene_encabezado=True, unidad_acel='g'):
    """Lector generico para un archivo de texto de dos columnas
    (tiempo, aceleracion). Ajusta columna_tiempo/columna_acel y el
    delimitador segun el archivo real que descargues.
    unidad_acel: 'g' o 'cm/s2' -- convierte internamente todo a 'g'."""
    datos = np.genfromtxt(ruta, delimiter=delimitador,
                           skip_header=1 if tiene_encabezado else 0)
    tiempo = datos[:, columna_tiempo]
    acel = datos[:, columna_acel]
    if unidad_acel == 'cm/s2':
        acel = acel / 981.0   # cm/s2 -> g
    dt = tiempo[1] - tiempo[0]
    return tiempo, acel, dt


# =============================================================================
# 2. ESPECTRO DE RESPUESTA (metodo de Newmark, aceleracion lineal)
# =============================================================================
def espectro_respuesta(acc_g, dt, periodos, zeta=0.05):
    """Calcula el espectro de pseudo-aceleracion (en g) de un registro,
    para cada periodo en 'periodos', mediante integracion de Newmark
    (beta=1/6, gamma=1/2 -- metodo de la aceleracion lineal, incondicional
    y suficientemente preciso para el paso de tiempo tipico de un
    acelerograma). acc_g: aceleracion del suelo en g. dt: paso de tiempo (s)."""
    acc = acc_g * 981.0   # trabajar en cm/s2 internamente
    # Newmark - aceleracion promedio constante (beta=1/4, gamma=1/2):
    # INCONDICIONALMENTE ESTABLE, a diferencia del metodo de aceleracion
    # lineal (beta=1/6), que se vuelve inestable cuando el paso de tiempo
    # del registro es grande respecto al periodo evaluado.
    beta, gamma = 1 / 4, 1 / 2
    Sa = np.zeros(len(periodos))

    for idx, T in enumerate(periodos):
        if T < 1e-6:
            Sa[idx] = np.max(np.abs(acc)) / 981.0   # T->0: Sa tiende al PGA
            continue

        wn = 2 * np.pi / T
        k = wn**2          # rigidez de un oscilador de masa unitaria
        c = 2 * zeta * wn  # amortiguamiento de un oscilador de masa unitaria

        a0 = 1 / (beta * dt**2) + gamma * c / (beta * dt)
        a1 = 1 / (beta * dt) + (gamma / beta - 1) * c
        a2 = (1 / (2 * beta) - 1) + dt * (gamma / (2 * beta) - 1) * c
        k_eff = k + a0

        u, v, a = 0.0, 0.0, -acc[0]   # condiciones iniciales
        u_max = 0.0
        for i in range(1, len(acc)):
            p_eff = -acc[i] + a0 * u + a1 * v + a2 * a
            u_nuevo = p_eff / k_eff
            v_nuevo = gamma / (beta * dt) * (u_nuevo - u) + (1 - gamma / beta) * v + \
                dt * (1 - gamma / (2 * beta)) * a
            a_nuevo = (u_nuevo - u) / (beta * dt**2) - v / (beta * dt) - \
                (1 / (2 * beta) - 1) * a
            u, v, a = u_nuevo, v_nuevo, a_nuevo
            u_max = max(u_max, abs(u))

        Sa[idx] = (wn**2 * u_max) / 981.0   # pseudo-aceleracion, de vuelta a g

    return Sa


# =============================================================================
# 3. ESPECTRO DE DISENO E.030 Y FACTOR DE ESCALA
# =============================================================================
def espectro_E030(periodos, Z, U, S, Tp, TL):
    """Espectro elastico de pseudo-aceleracion de la E.030 (en g),
    evaluado en cada periodo del arreglo 'periodos'."""
    Sa = np.zeros(len(periodos))
    for idx, T in enumerate(periodos):
        if T < Tp:
            C = 2.5
        elif T < TL:
            C = 2.5 * (Tp / T)
        else:
            C = 2.5 * (Tp * TL) / (T**2)
        Sa[idx] = Z * U * C * S
    return Sa


def factor_escala_registro(acc_g, dt, T1, Z, U, S, Tp, TL,
                            factor_rango=(0.2, 1.5), n_periodos=40):
    """Calcula el factor de escala tal que la media del espectro del
    registro, en el rango [factor_rango[0]*T1, factor_rango[1]*T1], iguale
    la media del espectro de diseno E.030 en ese mismo rango (criterio de
    escalado habitual de la norma)."""
    periodos = np.linspace(factor_rango[0] * T1, factor_rango[1] * T1, n_periodos)
    Sa_registro = espectro_respuesta(acc_g, dt, periodos)
    Sa_diseno = espectro_E030(periodos, Z, U, S, Tp, TL)

    media_registro = np.mean(Sa_registro)
    media_diseno = np.mean(Sa_diseno)
    factor = media_diseno / media_registro
    return factor, periodos, Sa_registro, Sa_diseno


# =============================================================================
# PRUEBA 1: con el registro REAL de REDACIS-CISMID, si esta en la misma
# carpeta (sismo del 23/06/2001, estacion MOQ001, Moquegua). Si no existe
# el archivo, esta prueba se salta y solo corre la sintetica de abajo.
# =============================================================================
import os

if __name__ == '__main__':
    archivo_real = 'MOQ001_20010623153300.txt'
    if os.path.exists(archivo_real):
        print("=" * 70)
        print(f"PRUEBA CON REGISTRO REAL: {archivo_real}")
        print("=" * 70)
        t, acc_g, dt = leer_cismid_redacis(archivo_real, componente='EW')
        print(f"Puntos leidos: {len(t)}  |  dt = {dt} s  |  duracion = {t[-1]:.1f} s")
        print(f"PGA (EW): {np.max(np.abs(acc_g)):.4f} g "
              f"({np.max(np.abs(acc_g))*981:.2f} cm/s2 -- debe coincidir con "
              f"el valor 'MAXIMUM ACCELERATION' del archivo)")

        factor, per, Sa_reg, Sa_e030 = factor_escala_registro(
            acc_g, dt, T1=0.396, Z=0.45, U=1.0, S=1.0, Tp=0.40, TL=2.50)
        print(f"Factor de escala (T1=0.396s, Zona 4/S1): {factor:.3f}")
        print()

# =============================================================================
# PRUEBA 2: con un acelerograma SINTETICO (solo para verificar que el codigo
# calcula bien -- NO es un registro real, no usar para el IDA final)
# =============================================================================
    np.random.seed(1)
    dt = 0.02
    n = 1500
    t = np.arange(n) * dt

    # Ruido blanco filtrado (Kanai-Tajimi simplificado) + envolvente, solo
    # para tener "algo" con forma de sismo y probar las funciones
    ruido = np.random.normal(0, 1, n)
    wg, zg = 2 * np.pi / 0.5, 0.6
    acc_filtrada = np.zeros(n)
    for i in range(2, n):
        acc_filtrada[i] = (ruido[i] - (-2 * zg * wg) * acc_filtrada[i - 1] * dt
                            - (wg**2) * acc_filtrada[i - 2] * dt**2) / 1.0
    envolvente = np.minimum(t / 2.0, 1.0) * np.exp(-0.15 * np.maximum(t - 6, 0))
    acc_sintetica_g = 0.15 * acc_filtrada / np.max(np.abs(acc_filtrada)) * envolvente

    print("=" * 70)
    print("PRUEBA: espectro de respuesta de un acelerograma SINTETICO")
    print("=" * 70)
    periodos_prueba = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0])
    Sa = espectro_respuesta(acc_sintetica_g, dt, periodos_prueba)
    for T, S in zip(periodos_prueba, Sa):
        print(f"  T = {T:.2f} s  ->  Sa = {S:.3f} g")
    print(f"  PGA del registro sintetico: {np.max(np.abs(acc_sintetica_g)):.3f} g")

    print()
    print("=" * 70)
    print("PRUEBA: factor de escala contra el espectro E.030 (Zona 4, S1)")
    print("=" * 70)
    T1 = 0.396   # el periodo del edificio de referencia
    factor, per, Sa_reg, Sa_e030 = factor_escala_registro(
        acc_sintetica_g, dt, T1, Z=0.45, U=1.0, S=1.0, Tp=0.40, TL=2.50)
    print(f"  Factor de escala calculado: {factor:.3f}")
    print(f"  Media Sa registro (sin escalar) en [0.2T1,1.5T1]: {np.mean(Sa_reg):.3f} g")
    print(f"  Media Sa diseno E.030 en [0.2T1,1.5T1]:           {np.mean(Sa_e030):.3f} g")
    print(f"  Verificacion: media(Sa_reg)*factor = {np.mean(Sa_reg)*factor:.3f} g "
          f"(debe igualar la media de diseno)")
