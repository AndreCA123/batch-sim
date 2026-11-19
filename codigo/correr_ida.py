# -*- coding: utf-8 -*-
"""
correr_ida.py
-------------
Análisis Dinámico Incremental (IDA; Vamvatsikos y Cornell 2002) de las 150
edificaciones de la población v3 con el pórtico no lineal IMK de
modelo_no_lineal.py. Retomable: salta los pares (edificio, registro) ya hechos.

Flujo por edificación
  1. Rehace el diseño (diseno_integrado) y construye el pórtico 2D IMK.
  2. T1 = primer periodo del modelo no lineal (T1_NL). La medida de intensidad
     es Sa(T1_NL, 5 %) del registro escalado.
  3. Para cada componente horizontal de cada registro "OK" de
     datos/inventario_registros.csv:
       - procesa el registro (línea base, filtro Butterworth 0.1-25 Hz,
         recorte por intensidad de Arias 0.1 %-99.9 %),
       - corre el IDA con el algoritmo hunt & fill: sube Sa con pasos
         crecientes hasta el colapso, refina el umbral de colapso por
         bisección y rellena los huecos más grandes de la curva.

Criterio de colapso: deriva de entrepiso máxima >= 10 % o falta de
convergencia numérica (FEMA P695 / Vamvatsikos y Cornell 2002).

Amortiguamiento: Rayleigh 5 % en los modos 1 y 3, con el coeficiente βK de
los elementos elásticos multiplicado por (n+1)/n, porque su rigidez ya está
multiplicada por ese factor (Ibarra y Krawinkler 2005; Zareian y Medina 2010).

Salida (en ../datos):
  ida_partes/ida_{id}_{registro}.csv   un archivo por par (edificio, componente),
                                       escrito apenas termina; sirve para retomar
  ida_150_v3.csv                       unión de todas las partes (paso final),
                                       una fila por análisis dinámico:
      id, n_pisos, T1_NL, registro, sismo_num, Mw, mecanismo, Sa_T1_g,
      factor_escala, IDR_max, piso_IDR_max, deriva_techo_max, colapso, motivo,
      segundos
  ida_errores.csv                      componentes que fallaron (la corrida sigue)

Cada tarea en paralelo es un par (edificio, componente); cada proceso guarda en
caché el diseño y el T1 de los edificios que ya preparó.

Uso en Colab o local:
    # cd codigo
    !python correr_ida.py --prueba              # 1 edificio, 2 registros
    !python correr_ida.py --workers 2           # toda la población
    !python correr_ida.py --ids 1 2 3 --workers 2
    !python correr_ida.py --unir                # solo une las partes en ida_150_v3.csv
"""
import os
import io
import csv
import glob
import time
import argparse
import traceback
import contextlib
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from scipy.signal import butter, filtfilt

import openseespy.opensees as ops

from rutas import DATOS, DATOS_IDA_PARTES, REGISTROS

CARPETA_REG = REGISTROS
F_INVENTARIO = DATOS / 'inventario_registros.csv'
F_INVENTARIO_NGA = DATOS / 'inventario_nga_sub.csv'   # auditar_nga_sub.py
CARPETA_NGA = REGISTROS / 'nga_sub'
F_OUT = DATOS / 'ida_150_v3.csv'
CARPETA_PARTES = DATOS_IDA_PARTES
F_ERRORES = DATOS / 'ida_errores.csv'

G = 9.81
# --- procesamiento de registros ---
F_BAJA, F_ALTA, ORDEN_FILTRO = 0.10, 25.0, 4
COLA_LIBRE = 3.0             # s de vibración libre al final
# --- IDA (hunt & fill) ---
SA_INICIAL = 0.10            # g
PASO_INICIAL = 0.15          # g
INCREMENTO_PASO = 0.05       # g, el paso crece en cada corrida del "hunt"
SA_MAXIMO = 5.0              # g, si no colapsa antes se detiene
N_BISECCION = 3
N_MAX_CORRIDAS = 15
IDR_COLAPSO = 0.10
# --- análisis transitorio ---
DT_ANALISIS = 0.01           # s (el registro se interpola si su dt es menor)
TOL_TEST = 1e-5              # tolerancia del test NormDispIncr
N_RIG = 10.0                 # mismo valor que en construir_portico_2d
# Tope de tiempo de reloj por corrida (s). Una corrida que no converge puede quedar
# en dt muy pequeños por horas (05/11: dos corridas colgaron el IDA toda la noche).
# Se pasa por variable de entorno para que llegue a los procesos hijos (fork y spawn).
T_MAX_CORRIDA = float(os.environ.get('IDA_T_MAX_CORRIDA', '900'))

CAMPOS = ['id', 'n_pisos', 'T1_NL', 'registro', 'sismo_num', 'Mw', 'mecanismo',
          'Sa_T1_g', 'factor_escala', 'IDR_max', 'piso_IDR_max', 'deriva_techo_max',
          'colapso', 'motivo', 'segundos']

REGISTROS = []               # se llena en main(); a los procesos hijos se les pasa con
                             # _iniciar_proceso (initializer), así funciona en Windows
                             # (spawn) y en Linux/Colab (fork)


_EDIFICIOS = {}              # caché por proceso: id -> (construir, T1)


def _iniciar_proceso(registros):
    """Initializer de cada proceso hijo: recibe los registros del padre."""
    global REGISTROS
    REGISTROS = registros


def _silenciar_opensees():
    """Manda los avisos de OpenSees (no convergencia, etc.) a un archivo
    descartable para que la bitácora muestre solo los mensajes del script."""
    ops.logFile(os.devnull, '-noEcho')


# =============================================================================
# 1. REGISTROS
# =============================================================================
def procesar_registro(acc_g, dt):
    """Línea base + taper + Butterworth pasa banda (fase cero) + recorte."""
    a = np.asarray(acc_g, dtype=float)
    t = np.arange(len(a)) * dt
    a = a - np.polyval(np.polyfit(t, a, 1), t)               # media y tendencia
    n_tap = max(int(0.05 * len(a)), 1)                       # taper coseno 5 %
    ventana = np.ones(len(a))
    rampa = 0.5 * (1 - np.cos(np.pi * np.arange(n_tap) / n_tap))
    ventana[:n_tap], ventana[-n_tap:] = rampa, rampa[::-1]
    a = a * ventana
    nyq = 0.5 / dt
    b, c = butter(ORDEN_FILTRO, [F_BAJA / nyq, min(F_ALTA / nyq, 0.99)], btype='band')
    a = filtfilt(b, c, a)
    ia = np.cumsum(a ** 2)
    ia = ia / ia[-1]
    i0 = max(int(np.argmax(ia >= 0.001)) - int(1.0 / dt), 0)
    i1 = int(np.argmax(ia >= 0.999))
    return np.concatenate([a[i0:i1 + 1], np.zeros(int(COLA_LIBRE / dt))])


def sa_periodo(acc_g, dt, T, zeta=0.05):
    """Pseudo-aceleración espectral (g) para un periodo, Newmark promedio."""
    from registros_sismicos import espectro_respuesta
    return float(espectro_respuesta(acc_g, dt, np.array([T]), zeta)[0])


def cargar_registros(f_inventario=F_INVENTARIO, max_registros=None):
    """Lee los registros 'OK' del inventario; cada componente horizontal es un
    registro independiente para el pórtico 2D."""
    from registros_sismicos import leer_cismid_redacis
    regs = []
    with open(f_inventario, newline='', encoding='utf-8') as f:
        filas = [r for r in csv.DictReader(f) if r.get('estado') == 'OK']
    for r in filas:
        ruta = os.path.join(CARPETA_REG, r['archivo'])
        for comp in ('EW', 'NS'):
            _, acc_g, dt = leer_cismid_redacis(ruta, componente=comp)
            regs.append(dict(nombre=f"{r['archivo'][:-4]}:{comp}",
                             acc=procesar_registro(acc_g, dt), dt=dt,
                             sismo_num=r.get('sismo_num', ''), Mw=r.get('Mw', ''),
                             mecanismo=r.get('mecanismo', '')))
    if max_registros:
        regs = regs[:max_registros]
    return regs


def cargar_registros_nga(f_inventario=F_INVENTARIO_NGA, max_registros=None):
    """Registros NGA-Sub (.AT2) con estado 'OK' en inventario_nga_sub.csv; una
    componente horizontal por fila. Nombre: NGA<RSN>_<sismo>:<H1|H2>."""
    from registros_sismicos import leer_peer_at2
    regs = []
    with open(f_inventario, newline='', encoding='utf-8') as f:
        filas = [r for r in csv.DictReader(f) if r.get('estado') == 'OK']
    for r in filas:
        _, acc_g, dt = leer_peer_at2(os.path.join(CARPETA_NGA, r['archivo']))
        regs.append(dict(nombre=f"NGA{r['RSN']}_{r['sismo']}:{r['comp']}",
                         acc=procesar_registro(acc_g, dt), dt=dt,
                         sismo_num=f"NGA_{r['sismo']}", Mw=str(r['Mw']),
                         mecanismo='Interfase'))
    if max_registros:
        regs = regs[:max_registros]
    return regs


# =============================================================================
# 2. MODELO
# =============================================================================
def preparar_edificio(e):
    """Devuelve una función que construye el pórtico de la edificación e
    (desde cero, con gravedad aplicada y amortiguamiento asignado)."""
    import poblacion_lhs
    from poblacion_lhs import CM_LOSA, CM_ACABADOS, CM_TABIQUES, R_SISMICO, U_USO, H_PISO
    from parametros_e030_2026 import ESTRATOS_SITIO_2026
    poblacion_lhs.ESTRATOS_SITIO = ESTRATOS_SITIO_2026
    from diseno_integrado import disenar_edificio_completo
    import modelo_no_lineal as mnl

    sitio = dict(fc=e['fc'], CV=e['CV'], CM_losa=CM_LOSA, CM_acabados=CM_ACABADOS,
                 CM_tabiques=CM_TABIQUES, Z=e['Z'], U=U_USO, S=e['S'],
                 Tp=e['Tp'], TL=e['TL'], R=R_SISMICO)
    geo = dict(n_pisos=e['n_pisos'], h_piso=H_PISO, n_vanos_x=e['n_vanos_x'],
               lx=float(e['lx']), n_vanos_y=e['n_vanos_y'], ly=float(e['ly']))
    with contextlib.redirect_stdout(io.StringIO()):
        ed = disenar_edificio_completo(geo, sitio)

    def construir():
        mod = mnl.construir_portico_2d(ed, geo, sitio, n_rig=N_RIG)
        T, _, _ = mnl.modal(mod, 3)
        aplicar_rayleigh_zm(mod, T)
        mod['T'] = T
        return mod
    return construir


def aplicar_rayleigh_zm(mod, T, xi=0.05, modos=(1, 3), n=N_RIG):
    """Rayleigh en los modos dados; βK inicial en los elementos elásticos
    multiplicado por (n+1)/n (Zareian y Medina 2010)."""
    m2 = min(modos[1], len(T))
    w1, w2 = 2 * np.pi / T[modos[0] - 1], 2 * np.pi / T[m2 - 1]
    a0 = xi * 2 * w1 * w2 / (w1 + w2)
    a1 = xi * 2 / (w1 + w2)
    ops.region(1, '-ele', *mod['elasticos'], '-rayleigh', 0.0, 0.0, a1 * (n + 1) / n, 0.0)
    nds = [nd for niv in mod['nudos_nivel'][1:] for nd in niv]
    ops.region(2, '-node', *nds, '-rayleigh', a0, 0.0, 0.0, 0.0)


# =============================================================================
# 3. HISTORIA DE RESPUESTA
# =============================================================================
ALGORITMOS = [('Newton',), ('KrylovNewton',), ('NewtonLineSearch', '-type', 'Bisection'),
              ('ModifiedNewton', '-initial')]


def correr_historia(mod, acc_g, dt_reg, factor):
    """Historia de respuesta con el registro multiplicado por 'factor'.
    Devuelve IDR máximo, piso, deriva de techo máxima y si colapsó."""
    ops.wipeAnalysis()
    ops.timeSeries('Path', 100, '-dt', dt_reg, '-values', *(acc_g * factor * G))
    ops.pattern('UniformExcitation', 100, 1, '-accel', 100)
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.test('NormDispIncr', TOL_TEST, 50)
    ops.algorithm('Newton')
    ops.integrator('Newmark', 0.5, 0.25)
    ops.analysis('Transient')

    nudos, hp, H = mod['nudos_piso'], mod['h_piso'], mod['altura']
    u0 = np.array([ops.nodeDisp(nd, 1) for nd in nudos])
    t_fin = len(acc_g) * dt_reg
    idr_max = np.zeros(len(nudos) - 1)
    techo_max, motivo = 0.0, ''
    t_reloj = time.time()
    while ops.getTime() < t_fin:
        if T_MAX_CORRIDA > 0 and time.time() - t_reloj > T_MAX_CORRIDA:
            motivo = 'tiempo_excedido'        # no se cuenta como colapso
            break
        ok = ops.analyze(1, DT_ANALISIS)
        if ok != 0:                                   # rescate: dt menor y otros algoritmos
            for alg in ALGORITMOS:
                ops.algorithm(*alg)
                for div in (2, 5, 10):
                    ok = ops.analyze(div, DT_ANALISIS / div)
                    if ok == 0:
                        break
                if ok == 0:
                    break
            ops.algorithm('Newton')
            if ok != 0:
                motivo = 'no_converge'
                break
        u = np.array([ops.nodeDisp(nd, 1) for nd in nudos]) - u0
        idr_max = np.maximum(idr_max, np.abs(np.diff(u)) / hp)
        techo_max = max(techo_max, abs(u[-1]) / H)
        if idr_max.max() >= IDR_COLAPSO:
            motivo = 'IDR>=10%'
            break
    ops.remove('loadPattern', 100)
    ops.remove('timeSeries', 100)
    return dict(IDR_max=float(idr_max.max()), piso_IDR_max=int(np.argmax(idr_max)) + 1,
                deriva_techo_max=float(techo_max),
                colapso=motivo in ('no_converge', 'IDR>=10%'), motivo=motivo)


def correr_historia_con_tope(mod, acc_g, dt_reg, factor):
    """Corre la historia en un proceso hijo (fork) y lo mata si pasa T_MAX_CORRIDA.
    Necesario porque OpenSees puede quedarse colgado DENTRO de ops.analyze (C++),
    donde el control de tiempo en Python nunca llega a ejecutarse (19/11: los 15
    workers quedaron atrapados 6 h). En Windows (sin fork) corre en el mismo proceso."""
    if T_MAX_CORRIDA <= 0 or not hasattr(os, 'fork'):
        return correr_historia(mod, acc_g, dt_reg, factor)
    import pickle
    import select
    import signal
    r_fd, w_fd = os.pipe()
    pid = os.fork()
    if pid == 0:                                   # hijo: corre y devuelve por el pipe
        os.close(r_fd)
        try:
            datos = pickle.dumps(correr_historia(mod, acc_g, dt_reg, factor))
        except BaseException as ex:  # noqa: BLE001
            datos = pickle.dumps({'_error': f'{type(ex).__name__}: {ex}'})
        try:
            with os.fdopen(w_fd, 'wb') as w:
                w.write(datos)
        finally:
            os._exit(0)
    os.close(w_fd)
    trozos, limite = [], time.time() + T_MAX_CORRIDA + 30
    with os.fdopen(r_fd, 'rb') as r:
        while True:
            resto = limite - time.time()
            if resto <= 0:
                break
            listo, _, _ = select.select([r], [], [], min(resto, 5.0))
            if listo:
                t = r.read1(1 << 16) if hasattr(r, 'read1') else os.read(r.fileno(), 1 << 16)
                if not t:
                    break                          # el hijo terminó de escribir
                trozos.append(t)
    if time.time() >= limite:
        try:
            os.kill(pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        os.waitpid(pid, 0)
        return dict(IDR_max=float('nan'), piso_IDR_max=0, deriva_techo_max=float('nan'),
                    colapso=False, motivo='tiempo_excedido')
    os.waitpid(pid, 0)
    res = pickle.loads(b''.join(trozos)) if trozos else {'_error': 'el proceso hijo no devolvió resultado'}
    if '_error' in res:
        raise RuntimeError(res['_error'])
    return res


# =============================================================================
# 4. IDA (hunt & fill)
# =============================================================================
def ida_registro(construir, reg, T1, sf_max=10.0):
    """IDA de un edificio con un registro. Devuelve una lista de corridas."""
    sa_reg = sa_periodo(reg['acc'], reg['dt'], T1)
    corridas = {}
    cortado = []                              # se llena si una corrida pasa el tope de tiempo

    def correr(sa):
        factor = sa / sa_reg
        if sf_max is not None and factor > sf_max:
            print(f"      [sf-max] {reg['nombre']} Sa={sa:.4f}g factor={factor:.2f} > {sf_max} omitido", flush=True)
            return None
        t0 = time.time()
        mod = construir()
        r = correr_historia_con_tope(mod, reg['acc'], reg['dt'], factor)
        r.update(Sa_T1_g=round(sa, 4), factor_escala=round(factor, 3),
                 segundos=round(time.time() - t0, 1))
        if r['motivo'] == 'tiempo_excedido':  # no se guarda (IDR desconocido) y se corta el registro
            print(f"      [tiempo] {reg['nombre']} Sa={sa:.4f}g > {T_MAX_CORRIDA:.0f} s, "
                  f"se detiene este registro", flush=True)
            cortado.append(sa)
            return None
        corridas[sa] = r
        return r['colapso']

    # hunt: pasos crecientes hasta el colapso
    sa, paso, sa_ok, sa_col = SA_INICIAL, PASO_INICIAL, 0.0, None
    while len(corridas) < N_MAX_CORRIDAS and sa <= SA_MAXIMO:
        factor = sa / sa_reg
        if sf_max is not None and factor > sf_max:
            print(f"      [sf-max] {reg['nombre']} Sa={sa:.4f}g factor={factor:.2f} > {sf_max} omitido", flush=True)
            break
        res = correr(sa)
        if res is None:
            break
        if res:
            sa_col = sa
            break
        sa_ok, sa, paso = sa, sa + paso, paso + INCREMENTO_PASO
    # bisección del umbral de colapso
    if sa_col is not None:
        for _ in range(N_BISECCION):
            if len(corridas) >= N_MAX_CORRIDAS:
                break
            medio = 0.5 * (sa_ok + sa_col)
            if sf_max is not None and (medio / sa_reg) > sf_max:
                print(f"      [sf-max] {reg['nombre']} Sa={medio:.4f}g factor={medio / sa_reg:.2f} > {sf_max} omitido", flush=True)
                break
            res = correr(medio)
            if res is None:
                break
            if res:
                sa_col = medio
            else:
                sa_ok = medio
    # fill: rellenar los huecos más grandes por debajo del colapso
    while len(corridas) < N_MAX_CORRIDAS and not cortado:
        niveles = [0.0] + sorted(s for s, r in corridas.items() if not r['colapso'])
        huecos = np.diff(niveles)
        if len(huecos) == 0 or huecos.max() < 0.05:
            break
        k = int(np.argmax(huecos))
        sa_fill = 0.5 * (niveles[k] + niveles[k + 1])
        if sf_max is not None and (sa_fill / sa_reg) > sf_max:
            print(f"      [sf-max] {reg['nombre']} Sa={sa_fill:.4f}g factor={sa_fill / sa_reg:.2f} > {sf_max} omitido", flush=True)
            break
        correr(sa_fill)
    return [corridas[s] for s in sorted(corridas)]


def _ruta_parte(id_e, nombre):
    return os.path.join(CARPETA_PARTES, f"ida_{id_e}_{nombre.replace(':', '_')}.csv")


def _escribir_parte(id_e, nombre, filas):
    """Escribe la parte en un temporal y la renombra: una parte a medias nunca
    queda con el nombre definitivo (y no cuenta como hecha al retomar)."""
    ruta = _ruta_parte(id_e, nombre)
    tmp = ruta + '.tmp'
    with open(tmp, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, CAMPOS)
        w.writeheader()
        w.writerows(filas)
    os.replace(tmp, ruta)


def procesar_componente(e, nombre, sf_max=10.0):
    """Trabajo de un proceso: IDA de una edificación con un componente.
    Devuelve (id, registro, filas, error); error es None si todo salió bien."""
    _silenciar_opensees()
    try:
        if e['id'] not in _EDIFICIOS:
            construir = preparar_edificio(e)
            _EDIFICIOS[e['id']] = (construir, float(construir()['T'][0]))
        construir, T1 = _EDIFICIOS[e['id']]
        reg = next(r for r in REGISTROS if r['nombre'] == nombre)
        filas = [dict(id=e['id'], n_pisos=e['n_pisos'], T1_NL=round(T1, 4),
                      registro=nombre, sismo_num=reg['sismo_num'],
                      Mw=reg['Mw'], mecanismo=reg['mecanismo'],
                      Sa_T1_g=r['Sa_T1_g'], factor_escala=r['factor_escala'],
                      IDR_max=round(r['IDR_max'], 5), piso_IDR_max=r['piso_IDR_max'],
                      deriva_techo_max=round(r['deriva_techo_max'], 5),
                      colapso=r['colapso'], motivo=r['motivo'], segundos=r['segundos'])
                 for r in ida_registro(construir, reg, T1, sf_max=sf_max)]
        if filas:
            _escribir_parte(e['id'], nombre, filas)
        return e['id'], nombre, filas, None
    except Exception as ex:  # noqa: BLE001 -- se anota y la corrida sigue
        lugar = traceback.extract_tb(ex.__traceback__)[-1]
        return e['id'], nombre, [], (f'{type(ex).__name__}: {ex} '
                                     f'({os.path.basename(lugar.filename)}:{lugar.lineno})')


# =============================================================================
# 5. PROGRAMA PRINCIPAL
# =============================================================================
def _leer_partes():
    filas = []
    for ruta in sorted(glob.glob(os.path.join(CARPETA_PARTES, 'ida_*.csv'))):
        try:
            with open(ruta, newline='', encoding='utf-8') as f:
                filas.extend(csv.DictReader(f))
        except Exception:
            pass
    return filas


def _hechos():
    return {(int(r['id']), r['registro']) for r in _leer_partes()}


def unir_partes():
    """Paso final: une todas las partes en ida_150_v3.csv (lo reescribe)."""
    filas = sorted(_leer_partes(), key=lambda r: (int(r['id']), r['registro'],
                                                  float(r['Sa_T1_g'])))
    with open(F_OUT, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, CAMPOS)
        w.writeheader()
        w.writerows(filas)
    n_comp = len({(r['id'], r['registro']) for r in filas})
    print(f'\n{len(filas)} análisis de {n_comp} componentes unidos -> {F_OUT}')


def _anotar_error(id_e, nombre, error):
    nuevo = not os.path.exists(F_ERRORES)
    with open(F_ERRORES, 'a', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, ['fecha', 'id', 'registro', 'error'])
        if nuevo:
            w.writeheader()
        w.writerow(dict(fecha=time.strftime('%Y-%m-%d %H:%M:%S'), id=id_e,
                        registro=nombre, error=error))


def main(ids=None, workers=1, max_registros=None, n_total=150, semilla=42, sf_max=10.0,
         fuente='cismid'):
    global REGISTROS
    from poblacion_lhs import generar_poblacion
    import poblacion_lhs
    from parametros_e030_2026 import ESTRATOS_SITIO_2026
    poblacion_lhs.ESTRATOS_SITIO = ESTRATOS_SITIO_2026

    REGISTROS = (cargar_registros_nga(max_registros=max_registros) if fuente == 'nga_sub'
                 else cargar_registros(max_registros=max_registros))
    nombres = [r['nombre'] for r in REGISTROS]
    print(f'{len(REGISTROS)} registros (componentes) cargados')
    for r in REGISTROS:
        print(f"   {r['nombre']:40s} Mw={r['Mw']:4s} {r['mecanismo']:10s} "
              f"PGA={np.max(np.abs(r['acc'])):.3f} g  dur={len(r['acc']) * r['dt']:.0f} s")

    pob = [e for e in generar_poblacion(n_total, semilla) if ids is None or e['id'] in ids]
    os.makedirs(CARPETA_PARTES, exist_ok=True)
    hechos = _hechos()
    # edificios más altos primero: son los que más tardan
    tareas = [(e, n) for e in sorted(pob, key=lambda e: -e['n_pisos'])
              for n in nombres if (e['id'], n) not in hechos]
    print(f"{len(tareas)} componentes pendientes en "
          f"{len({e['id'] for e, _ in tareas})} edificaciones (sf_max={sf_max})\n")

    t0 = time.time()
    n_hechos, n_errores = 0, 0

    def registrar(id_e, nombre, filas, error):
        nonlocal n_hechos, n_errores
        n_hechos += 1
        if error:
            n_errores += 1
            _anotar_error(id_e, nombre, error)
            print(f'[{id_e:3d}] {nombre:28s} ERROR: {error}', flush=True)
        else:
            col = [fl['Sa_T1_g'] for fl in filas if fl['colapso']]
            nc = [fl['Sa_T1_g'] for fl in filas if fl['motivo'] == 'no_converge']
            print(f"[{id_e:3d}] {nombre:28s} {len(filas):2d} corridas | "
                  f"Sa_col = {min(col) if col else '-':>6} g | "
                  f"no_conv = {nc if nc else '-'} | "
                  f"{sum(fl['segundos'] for fl in filas):5.0f} s", flush=True)
        print(f'      {n_hechos}/{len(tareas)} componentes, {n_errores} errores, '
              f'{time.time() - t0:7.0f} s acumulados', flush=True)

    if workers <= 1:
        for e, n in tareas:
            registrar(*procesar_componente(e, n, sf_max))
    else:
        with ProcessPoolExecutor(max_workers=workers, initializer=_iniciar_proceso,
                                 initargs=(REGISTROS,)) as ex:
            futuros = [ex.submit(procesar_componente, e, n, sf_max) for e, n in tareas]
            for fut in as_completed(futuros):
                try:
                    registrar(*fut.result())
                except Exception as exc:  # noqa: BLE001 -- p. ej. un proceso que murió
                    print(f'ERROR al recibir un resultado: {type(exc).__name__}: {exc}',
                          flush=True)
    unir_partes()


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='IDA de la población v3')
    ap.add_argument('--ids', type=int, nargs='*', help='solo estas edificaciones')
    ap.add_argument('--ids-archivo', type=str, help='archivo con ids separados por espacio o salto de línea')
    ap.add_argument('--sf-max', type=float, default=10.0, help='factor de escala máximo permitido (default: 10)')
    ap.add_argument('--workers', type=int, default=1, help='procesos en paralelo')
    ap.add_argument('--max-registros', type=int, help='usar solo los primeros N registros')
    ap.add_argument('--prueba', action='store_true', help='edificio 1 con 2 registros')
    ap.add_argument('--unir', action='store_true', help='solo unir las partes en ida_150_v3.csv')
    ap.add_argument('--t-max-corrida', type=float, default=None,
                    help='tope de tiempo de reloj por corrida en s (default 900; 0 = sin tope)')
    ap.add_argument('--registros', choices=['cismid', 'nga_sub'], default='cismid',
                    help='conjunto de registros: cismid (inventario_registros.csv) o nga_sub (inventario_nga_sub.csv)')
    a = ap.parse_args()
    if a.t_max_corrida is not None:
        os.environ['IDA_T_MAX_CORRIDA'] = str(a.t_max_corrida)
        T_MAX_CORRIDA = a.t_max_corrida
    if a.unir:
        unir_partes()
    elif a.prueba:
        main(ids=[1], workers=1, max_registros=2, sf_max=a.sf_max)
    else:
        ids = a.ids
        if a.ids_archivo:
            with open(a.ids_archivo, 'r', encoding='utf-8') as f:
                ids = [int(tok) for tok in f.read().split()]
        main(ids=ids, workers=a.workers, max_registros=a.max_registros, sf_max=a.sf_max,
             fuente=a.registros)

