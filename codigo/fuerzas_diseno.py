# -*- coding: utf-8 -*-
"""
fuerzas_diseno.py
-----------------
Puente entre el modelo elástico (modelo_parametrizado.py) y el diseño del
refuerzo (diseno_refuerzo_e060.py).

1. Construye el modelo con construir_y_analizar_edificio() (mismos nudos,
   elementos, masas y periodos que usa el diseño de secciones).
2. Retira la gravedad combinada y corre casos de carga:
     CM, CV, SX(+e), SX(-e), SY(+e), SY(-e), SX0, SY0 (sin excentricidad)
   Sismo estático E.030-2026 con la misma distribución en altura del modelo y
   excentricidad accidental 0.05·L (art. 37).
3. Acción bidireccional E.030-2026 (art. 28.1, 28.3 y 33.3): 100% en una
   dirección (con su excentricidad accidental) + 30% en la perpendicular
   (sin excentricidad), con todas las combinaciones de signo; la envolvente
   equivale a sumar valores absolutos |S100| + 0.3|S30|.
4. Superpone combinaciones E.060 (art. 9.2):
     1.4CM+1.7CV ; 1.25(CM+CV)±S ; 0.9CM±S
5. Arma los datos por elemento y llama a disenar_edificio().

Unidades: m, tonf.
Convención de signos (verificada con OpenSeesPy 3.8, localForce):
  vigas   : M_sagging(i) = My1 ; M_sagging(j) = -My2 ; V hacia arriba = Vz1
  columnas: Pu (compresión +) = N1 ; acción de pórtico X -> My ; pórtico Y -> Mz
"""
import numpy as np
import openseespy.opensees as ops

from modelo_parametrizado import construir_y_analizar_edificio, C_factor
from diseno_refuerzo_e060 import disenar_edificio, resumen

GAMMA_C = 2.4


def _tag(i, j, k, nx, ny):
    return 1 + i + j * (nx + 1) + k * (nx + 1) * (ny + 1)


def _elementos(n_pisos, nx, ny):
    """Reproduce la numeración de elementos de modelo_parametrizado.py."""
    cols, vx, vy = [], [], []
    e = 1
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(n_pisos):
                cols.append(dict(tag=e, ni=_tag(i, j, k, nx, ny),
                                 nj=_tag(i, j, k + 1, nx, ny), i=i, j=j, piso=k + 1))
                e += 1
    for k in range(1, n_pisos + 1):
        for j in range(ny + 1):
            for i in range(nx):
                vx.append(dict(tag=e, ni=_tag(i, j, k, nx, ny),
                               nj=_tag(i + 1, j, k, nx, ny), i=i, j=j, piso=k))
                e += 1
    for k in range(1, n_pisos + 1):
        for i in range(nx + 1):
            for j in range(ny):
                vy.append(dict(tag=e, ni=_tag(i, j, k, nx, ny),
                               nj=_tag(i, j + 1, k, nx, ny), i=i, j=j, piso=k))
                e += 1
    return cols, vx, vy


def _correr_caso(cargar):
    """Aplica un caso de carga, analiza, guarda localForce y deja el modelo limpio."""
    ops.timeSeries('Linear', 100)
    ops.pattern('Plain', 100, 100)
    cargar()
    ops.constraints('Transformation')      # requerido por rigidDiaphragm
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')
    if ops.analyze(1) != 0:
        raise RuntimeError('Falló el análisis de un caso de carga')
    f = {e: np.array(ops.eleResponse(e, 'localForce')) for e in ops.getEleTags()}
    ops.remove('loadPattern', 100)
    ops.remove('timeSeries', 100)
    ops.wipeAnalysis()
    ops.reset()
    return f


def obtener_fuerzas(p):
    """p: dict con los mismos argumentos de construir_y_analizar_edificio."""
    res = construir_y_analizar_edificio(**p)
    n, nx, ny = p['n_pisos'], p['n_vanos_x'], p['n_vanos_y']
    lx, ly, hp = p['lx'], p['ly'], p['h_piso']
    Lx, Ly = nx * lx, ny * ly
    A_v = p['b_viga'] * p['h_viga']
    CMne = p['CM_losa'] + p['CM_acabados'] + p['CM_tabiques']

    # retirar la gravedad combinada que dejó el modelo
    ops.remove('loadPattern', 1)
    ops.wipeAnalysis()
    ops.reset()
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')

    cols, vx, vy = _elementos(n, nx, ny)

    def trib(v):
        return ly if 0 < v['j'] < ny else ly / 2.0

    pp = p.get('incluir_peso_propio', True)
    A_c = p['b_col'] * p['h_col']

    def carga_CM():
        for v in vx:
            ops.eleLoad('-ele', v['tag'], '-type', '-beamUniform', 0.0,
                        -(CMne * trib(v) + A_v * GAMMA_C))
        if pp:   # mismo criterio que modelo_parametrizado.py
            for v in vy:
                ops.eleLoad('-ele', v['tag'], '-type', '-beamUniform', 0.0,
                            -A_v * GAMMA_C)
            for c in cols:
                ops.load(c['nj'], 0.0, 0.0, -A_c * GAMMA_C * hp, 0.0, 0.0, 0.0)

    def carga_CV():
        for v in vx:
            ops.eleLoad('-ele', v['tag'], '-type', '-beamUniform', 0.0,
                        -p['CV'] * trib(v))

    # fuerzas sísmicas estáticas (misma expresión que el modelo)
    W = res['peso_sismico_total']
    Pks = res['peso_sismico_por_nivel']
    alturas = [k * hp for k in range(1, n + 1)]

    def fuerzas(T):
        C = max(C_factor(T, p['Tp'], p['TL']), 0.11 * p['R'] / (p['Z'] * p['U'] * p['S']))
        V = p['Z'] * p['U'] * C * p['S'] / p['R'] * W
        kx = 1.0 if T <= 0.5 else (0.75 + 0.5 * T if T <= 2.5 else 2.0)
        den = sum(P * h ** kx for P, h in zip(Pks, alturas))
        return [V * P * h ** kx / den for P, h in zip(Pks, alturas)]

    FX, FY = fuerzas(res['periodos'][0]), fuerzas(res['periodos'][1])

    def carga_sismo(dirc, signo_e):   # signo_e = 0 -> sin excentricidad
        def f():
            for k in range(1, n + 1):
                m = 9000 + k
                if dirc == 'X':
                    ops.load(m, FX[k - 1], 0, 0, 0, 0, signo_e * 0.05 * Ly * FX[k - 1])
                else:
                    ops.load(m, 0, FY[k - 1], 0, 0, 0, signo_e * 0.05 * Lx * FY[k - 1])
        return f

    casos = {'CM': _correr_caso(carga_CM), 'CV': _correr_caso(carga_CV)}
    for d in ('X', 'Y'):
        for s, nom in ((1, '+e'), (-1, '-e'), (0, '0')):
            casos['S' + d + nom] = _correr_caso(carga_sismo(d, s))

    w = {v['tag']: (CMne * trib(v) + A_v * GAMMA_C, p['CV'] * trib(v)) for v in vx}
    w.update({v['tag']: (A_v * GAMMA_C if pp else 0.0, 0.0) for v in vy})
    return res, casos, cols, vx, vy, w


def combinaciones():
    """Lista de (fCM, fCV, [(caso_sismo, factor), ...]) según E.060 9.2 con la
    acción bidireccional 100% + 30% de la E.030-2026."""
    cb = [(1.4, 1.7, [])]
    for d, o in (('X', 'Y'), ('Y', 'X')):
        for e in ('+e', '-e'):
            for f100 in (1.0, -1.0):
                for f30 in (0.3, -0.3):
                    S = [('S' + d + e, f100), ('S' + o + '0', f30)]
                    cb.append((1.25, 1.25, S))
                    cb.append((0.9, 0.0, S))
    return cb


def preparar_datos_diseno(p):
    res, casos, cols, vx, vy, w = obtener_fuerzas(p)
    cbs = combinaciones()

    def comb(e, c):
        fCM, fCV, sismo = c
        f = fCM * casos['CM'][e] + fCV * casos['CV'][e]
        for caso, fs in sismo:
            f = f + fs * casos[caso][e]
        return f

    vigas = []
    for lista, dirc, L, ancho_col in ((vx, 'x', p['lx'], p['b_col']),
                                      (vy, 'y', p['ly'], p['h_col'])):
        # ancho_col = dimensión de la columna paralela a la viga (21.7.2.4)
        for v in lista:
            e = v['tag']
            Mneg, Mpos = 0.0, 0.0
            xs = np.linspace(0, L, 21)
            for c in cbs:
                f = comb(e, c)
                wu = c[0] * w[e][0] + c[1] * w[e][1]
                Ms = f[4] + f[2] * xs - wu * xs ** 2 / 2       # sagging +
                Mneg = max(Mneg, -Ms.min())
                Mpos = max(Mpos, Ms.max())
            vigas.append(dict(id=f"V{dirc.upper()}{e}", ni=v['ni'], nj=v['nj'], dir=dirc,
                              b=p['b_viga'], h=p['h_viga'], ln=L - ancho_col,
                              wu=1.25 * sum(w[e]), Mu_neg=Mneg, Mu_pos=Mpos,
                              c_par=ancho_col, piso=v['piso']))

    columnas = []
    for c0 in cols:
        e = c0['tag']
        combos = []
        for c in cbs:
            f = comb(e, c)
            combos.append((f[0], max(abs(f[4]), abs(f[10])), max(abs(f[5]), abs(f[11]))))
        Pg = casos['CM'][e][0] + 0.25 * casos['CV'][e][0]
        columnas.append(dict(id=f"C{e}", ni=c0['ni'], nj=c0['nj'],
                             bx=p['b_col'], by=p['h_col'], ln=p['h_piso'] - p['h_viga'],
                             P_grav=Pg, combos=combos, piso=c0['piso']))
    return res, vigas, columnas


def disenar_refuerzo_edificio(p, fy=4200.0, verbose=True):
    """Punto de entrada: parámetros del edificio -> refuerzo E.060."""
    res, vigas, columnas = preparar_datos_diseno(p)
    rv, rc = disenar_edificio(vigas, columnas, fc=p['fc'], fy=fy)
    piso = {d['id']: d['piso'] for d in vigas + columnas}
    for r in list(rv.values()) + list(rc.values()):
        r['piso'] = piso[r['id']]
    ok = all(r['ok'] for r in list(rv.values()) + list(rc.values()))
    if verbose:
        print(f"Periodos: {[round(t, 3) for t in res['periodos']]} s | "
              f"vigas {len(rv)}, columnas {len(rc)} | todo OK: {ok}")
    return {'modelo': res, 'vigas': rv, 'columnas': rc, 'ok': ok}
