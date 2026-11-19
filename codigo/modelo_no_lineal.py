# -*- coding: utf-8 -*-
"""
modelo_no_lineal.py
-------------------
Pórtico plano no lineal (OpenSeesPy, 2D) de una edificación diseñada con
diseno_integrado.disenar_edificio_completo(), para pushover e IDA.

Modelo (FEMA P695; Haselton et al. 2008; Zareian y Medina 2010):
  * Pórtico interior en la dirección X (línea j = 1), con el refuerzo diseñado
    elemento por elemento.
  * Vigas y columnas: elemento elástico + resortes rotacionales IMKPeakOriented
    en ambos extremos (parametros_imk.py). El EI del elemento se multiplica por
    (n+1)/n y el resorte tiene K0 = (n+1)·6EI/L (Ibarra y Krawinkler 2005).
  * Resortes de vigas asimétricos: rama (+) = momento positivo (acero inferior)
    en el extremo izquierdo; en el derecho se invierte (convención verificada
    en OpenSeesPy: hogging -> deformación (-) a la izquierda, (+) a la derecha).
  * Sin restricción de diafragma entre nudos del pórtico: la rigidez axial de
    las vigas lo cumple, y así se evita encadenar restricciones (un nudo a la
    vez maestro y esclavo), que el manejador Transformation no resuelve bien
    y que producía deriva lateral espuria por gravedad.
  * Masa sísmica por nivel = W_k/(ny+1): con nudos rígidos en planta y pórticos
    X idénticos, cada pórtico toma la misma fracción de la masa.
  * Gravedad (1.0 CM + 0.25 CV) tributaria del pórtico sobre las vigas y peso
    propio de columnas en los nudos; la columna leaning (rigidez a flexión
    despreciable, transformación P-Delta) recibe la diferencia entre la
    gravedad asociada a la masa del pórtico y la que el pórtico carga
    directamente (si esa diferencia es negativa se toma cero: el pórtico
    interior ya carga más gravedad que su fracción de masa, lo que es
    conservador para P-Delta).
  * Columnas con geomTransf PDelta.
  * Amortiguamiento de Rayleigh 5 % en los modos 1 y 3, solo en los elementos
    elásticos y proporcional a la rigidez inicial (Zareian y Medina 2010).

Verificación (pushover de las 150 edificaciones, 2026-10-09):
  * T1_NL / T1_elástico ≈ 2.2 en toda la población: la rigidez efectiva de
    Haselton (EIy/EIg) queda en su límite inferior 0.2 (vigas con ν = 0 y
    columnas con ν bajo).
  * Ω = 1.5-4.8 (mediana 1.8) y μT = 3.6-12.7 (mediana 8.4).
  * En ~la mitad de los edificios fluyen columnas de entrepisos intermedios a
    derivas grandes: con ΣMnc/ΣMnv = 1.2-1.3 en el diseño, el endurecimiento
    de las vigas (Mc/My ≈ 1.23) y My_col ≈ 0.93 Mn llevan la razón real a ~1.0.

Unidades: tonf, m, s.
"""
import math
import numpy as np
import openseespy.opensees as ops

from parametros_imk import imk_viga, imk_columna, material_imk

G = 9.81
GAMMA_C = 2.4


def _ijk(tag, nx, ny):
    """Inverso de la numeración de nudos de modelo_parametrizado.py."""
    t = tag - 1
    por_nivel = (nx + 1) * (ny + 1)
    k, r = divmod(t, por_nivel)
    j, i = divmod(r, nx + 1)
    return i, j, k


def construir_portico_2d(ed, geo, sitio, linea=1, n_rig=10.0, cv_factor=0.25,
                         verbose=False):
    """
    ed  : salida de disenar_edificio_completo(geo, sitio)
    geo : n_pisos, h_piso, n_vanos_x, lx, n_vanos_y, ly
    Devuelve un diccionario con la información del modelo construido en el
    dominio de OpenSees (nudos de control, masas, cargas, elementos).
    """
    n, hp = geo['n_pisos'], geo['h_piso']
    nx, ny, lx, ly = geo['n_vanos_x'], geo['n_vanos_y'], geo['lx'], geo['ly']
    b_v, h_v, b_c, h_c = ed['b_viga'], ed['h_viga'], ed['b_col'], ed['h_col']
    E = 15000.0 * math.sqrt(sitio['fc']) * 10.0          # tonf/m2

    # --- elementos del pórtico (línea j = linea) ---------------------------
    vigas, cols = {}, {}
    for v in ed['vigas'].values():
        if v['dir'] != 'x':
            continue
        i, j, k = _ijk(v['ni'], nx, ny)
        if j == linea:
            vigas[(i, k)] = v
    for c in ed['columnas'].values():
        i, j, k = _ijk(c['ni'], nx, ny)
        if j == linea:
            cols[(i, k)] = c                 # columna del entrepiso k -> k+1

    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    ops.geomTransf('PDelta', 1)              # columnas
    ops.geomTransf('Linear', 2)              # vigas
    ops.geomTransf('PDelta', 3)              # leaning

    def nudo(i, k):
        return 10000 + 100 * k + i

    for k in range(n + 1):
        for i in range(nx + 1):
            ops.node(nudo(i, k), i * lx, k * hp)
    for i in range(nx + 1):
        ops.fix(nudo(i, 0), 1, 1, 1)

    tag_mat, tag_ele, tag_nodo = 1, 1, 50000
    resortes = []                           # (ele, tipo, (i,k), extremo)
    elasticos = []

    def nodo_extremo(base_node):
        nonlocal tag_nodo
        x, y = ops.nodeCoord(base_node)
        ops.node(tag_nodo, x, y)
        ops.equalDOF(base_node, tag_nodo, 1, 2)
        tag_nodo += 1
        return tag_nodo - 1

    def resorte(nodo_j, nodo_e, p_pos, p_neg, L, tipo, pos, ext):
        nonlocal tag_mat, tag_ele
        fac = material_imk(ops, tag_mat, p_pos, p_neg, L, n=n_rig)
        ops.element('zeroLength', tag_ele, nodo_j, nodo_e, '-mat', tag_mat, '-dir', 6)
        resortes.append((tag_ele, tipo, pos, ext))
        tag_mat += 1
        tag_ele += 1
        return fac

    # --- columnas -----------------------------------------------------------
    for (i, k), c in sorted(cols.items()):
        p = imk_columna(c)['x']                          # flexión en el plano X
        nb, nt = nudo(i, k), nudo(i, k + 1)
        eb, et = nodo_extremo(nb), nodo_extremo(nt)
        fac = resorte(nb, eb, p, p, hp, 'columna', (i, k), 'inf')
        resorte(nt, et, p, p, hp, 'columna', (i, k), 'sup')
        I_mod = p['EIy'] / E * fac                       # peralte b_col en X
        ops.element('elasticBeamColumn', tag_ele, eb, et, b_c * h_c, E, I_mod, 1)
        elasticos.append(tag_ele)
        tag_ele += 1

    # --- vigas --------------------------------------------------------------
    for (i, k), v in sorted(vigas.items()):
        p = imk_viga(v)
        nl, nr = nudo(i, k), nudo(i + 1, k)
        el, er = nodo_extremo(nl), nodo_extremo(nr)
        fac = resorte(nl, el, p['pos'], p['neg'], lx, 'viga', (i, k), 'izq')
        resorte(nr, er, p['neg'], p['pos'], lx, 'viga', (i, k), 'der')
        I_mod = p['neg']['EIy'] / E * fac
        ops.element('elasticBeamColumn', tag_ele, el, er, b_v * h_v, E, I_mod, 2)
        elasticos.append(tag_ele)
        vigas[(i, k)]['_ele'] = tag_ele
        tag_ele += 1

    # --- columna leaning ------------------------------------------------------
    x_lean = (nx + 1) * lx
    lean = {}
    for k in range(n + 1):
        lean[k] = 20000 + k
        ops.node(lean[k], x_lean, k * hp)
    ops.fix(lean[0], 1, 1, 0)
    for k in range(1, n + 1):
        ops.equalDOF(nudo(0, k), lean[k], 1)
    for k in range(n):
        ops.element('elasticBeamColumn', tag_ele, lean[k], lean[k + 1], 10.0, E, 1e-6, 3)
        tag_ele += 1

    # --- masas y gravedad -------------------------------------------------------
    W_niv = ed['modelo']['peso_sismico_por_nivel']
    frac = 1.0 / (ny + 1)
    CMne = sitio['CM_losa'] + sitio['CM_acabados'] + sitio['CM_tabiques']
    w_viga = (CMne + cv_factor * sitio['CV']) * ly + b_v * h_v * GAMMA_C
    P_col = b_c * h_c * GAMMA_C * hp

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    masas, P_lean = [], []
    for k in range(1, n + 1):
        m = W_niv[k - 1] * frac / G
        for i in range(nx + 1):                      # masa repartida en el nivel
            ops.mass(nudo(i, k), m / (nx + 1), 1e-9, 1e-9)
        masas.append(m)
        directo = w_viga * nx * lx + P_col * (nx + 1) * (1.0 if k < n else 0.5)
        pl = max(W_niv[k - 1] * frac - directo, 0.0)
        P_lean.append(pl)
        if pl > 0:
            ops.load(lean[k], 0.0, -pl, 0.0)
        for i in range(nx + 1):
            ops.load(nudo(i, k), 0.0, -P_col, 0.0)
    for v in vigas.values():
        ops.eleLoad('-ele', v['_ele'], '-type', '-beamUniform', -w_viga)

    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.test('NormDispIncr', 1e-8, 50)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.1)
    ops.analysis('Static')
    ok = ops.analyze(10)
    ops.loadConst('-time', 0.0)
    ops.wipeAnalysis()

    return {'nudo_control': nudo(0, n), 'nudos_piso': [nudo(0, k) for k in range(n + 1)],
            'nudos_nivel': [[nudo(i, k) for i in range(nx + 1)] for k in range(n + 1)],
            'masas': masas, 'P_lean': P_lean, 'resortes': resortes,
            'elasticos': elasticos, 'altura': n * hp, 'h_piso': hp, 'n_pisos': n,
            'gravedad_ok': ok == 0, 'W_portico': sum(W_niv) * frac,
            'V_diseno_portico': ed['modelo']['Vx'] * frac}


def modal(mod, n_modos=3):
    ops.wipeAnalysis()                 # el eigen no admite otro sistema activo
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    n_modos = min(n_modos, mod['n_pisos'])
    vals = ops.eigen(n_modos)
    T = [2 * math.pi / math.sqrt(v) for v in vals]
    phi = [ops.nodeEigenvector(nd, 1, 1) for nd in mod['nudos_piso'][1:]]
    if phi[-1] < 0:
        phi = [-x for x in phi]
    return T, vals, phi


def aplicar_rayleigh(mod, xi=0.05, modos=(1, 3)):
    """Rayleigh en los modos dados; β·K_inicial solo en elementos elásticos."""
    _, vals, _ = modal(mod, max(modos))
    m2 = min(modos[1], len(vals))
    w1, w2 = math.sqrt(vals[modos[0] - 1]), math.sqrt(vals[m2 - 1])
    a0 = xi * 2 * w1 * w2 / (w1 + w2)
    a1 = xi * 2 / (w1 + w2)
    ops.region(1, '-ele', *mod['elasticos'], '-rayleigh', 0.0, 0.0, a1, 0.0)
    nds = [nd for niv in mod['nudos_nivel'][1:] for nd in niv]
    ops.region(2, '-node', *nds, '-rayleigh', a0, 0.0, 0.0, 0.0)
    return a0, a1


def pushover(mod, deriva_techo_max=0.06, paso=0.0005, caida=0.2, verbose=False):
    """Pushover con patrón proporcional a masa x primer modo (FEMA P695).
    Termina al alcanzar deriva_techo_max o cuando la resistencia cae a
    (1 - caida)·Vmax después del pico."""
    T, _, phi = modal(mod, 3)
    F = [m * f for m, f in zip(mod['masas'], phi)]
    F = [f / sum(F) for f in F]
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    for nds, f in zip(mod['nudos_nivel'][1:], F):
        for nd in nds:
            ops.load(nd, f / len(nds), 0.0, 0.0)

    H, ctrl = mod['altura'], mod['nudo_control']
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.test('NormDispIncr', 1e-6, 100)
    ops.algorithm('Newton')
    ops.analysis('Static')

    algoritmos = [('Newton',), ('KrylovNewton',), ('NewtonLineSearch', '-type', 'Bisection'),
                  ('ModifiedNewton', '-initial')]
    du = paso * H
    u0 = [ops.nodeDisp(nd, 1) for nd in mod['nudos_piso']]
    datos = [(0.0, 0.0, [0.0] * mod['n_pisos'])]
    Vmax = 0.0
    while datos[-1][0] < deriva_techo_max:
        ok = -1
        for alg in algoritmos:
            ops.algorithm(*alg)
            for sub in (1, 4, 16):
                ops.integrator('DisplacementControl', ctrl, 1, du / sub)
                ok = ops.analyze(sub)
                if ok == 0:
                    break
            if ok == 0:
                break
        if ok != 0:
            if verbose:
                print('  pushover: no converge en deriva', datos[-1][0])
            break
        ops.algorithm('Newton')
        # cortante basal = suma de fuerzas laterales aplicadas = λ·ΣF (ΣF = 1).
        # (nodeReaction no recoge la fuerza que pasa por las restricciones equalDOF)
        V = ops.getLoadFactor(2)
        u = [ops.nodeDisp(nd, 1) - a for nd, a in zip(mod['nudos_piso'], u0)]
        idr = [(u[k + 1] - u[k]) / mod['h_piso'] for k in range(mod['n_pisos'])]
        datos.append((u[-1] / H, V, idr))
        Vmax = max(Vmax, V)
        if V < (1 - caida) * Vmax and Vmax > 0:
            break

    drift = np.array([d[0] for d in datos])
    V = np.array([d[1] for d in datos])
    IDR = np.array([d[2] for d in datos])
    i_max = int(np.argmax(V))
    # rigidez efectiva secante a 0.6 Vmax y desplazamiento de fluencia efectivo
    i60 = int(np.argmax(V >= 0.6 * V[i_max]))
    Ke = V[i60] / drift[i60] if drift[i60] > 0 else np.nan
    dy = V[i_max] / Ke
    post = np.where((np.arange(len(V)) > i_max) & (V <= 0.8 * V[i_max]))[0]
    du_80 = drift[post[0]] if len(post) else drift[-1]
    rot = {}
    for ele, tipo, pos, ext in mod['resortes']:
        rot[(ele, tipo, pos, ext)] = abs(ops.eleResponse(ele, 'deformation')[0])
    return {'T': T, 'drift_techo': drift, 'V': V, 'IDR': IDR, 'Vmax': V[i_max],
            'drift_Vmax': drift[i_max], 'drift_y_eff': dy, 'drift_u': du_80,
            'mu_T': du_80 / dy, 'Omega': V[i_max] / mod['V_diseno_portico'],
            'C_max': V[i_max] / mod['W_portico'], 'rot_final': rot,
            'completo': len(post) > 0 or drift[-1] >= deriva_techo_max - 1e-9}
