# -*- coding: utf-8 -*-
"""
diseno_refuerzo_e060.py
-----------------------
Diseño automatizado del refuerzo de vigas y columnas de pórticos de concreto
armado según la Norma E.060, capítulo 21 (texto de referencia: Propuesta de
Norma E.060 Concreto Armado 2019, SENCICO), artículos 21.5 (vigas), 21.6
(columnas) y 21.7 (nudos) para pórticos.

v3 (2026-10-03): límites de confinamiento según 21.5.3.2 y 21.6.4.5/21.6.4.7,
cuantía máxima con εt >= 0.005 (21.5.2.1), Ash con 21-3c para Pu > 0.3f'cAg,
separación de ramas <= 350 mm (21.6.4.6) y diámetro de barra de viga limitado
por la dimensión de la columna (21.7.2.4: c >= 20 db).

Entrada : secciones (de diseno_automatizado_v2.py) + envolventes de fuerzas
          del análisis elástico con las combinaciones E.060.
Salida  : acero longitudinal, estribos y banderas de verificación por elemento,
          listos para parametros_imk.py.

Unidades de ENTRADA/SALIDA: m, tonf, tonf·m.   Internas: cm, kgf.

Simplificaciones (declararlas en la tesis):
  * Vigas: acero superior uniforme (máx. de los momentos negativos en caras) y
    acero inferior uniforme, una sola capa.
  * Columnas: armado simétrico, flexocompresión uniaxial por dirección
    (sin Bresler), barras de un solo diámetro.
  * Columna fuerte–viga débil por nudo y dirección (E.060 21.6.2.2).
"""
import math
import numpy as np

# ----------------------------------------------------------------------------
# Datos generales
# ----------------------------------------------------------------------------
BARRAS = {  # nombre: (diámetro cm, área cm2)
    '3/8': (0.953, 0.71), '1/2': (1.270, 1.29), '5/8': (1.588, 1.99),
    '3/4': (1.905, 2.84), '1': (2.540, 5.10),
}
LONG_VIGA = ['5/8', '3/4', '1']
LONG_COL = ['5/8', '3/4', '1']
ESTRIBO = '3/8'
RECUB = 4.0            # cm, recubrimiento libre
ES = 2.0e6             # kg/cm2
EPS_CU = 0.003
TONF = 1000.0          # kgf
PHI_V = 0.85


def beta1(fc):
    return 0.85 if fc <= 280 else max(0.65, 0.85 - 0.05 * (fc - 280) / 70)


def Ec(fc):
    return 15000.0 * math.sqrt(fc)          # E.060 8.5.2


def redondear_s(s, paso=2.5, s_min=5.0):
    return max(s_min, math.floor(s / paso) * paso)


# ----------------------------------------------------------------------------
# VIGAS
# ----------------------------------------------------------------------------
def as_flexion(Mu, b, d, fc, fy, phi=0.9):
    """Acero requerido (cm2) para Mu (kgf·cm), sección rectangular simple."""
    if Mu <= 0:
        return 0.0
    Rn = Mu / (phi * b * d ** 2)
    disc = 1 - 2 * Rn / (0.85 * fc)
    if disc <= 0:
        return float('inf')
    return 0.85 * fc / fy * (1 - math.sqrt(disc)) * b * d


def mn_viga(As, b, d, fc, fy):
    a = As * fy / (0.85 * fc * b)
    return As * fy * (d - a / 2)            # kgf·cm


def elegir_barras_viga(As_req, b, db_max=None):
    """Una capa, n>=2. Devuelve la opción con menor área que cumple.
    db_max (cm): diámetro máximo permitido por 21.7.2.4 (c >= 20 db)."""
    db_e = BARRAS[ESTRIBO][0]
    ancho = b - 2 * (RECUB + db_e)
    mejor = None
    for nom in LONG_VIGA:
        db, ab = BARRAS[nom]
        if db_max is not None and db > db_max + 1e-9:
            continue
        n = max(2, math.ceil(As_req / ab))
        libre = (ancho - n * db) / (n - 1)
        if libre < max(2.5, db):
            continue
        if mejor is None or n * ab < mejor['As']:
            mejor = {'n': n, 'barra': nom, 'As': n * ab, 'db': db}
    return mejor


def disenar_viga(v, fc=210.0, fy=4200.0):
    """
    v: dict con
      id, ni, nj, dir ('x'|'y'), b, h, ln (luz libre, m),
      wu (tonf/m, 1.25(CM+CV)), Mu_neg (tonf·m, máx. en caras),
      Mu_pos (tonf·m, máx. positivo en caras o centro),
      c_par (m, opcional): dimensión de la columna paralela a la viga (21.7.2.4)
    """
    b, h = v['b'] * 100, v['h'] * 100
    db_e = BARRAS[ESTRIBO][0]
    d = h - RECUB - db_e - 1.27            # asume barra ~1" en la estimación
    dp = h - d
    banderas = []

    As_min = max(0.7 * math.sqrt(fc) / fy * b * d, 0)        # E.060 10.5.2
    rho_t = 0.85 * beta1(fc) * fc / fy * 0.003 / (0.003 + 0.005)   # εt = 0.005
    As_max = min(rho_t, 0.025) * b * d                        # 21.5.2.1

    As_neg = max(as_flexion(v['Mu_neg'] * TONF * 100, b, d, fc, fy), As_min)
    As_pos = max(as_flexion(v['Mu_pos'] * TONF * 100, b, d, fc, fy), As_min)

    db_max = v['c_par'] * 100 / 20 if v.get('c_par') else None      # 21.7.2.4
    sup = elegir_barras_viga(As_neg, b, db_max)
    if sup is None or sup['As'] > As_max:
        return {'id': v['id'], 'ok': False, 'banderas': ['ACERO_SUP_EXCEDE: agrandar sección']}

    # 21.5.2.2: Mn+ en la cara >= 1/2 Mn-
    Mn_neg = mn_viga(sup['As'], b, d, fc, fy)
    As_pos = max(As_pos, as_flexion(0.5 * Mn_neg * 0.9, b, d, fc, fy))
    inf = elegir_barras_viga(As_pos, b, db_max)
    if inf is None or inf['As'] > As_max:
        return {'id': v['id'], 'ok': False, 'banderas': ['ACERO_INF_EXCEDE: agrandar sección']}
    Mn_pos = mn_viga(inf['As'], b, d, fc, fy)

    # --- Cortante por capacidad (21.5.4) -----------------------------------
    ln = v['ln'] * 100
    Mpr_neg = mn_viga(sup['As'], b, d, fc, 1.25 * fy)
    Mpr_pos = mn_viga(inf['As'], b, d, fc, 1.25 * fy)
    V_sis = (Mpr_neg + Mpr_pos) / ln
    Ve = V_sis + v['wu'] * TONF / 100 * ln / 2               # wu: kgf/cm
    Vc = 0.0 if V_sis >= 0.5 * Ve else 0.53 * math.sqrt(fc) * b * d
    Vs = max(Ve / PHI_V - Vc, 0.0)
    if Vs > 2.1 * math.sqrt(fc) * b * d:
        banderas.append('CORTANTE_EXCEDE: agrandar sección')
    Av = 2 * BARRAS[ESTRIBO][1]
    s_v = Av * fy * d / Vs if Vs > 0 else 1e9
    db_min = min(sup['db'], inf['db'])
    s_conf = redondear_s(min(d / 4, 6 * db_min, 15.0, s_v))          # 21.5.3.2
    s_resto = redondear_s(min(d / 2, s_v))                           # 21.5.3.4

    return {
        'id': v['id'], 'tipo': 'viga', 'ok': not banderas, 'banderas': banderas,
        'dir': v['dir'], 'ni': v['ni'], 'nj': v['nj'],
        'b': v['b'], 'h': v['h'], 'd': d / 100, 'dp': dp / 100, 'ln': v['ln'],
        'fc': fc, 'fy': fy,
        'sup': sup, 'inf': inf,
        'Mn_neg': Mn_neg / TONF / 100, 'Mn_pos': Mn_pos / TONF / 100,   # tonf·m
        'Mpr_neg': Mpr_neg / TONF / 100, 'Mpr_pos': Mpr_pos / TONF / 100,
        'Ve': Ve / TONF,
        'estribos': {'barra': ESTRIBO, 'ramas': 2, 's_conf': s_conf / 100, 'primero': 0.05,
                     'L_conf': 2 * h / 100, 's_resto': s_resto / 100},
    }


# ----------------------------------------------------------------------------
# COLUMNAS: diagrama de interacción
# ----------------------------------------------------------------------------
def capas_columna(arm, dim_par, dim_perp, direc):
    """
    Capas de acero para flexión con peralte = dim_par.
    arm: {'barra','nbx','nby'}; nbx = barras en cara paralela a X (largo bx),
         nby = barras en cara paralela a Y.
    direc 'x' -> peralte bx: capas extremas con nby barras, intermedias 2 barras.
    Devuelve lista de (y desde fibra comprimida, área).
    """
    db, ab = BARRAS[arm['barra']]
    dp = RECUB + BARRAS[ESTRIBO][0] + db / 2
    n_ext, n_lado = (arm['nby'], arm['nbx']) if direc == 'x' else (arm['nbx'], arm['nby'])
    ys = np.linspace(dp, dim_par - dp, n_lado)
    capas = []
    for i, y in enumerate(ys):
        n = n_ext if i in (0, n_lado - 1) else 2
        capas.append((y, n * ab))
    return capas, dp


def diagrama_interaccion(bw, h, capas, fc, fy, n=120):
    """Curva nominal (Pn, Mn) en kgf, kgf·cm. Pn>0 compresión."""
    b1 = beta1(fc)
    Ast = sum(a for _, a in capas)
    P0 = 0.85 * fc * (bw * h - Ast) + fy * Ast
    pts = [(-fy * Ast, 0.0)]
    for c in np.geomspace(0.02 * h, 3.0 * h, n):
        a = min(b1 * c, h)
        Cc = 0.85 * fc * a * bw
        P, M = Cc, Cc * (h / 2 - a / 2)
        for y, As in capas:
            eps = EPS_CU * (c - y) / c
            fs = max(-fy, min(fy, ES * eps))
            if y < a:
                fs -= 0.85 * fc
            P += fs * As
            M += fs * As * (h / 2 - y)
        pts.append((P, M))
    pts.append((P0, 0.0))
    pts = np.array(sorted(pts, key=lambda p: p[0]))
    return pts, P0, Ast


def phi_flexocompresion(Pn, fc, Ag):
    """E.060 9.3.2.2 (estribos): 0.70; transición a 0.90 para carga axial baja."""
    if Pn <= 0:
        return 0.9
    return max(0.70, 0.9 - 1.4 * Pn / (fc * Ag))


def M_en_P(pts, P):
    """Momento de la curva para una carga axial P (interpolación)."""
    return float(np.interp(P, pts[:, 0], pts[:, 1], left=0.0, right=0.0))


def verificar_PM(bw, h, capas, fc, fy, combos):
    """combos: lista de (Pu kgf, Mu kgf·cm). Devuelve razón demanda/capacidad máx."""
    pts, P0, Ast = diagrama_interaccion(bw, h, capas, fc, fy)
    Ag = bw * h
    phis = np.array([phi_flexocompresion(p, fc, Ag) for p in pts[:, 0]])
    dpts = np.column_stack([phis * pts[:, 0], phis * pts[:, 1]])
    Pmax = 0.80 * 0.70 * P0                                  # E.060 10.3.6.2
    dc = 0.0
    for Pu, Mu in combos:
        if Pu > Pmax:
            return float('inf')
        Mcap = M_en_P(dpts, Pu)
        dc = max(dc, abs(Mu) / Mcap if Mcap > 0 else float('inf'))
    return dc


def candidatos_armado(bx, by):
    """Armados simétricos ordenados por área (bx, by en cm)."""
    db_e = BARRAS[ESTRIBO][0]
    out = []
    for nom in LONG_COL:
        db, ab = BARRAS[nom]
        for nbx in range(2, 9):
            for nby in range(2, 9):
                ok = True
                for dim, n in ((bx, nbx), (by, nby)):
                    disp = dim - 2 * (RECUB + db_e)
                    libre = (disp - n * db) / (n - 1)
                    ctc = libre + db
                    if libre < max(1.5 * db, 4.0) or ctc > 20.0:  # 7.6.3 y práctica
                        ok = False
                if ok:
                    ntot = 2 * nbx + 2 * nby - 4
                    out.append({'barra': nom, 'nbx': nbx, 'nby': nby,
                                'n': ntot, 'As': ntot * ab})
    return sorted(out, key=lambda a: (a['As'], a['n']))


def evaluar_armado(c, arm, fc, fy):
    """True si el armado resiste todas las combinaciones en ambas direcciones."""
    bx, by = c['bx'] * 100, c['by'] * 100
    rho = arm['As'] / (bx * by)
    if rho < 0.01 or rho > 0.06:                             # 21.6.3.1
        return False, rho
    for direc, dim_par, dim_perp, k in (('x', bx, by, 1), ('y', by, bx, 2)):
        capas, _ = capas_columna(arm, dim_par, dim_perp, direc)
        combos = [(cb[0] * TONF, cb[k] * TONF * 100) for cb in c['combos']]
        if verificar_PM(dim_perp, dim_par, capas, fc, fy, combos) > 1.0:
            return False, rho
    return True, rho


def mn_columna(c, arm, direc, fc, fy, P=None):
    """Mn nominal (tonf·m) en la dirección dada para la P más desfavorable."""
    bx, by = c['bx'] * 100, c['by'] * 100
    dim_par, dim_perp = (bx, by) if direc == 'x' else (by, bx)
    capas, _ = capas_columna(arm, dim_par, dim_perp, direc)
    pts, _, _ = diagrama_interaccion(dim_perp, dim_par, capas, fc, fy)
    Ps = [cb[0] * TONF for cb in c['combos']] if P is None else [P * TONF]
    return min(M_en_P(pts, p) for p in Ps) / TONF / 100


def disenar_columna_flexion(c, fc=210.0, fy=4200.0, desde=0):
    """
    c: dict con id, ni (nudo inferior), nj (nudo superior), bx, by, ln (m),
       P_grav (tonf, CM+0.25CV), combos [(Pu, Mux, Muy), ...] en tonf, tonf·m.
       Mux = momento por acción del pórtico en X (peralte bx).
    """
    cands = candidatos_armado(c['bx'] * 100, c['by'] * 100)
    for i, arm in enumerate(cands[desde:], start=desde):
        ok, rho = evaluar_armado(c, arm, fc, fy)
        if ok:
            return i, arm, rho, cands
    return None, None, None, cands


def estribos_columna(c, arm, fc, fy, Mpr_tope=None):
    """Mpr_tope: {'x': (M_extremo_inf, M_extremo_sup), 'y': (...)} en tonf·m,
    momento máximo que las vigas pueden transmitir a cada extremo de la
    columna (21.6.5.1). None -> cota superior Ve = 2·Mpr_col/ln."""
    bx, by = c['bx'] * 100, c['by'] * 100
    db = BARRAS[arm['barra']][0]
    ae = BARRAS[ESTRIBO][1]
    Ag = bx * by
    Ach = (bx - 2 * RECUB) * (by - 2 * RECUB)
    ln = c['ln'] * 100
    Pu_min = min(cb[0] for cb in c['combos']) * TONF
    res = {'barra': ESTRIBO}
    s_lim = []
    for direc, dim_par, dim_perp, nb_cara in (('x', bx, by, arm['nby']),
                                              ('y', by, bx, arm['nbx'])):
        # Cortante de diseño por capacidad: Ve = 2 Mpr / ln  (21.6.5.1, cota superior)
        capas, dp = capas_columna(arm, dim_par, dim_perp, direc)
        pts, _, _ = diagrama_interaccion(dim_perp, dim_par, capas, fc, 1.25 * fy)
        Mpr = max(M_en_P(pts, cb[0] * TONF) for cb in c['combos'])
        if Mpr_tope and Mpr_tope.get(direc):
            m_inf, m_sup = (min(Mpr, m * TONF * 100) for m in Mpr_tope[direc])
            Ve = (m_inf + m_sup) / ln                                # 21.6.5.1
        else:
            Ve = 2 * Mpr / ln
        # "en ningún caso Vu menor que el cortante del análisis": no se dispone
        # aquí del cortante elástico por elemento; Ve por capacidad lo supera
        # en pórticos diseñados con R = 8.
        d = dim_par - dp
        Vc = 0.0 if Pu_min < 0.05 * fc * Ag else \
            0.53 * math.sqrt(fc) * (1 + Pu_min / (140 * Ag)) * dim_perp * d
        Vs = max(Ve / PHI_V - Vc, 0.0)
        bc = dim_perp - 2 * RECUB                       # núcleo ⟂ a las ramas
        req = bc * fc / fy * max(0.3 * (Ag / Ach - 1), 0.09)   # Ash/s (21-3a, 21-3b)
        Pu_max = max(cb[0] for cb in c['combos']) * TONF
        if Pu_max > 0.3 * fc * Ag:                             # 21-3c
            nl = 2 * arm['nbx'] + 2 * arm['nby'] - 4
            kf = max(fc / 1785.0 + 0.6, 1.0)
            kn = nl / (nl - 2)
            req = max(req, 0.2 * kf * kn * Pu_max / (fy * Ach) * bc)
        ramas = max(2, math.ceil(bc / 35.0) + 1)        # 21.6.4.6: <= 350 mm
        while True:
            Av = ramas * ae
            s_ash = Av / req
            s_v = Av * fy * d / Vs if Vs > 0 else 1e9
            if min(s_ash, s_v) >= 7.5 or ramas >= max(nb_cara, 2):
                break
            ramas += 1
        res['ramas_' + direc] = ramas
        res['Ve_' + direc] = Ve / TONF
        s_lim.append((s_ash, s_v))
    s0 = redondear_s(min(6 * db, 10.0,
                         *[min(a, v) for a, v in s_lim]))           # 21.6.4.5
    Lo = max(ln / 6, max(bx, by), 50.0)                              # 21.6.4.3
    s_resto = redondear_s(min(10 * db, min(bx, by), 25.0,
                              *[v for _, v in s_lim]))               # 21.6.4.7
    res.update({'s0': s0 / 100, 'Lo': Lo / 100, 's_resto': s_resto / 100})
    return res


# ----------------------------------------------------------------------------
# EDIFICIO COMPLETO + COLUMNA FUERTE – VIGA DÉBIL
# ----------------------------------------------------------------------------
def disenar_edificio(vigas, columnas, fc=210.0, fy=4200.0, max_iter=40):
    """Diseña todas las vigas y columnas y aplica ΣMnc >= 1.2 ΣMnv por nudo."""
    res_v = {v['id']: disenar_viga(v, fc, fy) for v in vigas}
    estado = {}
    res_c = {}
    for c in columnas:
        i, arm, rho, cands = disenar_columna_flexion(c, fc, fy)
        estado[c['id']] = {'i': i, 'arm': arm, 'rho': rho, 'cands': cands, 'dat': c}

    # Mapeo nudo -> elementos
    nudo_cols, nudo_vigas = {}, {}
    for c in columnas:
        for n in (c['ni'], c['nj']):
            nudo_cols.setdefault(n, []).append(c['id'])
    for v in vigas:
        for n in (v['ni'], v['nj']):
            nudo_vigas.setdefault(n, []).append(v['id'])

    for _ in range(max_iter):
        cambios = False
        for nudo, vids in nudo_vigas.items():
            cids = nudo_cols.get(nudo, [])
            if len(cids) < 2:          # nudo de azotea: exento (21.6.2.1)
                continue
            if any(estado[cid]['arm'] is None for cid in cids):
                continue
            for direc in ('x', 'y'):
                if any(estado[cid]['arm'] is None for cid in cids):
                    break                # ya no hay armado posible: se marca abajo
                vs = [res_v[vid] for vid in vids
                      if res_v[vid].get('ok') and res_v[vid]['dir'] == direc]
                if not vs:
                    continue
                Mnv = (sum(r['Mn_neg'] for r in vs) if len(vs) == 1 else
                       max(vs[0]['Mn_neg'] + vs[1]['Mn_pos'],
                           vs[0]['Mn_pos'] + vs[1]['Mn_neg']))
                Mnc = sum(mn_columna(estado[cid]['dat'], estado[cid]['arm'],
                                     direc, fc, fy) for cid in cids)
                if Mnc < 1.2 * Mnv:
                    cambios = True
                    for cid in cids:     # subir armado de ambas columnas
                        e = estado[cid]
                        i, arm, rho, _ = disenar_columna_flexion(
                            e['dat'], fc, fy, desde=e['i'] + 1)
                        e.update({'i': i, 'arm': arm, 'rho': rho})
                        if arm is None:
                            e['falla'] = 'COLUMNA_FUERTE_VIGA_DEBIL: no se cumple con rho<=6%'
        if not cambios:
            break

    # Momentos que las vigas pueden transmitir a cada nudo, repartidos entre las
    # columnas que concurren (21.6.5.1)
    def mpr_vigas(nudo, direc):
        vs = [res_v[vid] for vid in nudo_vigas.get(nudo, [])
              if res_v[vid].get('ok') and res_v[vid]['dir'] == direc]
        if not vs:
            return None
        if len(vs) == 1:
            return max(vs[0]['Mpr_neg'], vs[0]['Mpr_pos'])
        return max(vs[0]['Mpr_neg'] + vs[1]['Mpr_pos'], vs[0]['Mpr_pos'] + vs[1]['Mpr_neg'])

    for cid, e in estado.items():
        c = e['dat']
        tope = {}
        for direc in ('x', 'y'):
            ext = []
            for nudo in (c['ni'], c['nj']):
                m = mpr_vigas(nudo, direc)
                n_c = len(nudo_cols.get(nudo, [])) or 1
                ext.append(None if m is None else m / n_c)
            # base empotrada (sin vigas): sin tope en ese extremo
            if all(x is not None for x in ext):
                tope[direc] = tuple(ext)
            elif ext[1] is not None:
                tope[direc] = (1e9, ext[1])
        if e['arm'] is None:
            res_c[cid] = {'id': cid, 'tipo': 'columna', 'ok': False,
                          'banderas': [e.get('falla', 'RHO>6% o no cumple: agrandar sección')]}
            continue
        est = estribos_columna(c, e['arm'], fc, fy, tope)
        res_c[cid] = {
            'id': cid, 'tipo': 'columna', 'ok': True, 'banderas': [],
            'ni': c['ni'], 'nj': c['nj'], 'bx': c['bx'], 'by': c['by'],
            'ln': c['ln'], 'P_grav': c['P_grav'], 'fc': fc, 'fy': fy,
            'armado': e['arm'], 'rho': e['rho'],
            'Mn_x': mn_columna(c, e['arm'], 'x', fc, fy, c['P_grav']),
            'Mn_y': mn_columna(c, e['arm'], 'y', fc, fy, c['P_grav']),
            'estribos': est,
        }
    return res_v, res_c


def resumen(res_v, res_c):
    print(f"{'Elem':<6}{'Sección':<10}{'Longitudinal':<26}{'Estribos Ø3/8':<28}{'Estado'}")
    for r in res_v.values():
        if not r['ok'] and 'sup' not in r:
            print(f"{r['id']:<6}{'':<10}{'':<26}{'':<28}{r['banderas']}")
            continue
        sec = f"{r['b']*100:.0f}x{r['h']*100:.0f}"
        lon = f"sup {r['sup']['n']}Ø{r['sup']['barra']} / inf {r['inf']['n']}Ø{r['inf']['barra']}"
        e = r['estribos']
        est = f"1@5, @{e['s_conf']*100:.1f} en {e['L_conf']:.2f}m, r@{e['s_resto']*100:.0f}"
        print(f"{r['id']:<6}{sec:<10}{lon:<26}{est:<28}{'OK' if r['ok'] else r['banderas']}")
    for r in res_c.values():
        if not r['ok']:
            print(f"{r['id']:<6}{'':<10}{'':<26}{'':<28}{r['banderas']}")
            continue
        a, e = r['armado'], r['estribos']
        sec = f"{r['bx']*100:.0f}x{r['by']*100:.0f}"
        lon = f"{a['n']}Ø{a['barra']} (ρ={r['rho']*100:.2f}%)"
        est = f"@{e['s0']*100:.1f} en {e['Lo']:.2f}m, r@{e['s_resto']*100:.0f} ({e['ramas_x']}x{e['ramas_y']} ramas)"
        print(f"{r['id']:<6}{sec:<10}{lon:<26}{est:<28}OK")
