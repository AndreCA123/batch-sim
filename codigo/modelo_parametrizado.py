"""
=============================================================================
 MODELO 3D PARAMETRIZADO - FUNCION REUTILIZABLE
=============================================================================
Version refactorizada de modelo_3d_diafragma_rigido_VALIDADO.py: la misma
logica (grilla 3D, diafragma rigido, masas, modal, sismo estatico, derivas)
pero convertida en una funcion que recibe la geometria, las secciones y el
sitio como parametros, en vez de tenerlos fijos en el codigo.

Esto es lo que permite, en el siguiente paso, llamarla una vez por cada
una de las 150 edificaciones de tu poblacion parametrica, cambiando solo
los argumentos de entrada.

CORRECCIONES 2026-10-02:
  1. Inercia de vigas (ver comentario en 'secciones').
  2. Peso propio de vigas y columnas en el peso sismico y en la gravedad
     (parametro incluir_peso_propio=True; con False se reproduce la version
     anterior solo con la correccion 1).
CORRECCION 2026-10-03:
  3. Momento torsor accidental (E.030-2026 art. 37) y deriva medida en las
     esquinas de la planta (parametro excentricidad_accidental=True). No
     cambia periodos ni cortante; solo las derivas.
  Para reproducir el cortante del paper (validacion): incluir_peso_propio=False,
  excentricidad_accidental=False.

Verificacion de regresion: al final del archivo se llama a la funcion con
los datos del edificio de Llanos & Delgadillo:
  version original (con error de inercia) : T = 0.3960/0.3564/0.2853 s, V = 38.76 tonf
  solo correccion 1                        : T = 0.3360/0.3019/0.2453 s, V = 38.76 tonf
  correcciones 1 + 2                       : T = 0.4364/0.3921/0.3186 s, V = 61.46 tonf
  correcciones 1 + 2 + 3 (vigente)         : igual T y V; deriva max = 0.00680
  (paper: T = 0.371/0.348/0.293 s)
=============================================================================
"""

import numpy as np
import openseespy.opensees as ops


def torsion_constant(lado1, lado2):
    """Constante de torsion de Saint-Venant, seccion rectangular (Roark)."""
    a, b = max(lado1, lado2), min(lado1, lado2)
    return a * b**3 * (1 / 3 - 0.21 * (b / a) * (1 - (b**4) / (12 * a**4)))


def C_factor(T, Tp, TL):
    if T < Tp:
        return 2.5
    elif T < TL:
        return 2.5 * (Tp / T)
    return 2.5 * (Tp * TL) / (T**2)


def construir_y_analizar_edificio(
    # --- geometria ---
    n_pisos, h_piso, n_vanos_x, lx, n_vanos_y, ly,
    # --- secciones (m) ---
    b_viga, h_viga, b_col, h_col,
    # --- material ---
    fc,
    # --- cargas de gravedad (tonf/m2) ---
    CV, CM_losa, CM_acabados, CM_tabiques,
    # --- sitio sismico (E.030) ---
    Z, U, S, Tp, TL, R,
    # --- opciones ---
    nu=0.20, gamma_c=2.4, g=9.81, factor_reduccion_CV=0.25,
    incluir_peso_propio=True, excentricidad_accidental=True,
):
    """
    Construye el modelo 3D con diafragma rigido y devuelve un diccionario
    con los resultados clave para verificar cumplimiento de la E.030:
    periodos, cortante basal y derivas de entrepiso por direccion.

    Todas las unidades de entrada: metros, tonf, tonf/m2, segundos.
    """
    # --- geometria derivada ---
    x_coords = [i * lx for i in range(n_vanos_x + 1)]
    y_coords = [i * ly for i in range(n_vanos_y + 1)]
    z_coords = [i * h_piso for i in range(n_pisos + 1)]
    Lx_total, Ly_total = x_coords[-1], y_coords[-1]
    xc, yc = Lx_total / 2.0, Ly_total / 2.0
    area_planta = Lx_total * Ly_total

    # --- material ---
    E = 15000.0 * np.sqrt(fc) * 10.0   # tonf/m2
    G = E / (2 * (1 + nu))

    # --- secciones ---
    A_v = b_viga * h_viga
    # CORRECCION (2026-10-02): con geomTransf vecxz=(0,0,1) el eje local z de
    # las vigas es VERTICAL, por lo que la flexion vertical (gravedad y accion
    # de portico) ocurre alrededor del eje local y. Por eso Iy debe ser la
    # inercia fuerte b*h^3/12. Antes estaba invertido (la viga trabajaba con
    # su inercia debil h*b^3/12, como si fuera "echada").
    Iy_v = b_viga * h_viga**3 / 12.0   # flexion vertical (fuerte)
    Iz_v = h_viga * b_viga**3 / 12.0   # flexion horizontal (debil)
    J_v = torsion_constant(b_viga, h_viga)

    A_c = b_col * h_col
    Iz_c = b_col * h_col**3 / 12.0
    Iy_c = h_col * b_col**3 / 12.0
    J_c = torsion_constant(b_col, h_col)

    CM_no_estructural = CM_losa + CM_acabados + CM_tabiques

    # =====================================================================
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)

    node_tags = {}
    tag = 1
    for k, z in enumerate(z_coords):
        for j, y in enumerate(y_coords):
            for i, x in enumerate(x_coords):
                ops.node(tag, x, y, z)
                node_tags[(i, j, k)] = tag
                tag += 1

    for i in range(n_vanos_x + 1):
        for j in range(n_vanos_y + 1):
            ops.fix(node_tags[(i, j, 0)], 1, 1, 1, 1, 1, 1)

    transf_col, transf_viga_x, transf_viga_y = 1, 2, 3
    ops.geomTransf('Linear', transf_col, 1, 0, 0)
    ops.geomTransf('Linear', transf_viga_x, 0, 0, 1)
    ops.geomTransf('Linear', transf_viga_y, 0, 0, 1)

    ele_tag = 1
    for i in range(n_vanos_x + 1):
        for j in range(n_vanos_y + 1):
            for k in range(n_pisos):
                ni, nj = node_tags[(i, j, k)], node_tags[(i, j, k + 1)]
                ops.element('elasticBeamColumn', ele_tag, ni, nj,
                            A_c, E, G, J_c, Iy_c, Iz_c, transf_col)
                ele_tag += 1

    vigas_x = []
    for k in range(1, n_pisos + 1):
        for j in range(n_vanos_y + 1):
            for i in range(n_vanos_x):
                ni, nj = node_tags[(i, j, k)], node_tags[(i + 1, j, k)]
                ops.element('elasticBeamColumn', ele_tag, ni, nj,
                            A_v, E, G, J_v, Iy_v, Iz_v, transf_viga_x)
                vigas_x.append(ele_tag)
                ele_tag += 1

    vigas_y = []
    for k in range(1, n_pisos + 1):
        for i in range(n_vanos_x + 1):
            for j in range(n_vanos_y):
                ni, nj = node_tags[(i, j, k)], node_tags[(i, j + 1, k)]
                ops.element('elasticBeamColumn', ele_tag, ni, nj,
                            A_v, E, G, J_v, Iy_v, Iz_v, transf_viga_y)
                vigas_y.append(ele_tag)
                ele_tag += 1

    master_nodes = {}
    for k in range(1, n_pisos + 1):
        mtag = 9000 + k
        ops.node(mtag, xc, yc, z_coords[k])
        ops.fix(mtag, 0, 0, 1, 1, 1, 0)
        master_nodes[k] = mtag
        nodos_del_piso = [node_tags[(i, j, k)]
                           for i in range(n_vanos_x + 1) for j in range(n_vanos_y + 1)]
        ops.rigidDiaphragm(3, mtag, *nodos_del_piso)

    # CORRECCION (2026-10-02): el peso sismico incluye el peso propio de
    # vigas (X e Y) del nivel y de media columna arriba + media abajo.
    n_columnas = (n_vanos_x + 1) * (n_vanos_y + 1)
    long_vigas = n_vanos_x * (n_vanos_y + 1) * lx + (n_vanos_x + 1) * n_vanos_y * ly
    peso_vigas_nivel = A_v * gamma_c * long_vigas
    peso_col_entrepiso = n_columnas * A_c * gamma_c * h_piso

    peso_sismico_por_nivel = []
    for k in range(1, n_pisos + 1):
        peso_CM = CM_no_estructural * area_planta
        if incluir_peso_propio:
            frac_col = 1.0 if k < n_pisos else 0.5
            peso_CM += peso_vigas_nivel + frac_col * peso_col_entrepiso
        peso_CV = CV * area_planta
        peso_sismico = peso_CM + factor_reduccion_CV * peso_CV
        peso_sismico_por_nivel.append(peso_sismico)
        m_trans = peso_sismico / g
        m_rot = m_trans * (Lx_total**2 + Ly_total**2) / 12.0
        ops.mass(master_nodes[k], m_trans, m_trans, 1e-9, 1e-9, 1e-9, m_rot)
    peso_sismico_total = sum(peso_sismico_por_nivel)

    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for k in range(1, n_pisos + 1):
        for j in range(n_vanos_y + 1):
            for i in range(n_vanos_x):
                idx = (k - 1) * (n_vanos_y + 1) * n_vanos_x + j * n_vanos_x + i
                ele = vigas_x[idx]
                trib = ly if 0 < j < n_vanos_y else ly / 2.0
                w_total = (CM_no_estructural + CV) * trib + A_v * gamma_c
                ops.eleLoad('-ele', ele, '-type', '-beamUniform', 0.0, -w_total)
    if incluir_peso_propio:
        # CORRECCION (2026-10-02): peso propio de vigas Y (la losa sigue
        # siendo unidireccional hacia las vigas X) y de columnas (concentrado
        # en el nudo superior de cada tramo de columna).
        for ele in vigas_y:
            ops.eleLoad('-ele', ele, '-type', '-beamUniform', 0.0, -A_v * gamma_c)
        for k in range(1, n_pisos + 1):
            for i in range(n_vanos_x + 1):
                for j in range(n_vanos_y + 1):
                    ops.load(node_tags[(i, j, k)], 0.0, 0.0, -A_c * gamma_c * h_piso,
                             0.0, 0.0, 0.0)

    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')
    ok_gravedad = ops.analyze(1)
    ops.loadConst('-time', 0.0)

    eigenvalues = ops.eigen(3)
    periodos = [2 * np.pi / np.sqrt(val) for val in eigenvalues]

    def analisis_estatico(direccion, T_dir, pattern_tag):
        C = max(C_factor(T_dir, Tp, TL), 0.11 * R / (Z * U * S))
        ZUCS_R = (Z * U * C * S) / R
        V_basal = ZUCS_R * peso_sismico_total

        k_exp = 1.0 if T_dir <= 0.5 else (0.75 + 0.5 * T_dir if T_dir <= 2.5 else 2.0)
        alturas = [k * h_piso for k in range(1, n_pisos + 1)]
        denom = sum(p * (h**k_exp) for p, h in zip(peso_sismico_por_nivel, alturas))
        fuerzas = [V_basal * (p * (h**k_exp)) / denom
                   for p, h in zip(peso_sismico_por_nivel, alturas)]

        # CORRECCION (2026-10-03): momento torsor accidental Mt = F·e, con
        # e = 0.05 x dimension perpendicular (E.030-2026 art. 37), y deriva
        # evaluada en las cuatro esquinas (maximo desplazamiento relativo).
        L_perp = Ly_total if direccion == 'X' else Lx_total
        e_acc = 0.05 * L_perp if excentricidad_accidental else 0.0

        ops.timeSeries('Linear', pattern_tag)
        ops.pattern('Plain', pattern_tag, pattern_tag)
        for k in range(1, n_pisos + 1):
            fx = fuerzas[k - 1] if direccion == 'X' else 0.0
            fy = fuerzas[k - 1] if direccion == 'Y' else 0.0
            ops.load(master_nodes[k], fx, fy, 0.0, 0.0, 0.0, fuerzas[k - 1] * e_acc)

        ops.integrator('LoadControl', 1.0)
        ops.analyze(1)

        dof = 1 if direccion == 'X' else 2
        factor_inelastico = 0.75 * R
        esquinas = [(0, 0), (n_vanos_x, 0), (0, n_vanos_y), (n_vanos_x, n_vanos_y)]
        derivas = []
        for k in range(1, n_pisos + 1):
            d_k = 0.0
            for (i, j) in esquinas:
                u_sup = ops.nodeDisp(node_tags[(i, j, k)], dof)
                u_inf = ops.nodeDisp(node_tags[(i, j, k - 1)], dof)
                d_k = max(d_k, abs(u_sup - u_inf))
            derivas.append(factor_inelastico * d_k / h_piso)

        ops.remove('loadPattern', pattern_tag)
        ops.loadConst('-time', 0.0)
        return fuerzas, V_basal, derivas

    T_x, T_y = periodos[0], periodos[1]
    fuerzas_x, Vx, derivas_x = analisis_estatico('X', T_x, 2)
    fuerzas_y, Vy, derivas_y = analisis_estatico('Y', T_y, 3)

    return {
        'ok_gravedad': ok_gravedad == 0,
        'periodos': periodos,
        'peso_sismico_total': peso_sismico_total,
        'peso_sismico_por_nivel': peso_sismico_por_nivel,
        'Vx': Vx, 'Vy': Vy,
        'derivas_x': derivas_x, 'derivas_y': derivas_y,
        'deriva_max': max(max(derivas_x), max(derivas_y)),
        'cumple_deriva': max(max(derivas_x), max(derivas_y)) <= 0.007,
    }


# =============================================================================
# PRUEBA DE REGRESION: mismos datos del edificio de Llanos & Delgadillo.
# Debe reproducir T = [0.4364, 0.3921, 0.3186] s y V = 61.46 tonf (version vigente).
# =============================================================================
if __name__ == '__main__':
    resultado = construir_y_analizar_edificio(
        n_pisos=4, h_piso=3.0, n_vanos_x=3, lx=3.3, n_vanos_y=4, ly=2.9,
        b_viga=0.30, h_viga=0.40, b_col=0.40, h_col=0.45,
        fc=280.0,
        CV=0.200, CM_losa=0.300, CM_acabados=0.100, CM_tabiques=0.150,
        Z=0.45, U=1.00, S=1.00, Tp=0.40, TL=2.50, R=8.0,
    )

    print("=" * 70)
    print("PRUEBA DE REGRESION (debe coincidir con el modelo validado)")
    print("=" * 70)
    print(f"Periodos     : {[round(t, 4) for t in resultado['periodos']]}")
    print(f"Esperado     : [0.4364, 0.3921, 0.3186]")
    print(f"Vx = Vy      : {resultado['Vx']:.2f} tonf  (esperado: 61.46 tonf)")
    print(f"Deriva max   : {resultado['deriva_max']:.5f}  (esperado: 0.00680)")
    print(f"Cumple E.030 : {resultado['cumple_deriva']}")

    coincide = (abs(resultado['periodos'][0] - 0.4364) < 0.001 and
                abs(resultado['Vx'] - 61.46) < 0.01)
    print("-" * 70)
    print("RESULTADO DE LA PRUEBA:", "OK - coincide con el modelo validado" if coincide
          else "ALERTA - revisar, no coincide")
    print("=" * 70)
