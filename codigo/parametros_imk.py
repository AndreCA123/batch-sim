# -*- coding: utf-8 -*-
"""
parametros_imk.py
-----------------
Convierte el refuerzo diseñado (diseno_refuerzo_e060.py) en parámetros del
modelo de plasticidad concentrada Ibarra-Medina-Krawinkler (IMK).

  * My  : Panagiotakos y Fardis (2001), fluencia del acero o no linealidad del
          concreto (la menor curvatura de fluencia).
  * Mc/My, θp, θpc, λ, EIy/EIg : ecuaciones simplificadas de Haselton et al.
          (2008, PEER 2007/03; 2016, ACI Struct. J. 113(6)).
  * Resorte en OpenSeesPy: IMKPeakOriented con la modificación de rigidez
          de Ibarra y Krawinkler (2005) (resorte n veces más rígido).

Unidades de salida: tonf, m.  (Rotaciones en rad.)
"""
import math
from diseno_refuerzo_e060 import BARRAS, ESTRIBO, RECUB, ES, Ec

KGCM2_A_MPA = 0.0980665
# Límites de EIy/EIg de Haselton et al. (2008): [0.2, 0.6]. EI_MIN se puede
# cambiar para el análisis de sensibilidad de rigidez (curvas_capacidad.py --ei-min).
EI_MIN, EI_MAX = 0.2, 0.6


# ----------------------------------------------------------------------------
# Momento de fluencia (Panagiotakos & Fardis 2001) — unidades kgf, cm
# ----------------------------------------------------------------------------
def my_panagiotakos(b, d, dp, rho, rhop, rhov, N, fc, fy):
    ec = Ec(fc)
    n = ES / ec
    dl = dp / d
    # (a) fluencia del acero en tracción
    A = rho + rhop + rhov + N / (b * d * fy)
    B = rho + rhop * dl + 0.5 * rhov * (1 + dl) + N / (b * d * fy)
    ky = math.sqrt(n ** 2 * A ** 2 + 2 * n * B) - n * A
    phi_s = fy / (ES * (1 - ky) * d)
    # (b) no linealidad del concreto en compresión (εc ≈ 1.8 f'c / Ec)
    eps_c = 1.8 * fc / ec
    A2 = rho + rhop + rhov - N / (eps_c * ES * b * d)
    B2 = rho + rhop * dl + 0.5 * rhov * (1 + dl)
    ky2 = math.sqrt(n ** 2 * A2 ** 2 + 2 * n * B2) - n * A2
    phi_c = eps_c / (ky2 * d)
    if phi_c < phi_s:
        phi_y, k = phi_c, ky2
    else:
        phi_y, k = phi_s, ky
    My = b * d ** 3 * phi_y * (
        ec * k ** 2 / 2 * (0.5 * (1 + dl) - k / 3)
        + ES / 2 * ((1 - k) * rho + (k - dl) * rhop + rhov / 6 * (1 - dl)) * (1 - dl))
    return My, phi_y          # kgf·cm, 1/cm


# ----------------------------------------------------------------------------
# Ecuaciones de Haselton et al.
# ----------------------------------------------------------------------------
def haselton(nu, rho_sh, s, d, fc, Ls, H, a_sl=1.0):
    fcm = fc * KGCM2_A_MPA
    theta_p = 0.13 * (1 + 0.55 * a_sl) * 0.13 ** nu * (0.02 + 40 * rho_sh) ** 0.65 \
        * 0.57 ** (0.01 * fcm)
    theta_pc = min(0.76 * 0.031 ** nu * (0.02 + 40 * rho_sh) ** 1.02, 0.10)
    lam = 170.7 * 0.27 ** nu * 0.10 ** (s / d)
    McMy = 1.25 * 0.89 ** nu * 0.91 ** (0.01 * fcm)
    EIratio = min(max(-0.07 + 0.59 * nu + 0.07 * Ls / H, EI_MIN), EI_MAX)
    return {'theta_p': theta_p, 'theta_pc': theta_pc, 'lambda': lam,
            'McMy': McMy, 'EIy_EIg': EIratio}


def _empaquetar(My_kgcm, h, Ig, fc, hs, Ls):
    EIy = hs['EIy_EIg'] * Ec(fc) * Ig                       # kgf·cm2
    My = My_kgcm / 1000 / 100                               # tonf·m
    EIy_t = EIy / 1000 / 1e4                                # tonf·m2
    theta_y = My * Ls / (3 * EIy_t)
    return dict(hs, My=My, EIy=EIy_t, theta_y=theta_y)


def imk_viga(r):
    """Parámetros (+) y (−) de una viga diseñada. ν = 0."""
    fc, fy = r['fc'], r['fy']
    b, h, d, dp = r['b'] * 100, r['h'] * 100, r['d'] * 100, r['dp'] * 100
    e = r['estribos']
    s = e['s_conf'] * 100
    rho_sh = e['ramas'] * BARRAS[ESTRIBO][1] / (s * b)
    Ls = r['ln'] / 2
    hs = haselton(0.0, rho_sh, s, d, fc, Ls, r['h'])
    Ig = b * h ** 3 / 12
    out = {}
    for signo, At, Ac in (('neg', r['sup']['As'], r['inf']['As']),
                          ('pos', r['inf']['As'], r['sup']['As'])):
        My, _ = my_panagiotakos(b, d, dp, At / (b * d), Ac / (b * d), 0.0, 0.0, fc, fy)
        out[signo] = _empaquetar(My, h, Ig, fc, hs, Ls)
    return out


def imk_columna(r):
    """Parámetros por dirección ('x','y') de una columna diseñada (ν con P_grav)."""
    fc, fy = r['fc'], r['fy']
    bx, by = r['bx'] * 100, r['by'] * 100
    a, e = r['armado'], r['estribos']
    db, ab = BARRAS[a['barra']]
    N = r['P_grav'] * 1000
    nu = N / (bx * by * fc)
    s = e['s0'] * 100
    Ls = r['ln'] / 2
    out = {}
    for direc, h, b, n_ext, n_lado, ramas in (
            ('x', bx, by, a['nby'], a['nbx'], e['ramas_x']),
            ('y', by, bx, a['nbx'], a['nby'], e['ramas_y'])):
        dp = RECUB + BARRAS[ESTRIBO][0] + db / 2
        d = h - dp
        rho = n_ext * ab / (b * d)
        rhov = 2 * (n_lado - 2) * ab / (b * d)
        My, _ = my_panagiotakos(b, d, dp, rho, rho, rhov, N, fc, fy)
        rho_sh = ramas * BARRAS[ESTRIBO][1] / (s * b)
        hs = haselton(nu, rho_sh, s, d, fc, Ls, h / 100)
        out[direc] = _empaquetar(My, h, b * h ** 3 / 12, fc, hs, Ls)
        out[direc]['nu'] = nu
    return out


# ----------------------------------------------------------------------------
# Resorte IMK en OpenSeesPy
# ----------------------------------------------------------------------------
def material_imk(ops, tag, p_pos, p_neg, L, n=10.0, res=0.1, theta_u=0.20, D=1.0):
    """
    Crea uniaxialMaterial IMKPeakOriented (OpenSees >= 3.5; reemplaza al
    obsoleto ModIMKPeakOriented) para un resorte rotacional en el extremo de un
    elemento elástico de longitud L (m). Modificación de rigidez de Ibarra y
    Krawinkler (2005): K0 = (n+1)·6EI/L y el EI del elemento elástico se
    multiplica por (n+1)/n (factor que devuelve la función).
    Λ (energía histerética, E_t = Λ·My) = λ·θy  (λ normalizado de Haselton).
    """
    EI = p_pos['EIy']
    K0 = (n + 1) * 6 * EI / L
    Lam = p_pos['lambda'] * p_pos['theta_y']
    ops.uniaxialMaterial('IMKPeakOriented', tag, K0,
                         p_pos['theta_p'], p_pos['theta_pc'], theta_u,
                         p_pos['My'], p_pos['McMy'], res,
                         p_neg['theta_p'], p_neg['theta_pc'], theta_u,
                         p_neg['My'], p_neg['McMy'], res,
                         Lam, Lam, 0.0, Lam,       # Λ_S, Λ_C, Λ_A, Λ_K
                         1.0, 1.0, 1.0, 1.0,       # c_S, c_C, c_A, c_K
                         D, D)
    return (n + 1) / n
