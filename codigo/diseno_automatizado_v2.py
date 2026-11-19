"""
=============================================================================
 DISENO AUTOMATIZADO v2 - ESCALAMIENTO DIRECTO POR RIGIDEZ
=============================================================================
Reemplaza la estrategia anterior ("sumar 5 cm a ciegas cada vez") por un
calculo directo de cuanto hay que agrandar (o achicar) las secciones.

LA IDEA (en physical terms):
  La deriva de un portico es, en primera aproximacion, inversamente
  proporcional a la rigidez lateral del sistema. La rigidez de una
  seccion rectangular es proporcional al momento de inercia I = b*h^3/12.
  Si escalamos AMBAS dimensiones (b y h) de una seccion por un mismo
  factor de escala "s" (para mantener columnas proporcionadas y no
  volverlas listones angostos), el area escala como s^2 y el momento de
  inercia escala como s^4 (porque tanto b como h^3 crecen).

  Entonces, si la deriva actual excede el limite por un factor
  "factor = deriva_actual / deriva_objetivo", necesitamos multiplicar la
  rigidez por ese mismo factor, lo que significa escalar las dimensiones
  por s = factor^(1/4).

  Esto reemplaza el tanteo ciego por un salto calculado: en vez de sumar
  5 cm y volver a intentar cinco veces, se calcula de una vez cuanto hay
  que crecer (o encoger, si el diseno resulto sobrado) y se corrige en
  1 a 3 iteraciones en vez de 15.

  Se aplica el MISMO factor a vigas y columnas (aproximacion practica:
  endurece el portico completo de forma proporcional). No es exacto al
  100% porque la relacion deriva-rigidez real no es perfectamente lineal
  (por eso igual se itera y se verifica con el modelo completo cada vez,
  en vez de confiar ciegamente en un solo calculo).
=============================================================================
"""

from modelo_parametrizado import construir_y_analizar_edificio


SECCION_MINIMA = 0.25   # m -- minimo practico de columna/viga (E.060)


def redondear_5cm(valor):
    """Redondea hacia arriba al multiplo de 0.05 m mas cercano."""
    n = round(valor / 0.05)
    if n * 0.05 < valor - 1e-9:
        n += 1
    return round(n * 0.05, 2)


def predimensionar(n_pisos, lx, ly):
    """Punto de partida razonable (no reemplaza un calculo formal de
    cargas axiales; el lazo de ajuste es el que garantiza cumplimiento)."""
    b_col = redondear_5cm(0.25 + 0.03 * n_pisos)
    h_col = redondear_5cm(b_col + 0.05)
    luz_mayor = max(lx, ly)
    h_viga = redondear_5cm(luz_mayor / 11.0)
    b_viga = redondear_5cm(h_viga * 0.55)
    return b_viga, h_viga, b_col, h_col


def escalar_seccion(b, h, factor_rigidez):
    """Escala b y h por el mismo factor de forma tal que el momento de
    inercia (b*h^3) cambie exactamente por 'factor_rigidez'."""
    s = factor_rigidez ** 0.25
    b_new = max(redondear_5cm(b * s), SECCION_MINIMA)
    h_new = max(redondear_5cm(h * s), SECCION_MINIMA)
    return b_new, h_new


def disenar_edificio_automatico(
    n_pisos, h_piso, n_vanos_x, lx, n_vanos_y, ly,
    fc, CV, CM_losa, CM_acabados, CM_tabiques,
    Z, U, S, Tp, TL, R,
    deriva_limite=0.007, deriva_objetivo=0.0065,
    max_iteraciones=8, verbose=True,
):
    """
    deriva_limite: el limite normativo real (E.030) contra el que se
        verifica el cumplimiento.
    deriva_objetivo: un poco por debajo del limite (deja margen de
        seguridad para que el redondeo a 5 cm no empuje el diseno de
        vuelta fuera de norma en el ultimo paso).
    """
    b_viga, h_viga, b_col, h_col = predimensionar(n_pisos, lx, ly)
    historial = []

    for iteracion in range(1, max_iteraciones + 1):
        resultado = construir_y_analizar_edificio(
            n_pisos=n_pisos, h_piso=h_piso,
            n_vanos_x=n_vanos_x, lx=lx, n_vanos_y=n_vanos_y, ly=ly,
            b_viga=b_viga, h_viga=h_viga, b_col=b_col, h_col=h_col,
            fc=fc, CV=CV, CM_losa=CM_losa, CM_acabados=CM_acabados,
            CM_tabiques=CM_tabiques,
            Z=Z, U=U, S=S, Tp=Tp, TL=TL, R=R,
        )
        historial.append({
            'iteracion': iteracion,
            'b_viga': b_viga, 'h_viga': h_viga, 'b_col': b_col, 'h_col': h_col,
            'deriva_max': resultado['deriva_max'],
            'cumple': resultado['deriva_max'] <= deriva_limite,
        })
        if verbose:
            print(f"  Iter {iteracion}: columna {b_col:.2f}x{h_col:.2f} m, "
                  f"viga {b_viga:.2f}x{h_viga:.2f} m -> "
                  f"deriva max = {resultado['deriva_max']:.5f} "
                  f"({'CUMPLE' if resultado['deriva_max'] <= deriva_limite else 'no cumple'})")

        if resultado['deriva_max'] <= deriva_limite:
            return {
                'convergio': True, 'iteraciones': iteracion,
                'b_viga': b_viga, 'h_viga': h_viga, 'b_col': b_col, 'h_col': h_col,
                'resultado_final': resultado, 'historial': historial,
            }

        factor_rigidez = resultado['deriva_max'] / deriva_objetivo
        b_viga, h_viga = escalar_seccion(b_viga, h_viga, factor_rigidez)
        b_col, h_col = escalar_seccion(b_col, h_col, factor_rigidez)

    return {
        'convergio': False, 'iteraciones': max_iteraciones,
        'b_viga': b_viga, 'h_viga': h_viga, 'b_col': b_col, 'h_col': h_col,
        'resultado_final': resultado, 'historial': historial,
    }


# =============================================================================
# PRUEBAS (mismos dos casos que la version anterior, para comparar)
# =============================================================================
if __name__ == '__main__':
    print("=" * 74)
    print("CASO 1: edificio de Llanos & Delgadillo (4 pisos)")
    print("=" * 74)
    caso1 = disenar_edificio_automatico(
        n_pisos=4, h_piso=3.0, n_vanos_x=3, lx=3.3, n_vanos_y=4, ly=2.9,
        fc=280.0, CV=0.200, CM_losa=0.300, CM_acabados=0.100, CM_tabiques=0.150,
        Z=0.45, U=1.00, S=1.00, Tp=0.40, TL=2.50, R=8.0,
    )
    print(f"-> Convergio: {caso1['convergio']} en {caso1['iteraciones']} iteracion(es)")
    print(f"-> Seccion final: columna {caso1['b_col']:.2f}x{caso1['h_col']:.2f} m "
          f"(paper: 0.40x0.45 m), viga {caso1['b_viga']:.2f}x{caso1['h_viga']:.2f} m "
          f"(paper: 0.30x0.40 m)")

    print()
    print("=" * 74)
    print("CASO 2: edificio de 8 pisos que antes NO convergia en 15 iteraciones")
    print("=" * 74)
    caso2 = disenar_edificio_automatico(
        n_pisos=8, h_piso=3.0, n_vanos_x=3, lx=5.0, n_vanos_y=3, ly=5.0,
        fc=210.0, CV=0.200, CM_losa=0.300, CM_acabados=0.100, CM_tabiques=0.150,
        Z=0.45, U=1.00, S=1.05, Tp=0.60, TL=2.00, R=8.0,
    )
    print(f"-> Convergio: {caso2['convergio']} en {caso2['iteraciones']} iteracion(es)")
    print(f"-> Seccion final: columna {caso2['b_col']:.2f}x{caso2['h_col']:.2f} m, "
          f"viga {caso2['b_viga']:.2f}x{caso2['h_viga']:.2f} m")
