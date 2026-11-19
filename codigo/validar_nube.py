"""Recalcula un componente ya hecho en la laptop y lo compara con su parte guardada.
Uso: python codigo/validar_nube.py 54 NGA6002707_Coastal.Ecuador:H1"""
import csv, sys, os
import correr_ida as c
import poblacion_lhs
from parametros_e030_2026 import ESTRATOS_SITIO_2026
poblacion_lhs.ESTRATOS_SITIO = ESTRATOS_SITIO_2026

id_e, nombre = int(sys.argv[1]), sys.argv[2]
c.REGISTROS = [r for r in c.cargar_registros_nga() if r['nombre'] == nombre]
e = next(x for x in poblacion_lhs.generar_poblacion(150, 42) if x['id'] == id_e)
ref_f = os.path.join(c.CARPETA_PARTES, f"ida_{id_e}_{nombre.replace(':', '_')}.csv")
ref = {round(float(r['Sa_T1_g']), 4): r for r in csv.DictReader(open(ref_f, encoding='utf-8'))}
_, _, filas, err = c.procesar_componente(e, nombre)
assert err is None, err
dif = 0.0
for f in filas:
    r = ref.get(round(f['Sa_T1_g'], 4))
    if r is None:
        print('Sa nueva sin referencia', f['Sa_T1_g']); dif = 1e9; continue
    d = abs(f['IDR_max'] - float(r['IDR_max'])) / max(float(r['IDR_max']), 1e-9)
    dif = max(dif, d)
    print(f"Sa={f['Sa_T1_g']:.4f}  IDR nube={f['IDR_max']:.6f}  laptop={float(r['IDR_max']):.6f}  colapso {f['colapso']}/{r['colapso']}")
print(f'{len(filas)} corridas (laptop {len(ref)}); diferencia relativa máx. de IDR = {dif:.2e}')
sys.exit(0 if dif < 1e-3 and len(filas) == len(ref) else 1)
