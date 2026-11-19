"""Cuenta los componentes NGA-Sub pendientes por edificación (sin cargar registros).
Uso: python codigo/plan_nube.py "14 54 58" -> imprime JSON con los ids que tienen pendientes."""
import csv, glob, json, os, sys
from rutas import DATOS, DATOS_IDA_PARTES

arg = sys.argv[1].strip() if len(sys.argv) > 1 else ''
ids = [int(t) for t in arg.split()] if arg else \
      [int(t) for t in open(DATOS / 'ids_16.txt', encoding='utf-8').read().split()]
with open(DATOS / 'inventario_nga_sub.csv', newline='', encoding='utf-8') as f:
    comps = [f"NGA{r['RSN']}_{r['sismo']}:{r['comp']}" for r in csv.DictReader(f) if r.get('estado') == 'OK']
hechos = set()
for p in glob.glob(str(DATOS_IDA_PARTES / 'ida_*.csv')):
    with open(p, newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            hechos.add((int(r['id']), r['registro']))
pend = {i: sum((i, c) not in hechos for c in comps) for i in ids}
total = sum(pend.values())
print(f'{total} componentes pendientes: ' + ', '.join(f'{i}:{n}' for i, n in pend.items()), file=sys.stderr)
print(json.dumps([i for i in ids if pend[i] > 0]))
