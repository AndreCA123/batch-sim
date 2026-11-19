# IDA NGA-Sub en la nube (tesis de fragilidad sísmica, Arequipa)

Análisis dinámico incremental (OpenSeesPy) de 16 edificaciones de concreto armado con 40 componentes
de registros de interfase NGA-Sub, repartido en máquinas de GitHub Actions (una por edificación).

- `codigo/`: solo los módulos que necesita `correr_ida.py`.
- `datos/ida_partes/`: un CSV por (edificación, componente); el IDA salta los que ya existen.
- `registros_cifrados/`: registros NGA-Sub cifrados (la clave es un secreto del repositorio).
- Lanzar: pestaña **Actions → IDA NGA-Sub en la nube → Run workflow**. Si una máquina llega a 5 h 30 min,
  se guarda lo hecho y se lanza otra ronda sola (máximo 6 rondas).
