# batch-sim

Corridas numéricas por lotes en GitHub Actions.

- `codigo/`: scripts de cálculo.
- `datos/`: entradas y resultados parciales (un CSV por caso).
- `registros_cifrados/`: entradas cifradas.
- Lanzar: **Actions → Lote → Run workflow**. Si una máquina llega a 5 h 30 min, se guarda lo hecho y se lanza otra ronda sola (máx. 6).
