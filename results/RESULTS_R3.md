# Evidencia para el 3.er informe (segunda revisión) — puntos 4.1 y 4.2a

Scripts: `identif_proof.py`, `collide_check.py` (identificabilidad conjunta, 4.1);
`boot_calibration.py` + `boot_calibration_D.py` (calibración del bootstrap, 4.2a).
Seeds fijas (20260704 / 7). numpy 2.4.6, scipy 1.17.1.

## 4.1 — Identificabilidad conjunta de (k, phi): mecanismo de prueba + evidencia

Densidad circular de la mantisa: h(x) = c * r^{frac(x-phi)}, r = 10^{-1/k} in (0,1).
Tiene UNA sola discontinuidad, en x=phi, que parte [0,1) en dos ramas contiguas:
rama alta [phi,1) con densidad B*r^x, rama baja [0,phi) con densidad B*r*r^x.

| Ingrediente (identif_proof.py) | Resultado | Estado |
|---|---|---|
| (1) P(d+1)/P(d) invariante en phi para par de la misma rama | spread 8.5e-15 | OK |
| (2) R_d(beta)=I(d+1)/I(d) estrictamente monotona en beta (=> k identificado) | monotona para todo d | OK |
| (3) masa de la celda del salto estrictamente monotona en phi (=> phi identificado) | monotona | OK |
| (4) Jacobiano de (k,phi)->P de rango 2 en todo [1,K]x[0,1) | 2.o sing. val. 1.9e-3 (K=8), 5.0e-5 (K=50) > 0 | OK |

Busqueda global de colisiones (collide_check.py), min ||P-P'|| entre parametros separados:

| K | sep >= | min residual | ubicacion |
|---|---|---|---|
| 8  | 0.05 | 1.16e-4 | k~7.94<->8.00 (dk=0.06), **dphi=0.00** |
| 8  | 0.20 | 7.35e-4 | k~7.74<->8.00 (dk=0.26), **dphi=0.00** |
| 50 | 0.05 | 1.07e-5 | k~49.85<->50.00 (dk=0.15), **dphi=0.00** |
| 50 | 0.20 | 4.22e-5 | k~49.17<->50.00 (dk=0.83), **dphi=0.00** |

**Lectura.** En TODA cuasi-degeneracion dphi=0.00 (la fase nunca se confunde) y k esta en el
tope del rango (k->K), la cola de identificacion debil hacia Benford, ya declarada como
regimen saturado. No hay colision a k moderado. Conclusion: la inyectividad conjunta es
cierta en [1,K]x[0,1) para K finito y demostrable por el argumento de rama unica (k por un
cociente de rama estrictamente monotono; phi por la celda del salto). El unico mal
condicionamiento es en la direccion de k cuando k->K (limite Benford).

## 4.2a — Calibracion del test de equivalencia bootstrap (resolucion de paper)

`boot_calibration.py`: N=2000, grilla (k,phi)=200x120 (fina, para no sesgar psi_corr),
R=400 poblaciones, B=400 (decision) / 600 (cobertura), delta=0.10, alpha=0.05, seed 20260704.
Fit vectorizado por la identidad chi^2(res,P)=sum_d res_d^2/P_d - 1 (res y P suman 1).

| Condicion | Cantidad | Resultado | Esperado |
|---|---|---|---|
| (A) Nivel en psi_corr=delta (frontera H0), 3 formas base | P(consistente) | 0.025 ; 0.013 ; 0.018 | <= alpha=0.05 (conservador) |
| (B) H0 interior (psi_corr=0.149 > delta) | P(consistente) | 0.000 | ~ 0 |
| (C) Potencia (ley genuina, psi_corr~0), 2 casos | P(consistente) | 1.000 ; 1.000 | ~ 1 |
| (D) Frontera k=1 (ley de potencia) | P(consistente) | 1.000 | consistente |
| (D) Saturacion Benford (k=inf) | P(consistente) | 1.000 | consistente |
| (E) Cobertura de u_0.95 (nominal 0.95), 2 casos | P(u_0.95 >= psi_true) | 1.000 ; 1.000 | >= 0.95 |

**Lectura.** El bootstrap de equivalencia tiene nivel <= alpha (0.013-0.025 en la frontera,
conservador cerca de 0), potencia plena bajo H1, cobertura >= nominal, y decide bien en la
frontera k=1 y en la saturacion k=inf, justo los regimenes donde la asintotica chi-cuadrado
no aplica. Esto sustituye la asercion de la Prop. de equivalencia por evidencia. NOTA de
veracidad confirmada: con grilla gruesa (n_phi=16) el nivel salia 0.008 (mas conservador de
lo real); con la grilla fina sube a 0.013-0.025 (honesto), aun <= alpha. Por eso el estudio
usa grilla fina.
