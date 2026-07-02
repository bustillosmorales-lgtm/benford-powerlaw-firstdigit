# Recomputación de primer dígito sobre datos reales (2026-06-27)

Calculado con los scripts en `repro/` sobre los datos bajados en `repro2/`.
ψ = índice clásico (vs Benford); k̂ = exponente mín-chi-cuadrado; ψ_corr = índice corregido.
Umbrales: ψ<0.10 consistente, 0.10–0.30 borderline, ≥0.30 anómalo.

| Dataset (fuente real) | N | ψ real | veredicto real | paper decía | ¿coincide? |
|---|---|---|---|---|---|
| Pruitt/Laskowski 2016 retractado (cols continuas, Dryad 33f0n) | 10269 | 0.374 | anómalo | 0.49 F | sí (marca) |
| Gino moral-virtue (OSF, Likert 1-7) | 66844 | 0.138 | borderline ARTIFACTO | 0.231 F | NO aplica (Likert) |
| Election returns (county votes 2024, Dataverse VOQCHQ) | 20272 | 0.0007 | consistente | 1.84 F | **NO — es limpio** |
| Financial statements (SEC EDGAR 2026Q1, value) | 3.24M | 0.0019 | consistente | 2.23 F | **NO — es limpio** |
| Astronómico (HYG v41 distancias) | 109400 | 0.018 | consistente | 0.02 C | sí |
| World Bank población 2022 | 265 | 0.016 | consistente | 0.05 C | sí (cerca) |
| Constantes físicas (CODATA 2022) | 355 | 0.036 | consistente | 0.01 C | sí (cerca) |
| Sísmico (USGS magnitudes) | 500 | 7.47 | anómalo (acotado) | 3.00 FP | bounded-support |
| Fibonacci | 500 | 0.0003 | consistente | 0.02 C | sí |
| Stapel | -- | sin datos | -- | 6.47 F | no existe |

## Conclusión honesta
- Los casos LIMPIOS reproducen como limpios (astro, WB, constantes, Fibonacci, y además election y financial). El test no da falsas alarmas en datos auténticos.
- Pruitt (datos reales retractados) SÍ marca anómalo: un positivo genuino.
- DOS de los cuatro "fraudes" del paper (election, financial) son en realidad LIMPIOS en los datos públicos reales -> la sensibilidad 100% NO se sostiene.
- Gino es Likert: el primer dígito no aplica (artefacto), igual que el sísmico (soporte acotado).
- Stapel no tiene datos.

La validación forense, reconstruida sobre datos reales, no se sostiene como está. El paper publicable es el metodológico: teoría + demostración del falso positivo + comportamiento reproducible en datos limpios + Pruitt como ejemplo genuino + límites de alcance explícitos (acotado / Likert).
