# Bootstrap de calibración de la Tabla 4 (2026-07-02)

Script: `bootstrap.py` (B=3000, df=7, φ=0 por normalización, semilla 20260702).
Atiende 3.5 (reportar N·ψ_corr, gl, p-valores), 3.7 (bootstrap completo, incl. k̂=∞) y expone 3.6.

| Dataset | N | k̂ | ψ_corr | N·ψ_corr | p_boot | p_asint (χ²₇) | veredicto efecto | veredicto p |
|---|---|---|---|---|---|---|---|---|
| Sim power law k=2 | 2000 | 2.01 | 0.004 | 8.0 | 0.342 | 0.33 | consistente | consistente |
| Sim uniforme | 2000 | ∞ | 0.410 | 820.0 | 0.000 | ~0 | anómalo | anómalo |
| Election returns | 20272 | ∞ | 0.001 | 20.3 | 0.008 | 0.005 | consistente | **RECHAZA** |
| Financial (SEC) | 3.24M | ∞ | 0.002 | 6473 | 0.000 | 0 | consistente | **RECHAZA** |
| Physical constants | 355 | 7.71 | 0.029 | 10.3 | 0.173 | 0.17 | consistente | consistente |
| World Bank pob. | 265 | 40.0 | 0.016 | 4.2 | 0.788 | 0.75 | consistente | consistente |
| Astronómico HYG | 109400 | 8.57 | 0.012 | 1312.8 | 0.000 | ~0 | consistente | **RECHAZA** |
| Fibonacci | 500 | ∞ | 0.001 | 0.5 | 1.000 | 1.0 | consistente | consistente |
| Pruitt (retractado) | 10269 | ∞ | 0.386 | 3963.8 | 0.000 | 0 | anómalo | anómalo |
| Sísmico (acotado) | 500 | ∞ | 7.526 | 3763 | 0.000 | 0 | anómalo | anómalo |

El nivel del ruido nulo (q95 de N·ψ_corr bajo el nulo ajustado) ronda **14–15** en todas las filas, coherente con χ²₇ (indep. de N). El bootstrap y la asintótica coinciden: la calibración es correcta.

## Hallazgo central (esto decide el framing del paper)

La calibración **rechaza** justamente los auténticos de N grande — Election (N=20k), Astronómico (N=109k) y Financial (N=3.2M) — pese a tener ψ_corr diminuto (0.001–0.012). No es un error del bootstrap: es la propiedad estándar de todo test de bondad de ajuste. **Ningún dato real es EXACTAMENTE una ley de potencia**, así que con N grande cualquier desviación mínima excede el umbral χ²₇ y el test rechaza (potencia → 1).

Consecuencia directa para los reparos:
- **3.6 confirmado:** etiquetar Astronómico/Financial como "consistent" con un corte fijo de ψ_corr es incompatible con la calibración por p-valor. El evaluador tiene razón.
- **Resolución honesta:** ψ_corr debe leerse como **tamaño de efecto** (pequeño = cerca del baseline de ley de potencia), y hay que decir explícitamente que la significancia formal a N grande rechaza casi seguro (propiedad conocida, no defecto). Los veredictos son de tamaño de efecto, no de p-valor. Esto es exactamente lo que pide el evaluador (ψ leído relativo a 1/N).
- **3.9:** la frase "clasifica correctamente todos los auténticos como consistentes" NO se sostiene en sentido de p-valor; SÍ se sostiene en sentido de tamaño de efecto (ψ_corr chico y mucho menor que el clásico). Hay que reformular a lenguaje de tamaño de efecto.

## Lo que esto NO rompe

- La separación por tamaño de efecto sigue siendo real y fuerte: auténticos de ley de potencia dan ψ_corr ~0.001–0.03; uniforme/Pruitt/sísmico dan ψ_corr ~0.4–7.5 (dos-tres órdenes de magnitud). El corrected index SÍ discrimina.
- El bootstrap maneja las filas k̂=∞ (nulo = Benford) sin problema, que era el reparo 3.7.
