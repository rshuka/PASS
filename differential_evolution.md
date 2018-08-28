# Differential Evolution

Differential Evolution wurde 1997 in [1] vorgestellt.

## Initialisierung

Die Implementierung von SA in [2] kann über folgende Parameter konfiguriert werden:

Name|Beschreibung|Standardwert in [2]|Vorgeschlagener Wert in [3]|in PASS benutzter Wert
----|------------|-------------------|---------------------------|----------------------
_NP_|Populationsgröße|-|min(10 * number of parameters, 40)|min(10 * number of parameters, 40)
_gen_|number of generations|1|-|?
_F_|weight coefficient|0.8|0.8; "It has been found recently that selecting F from the interval [0.5, 1.0] randomly for each generation or for each difference vector, a technique called dither, improves convergence behaviour significantly, especially for noisy objective functions."|0.8
_CR_|crossover probability|0.9|0.9 für nicht separierbare Funktionen, 0.2 für separierbare Funktionen|0.9
_variant_|mutation variant|`DE/rand/1/exp`|"We mostly use `DE/rand/1/..`. or `DE/best/1/...`. The crossover method is not so important although Ken Price claims that binomial is never worse than exponential."|?
_ftol_|von pagmo2 definiertes Abbruchkriterium: Die Optimierung wird beendet, falls die Differenz zwischen dem besten und dem schlechtesten objective value kleiner ist als _ftol_.|1e-6|-|0
_xtol_|von pagmo2 definiertes Abbruchkriterium: Die Optimierung wird beendet, falls die Distanz zwischen dem besten und dem schlechtesten Agenten kleiner ist als _xtol_.|1e-6|-|0

Der pagmo2 Quellcode enthält Kommentare zu einigen _variant_-Werten:

* `DE/best/1/exp`
  The oldest DE variant but still not bad. However, we have found several optimization problems where misconvergence occurs.
* `DE/rand/1/exp`
  This is one of my favourite strategies. It works especially well when the "gbIter[]"-schemes experience misconvergence. Try e.g. F=0.7 and CR=0.5 as a first guess.
* `DE/rand-to-best/1/exp`
  This variant seems to be one of the best strategies. Try F=0.85 and CR=1. If you get misconvergence try to increase NP. If this doesn't help you should play around with all three control variables.
* `DE/best/2/exp`
  is another powerful variant worth trying
* `DE/rand/2/exp`
  seems to be a robust optimizer for many functions
* `.../bin`: Essentially same strategies but BINOMIAL CROSSOVER

## Pseudocode

```
let f(x) : ℝ^n → ℝ be the optimisation function, with lb_i ≤ x_i ≤ ub_i for i = 1, ..., n.
let P_old ∈ ℝ^{n×NP} = NP uniformly distributed random agents inside the problem boundaries.
let P_new ∈ ℝ^{n×NP} = NP uninitialized agents
let best_global = the global best agent in P_old
repeat _gen_ times:
  best_iter = best_global
  repeat for i = 1 to NP:
    r ∈ ℕ^5 = 5 distinct random indexes from range [1, NP]
    tmp = P_old_i
    j = uniformly distributed random index from range [1, n]
    if variant ∈ {`DE/best/1/exp`, `DE/rand/1/exp`, `DE/rand-to-best/1/exp`, `DE/best/2/exp`, `DE/rand/2/exp`}:
      repeat n times:
                ⎧best_iter_j + F * (P_old_{r_1,j} - P_old_{r_2,j})                                   , if variant == `DE/best/1/exp`
                ⎪P_old_{r_0,j} + F * (P_old_{r_1,j} - P_old_{r_2,j})                                 , if variant == `DE/rand/1/exp`
        tmp_j = ⎨tmp_j + F * (best_iter_j - tmp_j) + F * (P_old_{r_0,j} - P_old_{r_1,j})             , if variant == `DE/rand-to-best/1/exp`
                ⎪best_iter + (P_old_{r_0,j} + P_old_{r_1,j} - P_old_{r_2,j} - P_old_{r_3,j}) * F     , if variant == `DE/best/2/exp`
                ⎩P_old_{r_4,j} + (P_old_{r_0,j} + P_old_{r_1,j} - P_old_{r_2,j} - P_old_{r_3,j}) * F , if variant == `DE/rand/2/exp`
        j = (j + 1) mod n
        exit loop with probability CR
    else:
      repeat for L = 1 to n:
        with probability CR, or if L == n:
                  ⎧P_old_{r_0,j} + F * (P_old_{r_1,j} - P_old_{r_2,j})                                 , if variant == `DE/rand/1/bin`
                  ⎪best_iter + F * (P_old_{r_1,j} - P_old_{r_2,j})                                     , if variant == `DE/best/1/bin`
          tmp_j = ⎨tmp_j + F * (best_iter_j - tmp_j) + F * (P_old_{r_0,j} - P_old_{r_1,j})             , if variant == `DE/rand-to-best/1/bin`
                  ⎪best_iter_j + (P_old_{r_0,j} + P_old_{r_1,j} - P_old_{r_2,j} - P_old_{r_3,j}) * F   , if variant == `DE/best/2/bin`
                  ⎩P_old_{r_4,j} + (P_old_{r_0,j} + P_old_{r_1,j} - P_old_{r_2,j} - P_old_{r_3,j}) * F , if variant == `DE/rand/2/bin`
        j = (j + 1) mod n
    repeat for d = 1 to n:
      if not lb_d ≤ tmp_d ≤ ub_d:
        tmp_d = a uniformly distributed random number from [lb_d, ub_d]
    if f(tmp) ≤ f(P_old_i):
      P_new_i = tmp
      if f(tmp) ≤ f(best_global):
        best_global = tmp
    else:
      P_new_i = P_old_i
  swap P_old and P_new
  let worst = the agent in P_new with the highest objective value
  if |worst_1 - best_global_1| + ... + |worst_n - best_global_n| < xtol:
    terminate
  if |f(worst) - f(best_global)| < ftol:
    terminate
```

Die Wahl der zufälligen Indices ist in pagmo2 mit Durstenfelds Algorithmus implementiert:
```
r = (0, 0, 0, 0, 0)
idxs ∈ ℕ^NP = (0, 1, 2, ..., NP - 1)
repeat for j = 1 to 5:
  idx = a uniformly distributed random number from [0, NP - 1 - j]
  r_j = idxs_idx
  swap idxs_idx and idxs_{NP - 1 - j}
```

## Zusammenfassung

[1]: https://link.springer.com/article/10.1023/A:1008202821328
[2]: https://esa.github.io/pagmo2/docs/cpp/algorithms/de.html
[3]: http://www1.icsi.berkeley.edu/~storn/code.html#prac
