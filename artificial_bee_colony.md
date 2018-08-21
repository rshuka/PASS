# Artificial Bee Colony

ABC wurde erstmals 2005 _(laut eigener Aussage auf der offiziellen Webseite; aber das Paper ist erst von 2006 ...)_ in [1] veröffentlicht, die hier besprochene Version [2] ist von 2010. [1] stellt den Algorithmus folgendermaßen vor:

> In this work, ABC algorithm is used for optimizing multivariable functions and the results produced by ABC, Genetic Algorithm (GA), Particle Swarm Algorithm (PSO) and Particle Swarm Inspired Evolutionary Algorithm (PS-EA) have been compared. The results showed that ABC outperforms the other algorithms.

## Initialisierung

pagmo2 [3] implementiert exakt die Version des artifical bee colony Optimierers aus [2], Algorithm 2. Der Algorithmus kann über folgende Parameter konfiguriert werden.

Name|Beschreibung|Standardwert in [3]
----|------------|-------------------
_SN_|Populationsgröße|-
_limit_|Anzahl an Fehlversuchen, bevor ein Agent von "beschäftigt" auf "erkunden" wechselt|20
_MCN_|Iterationsanzahl|1

## Pseudocode

```
let f : ℝ^n -> ℝ be the optimisation function, lb ∈ ℝ^n its lower bounds, ub ∈ ℝ^n its upper bounds.
num_eval = 0
for s = 1, ..., SN:
  X(s) <- random solution by Eq. 1
  trial(s) <- 0
  num_eval <- num_eval + 1
repeat MCN times:
  // Employed bees phase
  mi <- {s : trial(s) = max(trial)}
  for s = 1, ..., SN:
    if trial(s) < limit or s != mi:
      x' <- a new solution produced by Eq. 2
      num_eval <- num_eval + 1
      if f(x') < f(X(s)):
        X(s) <- x'
        trial(s) <- 0
      else:
        trial(s) <- trial(s) + 1
  Memorize the best solution found so far.
  // Scout bee phase
  if trial(mi) >= limit:
    X(mi) <- random solution by Eq. 1
    f_mi <- f(X(mi))
    num_eval <- num_eval + 1
    trial(mi) <- 0
  Calculate probability values p_i for the solutions using fitness values by Eqs. 3 and 4
  // Onlooker bees phase
  s <- 1
  t <- 1
  while t <= SN:
    r <- rand(0, 1)
    if r < p(s)
      t <- t + 1
      x' <- a new solution produced by Eq. 2
      num_eval <- num_eval + 1
      if f(x') < f(X(s)):
        X(s) <- x'
        trial(s) <- 0
      else:
        trial(s) <- trial(s) + 1
    s <- (s mod SN) + 1
  Memorize the best solution found so far.
```

Gleichungen:
```
x_i = lb_i + rand(0,1) * (ub_i − lb_i) for i = 1, ..., n                     (1)

x' = X(s) + e_d * ϕ * (X(s)_d - X(k)_d) , d ∈ rand(1, n), ϕ ∈ rand(-1, 1),   (2)
                                          e_d is the dth unit vector,
                                          X(k) is a randomly chosen agent

                  ⎛ SN         ⎞-1
p_m = fit(X(s)) * ⎜ ∑  fit(x_i)⎟                                           (3)
                  ⎝i=1         ⎠

            ⎧ 1 / (1+f(X(s))) , falls f(X(s)) ≥ 0
fit(X(s)) = ⎨                                                               (4)
            ⎩ 1 + |f(X(s))|   , falls f(X(s)) < 0
```

## Zusammenfassung

Die erste Phase _employed bees_ weist Ähnlichkeiten zum PSO auf: Jeder Agent bewegt sich auf genau einen zufälligen anderen Agenten zu. Anders als beim PSO bewegen sich ABC Agenten jedoch nur in einer einzigen Dimension; außerdem bleiben sie stehen, anstatt zu einer Position mit einem schlechteren objective value zu wechseln.

Die zweite Phase _scout bees_ implementiert ein Neustart-Kriterium, das dann einsetzt, wenn sich ein Agent für mindestens `limit` Evaluationen nicht bewegt hat. Damit soll der _exploration_-Aspekt der Optimierung gewährleistet werden. Pro Iteration kann nur ein einzier Agent neu gestartet werden.

Die dritte Phase _onlooker bees_ dient dazu, eine gute _exploitation_ zu erreichen, indem diejenigen Agenten die meisten Evaluationen erhalten, die global betrachtet an den besten bisher gefundenen Positionen stehen. Da in dieser Phase noch einmal `SN` Evaluationen durchgeführt werden, ist die Population quasi 2 * `SN` groß, wobei nur die Hälfte dieser Population mit einem Gedächtnis versehen ist.

[1]: https://link.springer.com/article/10.1007/s10898-007-9149-x
[2]: https://abc.erciyes.edu.tr/pub/NevImpOfABC.pdf
[3]: https://esa.github.io/pagmo2/docs/cpp/algorithms/bee_colony.html
