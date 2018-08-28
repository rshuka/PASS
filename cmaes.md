# CMA-ES

## Initialisierung

Die Implementierung von SA in [2] kann über folgende Parameter konfiguriert werden:

Name|Beschreibung|Standardwert in [2]
----|------------|-------------------
_lam_|Populationsgröße|-
_gen_|number of generations.|1
_cc_|backward time horizon for the evolution path|(4 + μ_eff / n) / (n + 4 + 2 * μ_eff / n)
_cs_|makes partly up for the small variance loss in case the indicator is zero|(μ_eff + 2) / (n + μ_eff + 5)
_c1_|learning rate for the rank-one update of the covariance matrix|2 / ((n + 1.3)^2 + μ_eff)
_cμ_|learning rate for the rank- μ update of the covariance matrix|2 * (μ_eff - 2 + 1 / μ_eff) / ((n + 2)^2 + μ_eff)
_σ0_|initial step-size
_ftol_|stopping criteria on the x tolerance|1e-6
_xtol_|stopping criteria on the f tolerance|1e-6
_force_bounds_|when true the box bounds are enforced. The fitness will never be called outside the bounds but the covariance matrix adaptation mechanism will worsen|false

## Pseudocode

```
let f(x) : ℝ^n → ℝ be the optimisation function, with lb_i ≤ x_i ≤ ub_i for i = 1, ..., n.
let pop ∈ ℝ^{n×lam} = lam uniformly distributed random agents inside the problem boundaries.

mean ∈ ℝ^n = agent in pop with best objective value
μ = ⌊lam / 2⌋
weights ∈ ℝ^μ = (log(μ + 0.5), ..., log(μ + 0.5)) - (log(1), ..., log(μ))
weights = weights * (1 / ∑weights)
μ_eff = 1 / ||weights||^2
damps = 1 + 2 * max(0, √((μ_eff - 1) / (n + 1)) - 1) + cs
χ_n = √n * (1 - 1 / (4 * n) + 1 / (21 * n^2))
σ = σ0
B ∈ ℝ^{n×n} = identity matrix
              ⎛ max(ub_1 - lb_1), 1e-6)                0            ⎞
D ∈ ℝ^{n×n} = ⎜                          ⋱                          ⎟
              ⎝             0               max(ub_n - lb_n), 1e-6) ⎠
C ∈ ℝ^{n×n} = D * D
                     ⎛ 1 / D_{1,1}          0       ⎞
invsqrtC ∈ ℝ^{n×n} = ⎜ 	            ⋱               ⎟
                     ⎝     0             1 / D_{n,n}⎠
pc ∈ ℝ^n = (0, ..., 0)
ps ∈ ℝ^n = (0, ..., 0)
counteval = 0
eigeneval = 0

repeat gen times:
  newpop ∈ ℝ^{n×lam} = lam normal distributed random agents inside the problem boundaries.
  for each row r in newpop:
    r = σ * B * D * r + mean
  if ||σ * B * D * (last row of newpop)|| < xtol:
    terminate
  if |f(worst agent in pop) - f(best agent in pop)| < ftol:
    terminate
  // We fix the bounds. Note that this screws up the whole covariance matrix machinery and worsen performances considerably.
  if force_bounds:
    for each index (i,j) in newpop:
      if newpop_{i,j} < lb_i:
        newpop_{i,j} = lb_i
      else if newpop_{i,j} > ub_i:
        newpop_{i,j} = ub_i
  pop = newpop
  counteval = counteval + lam
  // 4 - We extract the elite from this generation.
  elite ∈ ℝ^{n×μ} = (best agend in pop, 2nd best agent in pop, ..., μth best agent in pop)
  // 5 - Compute the new mean of the elite storing the old one
  meanold = mean
         μ
  mean = ∑ (row i of elite) * weights_i
        i=1
  // 6 - Update evolution paths
  ps = (1 - cs) * ps + √(cs * (2 - cs) * μ_eff) * invsqrtC * (mean - meanold) / σ
         ⎧1 , if (||p||^2 / n / (1 - (1 - cs)^(2. * counteval / lam))) < (2 + 4 / (n + 1))
  hsig = ⎨
         ⎩0 , else
  pc = (1 - cc) * pc + hsig * √(cc * (2 - cc) * μ_eff) * (mean - meanold) / σ
  // 7 - Adapt Covariance Matrix
  C_old = C
      μ
  C = ∑ ((row i of elite) - meanold) * ((row i of elite) - meanold)^T * weights_i
     i=1
  C = C / σ^2
  C = (1 - c1 - cμ) * C_old + cμ * C + c1 * ((pc * pc^T) + (1 - hsig) * cc * (2 - cc) * C_old)
  // 8 - Adapt σ
  σ = σ * exp(min(0.6, (cs / damps) * (||ps|| / χ_n - 1)))
  // 9 - Perform eigen-decomposition of C
  if counteval - eigeneval > lam / (c1 + cμ) / n / 10:
      eigeneval = counteval
      C = (C + C^T) * 1/2
      if eigen decomposition of C is successful:
        B = eigen vectors of C

            ⎛√(max(1e-20, 1st eigen value of C))                    0                   ⎞
        D = ⎜                                     ⋱                                     ⎟
            ⎝               0                        √(max(1e-20, nth eigen value of C))⎠

                       ⎛1 / D_{1,1}         0      ⎞
        invsqrtC = B * ⎜             ⋱             ⎟ * B^T
                       ⎝    0           1 / D_{n,n}⎠
```

[2]: https://esa.github.io/pagmo2/docs/cpp/algorithms/cmaes.html
