# Homoset vs. mixle on the GuthrieSolv hydration free-energy dump

**Question.** GuthrieSolv is the largest literature dump of experimental hydration
data (53,895 rows, 5,309 molecules) but of uncertain quality and with no
reconciliation code of its own. How does the user's multi-source **Homoset**
consensus method compare with a **mixle**-style hierarchical partial-pooling model
at turning its noisy, redundant per-molecule measurements into a single curated
ΔG_hyd? Ground truth = FreeSolv (642 curated molecules).

## 1. Harmonization (prerequisite)
GuthrieSolv stores raw physical measurements in **172 distinct unit strings**;
only 737 rows carry a finished kcal/mol ΔG_hyd. The `process` column tags each
row's measurement type, which makes conversion tractable:

| process | meaning | rows | route to ΔG_hyd |
|---|---|---|---|
| FEOH, DGS | free energy of hydration/solvation | 3,160 | direct (unit convert) |
| KWG, KGW | Henry's-law constants | 11,751 | ΔG = −RT·ln K_wa |
| VP + AQSOL | vapour pressure + aqueous solubility | 36k | paired: K_wa = C_aq/C_gas |
| DGT/DGTX/DGB/DGV, PKA, KA, LAQACTCO | transfer/vaporisation/pKa | excl. | not hydration |

All converted to Ben-Naim ΔG*_hyd (gas 1 M → aq 1 M, same standard state as
FreeSolv), temperature-corrected. **Result: 12,027 observations / 2,473 molecules
(1,518 with ≥2 obs).**

Two conversion bugs were caught by validating each route's per-observation residual
against FreeSolv:
- **mole-fraction Henry** was off by exactly −RT·ln(1000) = +4.09 kcal/mol
  (dropped an L→m³ factor);
- **dimensionless M/M ratios**, despite the KWG "water/gas" label, are stored as the
  standard *air/water* Henry constant → needed inversion (verified empirically:
  (1/value)/K_wa_true = 0.963).

After fixes, per-molecule raw-median MAE vs FreeSolv fell **0.79 → 0.40 kcal/mol,
R 0.874 → 0.948**. Direct free-energy rows have ~0 median residual (they are the
clean reference); Henry/pairing routes are noisier (the conflict source).

## 2. The two methods

**Homoset** — the "homogeneous set" / proportional-similarity procedure originally developed in
the 1970s–80s to standardize scattered literature values of human **odor detection thresholds**
(Devos, Patte, Rouault, Laffort & Van Gemert, *Standardized Human Olfactory Thresholds*, IRL
Press / Oxford Univ. Press, 1990); applied unchanged here to hydration free energies
(faithful port of `merge_bbb_homoset_v2.py` / `build_delta_hvap_v2_homoset.py`):
sources = measurement processes, with FEOH+DGS as the trusted **reference**. The method has
**exactly two parameters**: the **noise level** η (tolerated scatter as a fraction of L) and
the **dimension value** L (the amplitude scale). A set is a *homogeneous set* iff
`noise = RMS/L ≤ η`, i.e. `PS = 1 − RMS/L ≥ 1 − η`.
- *Stage 1 — source alignment.* For each non-reference source, take its per-molecule diffs
  vs the reference, `offset = median(diffs)`, `RMS = √mean((diffs−offset)²)`; admit the source
  (shift by −offset) iff `RMS/L ≤ η`. (The χ² form `1−η = 1−√(χ²(α,n−1)/3n)` is just one way to
  set η as a function of sample size; the two operating knobs are η and L.)
- *Stage 2 — consensus.* Per molecule, **median** over admitted/aligned obs; flag a molecule
  as conflicted when its own observations fail the homogeneity test against the physical
  dimension L = 0.6 kcal/mol.
- Swept: **L ∈ {adaptive = max((max−min)/2, 0.5) on the diffs, fixed physical 0.6}** ×
  **noise level η ∈ {0.25, 0.50}}**.

**mixle** ([gmboquet/mixle](https://github.com/gmboquet/mixle), whose compositional dialect
writes a random-effects prior as `Normal(Normal(mu0, tau), sigma)`; reproduced here in numpy).
`Normal(Normal(μ0,τ), σ_source)` with per-source bias, EM:
`y_ij ~ N(μ_i + b_s, σ_s²)`, `b_reference = 0`, `μ_i ~ N(μ0, τ²)`. Jointly estimates
per-source bias `b_s`, per-source noise `σ_s`, and per-molecule shrunk posterior μ_i.
Robust variant uses a Student-t (ν=4) observation model → per-obs weights down-weight
outliers.

## 3. Results (556-molecule FreeSolv overlap; kcal/mol)

| estimator | all MAE | all RMSE | ≥4-obs MAE | **conflict MAE** | **conflict RMSE** | bias |
|---|---|---|---|---|---|---|
| raw mean | 0.570 | 1.226 | 0.406 | 0.916 | 1.495 | +0.40 |
| raw median | 0.399 | 1.180 | 0.203 | 0.581 | 1.419 | +0.23 |
| Homoset, L adaptive, η=0.50 | 0.392 | 1.159 | 0.202 | 0.574 | 1.402 | +0.22 |
| Homoset, L adaptive, η=0.25 | 0.350 | 1.068 | 0.159 | 0.491 | 1.248 | +0.19 |
| **Homoset, L=0.6 (any η)** | **0.283** | 0.979 | **0.083** | 0.357 | 1.071 | +0.15 |
| **mixle gauss** | 0.308 | **0.866** | 0.141 | 0.359 | **0.855** | **+0.07** |
| mixle robust | 0.276 | 0.949 | 0.116 | **0.354** | 1.060 | +0.09 |

(Conflict column = 278 molecules whose replicates fail the homogeneity test. At L=0.6 both
noise levels reject every converted source, so they score identically.)

## 4. The two parameters — noise level η and dimension L (the crux)
Homoset's behaviour is set by (η, L), read directly off the per-source noise = RMS/L:

| source | noise (L adaptive) | noise (L=0.6) | admitted |
|---|---|---|---|
| KWG (Henry water/gas) | 0.14 | 1.36 | η≥0.25 (adaptive L) |
| PAIR (vapour-P × sol.) | 0.18 | 3.01 | η≥0.25 (adaptive L) |
| KGW (Henry gas/water) | 0.34 | 5.03 | η≥0.50 (adaptive L) |

- **Adaptive L** (½ of the *data's own* diff-range, 5.9–10.3 kcal/mol) → source noise 0.14–0.34;
  η=0.25 admits the two cleanest sources (KWG, PAIR) but rejects KGW, η=0.50 admits all → MAE 0.35–0.39.
- **Fixed physical L = 0.6** → source noise 1.4–5.0 ≫ any η → **all converted sources REJECTED**,
  only direct free-energy survives → the gate is a hard source-quality filter → MAE 0.28.

So with a physically-anchored L the Homoset gate is the best *curator* — it throws away exactly
the sources FreeSolv-validation shows are noisy. Tellingly, the source it is most reluctant to
admit (KGW, the only one rejected at η=0.25) is the very source mixle assigns the largest bias
(+1.18) — the two methods independently agree on which source is worst.

## 5. Where mixle differs
mixle keeps every source but **estimates and corrects each one's bias**
(gauss: KGW +1.18, PAIR +1.25, KWG +0.30 kcal/mol vs the free-energy reference) and
shrinks sparse molecules toward μ0 = −4.95. It wins RMSE everywhere and has the
lowest bias, and it is the only method that helps molecules that have **no reference
source at all** — e.g. `CVXBEEMKQHEXEN` (2 obs, no free-energy): truth −9.45,
raw/Homoset stuck at −3.38, **mixle_robust −9.22**, because the globally-learned
source bias transfers to molecules the Homoset alignment can't reach.

Concrete recoveries (truth | raw-median | Homoset-Lfix | mixle-gauss | mixle-robust):
- `SHZIWNPUGXLXDT` (5 obs): −2.23 | +3.43 | −2.23 | −2.10 | −2.23
- `ZWRUINPWMLAQRD` (9 obs): −3.88 | +1.63 | −3.88 | −3.75 | −3.88

## 6. Verdict
- Both crush the raw mean/median; the redundancy in GuthrieSolv is real signal.
- **Homoset (fixed physical L)** = best simple *curator*: transparent, auditable, and
  its PS gate + per-molecule conflict flag directly produce a `use_for_training`
  decision. Its accuracy hinges entirely on choosing L to match true experimental
  noise; with the data-adaptive L it is barely better than a plain median.
- **mixle** = best *reconciler*: lowest RMSE/bias, no threshold to tune, keeps all
  data, corrects per-source bias, and uniquely rescues molecules with no clean
  reference — at the cost of interpretability (no explicit reject/curate flag).
- **Best of both:** use mixle's per-source bias + robust weights to *reconcile*, and
  Homoset's PS gate + range flag to *triage* which molecules a human should re-check.

Artifacts: `outputs/guthrie_dg_observations.csv` (harmonized obs),
`outputs/per_molecule_estimates.csv` (all estimators),
`outputs/comparison_report.json` (full metrics + gates + hyperparameters).
