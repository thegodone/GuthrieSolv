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

**Homoset** (faithful port of `merge_bbb_homoset_v2.py` / `build_delta_hvap_v2_homoset.py`):
sources = measurement processes, with FEOH+DGS as the trusted **reference**.
- *Stage 1 — source alignment via the PS gate.* For each non-reference source, take
  its per-molecule diffs vs the reference and compute
  `offset = median(diffs)`, `RMS = √mean((diffs−offset)²)`, **`PS = 1 − RMS/L`**;
  accept the source (and shift it by −offset) iff `PS ≥ ps_crit(n,α)`,
  `ps_crit = 1 − √(χ²(α,n−1)/3n)`, α = 0.01.
- *Stage 2 — consensus.* Per molecule, **median** over accepted/aligned obs; flag a
  molecule as conflicted if its range ≥ `outlier_threshold`.
- Swept (per your instruction): **L ∈ {adaptive = max((max−min)/2, 0.5) on the diffs,
  fixed physical 0.6}** × **threshold ∈ {2.0, 1.0} kcal/mol}**.

**mixle** (`Normal(Normal(μ0,τ), σ_source)` with per-source bias, EM):
`y_ij ~ N(μ_i + b_s, σ_s²)`, `b_reference = 0`, `μ_i ~ N(μ0, τ²)`. Jointly estimates
per-source bias `b_s`, per-source noise `σ_s`, and per-molecule shrunk posterior μ_i.
Robust variant uses a Student-t (ν=4) observation model → per-obs weights down-weight
outliers.

## 3. Results (556-molecule FreeSolv overlap; kcal/mol)

| estimator | all MAE | all RMSE | ≥4-obs MAE | **conflict MAE** | **conflict RMSE** | bias |
|---|---|---|---|---|---|---|
| raw mean | 0.570 | 1.226 | 0.406 | 1.134 | 1.698 | +0.40 |
| raw median | 0.399 | 1.180 | 0.203 | 0.706 | 1.611 | +0.23 |
| Homoset, L adaptive | 0.392 | 1.159 | 0.202 | 0.697 | 1.592 | +0.22 |
| **Homoset, L fixed 0.6** | **0.283** | 0.979 | **0.083** | 0.410 | 1.181 | +0.15 |
| **mixle gauss** | 0.308 | **0.866** | 0.141 | **0.382** | **0.875** | **+0.07** |
| mixle robust | 0.276 | 0.949 | 0.116 | 0.388 | 1.152 | +0.09 |

(The outlier threshold 2.0 vs 1.0 changes *only which molecules are flagged for
curation*, not the consensus value — so both give identical accuracy.)

## 4. What the noise scale L does (the crux)
Homoset's behaviour is entirely governed by L, exactly as you flagged:
- **Adaptive L** = half the *data's own* diff-range is large (5.9–10.3 kcal/mol here
  because a few molecules disagree wildly) → PS ≈ 0.86 > crit ≈ 0.47 → **every source
  passes** → the noisy Henry/pairing sources stay in → MAE 0.39.
- **Fixed physical L = 0.6** (FreeSolv's experimental noise) → RMS ≫ 0.6 → PS < 0 <
  crit → **all converted sources are REJECTED**, only direct free-energy survives → the
  PS gate becomes a hard source-quality filter → MAE 0.28.

So with a physically-anchored L the Homoset gate is the best *curator* — it throws
away exactly the sources that FreeSolv-validation shows are noisy.

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
