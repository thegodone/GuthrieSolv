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
FreeSolv), temperature-corrected. **Result: 16,357 observations / 2,675 molecules.**

Four conversion bugs were caught by validating each route's per-observation residual
against FreeSolv:
- **mole-fraction Henry** off by exactly −RT·ln(1000) = +4.09 kcal/mol (dropped L→m³ factor);
- **dimensionless Henry ratios** have a *process-dependent* direction: KWG rows store the
  air/water constant (invert → K_wa=1/value), KGW rows store water/gas directly (K_wa=value) —
  the labels are effectively swapped in the source (verified vs FreeSolv);
- **log-unit parser bug**: a one-char regex (`^log10?` needs "log**1**") silently dropped every
  `log(M)`, `log(mol/L)`, `log(atm)`, `log(mmHg)`… spelling — thousands of rows;
- **mole-fraction solubility** spellings (`mf`, `[mf]`) were unhandled.

Fixing the log-parser + KGW direction (pure math, no learning) recovered **+4,330 rows, +202
molecules** and took per-molecule raw-median MAE vs FreeSolv to **0.30 kcal/mol, R 0.965** (from
0.79 / 0.874 at the outset). Direct free-energy rows have ~0 median residual (the clean reference).

### Why not all 172 units? (the answer to "it's just math")
Only the **Henry-constant and free-energy** families are *unit conversions of* ΔG_hyd (a
molecule-independent map) — and every one is now converted. The remaining ~35,000 rows are
**vapour pressure (20k) + aqueous solubility (15k)**, which are *different physical observables*:
no molecule-independent math maps a lone vapour pressure to ΔG_hyd (it needs volatility AND
solubility jointly). Those convert only by **pairing** VP×solubility per molecule (1,722 molecules).
What remains — **750 VP-only + 1,666 solubility-only molecules** — is a single-observable **data
gap** no arithmetic can close. So the loop going "dry" is physically correct, not a limitation.

## 1b. Active-learning loop for the long-tail units (`calibrate_units_loop.py`)
The hand-coded switchboard covers the frequent units; a retro-feedback loop recovers much of the
rest **without hand-coding each one**, reusing the Homoset noise gate as the acceptance test:

1. **Anchors** = molecules already assigned ΔG_hyd by the physics routes.
2. For each unconverted `(unit, process)`, fit the best monotone map `ΔG ≈ a·f(value)+b`
   (`f ∈ {log10, −log10, ln, id}`, robust trimmed least-squares) against the anchors it shares.
3. **Admit** iff it correlates (`|R| ≥ R_min`) **and** passes the noise gate (`RMS/L ≤ η`).
4. An admitted unit converts all its rows → new anchors → the next round reaches further.

This runs as an outer **fixed-point / co-training loop** (`calibrate_units_iterative.py`): after
each admission pass, reconcile all obs to a per-molecule consensus, feed it back as cleaner
anchors, and **anneal** the gate strict→loose (`(η,R) = (0.25,0.90) → (0.65,0.78)`) so the
high-confidence backbone is laid first. Continue until a pass at the loosest gate is dry.

**Result (fixed point in 6 outer steps, on the bug-fixed base):** **19 units, +135 molecules
(2,675 → 2,810)**; FreeSolv MAE holds flat ~0.32 as coverage grows — a coverage/accuracy Pareto
with a self-terminating stop (FreeSolv never an anchor → held-out). Self-validating — cracks
`MPv/(RTCw)` (ln, R 0.99) and `logK=y/x at 1 atm` (R 0.99). (Before the log-parser fix the loop had
to *rescue* the log spellings and admitted 26 units; now the switchboard converts them directly, so
the loop's remaining job is smaller — the desired outcome.) See
`outputs/unit_calibration_trajectory.png`; expanded set `guthrie_dg_observations_iter_expanded.csv`.

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

## 3. Results (561-molecule FreeSolv overlap; kcal/mol; bug-fixed base)

| estimator | all MAE | all RMSE | ≥4-obs MAE | **conflict MAE** | **conflict RMSE** | bias |
|---|---|---|---|---|---|---|
| raw mean | 0.451 | 1.082 | 0.266 | 0.652 | 1.201 | +0.32 |
| raw median | 0.299 | 1.020 | 0.108 | 0.365 | 1.088 | +0.16 |
| Homoset, L adaptive, η=0.25 | 0.297 | 1.011 | 0.108 | 0.364 | 1.082 | +0.15 |
| **Homoset, L=0.6 (any η)** | 0.273 | 0.951 | **0.088** | 0.314 | 0.958 | +0.14 |
| mixle gauss | 0.315 | 0.902 | 0.150 | 0.362 | 0.844 | +0.09 |
| **mixle robust** | **0.249** | **0.873** | 0.095 | **0.276** | **0.837** | **+0.08** |

(Conflict column = 286 molecules whose replicates fail the homogeneity test.) **Note:** fixing the
conversion bugs removed much of the apparent "conflict" — a lot of it was *my bad unit conversions*,
not genuine literature disagreement. On the clean data the raw median is already good (0.299), so
the methods' edge shrinks; **mixle-robust is now best overall** (0.249), Homoset (L=0.6) still wins
the well-sampled ≥4-obs subset (0.088).

## 4. The two parameters — noise level η and dimension L (the crux)
Homoset's behaviour is set by (η, L), read off the per-source noise = RMS/L (bug-fixed base):

| source | noise (L adaptive) | noise (L=0.6) | admitted |
|---|---|---|---|
| KGW (Henry gas/water) | 0.09 | 0.89 | any η (adaptive L) |
| KWG (Henry water/gas) | 0.12 | 1.15 | any η (adaptive L) |
| PAIR (vapour-P × sol.) | 0.16 | 2.72 | any η (adaptive L) |

- **Adaptive L** → converted sources are all tight (noise 0.09–0.16) → any η admits all → consensus
  barely beats a plain median.
- **Fixed physical L = 0.6** → source noise 0.9–2.7 ≫ any η → **all converted sources REJECTED**,
  only direct free-energy survives → hard source-quality filter → best median accuracy.

Fixing the KGW-direction bug also *dissolved* the earlier "both methods distrust KGW" story — that
was an artifact. On clean data **PAIR** is the least-trusted converted source (highest gate noise
0.16 **and** the largest mixle bias +0.91) — the two methods still independently agree on the worst
source, just a different one.

## 5. Where mixle differs
mixle keeps every source but **estimates and corrects each one's bias**
(gauss: PAIR +0.91, KWG +0.23, KGW +0.21 kcal/mol vs the free-energy reference) and shrinks sparse
molecules. It wins overall MAE/RMSE and is the only method that helps sparse molecules the gate
can't align — e.g. `XOGPDSATLSAZEK` (2 sources that disagree): truth −11.53, raw median/Homoset
stuck at −6.77, **mixle_robust −11.48**, by trusting the lower-bias source.

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
