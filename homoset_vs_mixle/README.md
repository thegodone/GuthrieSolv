# Homoset vs. mixle — reconciling GuthrieSolv into ΔG_hyd

Turning the heterogeneous GuthrieSolv literature dump (53,895 records, 172 unit
strings, no reconciliation code) into a single trustworthy hydration free energy
per molecule, and comparing two curation philosophies against FreeSolv.

## Contents
| file | what |
|---|---|
| `harmonize_guthrie.py` | physics unit-conversion switchboard (172 units → Ben-Naim ΔG_hyd, kcal/mol), T-corrected, keyed by InChIKey |
| `calibrate_units_loop.py` | active-learning loop: calibrate long-tail units against anchors, admit via Homoset noise gate, feed back, loop until dry (single pass: 22 units) |
| `calibrate_units_iterative.py` | outer fixed-point / co-training wrapper: reconcile → re-anchor → anneal gate strict→loose → repeat until dry (26 units, +290 molecules, FreeSolv MAE 0.399→0.32) |
| `compare_methods.py` | Homoset PS-gate consensus (L sweep) vs mixle hierarchical partial-pooling (EM, ±robust) |
| `validate_vs_freesolv.py` | per-route + per-molecule accuracy gate against FreeSolv |
| `PAPER.tex` | arXiv-style preprint (compile with `pdflatex`/`tectonic`) |
| `REPORT.md` | full technical writeup |
| `report.html` | visual comparison |
| `outputs/` | harmonized observations, per-molecule estimates, metrics JSON |

## Reproduce
```bash
python harmonize_guthrie.py        # -> outputs/guthrie_dg_observations.csv  (12,027 obs / 2,473 mols)
python compare_methods.py          # -> outputs/comparison_report.json + per_molecule_estimates.csv
python validate_vs_freesolv.py     # sanity gate vs FreeSolv
```
Requires `rdkit`, `numpy`, `pandas`, `scipy`, and a local FreeSolv `database.json`.

## Headline (556-molecule FreeSolv overlap, kcal/mol)
| method | all MAE | all RMSE | conflict MAE | bias |
|---|---|---|---|---|
| raw median | 0.399 | 1.180 | 0.581 | +0.23 |
| Homoset (L=0.6, any η) | **0.283** | 0.979 | 0.357 | +0.15 |
| mixle (gaussian) | 0.308 | **0.866** | 0.359 | **+0.07** |
| mixle (robust-t) | 0.276 | 0.949 | **0.354** | +0.09 |

The **Homoset gate has two parameters** — the noise level η (tolerated RMS/L) and the dimension
value L; a source joins the homogeneous set iff `RMS/L ≤ η`. Homoset is the better transparent
**curator** (anchor L to real experimental noise, L=0.6, and it filters bad sources); mixle is the better
**reconciler** (lowest RMSE/bias, corrects per-source bias, rescues no-reference molecules).
Use both. See `REPORT.md` / `PAPER.tex` for details.

## References
- Homoset PS gate — the "homogeneous set" / proportional-similarity standardization procedure from
  Devos, Patte, Rouault, Laffort & Van Gemert, *Standardized Human Olfactory Thresholds* (IRL Press /
  Oxford Univ. Press, 1990), originally developed in the 1970s–80s to harmonize human odor detection
  thresholds; ported here from `merge_bbb_homoset_v2.py` / `build_delta_hvap_v2_homoset.py`.
- mixle hierarchical partial pooling — [gmboquet/mixle](https://github.com/gmboquet/mixle) (`Normal(Normal(mu0,tau),sigma)`), reproduced in numpy.
- Ground truth — [MobleyLab/FreeSolv](https://github.com/MobleyLab/FreeSolv).
