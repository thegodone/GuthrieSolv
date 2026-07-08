# Homoset vs. mixle — reconciling GuthrieSolv into ΔG_hyd

Turning the heterogeneous GuthrieSolv literature dump (53,895 records, 172 unit
strings, no reconciliation code) into a single trustworthy hydration free energy
per molecule, and comparing two curation philosophies against FreeSolv.

## Contents
| file | what |
|---|---|
| `harmonize_guthrie.py` | full unit-conversion switchboard (172 units → Ben-Naim ΔG_hyd, kcal/mol), T-corrected, keyed by InChIKey |
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
| raw median | 0.399 | 1.180 | 0.706 | +0.23 |
| Homoset (L fixed 0.6) | **0.283** | 0.979 | 0.410 | +0.15 |
| mixle (gaussian) | 0.308 | **0.866** | **0.382** | **+0.07** |
| mixle (robust-t) | 0.276 | 0.949 | 0.388 | +0.09 |

Homoset is the better transparent **curator** (its PS-gate noise scale L is the crux —
anchor it to real experimental noise and it filters bad sources); mixle is the better
**reconciler** (lowest RMSE/bias, corrects per-source bias, rescues no-reference molecules).
Use both. See `REPORT.md` / `PAPER.tex` for details.

## References
- Homoset PS gate — ported from the user's `merge_bbb_homoset_v2.py` / `build_delta_hvap_v2_homoset.py`.
- mixle hierarchical partial pooling — [gmboquet/mixle](https://github.com/gmboquet/mixle) (`Normal(Normal(mu0,tau),sigma)`), reproduced in numpy.
- Ground truth — [MobleyLab/FreeSolv](https://github.com/MobleyLab/FreeSolv).
