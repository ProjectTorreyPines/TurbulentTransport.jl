# QLNN_ukstep26_2 — STEP quasi-linear-weight NN bundle (DT-lumped)

TGLF-trained QL-weight heads for the UKAEA STEP design point, from the runTGLFdb
`ukstep26` corpus (production_run_2, `sat1_em_azf-1`: 69,273 input.gacode slices x 9 radii
0.1..0.9, ±1 major scans of RLTS_1 / RLTS_23 / RLNS_12, 11-run minor stencil incl. VEXB_SHEAR
and BETAE, ~5.3 M TGLF runs on OLCF Defiant, Aug 2026). Trained with TrainQLweightNN
`VERSION=ukstep26tglf_v2` (20-member ensembles, 500 epochs, residual MLP 4x64); see
`TrainQLweightNN/docs/QLNN_ukstep26_2_runbook.md` for the full provenance and scoring.

| file | head | outputs |
|---|---|---|
| `energy_regressor.bson` | QL energy weights | e/DT/imp x phi/apar (6) |
| `particle_regressor.bson` | QL particle weights | e/DT/imp x phi/apar/bpar (9) |
| `momentum_regressor.bson` | QL momentum weights | e/DT/imp x phi/apar (6) |
| `eigenvalue_regressor.bson` | gamma/ky, omega/ky | 2 |
| `stability_classifier.bson` | P(unstable) | 1 |
| `width_regressor.bson` | Gaussian width | 1 (37 inputs) |
| `momentum_sign` | `+1` — TJLF-native sign, as QLNN_d3d_1 | |
| `vexb_convention` | `legacy` — ig2it inputs generated before runTGLFdb 894e2ed | |
| `input.tglf.template` | the TGLF switches the database ran with + a real STEP case | |

## TGLF settings (UK STEP recommended)

`input.tglf.template` lists them. The STEP team's settings, chosen from comparisons with
nonlinear gyrokinetics, differ from the stock sat1 defaults in `NBASIS_MAX=10` (not 4),
`NBASIS_MIN=2`, `WIDTH=1.9`, `WIDTH_MIN=0.495`, `FILTER=2.0` and
`THETA_TRAPPED = 1 - min(0.7, max(0.3, sqrt(RMIN_LOC/RMAJ_LOC)))`; electromagnetic
(`USE_BPER/USE_BPAR`), `SAT_RULE=1` in GYRO units, `ALPHA_ZF=-1`, `KYGRID_MODEL=4`, `NKY=12`,
`NMODES=2`, `USE_AVE_ION_GRID=true`, `USE_MHD_RULE=false`. `run_qlnn` takes `SAT_RULE`,
`ALPHA_ZF` and `UNITS` from the input it is given, so pass them as in the template.

## Species: database unbundled, network lumped

The database ran TGLF with **NS=4 unbundled** species (e, D, T, lumped impurity; He ash
folded into the impurity charge- and Zeff-conserving). The **networks were trained on the
(e, DT, imp) lumping** of those inputs (D and T QL weights summed — exact, the saturation
rule's intensity factor is species-independent; inputs lumped per the IMAS
bulk/impurity convention). 35 inputs incl. `MASS_2` (density-weighted D/T mass, 1.00 pure D
.. 1.25 for 50/50 D-T — the only input separating the two populations), `ZS_3`, `MASS_3`.

`run_qlnn` / `qlnn_fluctuation_spectra` / `run_modeid_qlnn` lump automatically: an NS=4 or 5
`InputTJLF` with hydrogenic species 2 and 3 is converted with `qlnn_lump_dt` (a copy; the
caller's struct is untouched), an NS=3 input is used as is, and anything else errors. With
FUSE, `act.ActorTGLF.lump_ions = true` already yields the NS=3 layout.

## Scores (test split, robust metrics; TrainQLweightNN `scripts/eval_robust_metrics.sh`)

Spearman 0.74-0.86 on every apar/bpar channel, medAE/MAD < 1 on every channel, R² after a
-5 % |y| trim positive almost everywhere (v1 of this bundle fit only the tail envelope).
Eigenvalue gamma R² 0.945, width 0.981, stability accuracy 0.968 / F1 0.982. Momentum is
the weakest head (loss-space R² 0.41-0.63).

End-to-end (TrainQLweightNN `scripts/validate_qlnn_ukstep26_tjlf.jl`, 1500 database test shots,
`run_qlnn` vs `run_tjlf` on the NS=4 database layout with the settings above): Qe/Qi/Ge within a
factor 2 of TJLF for 69/62/58 % of shots (75/82/72 % of the shots with TJLF flux > 1), within 3x
for 80/75/74 %, Spearman 0.93/0.93/0.91, median bias -20 %; momentum within 2x for 33 %. The
largest fluxes are underpredicted (summed Qe ratio 0.53), consistent with the training tail cut.
