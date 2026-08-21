# QLNN_v39

v37 recipe (33 inputs incl. USE_AVE_ION_GRID, canon off, KY_SPECTRUM_MIN=0.15,
DERIV_KY_MIN=0) trained on the 54k joint d3d+mastu+nstx dataset
(qlnn_training_subset_merged_result_dict_dmn54p0_dict_test_54004.json). All six
heads retrained, ens20, SIZE_TAG=54k.

## Momentum sign convention (BREAKING vs v37/v18)

The momentum head was trained with MOMENTUM_SIGN_FLIP=true: the CGYRO momentum
QL targets were multiplied by -1 at prep, aligning the output with the TGLF
stress_par sign convention. The raw CGYRO DB (and every QLNN <= v38 momentum
head) uses the OPPOSITE sign — a 2026-08-20 per-fn meta-analysis showed the DB
momentum targets are internally consistent across all 175 batches but globally
anti-aligned with TGLF / sign(VPAR_SHEAR). The `momentum_sign` sidecar is set
to +1 accordingly (TurbulentTransport.loadqlnnbundle multiplies the momentum
NN output by it; the missing-file default of -1 is for CGYRO-convention
bundles like v37). Downstream consumers bypassing the sidecar must expect the
TGLF-aligned (inverted-vs-v37) Pi sign.
