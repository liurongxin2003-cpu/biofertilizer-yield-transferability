# Input formats and missing prerequisites

This document contains schema and path descriptions only, never scientific rows.
Every scientific entry requires authentic inputs; none runs end to end with public
files alone. config/input_contracts.json records the key prerequisites per script,
with explicit non-exhaustive scope.

## Core TSV files

Use UTF-8 with tab separation and the original header order. Keep IDs as strings;
do not coerce IDs into predictors or manufacture missing fields.

### Schutz_2018_phase0/analysis_matrix.tsv

Columns: observation_id, study_id, pair_id, logRR, crop, climate_zone, soil_texture, soil_organic_matter_pct, soil_ph_water, available_p_kg_ha, inoculant_function, n_applied_kg_ha, k_applied_kg_ha, p_applied_kg_ha, irrigation_status, fertilization_type, treatment_dry_yield_kg_ha, control_dry_yield_kg_ha, treatment_excel_row, control_excel_row.

### Schutz_2018_phase0/split_manifest.tsv

Columns: fold_id, design, fold_number, primary_or_secondary, held_out_level, held_out_function_code, train_n, test_n, purge_n, train_studies_n, test_studies_n, purge_studies_n, study_overlap_train_test, study_overlap_train_purge, train_study_ids, test_study_ids, purge_study_ids, train_ids, test_ids, purge_ids, train_ids_sha256, test_ids_sha256, purge_ids_sha256.

### Schutz_2018_phase0/data_dictionary.tsv

Columns: column_index, column, role, dtype, unit_or_code, source_field, missing_allowed, missing_token, missing_n, definition, future_handling.

Numeric missing values use the original NA definition. The matrix contains the
original logRR, source-yield/provenance and ID columns; only the frozen 12-field
predictor whitelist enters models. Split membership lists are exact original IDs,
not newly sampled splits; empty membership is NA.

## JSON/NPZ evidence formats

- explore.json: object with rows_data (authentic source-row objects) and schema
  (objects with field names), plus the original source reader's metadata.
- s1_probe.json: paragraphs list containing the authentic S1 reference text.
- delivery_tables.json: table-name to ordered row-object lists. Historical
  validators compare delivered TSVs against these private original values.
- frozen_fold_keys.json: frozen original fold definitions/member IDs from the
  prior phase; see its original reader rather than generating new assignments.
- Phase1 protocol.json: model/predictor/seed/bootstrap/decision definitions and
  authentic input fingerprints. config/phase1_protocol_reference.json preserves
  its original bytes as a reference; the expected evidence path is not populated.
- Phase1 run/fold_*.json: original per-fold preprocessing, train/test/purge identity,
  predictions and fit records. These actual values are withheld.
- Phase2 frozen_variants.json: objects containing authentic rows, split memberships,
  predictor groups and cohort metadata per frozen variant. No rows are bundled.
- bootstrap_draws.npz and oof_predictions.npz: original NumPy archives with named
  study multiplicities/draws and predictions. No substitute archives are provided.
- Figure plot-data directories: authentic TSV panel tables, QC JSON and sometimes
  saved geometry RDS/earlier images. Final scripts' read paths define exact fields;
  no draw data or image files are distributed.

## Authentic external_source_metadata.json

The schema specifies required keys and types without real values:
field_sites is a list of objects with year (integer), site (string), latitude and
longitude (finite numbers). mapped_predictors contains crop, climate_zone,
soil_texture, inoculant_function, irrigation_status and fertilization_type,
using the original labels/nominal codes. potassium_product_rate_kg_ha is a nonnegative
source product dose and potassium_k2o_mass_fraction is a fraction in [0,1].
The public algorithm retains the original multiplication by 0.8302 to obtain
elemental K. The file must be reconstructed from authentic source evidence;
the format checker is not a source-authenticity check.

## Specific blockers by entry

| Entry | Key reason it cannot run from public files alone |
| --- | --- |
| Schutz_2018_phase_minus2/scripts/audit.py | Authentic source/frozen inputs are withheld; no private-literal sidecar is required |
| Schutz_2018_phase_minus2/scripts/validate.py | Authentic source/frozen inputs are withheld; no private-literal sidecar is required |
| Schutz_2018_phase_minus1/scripts/audit.py | Authentic source/frozen inputs are withheld; no private-literal sidecar is required |
| Schutz_2018_phase_minus1/scripts/validate.py | Authentic source/frozen inputs are withheld; no private-literal sidecar is required |
| Schutz_2018_phase0/scripts/prepare.py | Authentic source/frozen inputs are withheld; no private-literal sidecar is required |
| Schutz_2018_phase0/scripts/validate.py | Authentic source/frozen inputs are withheld; no private-literal sidecar is required |
| Schutz_2018_phase1/scripts/freeze_protocol.py | Authentic source/frozen inputs are withheld; no private-literal sidecar is required |
| Schutz_2018_phase1/scripts/run_models.py | Authentic source/frozen inputs are withheld; no private-literal sidecar is required; exact recorded training dependencies also required |
| Schutz_2018_phase1/scripts/evaluate.py | Authentic source/frozen inputs are withheld; no private-literal sidecar is required |
| Schutz_2018_phase1/scripts/validate_outputs.py | Authentic source/frozen inputs are withheld; no private-literal sidecar is required |
| Schutz_2018_phase2/scripts/prepare.py | Authentic source/frozen inputs are withheld; no private-literal sidecar is required |
| Schutz_2018_phase2/scripts/run_models.py | Authentic source/frozen inputs are withheld; no private-literal sidecar is required; exact recorded training dependencies also required |
| Schutz_2018_phase2/scripts/evaluate.py | Authentic source/frozen inputs are withheld; no private-literal sidecar is required |
| Schutz_2018_phase2/scripts/validate_outputs.py | Authentic source/frozen inputs are withheld; no private-literal sidecar is required |
| Schutz_2018_external_Sato_Vermeersch/scripts/audit.py | Authentic source/frozen inputs are withheld; no private-literal sidecar is required |
| Schutz_2018_external_Vermeersch_validation/scripts/run_external.py | Authentic source/frozen inputs are withheld; no private-literal sidecar is required |
| Schutz_2018_figures/Figure1/figure1.R | Saved plotting tables/QC/caches and protected historical graphics are withheld; package/fonts/graphics devices also required |
| Schutz_2018_figures/Figure1/figure1_v2.R | Saved plotting tables/QC/caches and protected historical graphics are withheld; package/fonts/graphics devices also required |
| Schutz_2018_figures/Figure2/figure2_final.R | Saved plotting tables/QC/caches and protected historical graphics are withheld; package/fonts/graphics devices also required |
| Schutz_2018_figures/Figure3/figure3_v4.R | Saved plotting tables/QC/caches and protected historical graphics are withheld; package/fonts/graphics devices also required |
| Schutz_2018_figures/Figure3/figure3_final.R | Saved plotting tables/QC/caches and protected historical graphics are withheld; package/fonts/graphics devices also required |
| Schutz_2018_figures/Figure4/figure4_v3.R | Saved plotting tables/QC/caches and protected historical graphics are withheld; package/fonts/graphics devices also required |
| Schutz_2018_figures/Figure4/figure4_final.R | Saved plotting tables/QC/caches and protected historical graphics are withheld; package/fonts/graphics devices also required |
| Schutz_2018_figures/Figure5/figure5.R | Saved plotting tables/QC/caches and protected historical graphics are withheld; package/fonts/graphics devices also required |
| Schutz_2018_figures/Figure5/refine_weighted_voronoi.R | Saved plotting tables/QC/caches and protected historical graphics are withheld; package/fonts/graphics devices also required |
| Schutz_2018_figures/Figure6/figure6_draft.R | Saved plotting tables/QC/caches and protected historical graphics are withheld; package/fonts/graphics devices also required |
| Schutz_2018_figures/Figure6/figure6_v2.R | Saved plotting tables/QC/caches and protected historical graphics are withheld; package/fonts/graphics devices also required |
| Schutz_2018_figures/Figure6/figure6_v3.R | Saved plotting tables/QC/caches and protected historical graphics are withheld; package/fonts/graphics devices also required |
| Schutz_2018_figures/Figure6/figure6_final.R | Saved plotting tables/QC/caches and protected historical graphics are withheld; package/fonts/graphics devices also required |
| Schutz_2018_phase_minus2/scripts/explore.py | Source workbook not distributed; executing this historical script writes private explore.json containing raw rows |

Reconstruction scripts need upstream publisher attachments and prior phase evidence;
model scripts additionally need the exact frozen protocols/splits and recorded
scientific software; evaluators/validators need authentic fits/predictions/draws
and historical delivery/QC files. R scripts additionally need saved plotting
tables/caches and protected historic graphics. These limitations are not repaired
by hiding constants or by supplying dummy values.

Current environment discovery and exact per-entry results are recorded in the
external CODE_ONLY_FINAL_VERIFICATION.md, not asserted as universal requirements
being installed on every machine.
