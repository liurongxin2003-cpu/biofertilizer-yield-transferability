# Run guide: useful source, authentic inputs required

## 1. Read-only checks available without research data

verify_code_only.py checks the file allow-list, SHA-256, Python syntax and the
source identity map. preflight.py checks input presence and module discovery,
without importing research scripts or reading study observations. Its prerequisites
are key examples, not an exhaustive dynamic dependency resolver.

    python code/preflight.py --stage phase1 --require-inputs
    python code/preflight.py --inputs-root /private/phase-directories

The first command should stop with exit 2 when inputs are absent. The second
accepts a directory containing Schutz_2018_phase0, Schutz_2018_phase1, etc.
Presence of files is not proof of valid content or scientific reproducibility.

R syntax checks use parse(file=...) with Rscript --vanilla. Do not source() the
scientific scripts to test syntax. Most archived scripts execute at top level;
even --help is not a safe test for those scripts.

## 2. Original source acquisition and input reconstruction

Obtain the authentic publisher attachments listed in DATA_SOURCES.md. Do not
redistribute them as part of this code-only package. For the retained source adapter,
--sources-root must contain the pipeline/... paths in provenance/linked_sources.tsv:

    python code/prepare_linked_sources.py --sources-root /private/downloads --check-only

The adapter verifies original attachment sizes/SHA-256 and extracts S1 reference
paragraphs deterministically. --write-workspace creates a NEW directory outside
this package and writes the authentic attachments plus minimal s1_probe.json;
use only a private destination after author authorisation.

The restored historical pipeline/Schutz_2018_phase_minus2/scripts/explore.py reads
the workbook's active worksheet and writes evidence/explore.json containing
actual raw rows, field summaries and last-row previews. The authentic workbook's
active worksheet was checked structurally to be Sheet3; no cells were read for
that structural check. Before any future run, verify the accepted workbook hash
and active-sheet identity. This script was not executed during package preparation.

Then the original audit.py pairs source arms using explore.json and S1 reference
paragraphs. These are genuine missing inputs, not hidden method constants.

## 3. Scientific workflow order and non-runnable entries

Phase−2 extraction/pairing → Phase−1 cohort/availability audit → Phase0 matrix and
frozen split reconstruction → Phase1 fixed models and evaluation → Phase2
prespecified sensitivities and evaluation. External structural audit precedes
external validation. Details and schemas are in docs/INPUT_FORMATS.md.

run_models.py and run_external.py train models. evaluate.py and some legacy
validate scripts generate/reconstruct bootstrap draws. freeze_protocol.py and
prepare.py create evidence/protocol files; historical validate scripts expect
delivery tables, original reports or prior QC artefacts. None were run here.

Use a NEW private workspace with authentic frozen inputs for any separately
authorised scientific reproduction. Preserve recorded software versions and hashes.
Do not copy config/*reference.json over the workflow's evidence/protocol.json:
the Phase1 reference is archived byte-identical, but dynamic evidence/runtime files
and exact path placement are additional requirements.

## 4. External metadata: one genuine source input, not a literal bag

The minimally edited external structural audit reads
Schutz_2018_external_Sato_Vermeersch/evidence/external_source_metadata.json.
The exact format is config/external_source_metadata.schema.json. It contains
authentic site/year/GPS observations, six mapped source predictors, and K product
dose/composition. No values or substitute rows are distributed here.

    python code/validate_external_metadata.py --check /private/external_source_metadata.json

A format pass does not verify source authenticity, freeze identity or results.
The field ordering, grid lookup, strict/missing predictor mapping and original
K2O-to-elemental-K conversion remain in public source. The private author-side
check reconstructed this metadata from authentic archived source for code-equivalence
and format testing only; it is not part of this public directory.

Regenerating audit JSON can alter provenance hashes because a new metadata input
exists. Such outputs must remain separate from frozen original evidence; no
byte-identical audit regeneration or external training was claimed here.

## 5. R figure dependencies and Figure 3 template identity

Keep the final source choices: Figure1/figure1_v2.R, Figure2/figure2_final.R,
Figure3/figure3_final.R, Figure4/figure4_final.R, Figure5/figure5.R,
Figure6/figure6_final.R. Historical R scripts are retained only where required
as source/protected dependencies; they do not replace final figures.

Actual plotting inputs and protected draft/final image assets are absent. Some
scripts require both data tables and earlier graphics/QC merely to verify that
those earlier files were unchanged. Figure5 may generate/refine a geometry cache
if run; this was not done here. Figure6 reads historical QC source-path mappings
and protected draft/v2/v3 assets: supplying numerical tables alone is insufficient.
Cross-platform fonts/devices, RcppCGAL and other packages remain limitations.

Figure3_v4's four hard-coded response-span scalars now come from its existing
means table when run. Its drawing/composition slice is unchanged. The final wrapper
checks either the original template hash or the explicit original/final pair in
provenance/source_code_map.tsv. It does not pretend an edited hash equals the
original, mutate frozen QC, or weaken scientific data-file hash checks.
Source-text replacement compatibility was checked without plotting.

## 6. What the verification establishes

Syntax, file/source identities, preserved method/function expressions, standalone
tools, authentic metadata-format/equivalence checks and pure-program helper tests
are distinct from scientific replication. No fake study dataset, predicted values
or bootstrap outputs were constructed. Full training, statistical evaluation and
Figure1–6 generation are NOT VERIFIED.
