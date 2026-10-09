# Biofertilizer yield transferability — code-only

This repository preserves the research's Python analysis and statistical-validation
implementations and R scientific-visualization source. It is a code repository,
not a complete research-data archive. Research datasets, predictions, bootstrap
results and manuscript figure files are not distributed in this code-only version.

## Source preservation

Of the 29 scientific scripts in the earlier code-only candidate, 26 are now
byte-identical to the private archived originals. Three have narrowly documented
public changes: external source metadata are read from an authentic input file;
one Figure 3 caption obtains its result values from its plot input; the final
Figure 3 wrapper checks the changed template's explicit original/public identities.
No algorithm function, model hyperparameter, statistical formula or plotting
function was masked behind a private-literal lookup.

The original Phase−2 workbook-reading script explore.py was also restored because
the audit actually depends on its explore.json output. This makes 30 archived
scientific source files. The filename is historical: running it reads real source
data and writes private row/schema output. It was NOT executed here.

All original scientific sources and frozen results remain private and unchanged.
provenance/source_code_map.tsv distinguishes exact originals from public derivatives.
There is no EJA_PRIVATE_LITERALS dependency and no private-source restoration step.

## What is intentionally absent

No source attachments, analysis matrices, cohort/fold membership records, OOF
predictions, bootstrap draws/results, fitted-model outputs, supplementary data
tables, actual plot-data tables, image/PDF output, geometry caches or historical
report payloads are included. No dummy datasets are supplied.

Method constants, scientific variable names, units, model settings, seeds,
thresholds, plotting styles, explanatory comments, frozen cohort-size/QC contracts
and source-file fingerprints remain visible because they explain the implementation.
They are not individual study records or vectors of scientific results.
Some historical captions/comments describe the study context or conclusions;
this package does not claim to remove every previously disclosed fact.

## Safe independent checks

From the package root with Python 3.12 or later:

    python code/verify_code_only.py
    python code/preflight.py
    python code/preflight.py --require-inputs
    python code/validate_external_metadata.py --schema
    python code/prepare_linked_sources.py --help

The require-inputs command is expected to exit 2 when scientific inputs are absent.
These commands do not run research scripts, fit models, generate bootstrap samples
or create figures. Syntax verification and program checks are not scientific
replication. See RUN_GUIDE.md and docs/INPUT_FORMATS.md for real prerequisites.

## Run capability

The scientific workflows cannot run using this package alone because authentic
datasets, frozen protocols/memberships, predictions/archives, plotting tables and
protected historical assets are deliberately absent. Several exact dependencies
and environment details must also be restored. No training or figure regeneration
was performed; full end-to-end and cross-platform reproduction remain NOT VERIFIED.

Two independent source-preparation steps are supplied as code: obtaining/verifying
publisher attachments and extracting the S1 paragraph input, and the original
workbook reader producing explore.json. These write private data and must be used
only in a separate private workspace following the guide, never in the public tree.

## Publication history and licences

The previous public GitHub v1.0.0 contained the complete scientific reproducibility
materials. It was publicly available; the code-only change does not erase that
history. The default branch is rebuilt as a new code-only root history; the old
v1.0.0 release, uploaded archive and tag are withdrawn as part of this transition.
The complete scientific archive and original Git history are retained privately.
Previously downloaded copies, caches or access to unreferenced old objects may
persist. No new Release or uploaded ZIP accompanies this code-only transition.

Project-authored code and specified software files: MIT (2026 Rongxin Liu and
Yusheng Niu). Project-authored explanatory documents: separately scoped CC BY 4.0.
Upstream licences/attribution remain intact. See LICENCE_SCOPE.md,
THIRD_PARTY_NOTICES.md, DATA_SOURCES.md and CITATION.cff. No new DOI or
release version has been assigned to this code-only source snapshot.
