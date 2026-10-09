# Third-party notices and software requirements

No installed Python/R library, runtime, SDK, map dataset or source-paper attachment
is bundled. environment/third_party_notices/ preserves original R package
DESCRIPTION and licence/attribution records for actual plotting dependencies and
their recorded dependency closure. Necessary upstream names/contact/attribution
are preserved, not relabelled as project-owned material.

Python training dependencies and original direct pins are in
environment/requirements-training-recorded.txt. Reconstruction additionally
requires openpyxl and python-docx; their original exact reconstruction versions
were not frozen in the training protocol. These are separately installed packages.

environment/R_packages_recorded.tsv is a recorded metadata list, not a transitive
lock file or tested clean installation. RcppCGAL is referenced as a compiled
dependency but its version was not recorded in that original table. System
libraries, fonts and device/platform behaviour require separate verification.
No version or successful installation is invented here.
