# Nextclade Dataset for "Enterovirus D68" based on reference "AY426531.1"

| Key              | Value                                                                 |
|------------------|-----------------------------------------------------------------------|
| authors          | Nadia Neuner-Jehle, Emma B. Hodcroft                                  |
| name             | Enterovirus D68                                                       |
| reference        | AY426531.1                                                            |
| dataset path     | ...                                                                   |



The dataset represents the species *Enterovirus D68*, which is part of the *Enterovirus D* species group. It is based on the complete genome of the Fermon strain (AY426531.1), one of the earliest sequenced EV-D68 isolates. This dataset supports genogroup classification, mutation annotation, and quality control in clinical and surveillance contexts.

The example sequences include a representative subsample of EV-D68 genotypes across major genogroups A (A1–A2), B (B1–B3), C, and D.

## Scope of this dataset

This dataset is intended for the analysis of *Enterovirus D68* genomes. It is not suitable for other enterovirus types, which may be flagged as outliers or unassigned.

## Features

This dataset was generated using a custom workflow based on public EV-D68 sequences from GenBank and curated trees representing the known genogroup structure. Clade and genogroup naming is informed by literature (e.g. Midgley et al., Schuffenecker et al.) but not tied to a fixed lineage naming system.

Quality control thresholds were adjusted for the diversity and coding structure of enteroviruses. In particular:
- `frameShifts` are strictly penalized (scoreWeight: 100)
- `divergence.maxDivergence` set to 0.15
- Suitable for full genomes or full VP1 sequences

## What is a Nextclade dataset?

Read more about Nextclade datasets in the [Nextclade documentation](https://docs.nextstrain.org/projects/nextclade/en/stable/user/datasets.html).