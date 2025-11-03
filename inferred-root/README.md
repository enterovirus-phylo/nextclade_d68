# Inferred Ancestral Sequence for Enterovirus D68

This repository provides a reproducible workflow for generating a static inferred ancestral ("root") sequence for Enterovirus D68 (EV-D68), designed for use as a custom reference in [Nextclade](https://clades.nextstrain.org/) analyses.


## Overview

Phylogenetic analyses, such as those performed by Nextclade and Augur, benefit from a high-quality, dataset-representative reference sequence. This workflow infers such a static sequence by reconstructing the ancestral sequence from your entire  dataset, allowing for more accurate mutation and clade assignments.

**Workflow summary:**
1. **Phylogenetic Tree Construction:** All sequences in the provided dataset are aligned and a maximum-likelihood tree is built.
2. **Ancestral Sequence Inference:** The [Augur](https://github.com/nextstrain/augur) toolkit is used to infer the ancestral (root) sequence, labeled `NODE_0000000` in the output FASTA.
3. **Gap Correction:** The script [`fix_root_gaps.py`](scripts/fix_root_gaps.py) replaces any gaps (`-`) or ambiguous bases (`N`) in the inferred sequence with the corresponding reference nucleotides, ensuring a contiguous and biologically plausible inferred sequence.
4. **Export:** The cleaned ancestral sequence FASTA is available for use as a custom reference in Nextclade. The corresponding metadata is also provided and should be kept in sync.


## Getting Started
### Running the Workflow

To generate the inferred ancestral root:

```bash
snakemake --cores 9 all
```

This will:
- Build a phylogenetic tree and infer the ancestral root
- Extract and gap-correct the root sequence (`NODE_0000000`)
- Output a finalized FASTA and updated metadata

### Output Files

- `results/ancestral_sequences.fasta`: All inferred ancestral sequences (from Augur)
- `resources/inferred_root.fasta`: Gap-corrected ancestral sequence


## Updating the Ancestral Root or Metadata

If you update your dataset or rerun the workflow:
- **Be sure to also update the corresponding [metadata file](../resources/static_inferred_metadata.tsv)**  to match the new inferred sequence.
- Record the date of the latest update (see below).

**Latest static inferred sequence generated:**  
📅 2025-10-17


---
## FAQ

**Q: Why use a static inferred ancestor for Nextclade?**  
A: Using an inferred ancestral sequence representative of your dataset improves mutation calling and clade assignment, especially for highly variable viruses like EV-D68.

**Q: How is the inferred sequence cleaned?**  
A: Gaps or ambiguous bases in the inferred root are replaced positionally with nucleotides from the reference sequence using [`fix_root_gaps.py`](scripts/fix_root_gaps.py).

**Q: Should I also use this with other enterovirus datasets?**  
A: Yes! This approach is recommended for other enterovirus types (EV-A71, CVA16, etc.) — a Nextclade Dataset Template is available under [enterovirus-phylo/dataset-template-inferred-root](https://github.com/enterovirus-phylo/dataset-template-inferred-root).



## Author & Contact

- Maintainers: Nadia Neuner-Jehle, Alejandra González-Sánchez and Emma B. Hodcroft ([hodcroftlab](https://github.com/hodcroftlab))
- For questions or suggestions, please [open an issue](https://github.com/hodcroftlab/nextclade_d68/issues) or email: eve-group[at]swisstph.ch



