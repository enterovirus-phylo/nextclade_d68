# Inferred Ancestral Sequence for Enterovirus D68

This pipeline generates a static inferred ancestral sequence—a stable custom reference—for Enterovirus D68 (EV-D68), which is integrated into the main Nextclade dataset build.

High-quality clade and mutation assignment depend on a representative reference. This workflow infers an ancestral sequence and inserts it early into the dataset, allowing Nextclade to use the inferred ancestor as a consistent reference sequence.

## Key Concepts

- The top-level Snakefile can run a dedicated sub-workflow (`inferred-root`) to infer a static ancestral/root sequence.
- When `STATIC_ANCESTRAL_INFERENCE = True`, the pipeline:
  1. Runs the `inferred-root` sub-workflow (`inferred-root/Snakefile`) to infer a root sequence and write it as `resources/inferred-root.fasta`.
  2. Concatenates the inferred root FASTA with your input sequences to produce `results/sequences_with_ancestral.fasta`.
  3. Merges the project metadata with a small metadata TSV for the ancestral sequence to produce `results/metadata_with_ancestral.tsv`.
  4. Continues the standard downstream pipeline (`index → filter → align → tree → refine → ancestral → clades → export → assemble dataset`), now operating on the sequence set that includes the static inferred ancestor.

Because the inferred ancestral sequence is injected early, all downstream steps consistently include the static inferred ancestor.

## Files and Configuration

- **Parameters:** Adjust in the first 8 lines of the Snakefile.
- **Outgroup sequences:** Update the outgroup list and add FASTA files (one sequence per file) in `resources/outgroup/`, named by their accession IDs (e.g., `AY426531.1.fasta`).

## Important Snakefile Options

- `STATIC_ANCESTRAL_INFERENCE`: Enable or disable static ancestral insertion.
- `INFERRENCE_RERUN`: Trigger rerunning the inferred-root workflow.
- `INFERRED_ANCESTOR`: Path to the inferred root FASTA (default: `resources/inferred-root.fasta`).
- `ID_FIELD`: Metadata ID column used during augur merge (default: `"accession"`).

## Inputs

- Your input sequences and metadata.
- `include.txt` — list of sequences to include.

## Key Outputs

- `resources/inferred-root.fasta` — inferred root sequence produced by the workflow.
- `results/sequences_with_ancestral.fasta` — input sequences plus the ancestral sequence.
- `results/metadata_with_ancestral.tsv` — merged metadata including the ancestral sequence.

---
## Requirements 

- Snakemake
- Augur (nextstrain/augur)
- Nextclade v3 (nextclade3)
- jq (for dataset mutLabel merging)
- Standard Python packages used by scripts in the repo (Biopython, pandas, etc.)

## How the Static Inferred-Root Workflow Works

1. **Subsample & Align:**  
   A representative subset (up to `MAX_SEQS`) is selected via `augur filter`. Outgroup sequences from `resources/outgroup/` are always included. The selected sequences are aligned using MAFFT.

2. **Tree Construction:**  
   A maximum-likelihood tree is inferred with `augur tree`.

3. **Outgroup Rooting & Ancestral Sequence Extraction:**  
   The script [`pick_ancestral_sequence.py`](scripts/pick_ancestral_sequence.py) reroots the tree on the provided outgroup(s), identifies the MRCA of the ingroup (EV-D68 sequences), and extracts the ancestral sequence at this node. It replaces gaps with reference nucleotides to ensure a contiguous, biologically plausible sequence.

4. **Post-processing & Export:**  
   The cleaned ancestral FASTA is saved (`resources/inferred-root.fasta`) along with a matching metadata row (`resources/static_inferred_metadata.tsv`). These are injected back into the main pipeline to ensure consistent inclusion of the static inferred ancestor.

---


## Quick run
Run the inferred-root subworkflow:

```bash
snakemake --cores 9 all_sub
```


---
## Tips and Cautions
- Ensure the metadata file [static_inferred_metadata.tsv](../resources/static_inferred_metadata.tsv) is provided and kept up-to-date.
- The inferred root sequence header must match the metadata ID exactly to guarantee correct merging and labeling.

- After running static inference, verify `results/sequences_with_ancestral.fasta` and `results/metadata_with_ancestral.tsv` to confirm the ancestral sequence has been added and correctly labeled.

## Additional Resources

- A Nextclade Dataset Template using this inferred-root approach is available here: [enterovirus-phylo/dataset-template-inferred-root](https://github.com/enterovirus-phylo/dataset-template-inferred-root).

## Author & Contact

- Maintainers: Nadia Neuner-Jehle, Alejandra González-Sánchez and Emma B. Hodcroft ([eve-lab.org](https://eve-lab.org/))
- For questions or suggestions, please [open an issue](https://github.com/enterovirus-phylo/nextclade_d68/issues) or email: eve-group[at]swisstph.ch



