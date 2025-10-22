# Nextclade Workflow for Enterovirus D68

This repository contains a robust, reproducible workflow for building a custom [Nextclade](https://github.com/nextstrain/nextclade) dataset for Enterovirus D68 (EV-D68). It enables you to generate reference and annotation files, download and process sequence data, infer an ancestral root, and create all files needed for Nextclade analyses and visualization.

---

## Folder Structure

Follow the [Nextclade example workflow](https://github.com/nextstrain/nextclade_data/tree/master/docs/example-workflow) or use the structure below:

```bash
mkdir -p dataset data ingest resources results scripts
```

---

## Workflow Overview

This workflow is composed of several modular steps:

1. **Reference Generation**  
   Extracts relevant reference and annotation files from GenBank.
2. **Dataset Ingest**  
   Downloads and processes sequences and metadata from NCBI Virus.
3. **Phylogenetic Root Inference (optional)**  
   Infers a dataset-specific ancestral root sequence to use as a reference in Nextclade, improving mutation and clade assignments.
4. **Augur Phylogenetics & Nextclade Preparation**  
   Builds trees, prepares multiple sequence alignments, and generates all files required for Nextclade and Auspice.
5. **Visualization & Analysis**  
   Enables both command-line and web-based Nextclade analyses, including local dataset hosting.

---
## Setup Instructions

### 1. Generate Reference Files

Run the script to extract the reference FASTA and genome annotation from GenBank:

```bash
python3 scripts/generate_from_genbank.py --reference "AY426531.1" --output-dir dataset/
```

During the script execution, follow the prompts for CDS annotation selection.
   - `[0]`
   - `[product]` or `[leave empty for manual choice]` to select proteins.
   - `[2]`.

**Outputs:**
- `dataset/reference.fasta`
- `dataset/genome_annotation.gff3`

---

### 2. Configure `pathogen.json`

Edit `pathogen.json` to:
- Reference your generated files (`reference.fasta`, `genome_annotation.gff3`)
- Update metadata and QC settings as needed  
> [!WARNING]  
> If QC is not set, Nextclade will skip quality checks.

See the [Nextclade pathogen config documentation](https://docs.nextstrain.org/projects/nextclade/en/latest/user/input-files/05-pathogen-config.html) for details.

---

### 3. Prepare GenBank Reference

Copy your GenBank file to `resources/reference.gb`.  
Edit protein names and features as necessary for your use case.

---
### 4. Update the `Snakefile`

- Adjust the first 29 lines for correct file paths and parameters.
- Ensure required files are available:
  - `data/sequences.fasta`
  - `data/metadata.tsv`
  - `resources/auspice_config.json`

Sequences and metadata can be downloaded automatically via the ingest process (see below).

---

## Subprocesses

### Ingest

Automates downloading of EV-D68 sequences and metadata from NCBI Virus.  
See [ingest/README.md](ingest/README.md) for specifics.

**Required packages:**  
`csvtk, nextclade, tsv-utils, seqkit, zip, unzip, entrez-direct, ncbi-datasets-cli` (installable via conda-forge/bioconda)

---

### Inferred Root (Optional but Recommended)

The `inferred-root/` directory contains a self-contained pipeline to infer a dataset-specific ancestral root, which can be used as a reference for Nextclade. This enhances mutation and clade call accuracy for EV-D68 datasets.

- **See:** [`inferred-root/README.md`](inferred-root/README.md) for details.
- To enable, set `ANCESTRAL_ROOT_INFERRENCE = True` in your config and run with  
  `--config root_inference_confirmed=true`.
- Without confirmation, the workflow will halt and display an opt-in message.

> [!NOTE]  
> To skip the inferred root step, leave `ANCESTRAL_ROOT_INFERRENCE = False`.

### **Template for other enteroviruses:**  
If you want to apply this approach to other enterovirus types (e.g., EV-A71, CVA16), a [Nextclade Dataset Template for Inferred Root](https://github.com/enterovirus-phylo/dataset-template-inferred-root) is available and recommended for reuse.

---

## Running the Workflow

To generate the Auspice JSON and a Nextclade example dataset:

```bash
snakemake --cores 9 all --config root_inference_confirmed=true
```

This will:
- Build the reference tree and produce the Nextclade dataset in `dataset/`
- Run Nextclade on the example sequences in `out-dataset/sequences.fasta`
- Output results to `test_out/` (alignment, translations, summary TSV)


### Labeling Mutations of Interest
To label mutations of interest, execute the `mutLabels` rule as a standalone instance. They will be added to the `out-dataset/pathogen.json` file.

---
## Visualizing Your Custom Nextclade Dataset

To use the dataset in Nextclade Web, serve it locally:

```bash
serve --cors out-dataset -l 3000
```

Then open:

```
https://master.clades.nextstrain.org/?dataset-url=http://localhost:3000
```

- Click "Load example", then "Run"
- You may want to reduce "Max. nucleotide markers" to 500 under "Settings" → "Sequence view" to optimize performance

---

## Author & Contact
- Maintainers: Nadia Neuner-Jehle, Alejandra Gonzalez Sanchez and Emma B. Hodcroft ([hodcroftlab](https://github.com/hodcroftlab))
- For questions or suggestions, please [open an issue](https://github.com/hodcroftlab/nextclade_d68/issues) or email: eve-group[at]swisstph.ch

## Troubleshooting and Further Help

- For issues, see the [official Nextclade documentation](https://docs.nextstrain.org/projects/nextclade/en/stable/index.html#) or [open an issue](https://github.com/hodcroftlab/nextclade_d68/issues).
- For details on the inferred root workflow, see [`inferred-root/README.md`](inferred-root/README.md).
- For adapting to other enteroviruses, see the [dataset-template-inferred-root](https://github.com/enterovirus-phylo/dataset-template-inferred-root).

---

This guide provides a structured, scalable approach to building and using high-quality Nextclade datasets for EV-D68 — and can be adapted for other enterovirus types as well.

## Task List
- [x] Integrate ancestral inferred-root into workflow (https://github.com/hodcroftlab/nextclade_d68/pull/2)
- [x] Validate clade assignment of fragmented sequences in Nextclade (`testing/`)
- [x] Ensure novel recombinants get assigned to the root (issue https://github.com/hodcroftlab/nextclade_d68/issues/3) -> recombinant feature in testing; QC label
- [ ] Review and validate EV-D68 nomenclature, including robustness with recombinant sequences
- [ ] Integrate epitope mutation information as tree coloring and/or display in the Nextclade results table
