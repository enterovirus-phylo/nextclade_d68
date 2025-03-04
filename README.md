# Nextclade Setup for D68

## Folder Structure
First, create the necessary folder structure as shown in the [example workflow](https://github.com/nextstrain/nextclade_data/tree/master/docs/example-workflow):

```
dataset/
profiles/
resources/
rules/
scripts/
results/
```

You can create these directories using the following command:
```bash
mkdir dataset profiles resources rules scripts results
```

---
## Steps to Set Up The Workflow

### 1. Run `generate_from_genbank.py`
This script (located in `scripts/`) generates reference files from GenBank.

Run the following command:
```bash
python3 scripts/generate_from_genbank.py --reference "AY426531.1" --output-dir dataset/
```

During execution, you may be asked to provide CDS annotations. You can use the following codes to specify the CDS automatically:
   - `[0]`
   - `[product]` or `[leave empty for manual choice]` to select proteins.
   - `[2]`.

The script will generate:
- `dataset/reference.fasta`
- `dataset/genome_annotation.gff3`

---
### 2. Update `pathogen.json`
Modify `pathogen.json` to:
- Ensure file names match the generated reference files.
- Update attributes as needed.
- Adjust the Quality Control (QC) settings if necessary. If QC is not configured, Nextclade will not perform any checks.

For more details on configuration, refer to the [Nextclade documentation](https://docs.nextstrain.org/projects/nextclade/en/latest/user/input-files/05-pathogen-config.html).

---
### 3. Prepare `reference.gb`
- Copy the `reference.gb` file into the `resources/` directory.
- Modify protein names as needed to match your requirements.

---
### 4. Update the `Snakefile`
- Modify lines 1-13 to adjust paths and parameters.
- Ensure all necessary files for the Augur pipeline are present, including:
  - `tree.json`
  - `auspice_config.json`
- These files are essential for building the reference tree and running Nextclade.

---
This guide provides a structured workflow for setting up Nextclade for D68. If you encounter issues, refer to the official documentation or seek support from the Nextstrain community.