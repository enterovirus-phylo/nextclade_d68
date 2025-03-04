# Nextclade setup for D68
## First steps
Create a folder structure similar to this [example-workflow](https://github.com/nextstrain/nextclade_data/tree/master/docs/example-workflow).
- dataset/
- profiles/
- resources/
- rules/ 
- scripts/ 

```bash
mkdir dataset profiles resources rules scripts results
```

1. **Run `generate_from_genbank.py` Script:**  
   Execute the script (located in `scripts/`) to generate required reference files:
   ```bash
   python3 scripts/generate_from_genbank.py --reference "AY426531.1" --output-dir dataset/
   ```

   During execution, you may be asked to provide CDS annotations. You can use the following codes to specify the CDS automatically:
   - `[0]`
   - `[product]` or `[leave empty for manual choice]` to select proteins.
   - `[2]`.
   
   The generated files `reference.fasta` and `genome_annotation.gff3` will be saved in the dataset sub-directory.

2. **Update the files and attributes in `pathogen.json`**
   The file names should match the files in the directory. Updated the atrributes. 
   If needed, adapted the Quality control (qc) configuration. If not provided, Nextclade does not do any QC checks. For more details, please check [here](https://docs.nextstrain.org/projects/nextclade/en/latest/user/input-files/05-pathogen-config.html).

3. **Create/ update `reference.gb`**
   Copy the reference.gb file over to resources/. Adapt the protein names to your liking.


4. **Run the `Snakefile`**
   Lines 1-13: adapt the paths and parameters.
      Create/ copy over all the files needed for the augur pipeline to create the reference `tree.json` and run nextclade:
      e.g. `auspice_config.json`




