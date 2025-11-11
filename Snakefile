# Set the parameters
REFERENCE_ACCESSION =   "AY426531"
TAXON_ID =              42789
GENES =                 ["VP4", "VP2", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D"]
ALLOWED_DIVERGENCE =    "1000" # TODO: lower this threshold to exclude outliers
MIN_DATE =              "1990-01-01"
MIN_LENGTH =            "6000" # is 6000 for whole genome build on Nextstrain
MAX_SEQS =              "1200" #TODO: set to 10000 for testing
ROOTING =               "ancestral_sequence"  # mid_point, outgroup, reference, ancestral sequence
ID_FIELD=               "accession" # either accession or strain, used for meta-id-column in augur

# Set the paths
SEQUENCES =             "data/sequences.fasta"
METADATA =              "data/metadata.tsv"

GFF_PATH =              "dataset/genome_annotation.gff3" 
PATHOGEN_JSON =         "dataset/pathogen.json"
README_PATH =           "dataset/README.md"
CHANGELOG_PATH =        "dataset/CHANGELOG.md"
REFERENCE_PATH =        "dataset/reference.fasta"

GENBANK_PATH =          "resources/reference.gbk"
AUSPICE_CONFIG =        "resources/auspice_config.json"
EXCLUDE =               "resources/exclude.txt"
CLADES =                "resources/clades.tsv"
ACCESSION_STRAIN =      "resources/accession_strain.tsv"
INCLUDE_EXAMPLES =      "resources/include_examples.txt"
COLORS =                "resources/colors.tsv"
COLORS_SCHEMES =        "resources/color_schemes.tsv"
INFERRED_ANCESTOR =     "resources/inferred-root.fasta"

FETCH_SEQUENCES = True
STATIC_ANCESTRAL_INFERRENCE = True # whether to use the static inferred ancestral sequence
INFERRENCE_RERUN = False # whether to rerun the inference of the ancestral sequence worfkflow (inferred-root)

INFERRED_SEQ_PATH = "results/sequences_with_ancestral.fasta" if STATIC_ANCESTRAL_INFERRENCE else SEQUENCES
INFERRED_META_PATH = "results/metadata_with_ancestral.tsv" if STATIC_ANCESTRAL_INFERRENCE else "results/metadata.tsv"

include: "scripts/workflow_messages.snkm"

rule all:
    input:
        auspice = "results/auspice.json",
        augur_jsons = "test_out/",
        data = "dataset.zip",
        seqs = "results/example_sequences.fasta",
        json = "out-dataset/pathogen.json",
        **({"root": INFERRED_ANCESTOR} if STATIC_ANCESTRAL_INFERRENCE else {})

rule testing:
    input:
        "testing/EV-D68_fragments.fasta",
        "testing/EV-D68_recombinants.fasta"


if FETCH_SEQUENCES == True:
    rule fetch:
        input:
            dir = "ingest"
        output:
            sequences=SEQUENCES,
            metadata=METADATA
        threads: workflow.cores
        shell:
            """
            cd {input.dir} 
            snakemake --cores {threads} all
            cd ../
            """

rule curate:
    message:
        """
        Cleaning up metadata with augur merge & augur curate
        """
    input:
        meta=METADATA,  # Path to input metadata file
        strains = ACCESSION_STRAIN  # Strain - accession lookup table
    params:
        strain_id_field = ID_FIELD,
    output:
        metadata = "results/metadata.tsv",  # Final output file for publications metadata
    shell:
        """
        augur merge --metadata metadata={input.meta} strains={input.strains}\
            --metadata-id-columns {params.strain_id_field} \
            --output-metadata metadata.tmp
        augur curate normalize-strings \
            --metadata metadata.tmp \
            --id-column {params.strain_id_field} \
            --output-metadata {output.metadata}

        rm metadata.tmp
        """

rule add_reference_to_include:
    """
    Create an include file for augur filter
    """
    input:
        "resources/include.txt",
    output:
        "results/include.txt",
    shell:
        """
        cat {input} >> {output}
        echo "{REFERENCE_ACCESSION}" >> {output}
        echo ancestral_sequence >> {output}
        """

if STATIC_ANCESTRAL_INFERRENCE and INFERRENCE_RERUN:
    rule static_inferrence:
        message:
            """
            Running "inferred-root" snakefile for inference of the ancestral root. 
            This reference will be included in the Nextclade reference tree.
            WARNING: This will overwrite your {output.seq} & {output.meta} files!
            """
        input:
            dir = "inferred-root",
            dataset_path = "dataset",
            meta = rules.curate.output.metadata,
            seq = SEQUENCES,
            meta_ancestral = "resources/static_inferred_metadata.tsv",
            include = "results/include.txt"
        params:
            strain_id_field = ID_FIELD,
        output:
            inref = INFERRED_ANCESTOR,
            seq = INFERRED_SEQ_PATH,
            meta = INFERRED_META_PATH,
        threads: workflow.cores
        shell:
            r"""
            set -euo pipefail

            echo "Cleaning previous results..."
            rm -rf {input.dir}/results/* {input.dir}/resources/inferred-root.fasta

            echo "Running inferred-root workflow..."
            cd {input.dir}
            snakemake --cores {threads} all_sub
            cd - > /dev/null

            echo "Combining sequences with ancestral root..."
            cat {input.seq} {output.inref} > {output.seq}

            echo "Merging metadata..."
            augur merge \
                --metadata metadata={input.meta} ancestral={input.meta_ancestral} \
                --metadata-id-columns {params.strain_id_field} \
                --output-metadata {output.meta}

            echo "Static ancestral inference completed successfully!"
            """

rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering
        """
    input:
        sequences = INFERRED_SEQ_PATH,
    output:
        sequence_index = "results/sequence_index.tsv"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """

rule filter:
    """
    Exclude sequences from before {MIN_DATE} and subsample to {MAX_SEQS} sequences.
    Only take sequences longer than {MIN_LENGTH}
    """
    input:
        sequences = INFERRED_SEQ_PATH,
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = INFERRED_META_PATH,
        include = rules.add_reference_to_include.output,
    output:
        filtered_sequences = "results/filtered_sequences_raw.fasta",
        filtered_metadata = "results/filtered_metadata_raw.tsv",
    params: 
        min_date="" if MIN_DATE == "" else "--min-date " + MIN_DATE,
        min_length="" if MIN_LENGTH == "" else "--min-length " + MIN_LENGTH,
        max_seqs=MAX_SEQS,
        categories = "country year", #TODO: add subsampling per category?
        strain_id_field = ID_FIELD,
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            {params.min_length} \
            {params.min_date} \
            --include {input.include} \
            --subsample-max-sequences {params.max_seqs} \
            --output-sequences {output.filtered_sequences} \
            --output-metadata {output.filtered_metadata}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference} using Nextclade3.
        """
    input:
        sequences = rules.filter.output.filtered_sequences,
        reference = REFERENCE_PATH,
        annotation = GFF_PATH,
    output:
        alignment = "results/aligned.fasta",
        tsv = "results/nextclade.tsv",
    params:
        translation_template = lambda w: "results/translations/cds_{cds}.translation.fasta",
        penalty_gap_extend = 3, #make longer gaps more costly - default is 0
        penalty_gap_open = 16,  #make gaps more expensive relative to mismatches - default is 13
        penalty_gap_open_in_frame = 7, #make gaps more expensive relative to mismatches - default is 7
        penalty_gap_open_out_of_frame = 30, #make out of frame gaps more expensive - default is 8 # prev was 19
        kmer_length = 6, #reduce to find more matches - default is 10
        kmer_distance = 15, #reduce to try more seeds - default is 50
        min_match_length = 25, #reduce to keep more seeds - default is 40
        allowed_mismatches = 10, #increase to keep more seeds - default is 8
        min_length = 100, # min_length - default is 100
        gap_alignment_side = "right", # default is left
    shell:
        """
        nextclade3 run \
        -j {threads} \
        {input.sequences} \
        --input-ref {input.reference} \
        --input-annotation {input.annotation} \
        --alignment-preset high-diversity \
        --penalty-gap-open {params.penalty_gap_open} \
        --penalty-gap-extend {params.penalty_gap_extend} \
        --penalty-gap-open-in-frame {params.penalty_gap_open_in_frame} \
        --penalty-gap-open-out-of-frame {params.penalty_gap_open_out_of_frame} \
        --gap-alignment-side {params.gap_alignment_side} \
        --kmer-length {params.kmer_length} \
        --kmer-distance {params.kmer_distance} \
        --min-match-length {params.min_match_length} \
        --allowed-mismatches {params.allowed_mismatches} \
        --min-length {params.min_length} \
        --max-alignment-attempts 5 \
        --include-reference false \
        --output-tsv {output.tsv} \
        --output-translations {params.translation_template} \
        --output-fasta {output.alignment} 
        """
        #         --penalty-mismatch 



rule get_outliers:
    """
    Automatically identify sequences with >{ALLOWED_DIVERGENCE} substitutions
    (likely to be sequencing errors or low quality/misannotated sequences) and put them in outliers.txt
    """
    input:
        nextclade = rules.align.output.tsv,
    output:
        outliers = "results/outliers.txt",
        tmp = "tmp/outliers.txt",
    params:
        allowed_divergence = lambda w: ALLOWED_DIVERGENCE,
    shell:
        """
        tsv-filter -H -v --is-numeric totalSubstitutions {input.nextclade} \
        > {output.tmp}
        tsv-filter -H \
            --is-numeric totalSubstitutions \
            --gt totalSubstitutions:{params.allowed_divergence} \
            {input.nextclade} \
        | tail -n +2 >> {output.tmp}
        cat {output.tmp} \
        | tsv-select -H -f seqName \
        | tail -n +2 > {output.outliers}
        """


rule exclude:
    """
    Rule to allow for manual and automatic exclusion of sequences
    without triggering a new subsampling that could
    surface new bad sequences resulting in an infinite loop
    """
    input:
        sequences = rules.align.output.alignment,
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = INFERRED_META_PATH,
        exclude = EXCLUDE,
        outliers = rules.get_outliers.output.outliers,
        example = INCLUDE_EXAMPLES,

    params:
        strain_id_field = ID_FIELD,
    output:
        filtered_sequences = "results/filtered_aligned.fasta",
        filtered_metadata = "results/filtered_metadata.tsv",
        strains = "results/tree_strains.txt",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --exclude {input.exclude} {input.outliers} {input.example} \
            --output-sequences {output.filtered_sequences} \
            --output-metadata {output.filtered_metadata} \
            --output-strains {output.strains}
        """


rule tree:
    message:
        """
        Creating a maximum likelihood tree
        """
    input:
        alignment = rules.exclude.output.filtered_sequences,
    output:
        tree = "results/tree_raw.nwk",
    threads: 9
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --nthreads {threads}\
            --output {output.tree} \
        """

rule refine:
    input:
        tree=rules.tree.output.tree,
        alignment=rules.exclude.output.filtered_sequences,
    output:
        tree="results/tree.nwk",
        node_data="results/branch_lengths.json",
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --root {ROOTING} \
            --keep-polytomies \
            --divergence-unit mutations-per-site \
            --output-node-data {output.node_data} \
            --output-tree {output.tree}
        """


rule ancestral:
    input:
        tree=rules.refine.output.tree,
        alignment=rules.exclude.output.filtered_sequences,
        annotation=GENBANK_PATH,
    output:
        node_data="results/muts.json",
        ancestral_sequences="results/ancestral_sequences.fasta",
    params:
        translation_template=r"results/translations/cds_%GENE.translation.fasta",
        output_translation_template=r"results/translations/cds_%GENE.ancestral.fasta",
        genes=" ".join(GENES),
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --annotation {input.annotation} \
            --root-sequence {input.annotation} \
            --genes {params.genes} \
            --translations {params.translation_template} \
            --output-node-data {output.node_data} \
            --output-translations {params.output_translation_template}\
            --output-sequences {output.ancestral_sequences}
        """


rule clades:
    input:
        tree = rules.refine.output.tree,
        mutations = rules.ancestral.output.node_data,
        clades = CLADES
    output:
        json = "results/clades.json",
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.mutations} \
            --clades {input.clades} \
            --output-node-data {output.json}
        """

rule get_dates:
    """Create ordering for color assignment"""
    input:
        metadata = rules.exclude.output.filtered_metadata
    output:
        ordering = "results/color_ordering.tsv"
    run:
        import pandas as pd
        column = "date"
        meta = pd.read_csv(input.metadata, delimiter='\t')

        if column not in meta.columns:
            print(f"The column '{column}' does not exist in the file.")
            sys.exit(1)

        deflist = meta[column].dropna().tolist()
        # Store unique values (ordered)
        deflist = sorted(set(deflist))
        if "XXXX-XX-XX" in deflist:
            deflist.remove("XXXX-XX-XX")

        result_df = pd.DataFrame({
            'column': ['date'] * len(deflist),
            'value': deflist
        })

        result_df.to_csv(output.ordering, sep='\t', index=False, header=False)
        
rule epitopes:
    input:
        anc_seqs = rules.ancestral.output.node_data,
        tree = rules.refine.output.tree,
    output:
        node_data = "results/epitopes.json"
    params:
        translation = "results/translations/cds_VP1.ancestral.fasta",
        # Original VP1-relative epitope positions (amino acids)
        epitopes = {
            'BC': [90, 92, 95, 97, 98, 103], 
            'DE': [140, 141, 142, 143, 144, 145, 146, 147, 148],
            'CTERM': [280, 283, 284, 288, 290, 297, 299, 301, 304, 305, 306, 308]
        },
        min_count = 4,  # Minimum count to keep a specific epitope sequence
    run:
        import json
        from collections import defaultdict
        from augur.translate import safe_translate
        from Bio import Phylo
        from Bio import SeqIO
        import ipdb

        manyXList = ["XXXXXXXXXXXX", "KEXXXXXXXXXX", "KERANXXXXXXX", "KERXXXXXXXXX", "KERAXXXXXXXX"]
        
        # with open(input.anc_seqs) as fh:
        #     anc = json.load(fh)["nodes"]

        # Read translation files
        vp1_anc = SeqIO.to_dict(SeqIO.parse(params.translation, "fasta"))

        T = Phylo.read(input.tree, 'newick')
        for node in T.find_clades(order='preorder'):
            for child in node:
                child.parent = node

        nodes = {}
        epitope_counts = {epi: defaultdict(int) for epi in params.epitopes}
        

        for node in T.find_clades(order='preorder'):
            # ipdb.set_trace()
            n = node.name
            aa = vp1_anc[n].seq
            nodes[n] = {}
            for epi, pos in params.epitopes.items():
                pos = [p - 1 for p in pos]  # Convert to 0-based indexing
                nodes[n][epi] = "".join([aa[p] for p in pos])
                if epi == 'CTERM':
                    if nodes[n]['CTERM'] in manyXList:
                        nodes[n]['CTERM'] = "many x"
                    elif 'X' in nodes[n]['CTERM']:
                        nodes[n]['CTERM'] = nodes[node.parent.name]['CTERM']
                if not n.startswith('NODE_'):
                    epitope_counts[epi][nodes[n][epi]] += 1

        for node in nodes:
            for epi, seq in nodes[node].items():
                min_count2 = params.min_count if epi != "CTERM" else 6
                if epi == "CTERM" and seq in manyXList:
                    nodes[node][epi] = 'many X'
                elif epitope_counts[epi][seq] < min_count2:
                    nodes[node][epi] = 'other'

        with open(output.node_data, 'w') as fh:
            json.dump({"epitopes": params.epitopes, "nodes": nodes}, fh)

rule colors:
    """Assign colors based on ordering"""
    input:
        ordering=rules.get_dates.output.ordering,
        color_schemes=COLORS_SCHEMES,
        colors=COLORS,
    output:
        colors="results/colors_dates.tsv",
        final_colors="results/final_colors.tsv"
    shell:
        """
        python3 scripts/assign-colors.py \
            --ordering {input.ordering} \
            --color-schemes {input.color_schemes} \
            --output {output.colors}

        echo -e '\ndate\tXXXX-XX-XX\t#a6acaf' >> {output.colors}

        cat {output.colors} {input.colors} >> {output.final_colors}
        """

rule export: 
    input:
        tree = rules.refine.output.tree,
        metadata = rules.exclude.output.filtered_metadata,
        mutations = rules.ancestral.output.node_data,
        branch_lengths = rules.refine.output.node_data,
        clades = rules.clades.output.json, # dummy_clades if not set yet
        auspice_config = AUSPICE_CONFIG,
        colors = rules.colors.output.final_colors,
        epitopes = rules.epitopes.output.node_data,
    params:
        strain_id_field = ID_FIELD,
        fields="region country date",
    output:
        auspice = "results/auspice.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --auspice-config {input.auspice_config} \
            --color-by-metadata {params.fields} \
            --colors {input.colors} \
            --node-data {input.mutations} {input.branch_lengths} {input.clades} {input.epitopes}\
            --output {output.auspice}
        """

rule extract_clades_tsv:
    input:
        json=rules.clades.output.json,
    output:
        tsv = "results/clades_metadata.tsv"
    run:
        import json
        import csv

        with open(input.json) as f:
            data = json.load(f)

        nodes = data.get("nodes", {})

        with open(output.tsv, "w", newline="") as out_f:
            writer = csv.writer(out_f, delimiter="\t")
            writer.writerow(["accession", "clade"])

            for accession, values in nodes.items():
                clade = values.get("clade_membership", None)
                if clade:
                    writer.writerow([accession, clade])


rule subsample_example_sequences:
    input:
        all_sequences = INFERRED_SEQ_PATH,
        metadata = INFERRED_META_PATH,
        exclude = EXCLUDE,
        outliers = rules.get_outliers.output.outliers,
        incl_examples = INCLUDE_EXAMPLES,
        clades =  rules.extract_clades_tsv.output.tsv,
        tree_strains = "results/tree_strains.txt",  # strains in the tree
    output:
        example_sequences = "results/example_sequences.fasta",
    params:
        strain_id_field = ID_FIELD,
    shell:
        """
        augur merge \
            --metadata metadata={input.metadata} clades={input.clades} \
            --metadata-id-columns {params.strain_id_field} \
            --output-metadata metadata.tmp
        augur filter \
            --sequences {input.all_sequences} \
            --metadata metadata.tmp \
            --metadata-id-columns {params.strain_id_field} \
            --min-length 4000 \
            --include {input.incl_examples} \
            --exclude {input.exclude} {input.outliers} \
            --exclude-ambiguous-dates-by year \
            --min-date 2012 --group-by clade \
            --subsample-max-sequences 15  \
            --probabilistic-sampling \
            --output-sequences {output.example_sequences}
        rm metadata.tmp
        """
        # seqkit grep -v -f {input.tree_strains} {input.all_sequences} \
        # | seqkit sample -n 100 -s 41 > {output.example_sequences}

rule assemble_dataset:
    input:
        tree = rules.export.output.auspice,
        reference = REFERENCE_PATH,
        annotation = GFF_PATH,
        sequences = rules.subsample_example_sequences.output.example_sequences,
        pathogen = PATHOGEN_JSON,
        readme = README_PATH,
        changelog = CHANGELOG_PATH,
    params:
        pathogen = "out-dataset/pathogen.json",
    output:
        tree = "out-dataset/tree.json",
        reference = "out-dataset/reference.fasta",
        annotation = "out-dataset/genome_annotation.gff3",
        sequences = "out-dataset/sequences.fasta",
        readme = "out-dataset/README.md",
        changelog = "out-dataset/CHANGELOG.md",
        dataset_zip = "dataset.zip",
    shell:
        """
        cp {input.tree} {output.tree}
        cp {input.reference} {output.reference}
        cp {input.annotation} {output.annotation}
        cp {input.sequences} {output.sequences}
        cp {input.pathogen} {params.pathogen}
        cp {input.readme} {output.readme}
        cp {input.changelog} {output.changelog}
        zip -rj dataset.zip  out-dataset/*
        """


rule test:
    input:
        dataset = rules.assemble_dataset.output.dataset_zip,
        sequences = rules.assemble_dataset.output.sequences,
    output:
        output = directory("test_out"),
    shell:
        """
        nextclade3 run \
            --input-dataset {input.dataset} \
            --output-all {output.output} \
            {input.sequences}
        """

rule mutLabels:
    input:
        d = directory("test_out"),
        json = PATHOGEN_JSON,
    params:
        "results/virus_properties.json",
    output:
        "out-dataset/pathogen.json"
    shell:
        """
        python3 scripts/generate_virus_properties.py
        jq --slurpfile v {params} \
           '.mutLabels.nucMutLabelMap = $v[0].nucMutLabelMap |
            .mutLabels.nucMutLabelMapReverse = $v[0].nucMutLabelMapReverse' \
           {input.json} > {output}
        zip -rj dataset.zip  out-dataset/*
        """


rule fragment_testing:
    input:
        nextstrain = "testing/nextstrain_d68_vp1.tsv",
        sequences = "results/aligned.fasta",
    output:
        fragments = "testing/EV-D68_fragments.fasta"
    params:
        length = range(100, 3000, 100),  # lengths from 200 to 3000
        gene = ["VP1", "3D"]  # genes to sample from; atm only VP1 and 3D supported
    run:
        import os
        import random
        from Bio import SeqIO
        import pandas as pd

        # Read all sequences from the input file
        records = list(SeqIO.parse(input.sequences, "fasta"))
        os.makedirs(os.path.dirname(output.fragments), exist_ok=True)

        # filter records in nextstrain file
        ns_ids = list(pd.read_csv(input.nextstrain).accession)
        records = [r for r in records if r.id in ns_ids]

        with open(output.fragments, "w") as out_handle:
            for length in params.length:
                record = random.choice(records)
                seq_len = len(record.seq)
                if "VP1" in params.gene or "3D" in params.gene:
                    if "VP1" in params.gene: 
                        seq1 = record.seq[2389:3315]
                        l = len(seq1) - seq1.count("-") - seq1.count("N")
                        if l > length:
                            s = random.randint(0, l - length)
                            seq1 = seq1[s:s+length]
                            header = f"{record.id}_partial_{length}_VP1"
                            out_handle.write(f">{header}\n{seq1}\n")
                    if "3D" in params.gene:
                        seq2 = record.seq[5926:7296]
                        l = len(seq2)
                        if l > length:
                            s = random.randint(0, l - length)
                            seq2 = seq2[s:s+length]
                            header = f"{record.id}_partial_{length}_3D"
                            out_handle.write(f">{header}\n{seq2}\n")
                else: 
                    print(f"Gene {params.gene} not recognized.")
                        
                while seq_len < length:
                    record = random.choice(records)
                    seq_len = len(record.seq)
                start = random.randint(0, seq_len - length)
                fragment_seq = record.seq[start:start+length]
                header = f"{record.id}_partial_{length}"
                out_handle.write(f">{header}\n{fragment_seq}\n")


rule recombinant_testing:
    input:
        sequences = SEQUENCES,
        nextstrain = "testing/nextstrain_d68_vp1.tsv",
        clades = "results/clades_metadata.tsv",
        evD_seq = "testing/EV-D_sequence.fasta"
    output:
        recombinants = "testing/EV-D68_recombinants.fasta"
    params:
        inter_recombinants = 10,
        intra_recombinants = 10,
        min_length = 3500,
    run:
        import random
        from Bio import SeqIO
        import pandas as pd

        def eligible(records, ml):
            return [r for r in records if len(r) >= ml]

        # Load sequences and filter by Nextstrain IDs & min_length
        seqs = list(SeqIO.parse(input.sequences, "fasta"))
        ns_ids = list(pd.read_csv(input.nextstrain, sep="\t").accession)
        
        seqs = eligible([r for r in seqs if r.id in ns_ids], params.min_length)

        # Map clade assignments
        clade_map = pd.read_csv(input.clades, sep="\t").set_index("accession")["clade"].to_dict()
        clade2seqs = {}
        for r in seqs:
            clade = clade_map.get(r.id, "NA")
            clade2seqs.setdefault(clade, []).append(r)
        clades = [c for c in clade2seqs if c != "NA" and len(clade2seqs[c]) > 0]

        # EV-D sequences for intertypic recombination
        evd = eligible(list(SeqIO.parse(input.evD_seq, "fasta")), params.min_length)

        with open(output.recombinants, "w") as out:
            # Intra-typic: between clades
            for i in range(params.intra_recombinants):
                c1, c2 = random.sample(clades, 2)
                p1, p2 = random.choice(clade2seqs[c1]), random.choice(clade2seqs[c2])
                minlen = min(len(p1.seq), len(p2.seq))
                if minlen < params.min_length: continue
                x = random.randint(1, minlen-1)
                out.write(f">intra_{p1.id}_{c1}_{x}_{p2.id}_{c2}\n{p1.seq[:x]}{p2.seq[x:]}\n")

            # Inter-typic: D68 x EV-D
            for i in range(params.inter_recombinants):
                p1 = random.choice(seqs)
                p2 = random.choice(evd)
                minlen = min(len(p1.seq), len(p2.seq))
                if minlen < params.min_length: continue
                x = random.randint(1, minlen-1)
                out.write(f">inter_{p1.id}_D68_{x}_{p2.id}_D\n{p1.seq[:x]}{p2.seq[x:]}\n")


rule clean:
    shell:
        """
        rm ingest/data/* data/*
        rm -r results out-dataset test_out dataset.zip tmp
        rm resources/inferred-root.fasta inferred-root/resources/inferred-root.fasta
        rm -r inferred-root/results/*
        """