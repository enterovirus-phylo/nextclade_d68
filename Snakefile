REFERENCE_ACCESSION = "AY426531.1"
TAXON_ID = 42789
GENES = ["VP4", "VP2", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D"]
ALLOWED_DIVERGENCE = "500"
MIN_DATE = "1960-01-01"
MIN_LENGTH = "6500"
MAX_SEQS = "100"
GFF_PATH = "../dataset/genome_annotation.gff3"
GENBANK_PATH = "resources/ev_d68_reference_genome.gb"
REFERENCE_PATH = "../dataset/reference.fasta"
README_PATH = "../dataset/README.md"
CHANGELOG_PATH = "../dataset/CHANGELOG.md"
ROOTING = "mid_point"  # alternative root using outgroup, e.g. the reference "AY426531.1"
AUSPICE_CONFIG = "resources/auspice_config.json"
EXCLUDE = "resources/exclude.txt"
METADATA = "data/metadata.tsv"

FETCH_SEQUENCES = False

if FETCH_SEQUENCES == True:

    include: "rules/fetch_from_ncbi.smk"


rule add_reference_to_include:
    """
    Create an include file for augur filter
    """
    input:
        include="resources/include.txt",
    output:
        "results/include.txt",
    shell:
        """
        echo "{REFERENCE_ACCESSION}" >> results/include.txt
        """


rule filter:
    """
    Exclude sequences from before {MIN_DATE} and subsample to {MAX_SEQS} sequences.
    Only take sequences longer than {MIN_LENGTH}
    """
    input:
        sequences="data/sequences.fasta",
        metadata=METADATA,
        include=rules.add_reference_to_include.output,
    output:
        filtered_sequences="results/filtered_sequences_raw.fasta",
        filtered_metadata="results/filtered_metadata_raw.tsv",
    params:
        min_date="" if MIN_DATE == "" else "--min-date " + MIN_DATE,
        min_length="" if MIN_LENGTH == "" else "--min-length " + MIN_LENGTH,
        max_seqs=MAX_SEQS,
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            {params.min_length} \
            {params.min_date} \
            --include {input.include} \
            --subsample-max-sequences {params.max_seqs} \
            --output {output.filtered_sequences} \
            --output-metadata {output.filtered_metadata}
        """


rule align:
    input:
        sequences=rules.filter.output.filtered_sequences,
        reference=REFERENCE_PATH,
        annotation=GFF_PATH,
    output:
        alignment="results/aligned.fasta",
        tsv="results/nextclade.tsv",
    params:
        translation_template=lambda w: "results/translations/cds_{cds}.translation.fasta",
    shell:
        """
        nextclade3 run \
            {input.sequences} \
            --input-ref {input.reference} \
            --input-annotation {input.annotation} \
            --output-translations {params.translation_template} \
            --output-tsv {output.tsv} \
            --output-fasta {output.alignment}
        """


rule get_outliers:
    """
    Automatically identify sequences with >{ALLOWED_DIVERGENCE} substitutions
    (likely to be sequencing errors or low quality/misannotated sequences) and put them in outliers.txt
    """
    input:
        nextclade=rules.align.output.tsv,
    output:
        outliers="results/outliers.txt",
        tmp="tmp/outliers.txt",
    params:
        allowed_divergence=lambda w: ALLOWED_DIVERGENCE,
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
        sequences=rules.align.output.alignment,
        metadata=METADATA,
        exclude=EXCLUDE,
        outliers=rules.get_outliers.output.outliers,
    output:
        filtered_sequences="results/filtered_aligned.fasta",
        filtered_metadata="results/filtered_metadata.tsv",
        strains="results/tree_strains.txt",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} {input.outliers} \
            --output {output.filtered_sequences} \
            --output-metadata {output.filtered_metadata} \
            --output-strains {output.strains}
        """


rule tree:
    input:
        alignment=rules.exclude.output.filtered_sequences,
    output:
        tree="results/tree_raw.nwk",
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
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
            --divergence-units mutations \
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
            --output-translations {params.output_translation_template}
        """


rule dummy_clades:
    """
    Nextclade requires clade membership to be specified for each node
    in the tree. This rule creates a dummy clade membership for each node
    """
    input:
        rules.refine.output.node_data,
    output:
        "results/dummy_clades.json",
    shell:
        """
        jq '.nodes |= map_values({{"clade_membership": "dummy"}})' {input} > {output}
        """


rule export:
    input:
        tree=rules.refine.output.tree,
        metadata=rules.exclude.output.filtered_metadata,
        mutations=rules.ancestral.output.node_data,
        branch_lengths=rules.refine.output.node_data,
        clades=rules.dummy_clades.output,
        auspice_config=AUSPICE_CONFIG,
    output:
        auspice="results/auspice.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --auspice-config {input.auspice_config} \
            --node-data {input.mutations} {input.branch_lengths} {input.clades} \
            --output {output.auspice}
        """


rule subsample_example_sequences:
    input:
        all_sequences="data/sequences.fasta",
        tree_strains="results/tree_strains.txt",
    output:
        example_sequences="results/example_sequences.fasta",
    shell:
        """
        # Exclude tree sequences from all sequences
        seqkit grep -v -f {input.tree_strains} {input.all_sequences} \
        | seqkit sample -n 100 -s 42 > {output.example_sequences}
        """


rule assemble_dataset:
    input:
        tree="results/auspice.json",
        reference=REFERENCE_PATH,
        annotation=GFF_PATH,
        sequences="results/example_sequences.fasta",
        pathogen="resources/pathogen.json",
        readme=README_PATH,
        changelog=CHANGELOG_PATH,
    output:
        tree="dataset/tree.json",
        reference="dataset/reference.fasta",
        annotation="dataset/genome_annotation.gff3",
        sequences="dataset/sequences.fasta",
        pathogen="dataset/pathogen.json",
        readme="dataset/README.md",
        changelog="dataset/CHANGELOG.md",
        dataset_zip="dataset.zip",
    shell:
        """
        cp {input.tree} {output.tree}
        cp {input.reference} {output.reference}
        cp {input.annotation} {output.annotation}
        cp {input.sequences} {output.sequences}
        cp {input.pathogen} {output.pathogen}
        cp {input.readme} {output.readme}
        cp {input.changelog} {output.changelog}
        zip -rj dataset.zip  dataset/*
        """


rule test:
    input:
        dataset="dataset.zip",
        sequences="dataset/sequences.fasta",
    output:
        output=directory("test_out"),
    shell:
        """
        nextclade3 run \
            --input-dataset {input.dataset} \
            --output-all {output.output} \
            {input.sequences}
        """