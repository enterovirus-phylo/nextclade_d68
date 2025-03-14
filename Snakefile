# Set the parameters
REFERENCE_ACCESSION =   "AY426531.1"
TAXON_ID =              42789
GENES =                 ["VP4", "VP2", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D"]
ALLOWED_DIVERGENCE =    "1000" # was 
MIN_DATE =              "1950-01-01"
MIN_LENGTH =            "6000" # was 6000 for whole genome build on Nextstrain
MAX_SEQS =              "10000" #TODO: set to 10000 for testing
ROOTING =               "mid_point"  # alternative root using outgroup, e.g. the reference "AY426531.1"
ID_FIELD=               "accession" # either accession or strain, used for meta-id-column in augur

# Set the paths
GFF_PATH =              "dataset/genome_annotation.gff3"
PATHOGEN_JSON =         "dataset/pathogen.json"
GENBANK_PATH =          "resources/reference.gbk"
REFERENCE_PATH =        "dataset/reference.fasta"
README_PATH =           "dataset/README.md"
CHANGELOG_PATH =        "dataset/CHANGELOG.md"
AUSPICE_CONFIG =        "resources/auspice_config.json"
EXCLUDE =               "resources/exclude.txt"
SEQUENCES =             "data/sequences.fasta"
METADATA =              "data/metadata.tsv"
CLADES =                "resources/clades.tsv"

FETCH_SEQUENCES = False

rule all_augur:
    input: 
        augur_jsons = "results/auspice.json"

rule all_nextclade:
    input: 
        augur_jsons = "test_out/"


if FETCH_SEQUENCES == True:

    include: "ingest/Snakefile"


rule add_reference_to_include:
    """
    Create an include file for augur filter
    """
    input:
        include = "resources/include.txt",
    output:
        "results/include.txt",
    shell:
        """
        echo "{REFERENCE_ACCESSION}" >> results/include.txt
        """

rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering
        """
    input:
        sequences = SEQUENCES,
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
        sequences = SEQUENCES,
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = METADATA,
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
            --output {output.filtered_sequences} \
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
        #high-diversity 
        penalty_gap_extend = 1, #make longer gaps more costly - default is 0
        penalty_gap_open = 13,  #make gaps more expensive relative to mismatches - default is 13
        penalty_gap_open_in_frame = 18, #make gaps more expensive relative to mismatches - default is 7
        penalty_gap_open_out_of_frame = 23, #make out of frame gaps more expensive - default is 8 # prev was 19
        kmer_length = 6, #reduce to find more matches - default is 10
        kmer_distance = 25, #reduce to try more seeds - default is 50
        min_match_length = 30, #reduce to keep more seeds - default is 40
        allowed_mismatches = 15, #increase to keep more seeds - default is 8
        min_length = 30, # min_length - default is 100
        #cost of a mutation is 4
    shell:
        """
        nextclade3 run \
        -j {threads} \
        {input.sequences} \
        --input-ref {input.reference} \
        --input-annotation {input.annotation} \
        --penalty-gap-open {params.penalty_gap_open} \
        --penalty-gap-extend {params.penalty_gap_extend} \
        --penalty-gap-open-in-frame {params.penalty_gap_open_in_frame} \
        --penalty-gap-open-out-of-frame {params.penalty_gap_open_out_of_frame} \
        --kmer-length {params.kmer_length} \
        --kmer-distance {params.kmer_distance} \
        --min-match-length {params.min_match_length} \
        --allowed-mismatches {params.allowed_mismatches} \
        --min-length {params.min_length} \
        --include-reference false \
        --output-tsv {output.tsv} \
        --output-translations {params.translation_template} \
        --output-fasta {output.alignment} 
        """


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
        metadata = METADATA,
        exclude = EXCLUDE,
        outliers = rules.get_outliers.output.outliers,
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
            --exclude {input.exclude} {input.outliers} \
            --output {output.filtered_sequences} \
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

rule clades:
    input:
        tree = rules.refine.output.tree,
        mutations = rules.ancestral.output.node_data,
        clades = CLADES
    output:
        "results/clades.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.mutations} \
            --clades {input.clades} \
            --output-node-data {output}
        """

rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.exclude.output.filtered_metadata,
        mutations = rules.ancestral.output.node_data,
        branch_lengths = rules.refine.output.node_data,
        clades = rules.clades.output, # dummy_clades if not set yet
        auspice_config = AUSPICE_CONFIG,
    params:
        strain_id_field = ID_FIELD,
    output:
        auspice = "results/auspice.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --auspice-config {input.auspice_config} \
            --node-data {input.mutations} {input.branch_lengths} {input.clades} \
            --output {output.auspice}
        """


rule subsample_example_sequences:
    input:
        all_sequences = SEQUENCES,
        tree_strains = rules.exclude.output.strains,
    output:
        example_sequences = "results/example_sequences.fasta",
    shell:
        """
        # Exclude tree sequences from all sequences
        seqkit grep -v -f {input.tree_strains} {input.all_sequences} \
        | seqkit sample -n 100 -s 42 > {output.example_sequences}
        """


rule assemble_dataset:
    input:
        tree = rules.export.output.auspice,
        reference = REFERENCE_PATH,
        annotation = GFF_PATH,
        sequences = rules.subsample_example_sequences.output.example_sequences,
        pathogen = PATHOGEN_JSON,
        readme = README_PATH,
        changelog = CHANGELOG_PATH,
    output:
        tree = "out-dataset/tree.json",
        reference = "out-dataset/reference.fasta",
        annotation = "out-dataset/genome_annotation.gff3",
        sequences = "out-dataset/sequences.fasta",
        pathogen = "out-dataset/pathogen.json",
        readme = "out-dataset/README.md",
        changelog = "out-dataset/CHANGELOG.md",
        dataset_zip = "dataset.zip",
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

rule clean:
    shell:
        """
        rm -rf results
        rm -r dataset/tree.json
        """