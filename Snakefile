configfile: "config.yaml"


samples = {
    file.strip(".vcf"): os.path.join(root, file)
    for root, _, files in os.walk(config["input_path"])
    for file in files
    if file.endswith(".vcf")
}

print(samples)

rule all:
    input:
        expand("annotated/{sample}.annotated.vcf", sample=samples.keys()),


rule annotate_variants:
    input:
        calls=lambda wc: samples[wc.sample],
        cache="resources/vep/cache",
        plugins="resources/vep/plugins",
        fasta="refs/genome.fasta"
    output:
        calls="annotated/{sample}.annotated.vcf",
        stats=temp("annotated/{sample}.stats.html"),
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=config["annotations"]["vep"]["plugins"],
        extra="{}".format(
            config["annotations"]["vep"]["params"]
        ),
    log:
        "logs/vep/{sample}.annotate.log",
    threads: max(workflow.cores / len(samples), 1) if len(samples) else 1
    wrapper:
        "0.78.0/bio/vep/annotate"


rule filter_by_annotation:
    input:
        "get_annotated_bcf",
    output:
        "results/calls/{group}.{filter}.filtered_ann.bcf",
    log:
        "logs/filter-calls/annotation/{group}.{filter}.log",
    params:
        filter=lambda w: config["calling"]["filter"][w.filter],
    conda:
        "../envs/vembrane.yaml"
    shell:
        "vembrane filter {params.filter:q} {input} --output-fmt bcf --output {output} &> {log}"


rule get_genome:
    output:
        "refs/genome.fasta"
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    log:
        "logs/get_genome.log"
    cache: True  # save space and time with between workflow caching (see docs)
    wrapper:
        "0.78.0/bio/reference/ensembl-sequence"

rule get_vep_cache:
    output:
        directory("resources/vep/cache"),
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    log:
        "logs/vep/cache.log",
    cache: True  # save space and time with between workflow caching (see docs)
    wrapper:
        "0.78.0/bio/vep/cache"


rule get_vep_plugins:
    output:
        directory("resources/vep/plugins"),
    params:
        release=config["ref"]["release"],
    log:
        "logs/vep/plugins.log",
    wrapper:
        "0.78.0/bio/vep/plugins"
