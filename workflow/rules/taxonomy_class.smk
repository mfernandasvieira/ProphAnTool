rule taxonomy_classification:
    input:
        contigs = f"{sys.path[0]}/../../results/fasta_filter.fasta"

    output:
        f"{sys.path[0]}/../../results/taxonomy_class"

    conda:
        "../envs/phaGCN.yaml"

    shell:
        "python ./PhaBOX/PhaGCN_single.py -f {input.contigs} -out {output}"

