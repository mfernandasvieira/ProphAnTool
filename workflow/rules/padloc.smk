rule run_padloc:
    input:
        fasta = f"{sys.path[0]}/../../results/fasta_filter.fasta",

    output:
        f"{sys.path[0]}/../../results/fasta_filter_padloc.csv"

    conda:
        "../envs/padloc.yaml"

    shell:
        "padloc --fna {input.fasta}"

