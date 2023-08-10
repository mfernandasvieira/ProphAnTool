rule create_fasta:
    input:
        f"{sys.path[0]}/../../results/fasta_filter.fasta",

    output:
