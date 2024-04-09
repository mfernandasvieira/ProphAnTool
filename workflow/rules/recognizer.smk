from subprocess import Popen
import glob

rule recognizer:
    input:
        f"{sys.path[0]}/../../results/fasta_cds_recognizer.fasta"

    output:
        #[expand(f"{ZZ}/recognizer_output") for ZZ in get_zz_folders()]
        #dynamic('phastest_results/{folder_name_3}/CDS/recognizer_output_{fasta_number}/COG_report.tsv')
        directory(f"{sys.path[0]}/../../results/output_recognizer")

    threads:
        config["threads"]

    params:
        resources_directory = config["resources_directory"],
        recognizer_databases=','.join(config["recognizer_databases"]),

    conda:
        "../envs/recognizer.yaml"

    shell:
        "recognizer -f {input} -t {threads} -o {output} "
        "-rd {params.resources_directory} -dbs {params.recognizer_databases}"