from subprocess import Popen

rule recognizer:
    input:
        fasta_files = glob("phastest_results/ZZ_*/CDS/*.fasta")

    output:
        #[expand(f"{ZZ}/recognizer_output") for ZZ in get_zz_folders()]
        dynamic('phastest_results/{folder_name_3}/CDS/recognizer_output_{fasta_number}/reCOGnizer_results.xlsx')

    threads:
        config["threads"]

    params:
        resources_directory = config["resources_directory"],
        recognizer_databases=','.join(config["recognizer_databases"]),
        download_cdd_resources='' if not config['download_cdd_resources'] else ' -dr'

    run:
        for file in input.fasta_files:
            folder_name = os.path.basename(os.path.dirname(os.path.dirname(file)))
            fasta_number = int(re.search(r"(\d+)\.fasta$", file).group(1))
            recognizer_directory = f"phastest_results/{folder_name}/CDS/recognizer_output_{fasta_number}"
            if os.path.isdir(os.path.join(recognizer_directory)):
                print(f"{recognizer_directory} exists for {file}")
            else:
                print("Run recognizer with",file,"from f",folder_name)
                command = f'recognizer -f {file} -o {recognizer_directory} -rd /home/fernanda/resources_directory'
                result = Popen(command,shell=True).communicate()
                print(result)