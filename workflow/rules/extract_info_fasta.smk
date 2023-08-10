from Bio import SeqIO

rule process_fasta:
    input:
        prophages = f"{sys.path[0]}/../../results/prophages.json",
        fasta_files = glob("phastest_results/ZZ_*/phage_regions.fna")

    output:
        #[expand(f"{ZZ}/fasta_regions") for ZZ in get_zz_folders()]
        dynamic('phastest_results/{folder_name_1}/fasta_regions/{region}/{seq_rec_1}.fasta')

    run:
        with open(input.prophages,'r') as json_file:
            data = json.load(json_file)

        updated_data = {}

        for key, value in data.items():
            updated_value = {}
            for inner_key, inner_value in value.items():
                updated_value[inner_key] = {item: {} for item in inner_value}
            updated_data[key] = updated_value

        for file in input.fasta_files:
            folder_name = os.path.basename(os.path.dirname(file))
            with open(file) as f:
                record = SeqIO.parse(f,"fasta")
                for seq_rec in record:
                    file_count = 0
                    seqLen = len(seq_rec)
                    if 'incomplete' in updated_data[folder_name]:
                        if seq_rec.id in updated_data[folder_name]['incomplete']:
                            updated_data[folder_name]['incomplete'][seq_rec.id] = seqLen
                    if 'intact' in updated_data[folder_name]:
                        if seq_rec.id in updated_data[folder_name]['intact']:
                            updated_data[folder_name]['intact'][seq_rec.id] = seqLen
                    if 'questionable' in updated_data[folder_name]:
                        if seq_rec.id in updated_data[folder_name]['questionable']:
                            updated_data[folder_name]['questionable'][seq_rec.id] = seqLen
                    file_count = file_count + 1

                    if 'incomplete' in updated_data[folder_name] and seq_rec.id in updated_data[folder_name][
                        'incomplete']:
                        output_dir = f'phastest_results/{folder_name}/fasta_regions/incomplete'
                        if not os.path.exists(output_dir):
                            os.makedirs(output_dir, exist_ok=True)
                        with open(os.path.join(output_dir, f"{seq_rec.id}.fasta"),"w") as fho:
                            SeqIO.write(seq_rec,fho,"fasta")

                    if 'questionable' in updated_data[folder_name] and seq_rec.id in updated_data[folder_name][
                        'questionable']:
                        output_dir = f'phastest_results/{folder_name}/fasta_regions/questionable'
                        if not os.path.exists(output_dir):
                            os.makedirs(output_dir, exist_ok=True)
                        with open(os.path.join(output_dir, f"{seq_rec.id}.fasta"),"w") as fho:
                            SeqIO.write(seq_rec,fho,"fasta")

                    if 'intact' in updated_data[folder_name] and seq_rec.id in updated_data[folder_name][
                        'intact']:
                        output_dir = f'phastest_results/{folder_name}/fasta_regions/intact'
                        if not os.path.exists(output_dir):
                            os.makedirs(output_dir, exist_ok=True)
                        with open(os.path.join(output_dir,f"{seq_rec.id}.fasta"),"w") as fho:
                            SeqIO.write(seq_rec,fho,"fasta")
            if file_count == 0:
                raise Exception("No valid sequence in fasta file")

