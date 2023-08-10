from Bio import SeqIO

rule create_fasta:
    input:
        prophages = f"{sys.path[0]}/../../results/prophages.json",
        fasta_files = glob("phastest_results/ZZ_*/phage_regions.fna")

    output:
        f"{sys.path[0]}/../../results/dna_fasta.fasta"

    run:
        with open(input.prophages,'r') as json_file:
            data = json.load(json_file)

        for key, value in data.items():
            for inner_key, inner_value in value.items():
                value[inner_key] = list(inner_value)

        for file in input.fasta_files:
            folder_name = os.path.basename(os.path.dirname(file))
            with open(file) as f:
                fasta_sequences = SeqIO.parse(f, 'fasta')
                with open(output[0], 'a') as out_file:
                    for seq_record in fasta_sequences:
                        if 'intact' in data[folder_name]:
                            for i in data[folder_name]['intact']:
                                if seq_record.id == i:
                                    seq_record.id = f"{seq_record.id}_{folder_name}"
                                    seq_record.description = ' '.join(seq_record.description.split()[1:])
                                    SeqIO.write(seq_record,out_file,'fasta')