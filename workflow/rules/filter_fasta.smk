from Bio import SeqIO
rule extract_att:
    input:
        fasta_files= glob("phastest_results/ZZ_*/phage_regions.fna"),
        prophages= f"{sys.path[0]}/../../results/prophages.json",
        detail_files = glob("phastest_results/ZZ_*/detail.txt")

    output:
        f"{sys.path[0]}/../../results/df_att.csv",
        f"{sys.path[0]}/../../results/fasta_filter.fasta"

    run:
        with open(input.prophages,'r') as json_file:
            data = json.load(json_file)

        for key, value in data.items():
            for inner_key, inner_value in value.items():
                value[inner_key] = list(inner_value)

        data_for_df = []
        for file in input.detail_files:
            folder_name = os.path.basename(os.path.dirname(file))
            with open(file) as f:
                cds_df = pd.read_table(f,header=[1],sep=r"\s{3,}",engine="python")
                m = cds_df.CDS_POSITION.str.contains('##').cumsum()
                d = {f'df{i}': g for i, g in cds_df.groupby(m)}
                for k, v in d.items():
                    s = v.iloc[[0]]
                    region = re.findall(r'\d+',str(s['CDS_POSITION']))
                    if 'intact' in data[folder_name]:
                        if len(region) >= 2 and region[1] in data[folder_name]['intact']:
                            result = re.search('####(.*)',str(s['CDS_POSITION']))
                            region_name = result.group(1).replace(" ","")
                            r = re.findall(r'\d+',region_name)
                            if r[0] == region[1]:
                                df = v.iloc[1:]
                                mask = df['BLAST_HIT'] == 'attL'
                                prophage_pro_seq_values = df.loc[mask, 'PID']
                                for value in prophage_pro_seq_values:
                                    sequence_start = value
                                    data_for_df.append([folder_name, region[1], sequence_start])
        df = pd.DataFrame(data_for_df, columns=['Folder_name', 'Region', 'attR'])
        df.to_csv(output[0])
        for file in input.fasta_files:
            folder_name = os.path.basename(os.path.dirname(file))
            with open(file) as f:
                fasta_sequences = SeqIO.parse(f, 'fasta')
                with open(output[1], 'a') as out_file:
                    for seq_record in fasta_sequences:
                        if 'intact' in data[folder_name]:
                            for i in data[folder_name]['intact']:
                                if seq_record.id == i:
                                    matching_rows = df[
                                        (df['Region'] == seq_record.id) & (df['Folder_name'] == folder_name)]
                                    matching_rows = matching_rows.drop_duplicates(subset=['Folder_name',
                                                                                          'Region'],keep='first')
                                    if not matching_rows.empty:
                                        for index, row in matching_rows.iterrows():
                                            att_value = row['attR']
                                            start_index = -1
                                            end_index = -1
                                            first_occurrence = True

                                            if first_occurrence and att_value in seq_record.seq:
                                                start_index = seq_record.seq.find(att_value)
                                                first_occurrence = False

                                            if att_value in seq_record.seq:
                                                end_index = seq_record.seq.rfind(att_value)
                                                if end_index == start_index:
                                                    end_index = len(seq_record.seq)

                                        if start_index != -1:
                                            extracted_sequence = seq_record.seq[start_index:end_index]
                                        else:
                                            print('Error')
                                        seq_record.seq = extracted_sequence
                                        seq_record.id = f"{seq_record.id}_{folder_name}"
                                        seq_record.description = ' '.join(seq_record.description.split()[1:])
                                        SeqIO.write(seq_record,out_file,'fasta')
                                    else:
                                        seq_record.id = f"{seq_record.id}_{folder_name}"
                                        seq_record.description = ' '.join(seq_record.description.split()[1:])
                                        SeqIO.write(seq_record,out_file,'fasta')