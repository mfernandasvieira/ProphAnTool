
rule process_detail:
    input:
        prophages = f"{sys.path[0]}/../../results/prophages.json",
        fasta_files = glob("phastest_results/ZZ_*/detail.txt")

    output:
        #[expand(f"{ZZ}/CDS") for ZZ in get_zz_folders()]
        dynamic('phastest_results/{folder_name_2}/CDS/{seq_rec_2}.fasta')

    run:
        with open(input.prophages,'r') as json_file:
            data = json.load(json_file)

        for key, value in data.items():
            for inner_key, inner_value in value.items():
                value[inner_key] = list(inner_value)

        for file in input.fasta_files:
            folder_name = os.path.basename(os.path.dirname(file))
            with open(file) as f:
                cds_df = pd.read_table(f, header=[1], sep=r"\s{3,}", engine="python")
                #print(cds_df)
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
                                v = v.iloc[1:]
                                idx_start = v[v['BLAST_HIT'] == 'attL'].index
                                idx_end = v[v['BLAST_HIT'] == 'attR'].index
                                for i in idx_start:
                                    for j in idx_end:
                                        v = v.loc[i + 1:j - 1]
                                #output_dir = f'phastest_results/{folder_name}/%s.csv'
                                #v.to_csv(output_dir % str(r[0]))

                                output_dir = f'phastest_results/{folder_name}/CDS'
                                if not os.path.exists(output_dir):
                                    os.makedirs(output_dir,exist_ok=True)
                                with open(os.path.join(output_dir, f"{str(r[0])}.fasta"), "w") as fileOutput:
                                    for ind in v.index:
                                        fileOutput.write("> " + str(v['CDS_POSITION'][ind]) + ' ' +
                                                         str(v['BLAST_HIT'][ind]) + "\n")
                                        fileOutput.write(str(v['prophage_PRO_SEQ'][ind]) + "\n")
                                    fileOutput.close()