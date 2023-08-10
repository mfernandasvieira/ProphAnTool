rule create_fasta:
    input:
        cog_report = glob("phastest_results/ZZ_*/CDS/recognizer_output_*/COG_report.tsv"),

    output:
        f"{sys.path[0]}/../../results/cog_extraction.csv"

    run:
        for file in input.cog_report:
            folder_name = os.path.basename(os.path.dirname(file))
            df = pd.read_csv(file, sep='\t')
            keywords_to_exclude = ["phage", "integrase", "capsid", "terminase",
                                   "Replication, recombination and repair",
                                   "Mobilome: prophages, transposons", "Transcription"]
            filter_condition = df["Functional category"].str.contains('|'.join(keywords_to_exclude),
                case=False) | \
                               df["Protein description"].str.contains('|'.join(keywords_to_exclude),
                                   case=False) | \
                               df["DB description"].str.contains('|'.join(keywords_to_exclude),
                                   case=False)
            filtered_df = df[~filter_condition]
            if os.path.isfile(output[0]):
                cog = pd.read_csv(output[0], index_col=0)
                df_concat = pd.concat([cog, filtered_df],axis=0,ignore_index=True)
                df_concat.to_csv(output[0])
            else:
                filtered_df.to_csv(output[0])