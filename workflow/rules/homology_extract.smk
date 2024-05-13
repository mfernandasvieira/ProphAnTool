import pandas as pd
from glob import glob
import re

import warnings
warnings.filterwarnings("ignore")

rule homology_extract:
    input:
        glob("phastest_results/GCA_*/summary.txt")

    output:
        f"{sys.path[0]}/../../results/homology_extraction.csv"

    run:
        all_data = []
        regions = {}
        for file in input:
            print(file)
            if os.path.getsize(file) == 0:
                print(f"File '{file}' is empty (0KB).")
            else:
                folder_name = os.path.basename(os.path.dirname(file))
                regions[folder_name] = {}
                count_intact = 0
                count_incomplete = 0
                count_questionable = 0
                with open(file, 'r') as f:
                    summary_file = f.read()
                    search_words = ["genome", "sequence"]
                    list_of_words = summary_file.split()
                    last_occurrence = max((list_of_words.index(word) for word in search_words if
                                           word in list_of_words),default=None)
                    if last_occurrence is not None:
                        bacterial_gs = list_of_words[last_occurrence + 1]
                        df = pd.read_table(file,header=[22],sep=r"\s{2,}",engine="python")
                    else:
                        df = pd.read_table(file,header=[21],sep=r"\s{2,}",engine="python")
                    df = df.iloc[1:, :]
                    if not df.empty:
                        df['REGION_LENGTH'] = df['REGION_LENGTH'].str.replace('Kb', '').astype(float)
                        df = df[~(df['COMPLETENESS(score)'].astype(str).str.startswith('intact') & (df['REGION_LENGTH'] < 15))]

                        # Filter rows where COMPLETENESS(score) starts with 'intact'
                        intact_rows = df[df['COMPLETENESS(score)'].str.startswith('intact')]

                        # Split 'MOST_COMMON_PHAGE_NAME(hit_genes_count)' column based on commas
                        split_columns = intact_rows[
                            'MOST_COMMON_PHAGE_NAME(hit_genes_count)'].str.split(',',expand=True).apply(lambda x: x.map(lambda y: re.search(r'_(.*?)_(.*?)_NC_',y).group(2) if isinstance(y,str) else None))
                        selected_columns = intact_rows[
                            ['REGION', 'COMPLETENESS(score)', 'MOST_COMMON_PHAGE_NAME(hit_genes_count)']]

                        # Concatenate the selected columns with the split columns
                        combined_columns = pd.concat([selected_columns, split_columns],axis=1)

                        # Add a new column for folder name
                        combined_columns.insert(0,'FOLDER_NAME',folder_name)

                        all_data.append(combined_columns)

        combined_df = pd.concat(all_data, ignore_index=True)

        combined_df.to_csv(output[0])