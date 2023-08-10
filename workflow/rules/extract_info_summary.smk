import os.path
import json
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


rule process_summary:
    input:
        glob("phastest_results/ZZ_*/summary.txt")

    output:
        f"{sys.path[0]}/../../results/result_dict.csv",
        f"{sys.path[0]}/../../results/prophages_frequency.png",
        f"{sys.path[0]}/../../results/distribution_plots.png",
        f"{sys.path[0]}/../../results/prophages.json"

    run:
        for file in input:
            folder_name = os.path.basename(os.path.dirname(file))
            regions[folder_name] = {}
            count_intact = 0
            count_incomplete = 0
            count_questionable = 0
            with open(file, 'r') as f:
                summary_file = f.read()
                search_word = "sequence"
                list_of_words = summary_file.split()
                if search_word in list_of_words:
                    bacterial_gs = list_of_words[list_of_words.index(search_word) + 1]
                    df = pd.read_table(file,header=[22],sep=r"\s{2,}",engine="python")
                else:
                    df = pd.read_table(file,header=[21],sep=r"\s{2,}",engine="python")
                df = df.iloc[1:, :]
                df['REGION_LENGTH'] = df['REGION_LENGTH'].str.replace('Kb', '').astype(float)
                df = df[~(df['COMPLETENESS(score)'].str.startswith('intact') & (df['REGION_LENGTH'] < 15))]
                for item in df['COMPLETENESS(score)']:
                    if item.startswith('intact'):
                        found_intact_regions = df.loc[df['COMPLETENESS(score)'] == item, 'REGION'].values
                        if 'intact' in regions[folder_name]:
                            regions[folder_name]['intact'].update(found_intact_regions)
                        else:
                            regions[folder_name]['intact'] = {row for row in found_intact_regions}
                        count_intact += 1
                    elif item.startswith('questionable'):
                        found_questionable_regions = df.loc[df['COMPLETENESS(score)'] == item, 'REGION'].values
                        if 'questionable' in regions[folder_name]:
                            regions[folder_name]['questionable'].update(found_questionable_regions)
                        else:
                            regions[folder_name]['questionable'] = {row for row in found_questionable_regions}
                        count_questionable += 1
                    else:
                        found_incomp_regions = df.loc[df['COMPLETENESS(score)'] == item, 'REGION'].values
                        if 'incomplete' in regions[folder_name]:
                            regions[folder_name]['incomplete'].update(found_incomp_regions)
                        else:
                            regions[folder_name]['incomplete'] = {row for row in found_incomp_regions}
                        count_incomplete += 1


                data = {'Folder': folder_name, 'Intact': count_intact, 'Incomplete': count_incomplete,
                        'Questionable': count_questionable, 'Bact_genome_size': bacterial_gs[:-1]}

                output_dir = str(output)
                if os.path.isfile(output[0]):
                    df_exist = pd.read_csv(output[0],index_col=0)
                    df_to_save = pd.concat([pd.DataFrame([data]), df_exist],ignore_index=True)
                    df_to_save['Total_prophages'] = df_to_save['Intact'] + df_to_save['Incomplete'] + df_to_save[
                        'Questionable']
                    df_to_save.to_csv(output[0], index=True)
                else:
                    df_to_save = pd.DataFrame([data])
                    df_to_save.to_csv(output[0])

        sns.boxplot(data=df_to_save)
        plt.ylim(bottom=0)
        plt.savefig(output[1], dpi=300)

        # Create distribution plot
        intact = df_to_save['Intact'].to_list()
        incomplete = df_to_save['Incomplete'].to_list()
        questionable = df_to_save['Questionable'].to_list()
        bacterial_genome_size = df_to_save['Bact_genome_size'].to_list()
        df_to_save['Bact_genome_size'] = df_to_save['Bact_genome_size'].astype(float)
        df_to_save['Bact_genome_size'] = (df_to_save['Bact_genome_size'] / 1000000).round(1)
        df_to_save['Defective'] = df_to_save['Incomplete'] + df_to_save['Questionable']
        df_to_save = df_to_save.set_index('Folder')
        grouped_data = df_to_save.groupby('Bact_genome_size').mean()
        error = df_to_save.groupby('Bact_genome_size').sem()


        fig, axes = plt.subplots(nrows=1,ncols=3,figsize=(12, 4))
        axes[0].errorbar(grouped_data.index,grouped_data['Total_prophages'],yerr=error[
            'Total_prophages'],fmt='D',color="black",capsize=3,markersize=5,markeredgecolor='black',markeredgewidth=1)
        axes[0].set_xlabel('Bacterial genome size (Mb)')
        axes[0].set_ylabel('Average number of prophages')
        axes[0].spines[['right', 'top']].set_visible(False)

        axes[1].errorbar(grouped_data.index,grouped_data['Intact'],yerr=error[
            'Intact'],label='Intact',fmt='D',color="black",capsize=3,markersize=5,markeredgecolor='black',markeredgewidth=1)
        axes[1].set_xlabel('Bacterial genome size (Mb)')
        axes[1].set_ylabel('Average number of intact prophages')
        axes[1].spines[['right', 'top']].set_visible(False)

        axes[2].errorbar(grouped_data.index,grouped_data['Defective'],yerr=error[
            'Defective'],label='Defective',fmt='D',color="black",capsize=3,markersize=5,markeredgecolor='black',markeredgewidth=1)
        axes[2].set_xlabel('Bacterial genome size (Mb)')
        axes[2].set_ylabel('Average number of defective prophages')
        axes[2].spines[['right', 'top']].set_visible(False)

        plt.tight_layout()

        ymin = min(grouped_data[['Intact', 'Incomplete', 'Questionable']].min().values)
        ymax = max(grouped_data[['Intact', 'Incomplete', 'Questionable']].max().values)
        for ax in axes:
            ax.set_ylim(0,12)
        plt.savefig(output[2],dpi=300)

        for key, value in regions.items():
            for inner_key, inner_value in value.items():
                value[inner_key] = list(inner_value)

        with open(output[3], "w", encoding='utf-8') as json_file:
            json.dump(regions, json_file, ensure_ascii=False, indent=2)

