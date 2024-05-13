import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import warnings
warnings.filterwarnings("ignore")

rule create_tab_abricate:
    input:
        fasta = f"{sys.path[0]}/../../results/fasta_filter.fasta",

    output:
        tab="{database}.tab"

    conda:
        "../envs/abricate.yaml"

    params:
        database = config["abricate_databases"]

    shell:
        """
        for db in {params.database};
        do
            abricate {input.fasta} > results/$db.tab --db $db;
        done
        """

rule create_tab_summary:
    input:
        expand("results/{abricate_database}.tab", abricate_database = config['abricate_databases'])

    output:
        f"{sys.path[0]}/../../results/abricate_summary.tab"

    conda:
        "../envs/abricate.yaml"

    shell:
        "abricate --summary {input} > {output}"

rule create_plot_abricate:
    input:
        resfinder =  f"{sys.path[0]}/../../results/resfinder.tab",
        card = f"{sys.path[0]}/../../results/card.tab",
        vfdb = f"{sys.path[0]}/../../results/vfdb.tab"

    output:
        f"{sys.path[0]}/../../results/heatmap_virulence.png",
        f"{sys.path[0]}/../../results/dataframe_args.csv"

    run:
        card = pd.read_csv(input.card,sep='\t')
        resfinder = pd.read_csv(input.resfinder, sep='\t')
        vfdb = pd.read_csv(input.vfdb, sep='\t')
        #dataframes = [resfinder, card, vfdb]

        #modified_dataframes = []
        #for df in dataframes:
        #    #df['SEQUENCE'] = df['SEQUENCE'].str.replace(r'^\d+_','',regex=True)
        #    df = df.groupby(['SEQUENCE', 'GENE']).size().unstack(fill_value=0)
        #    df.reset_index(inplace=True)
        #    df.rename(columns={'SEQUENCE': 'Folder'},inplace=True)
        #    modified_dataframes.append(df)
        #merged_df = pd.merge(modified_dataframes[0], modified_dataframes[1], on='Folder', how='outer', suffixes=('_df1', '_df2'))
        #merged_df = pd.merge(merged_df,modified_dataframes[2],on='Folder', how='outer')

        merged_df = pd.concat([card, resfinder, vfdb])
        merged_df = merged_df.drop_duplicates(subset=['START', 'END'],keep='first')
        df = merged_df.groupby(['SEQUENCE', 'GENE']).size().unstack(fill_value=0)
        df.reset_index(inplace=True)
        df.rename(columns={'SEQUENCE': 'Folder'},inplace=True)

        df.iloc[:, 1:] = df.iloc[:, 1:].fillna(0).astype(int)
        df = df.set_index('Folder')
        df = df.astype(int)
        df.to_csv(output[1])
        #custom_labels = ["Strain{}".format(i) for i in range(1, len(df) + 1)]
        #df.index = custom_labels

        #np.random.seed(42)
        #matrix_data = np.random.randint(2,size=(5, 5))

        colors = ['#FAA0A0', '#9DD6AD']
        cmap = sns.color_palette(colors)

        vmin, vmax = 0.5, 1.5

        # Plot the heatmap
        plt.figure(figsize=(12, 10))  # Set the size of the figure
        sns.heatmap(df, cmap=cmap, annot=False,cbar=False,vmin=vmin,vmax=vmax, linecolor='white', linewidths=0.3)
        #sns.heatmap(df, cmap='flare',annot=False, fmt="d", linecolor='white', linewidths=0.3)

        #c_bar = ax.collections[0].colorbar
        #c_bar.set_ticks(list(set(merged_df.values.flatten())))

        legend_elements = [
            Patch(facecolor='#FAA0A0',edgecolor='black',label='Absence'),
            Patch(facecolor='#9DD6AD',edgecolor='black',label='Presence')
        ]
        plt.legend(handles=legend_elements,bbox_to_anchor=(1.05, 1),loc='upper left')
        plt.xlabel('Virulence genes')
        plt.ylabel('Intact prophage regions')
        plt.savefig(output[0], dpi=300, bbox_inches='tight')