import matplotlib
import pandas as pd

import warnings
warnings.filterwarnings("ignore")

matplotlib.use('Agg')  # Use the Agg backend
import matplotlib.pyplot as plt

rule taxonomy_classification:
    input:
        contigs = f"{sys.path[0]}/../../results/fasta_filter.fasta"

    output:
        f"{sys.path[0]}/../../results/taxonomy_class/out/phagcn_prediction.csv"

    conda:
        "../envs/phaGCN.yaml"

    shell:
        """
        bash ./workflow/scripts/install_requirements.sh 
        python PhaBOX/PhaGCN_single.py --contigs {input.contigs} --rootpth results/taxonomy_class --out out/ --threads 15
        """

rule taxonomy_plot:
    input:
        taxonomy_csv = rules.taxonomy_classification.output

    output:
        f"{sys.path[0]}/../../results/taxonomy_plot.png"

    run:
        print(input)

        df = pd.read_csv(str(input.taxonomy_csv))
        print(df)
        df['Pred'] = df['Pred'].apply(lambda x: "unknown" if "no_family" in x else x)

        taxonomy_counts = df['Pred'].value_counts()

        bbox_props = dict(boxstyle='square,pad=0.3',fc='w',ec='k',lw=0.72)
        kw = dict(xycoords='data',textcoords='data',arrowprops=dict(arrowstyle='-',color='gray'),zorder=0,va='center')

        fig1, ax1 = plt.subplots()

        # Pie chart without labels
        wedges, _, autotexts = ax1.pie(taxonomy_counts,labels=None,
            startangle=90,autopct='',textprops={'fontsize': 7})

        # Custom legend on the right side
        legend_labels = ['{} - {:.2f}%'.format(label,size) for label, size in
                         zip(taxonomy_counts.index,(taxonomy_counts / taxonomy_counts.sum()) * 100)]
        legend_colors = [p.get_facecolor() for p in wedges]
        legend_handles = [plt.Rectangle((0, 0),1,1,fc=color,edgecolor='black') for color in legend_colors]
        ax1.legend(legend_handles,legend_labels,loc='center left',bbox_to_anchor=(1, 0.5),fontsize='small')


        plt.title('Taxonomy Classification',loc='center',y=1.08)
        plt.savefig(output[0],dpi=300,bbox_inches='tight')


