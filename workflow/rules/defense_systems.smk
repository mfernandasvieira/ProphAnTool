import matplotlib
import pandas as pd
import seaborn as sns

matplotlib.use('Agg')  # Use the Agg backend
import matplotlib.pyplot as plt

rule defense_systems:
    input:
        f"{sys.path[0]}/../../results/fasta_filter.fasta"

    output:
        directory(f"{sys.path[0]}/../../results/padloc_results")

    conda:
        "../envs/padloc.yaml"

    shell:
        """
        padloc --db-update
        mkdir -p results/padloc_results
        padloc --fna {input} --outdir {output}
        """


rule defense_systems_plot:
    input:
        defense_folder = rules.defense_systems.output
        #f"{sys.path[0]}/../../results/padloc_results/fasta_filter.fasta_padloc.csv"

    output:
        f"{sys.path[0]}/../../results/defense_systems.png"

    run:
        for file in os.listdir(str(input.defense_folder)):
            if "fasta_filter.fasta_padloc.csv" in file:
                print(file)
                file_path = os.path.join(str(input.defense_folder),file)
                df = pd.read_csv(file_path)
                '''
                # Count occurrences of each system per unique seqid
                system_counts = df.groupby(['seqid', 'system']).size().reset_index(name='count')
                print(system_counts)

                # Pivot the table to have systems as columns
                pivot_table = system_counts.pivot(index='seqid',columns='system',values='count').fillna(0)

                # Create a horizontal bar plot
                plt.figure(figsize=(10, 6))
                sns.barplot(data=pivot_table,orient='h',palette='viridis',errorbar=None,estimator=sum)  # Set estimator to sum

                # Add labels and title
                plt.xlabel('Strains Encoding Antiphage System')
                plt.ylabel('Antiphage Systems')
                plt.title('Occurrences of Defense Mechanisms')
                plt.savefig(output[0], dpi=300, bbox_inches='tight')
                '''
                # Remove duplicate occurrences of seqid and system combinations
                unique_occurrences = df.drop_duplicates(subset=['seqid', 'system'])

                # Count occurrences of each system per unique seqid
                system_counts = unique_occurrences.groupby(['seqid', 'system']).size().reset_index(name='count')

                # Pivot the table to have systems as columns
                pivot_table = system_counts.pivot(index='seqid', columns='system', values='count').fillna(0)

                # Create a horizontal bar plot
                plt.figure(figsize=(10, 6))
                sns.barplot(data=pivot_table, orient='h', palette='viridis', errorbar=None, estimator=sum)  # Set estimator to sum

                # Add labels and title
                plt.xlabel('No. of prophages encoding defense systems')
                plt.ylabel('Defense Mechanisms')
                plt.title('Frequency of defense systems among intact prophages')
                plt.savefig(output[0], dpi=300, bbox_inches='tight')