import matplotlib
matplotlib.use('Agg')  # Use the Agg backend
import os.path
import json
import sys
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, RegularPolygon
from matplotlib.path import Path
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D
from glob import glob

import warnings
warnings.filterwarnings("ignore")


rule process_summary:
    input:
        glob("phastest_results/GCA_*/summary.txt")

    output:
        f"{sys.path[0]}/../../results/result_dict.csv",
        f"{sys.path[0]}/../../results/prophages_frequency.png",
        f"{sys.path[0]}/../../results/distribution_plots.png",
        f"{sys.path[0]}/../../results/prophages.json"

    run:
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
                    else:
                        count_intact = 0
                        count_questionable = 0
                        count_incomplete = 0

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
                        print(df_to_save)
                        df_to_save.to_csv(output[0])

        for index, row in df_to_save.iterrows():
            print("Folder:",row['Folder'],"Bact_genome_size:",row['Bact_genome_size'])
        df = df_to_save[df_to_save["Total_prophages"] > 0]
        sns.boxplot(data=df)
        plt.ylim(bottom=0)
        plt.ylabel('Number of prophages per genome')
        plt.savefig(output[1], dpi=300)

        # Create distribution plot
        df_to_save["Bact_genome_size"] = df_to_save["Bact_genome_size"].astype(int)
        print(df_to_save.dtypes)
        print(df_to_save)
        df_to_save = df_to_save[df_to_save["Bact_genome_size"] > 1000000]
        print(df_to_save)

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
        axes[0].set_xlabel('Bacterial genome size (Mb)')     #Plasmid sequence size
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

rule plot_radar_chart:
    input:
        f"{sys.path[0]}/../../results/result_dict.csv"

    output:
        f"{sys.path[0]}/../../results/radar_chart.png"

    run:
        def radar_factory(num_vars, frame='circle'):
            # calculate evenly-spaced axis angles
            theta = np.linspace(0,2 * np.pi,num_vars,endpoint=False)

            class RadarTransform(PolarAxes.PolarTransform):

                def transform_path_non_affine(self, path):
                    # Paths with non-unit interpolation steps correspond to gridlines,
                    # in which case we force interpolation (to defeat PolarTransform's
                    # autoconversion to circular arcs).
                    if path._interpolation_steps > 1:
                        path = path.interpolated(num_vars)
                    return Path(self.transform(path.vertices),path.codes)

            class RadarAxes(PolarAxes):

                name = 'radar'
                PolarTransform = RadarTransform

                def __init__(self, *args, **kwargs):
                    super().__init__(*args,**kwargs)
                    # rotate plot such that the first axis is at the top
                    self.set_theta_zero_location('N')

                def fill(self, *args, closed=True, **kwargs):
                    """Override fill so that line is closed by default"""
                    return super().fill(closed=closed,*args,**kwargs)

                def plot(self, *args, **kwargs):
                    """Override plot so that line is closed by default"""
                    lines = super().plot(*args,**kwargs)
                    for line in lines:
                        self._close_line(line)

                def _close_line(self, line):
                    x, y = line.get_data()
                    # FIXME: markers at x[0], y[0] get doubled-up
                    if x[0] != x[-1]:
                        x = np.append(x,x[0])
                        y = np.append(y,y[0])
                        line.set_data(x,y)

                def set_varlabels(self, labels):
                    self.set_thetagrids(np.degrees(theta),labels)

                def _gen_axes_patch(self):
                    # The Axes patch must be centered at (0.5, 0.5) and of radius 0.5
                    # in axes coordinates.
                    if frame == 'circle':
                        return Circle((0.5, 0.5),0.5)
                    elif frame == 'polygon':
                        return RegularPolygon((0.5, 0.5),num_vars,
                            radius=.5,edgecolor="k")
                    else:
                        raise ValueError("Unknown value for 'frame': %s" % frame)

                def _gen_axes_spines(self):
                    if frame == 'circle':
                        return super()._gen_axes_spines()
                    elif frame == 'polygon':
                        # spine_type must be 'left'/'right'/'top'/'bottom'/'circle'.
                        spine = Spine(axes=self,
                            spine_type='circle',
                            path=Path.unit_regular_polygon(num_vars))
                        # unit_regular_polygon gives a polygon of radius 1 centered at
                        # (0, 0) but we want a polygon of radius 0.5 centered at (0.5,
                        # 0.5) in axes coordinates.
                        spine.set_transform(Affine2D().scale(.5).translate(.5,.5)
                                            + self.transAxes)
                        return {'polar': spine}
                    else:
                        raise ValueError("Unknown value for 'frame': %s" % frame)

            register_projection(RadarAxes)
            return theta


        def plot_percentage_labels(ax, theta, data):
            for d in data:
                for i, y in enumerate(d):
                    ax.text(theta[i],y,f"{y:.0%}",va='center',ha='center',fontsize=10,weight='bold')


        def example_data(defective_more_than_one, intact_more_than_one, defective_prophages, intact_prophages, p_with_prophages):
            # Your data here
            data = [
                ['More than one defective', 'More than one intact', 'Defective prophages', 'Intact prophages',
                 'Total prophages'],
                ('Sample Data', [
                    [defective_more_than_one, intact_more_than_one, defective_prophages, intact_prophages, p_with_prophages]
                ])
            ]
            return data



        df = pd.read_csv(input[0],index_col=0)
        df = df[df["Bact_genome_size"] > 1000000]
        print(df)
        df['Defective'] = df['Incomplete'] + df['Questionable']
        df.to_csv('test_dataframe.csv')

        defective_more_than_one = len(df[df["Defective"] > 1]) / df.shape[0]
        print('total:',df.shape[0])
        print('defective_more_than_one', len(df[df["Defective"] > 1]))

        intact_more_than_one =len(df[df["Intact"] > 1]) / df.shape[0]
        print('intact_more_than_one',len(df[df["Intact"] > 1]))

        defective_prophages = len(df[df["Defective"] > 0]) / df.shape[0]
        print('defective_prophages', len(df[df["Defective"] > 0]))

        intact_prophages = len(df[df["Intact"] > 0]) / df.shape[0]
        print('intact_prophages',len(df[df["Intact"] > 0]))

        total_strains = len(df)
        strains_with_no_prophages = len(df[(df['Intact'] == 0) &
                                           (df['Incomplete'] == 0) &
                                           (df['Defective'] == 0)])
        p_with_prophages = ((total_strains - strains_with_no_prophages) / total_strains)

        if p_with_prophages == 0 and total_strains > 0:
            p_with_prophages = 1

        print(defective_more_than_one, intact_more_than_one, defective_prophages, intact_prophages, p_with_prophages)

        N = 5
        theta = radar_factory(N,frame='polygon')

        data = example_data(defective_more_than_one, intact_more_than_one, defective_prophages, intact_prophages, p_with_prophages)
        spoke_labels = data.pop(0)

        fig, ax = plt.subplots(figsize=(6, 6),subplot_kw=dict(projection='radar'))

        ax.set_rgrids([0.25, 0.5, 0.75, 1])
        ax.set_yticklabels(['{:.0f}%'.format(t * 100) for t in plt.gca().get_yticks()],color='grey')

        d = data[0][1][0]
        ax.plot(theta, d, color='#77DD77')
        ax.fill(theta, d, facecolor='#77DD77',alpha=0.45,label='_nolegend_')
        ax.set_varlabels(spoke_labels)
        ax.tick_params(axis='both', which='major', pad=25)
        plot_percentage_labels(ax, theta, data[0][1])
        plt.savefig(output[0], dpi=300, bbox_inches='tight')



