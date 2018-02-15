import argparse
import pandas as pd
import numpy as np
from bokeh.io import show
from bokeh.models import ColumnDataSource, HoverTool, Button, Range1d, TextInput
from bokeh.plotting import figure
from bokeh.layouts import layout
from bokeh.io import curdoc
import sys
import os
from time import localtime, strftime


def chrom_number(chrom):
    """
    Used to order chromosomes.
    """
    x = chrom.split('chr')[1]
    try:
        return int(x)
    except ValueError:
        if x == 'X':
            return 23
        elif x == 'Y':
            return 24
        else:
            return None


def genome_file_handling(chromosome_file):
    """
    Read through genome file to get lengths of chromosomes and get an ordered list of chromosome for x-axis
    """

    genome_lens = {}
    with open(chromosome_file, 'r') as genome:
        for entry in genome.readlines():
            chrom_name, chrom_len = entry.strip().split('\t')

            # Skip unwanted chr
            if any(s in chrom_name for s in ['chrM', 'chrEBV', 'random', 'chrUn']):
                continue

            genome_lens[chrom_name] = int(chrom_len)

    chrom_order = list(genome_lens.keys())
    chrom_order = sorted(chrom_order, key=lambda x: chrom_number(x))

    # Create translation key that can set the relative start point for a particular chromosome
    translate_key = {}
    current_len = 0
    for chrom in chrom_order:
        translate_key[chrom] = current_len
        current_len += genome_lens[chrom]

    return chrom_order, translate_key


def get_cluster_df(input_bed, translate_key, chrom_order):
    """
    Read through bed file to create dataframe with information.
    """
    cluster_regions = []

    with open(input_bed, 'r') as bed:
        for interval in bed.readlines():

            chrom, start, end, cluster_id, weight, direction = interval.strip().split('\t')
            if chrom in chrom_order:
                cluster_region = (int(cluster_id),
                                  str(chrom),
                                  int(start),
                                  int(end),
                                  int(start) + translate_key[chrom],
                                  int(end) + translate_key[chrom])

                cluster_regions.append(cluster_region)

    return pd.DataFrame(data=cluster_regions, columns=['Cluster', 'Chromosome', 'Start', 'End', 'xmin', 'xmax'])


def get_cluster_info(dataframe, phasing_threshold=100000, include_fragments=True, include_chromosomes=True,
                     include_phasing=True):
    """
    Get basic cluster information from dataframe
    """

    cluster_data = []
    # Loop over each cluster
    for cluster in dataframe.Cluster.unique():
        # Extract dataframe for cluster
        cluster_dataframe = dataframe[dataframe.Cluster == cluster]

        # Get number of fragments
        fragments = 0
        if include_fragments:
            fragments = len(cluster_dataframe)

        # Get number of chromosomes
        chromosomes = 0
        if include_chromosomes:
            chromosomes = len(cluster_dataframe.Chromosome.unique())

        # Get number of phased reads i.e pairs of neigbouring fragments within threshold
        phasing = 0
        if include_phasing:
            if fragments > 1:
                for chromosome in cluster_dataframe.Chromosome.unique():
                    fragment_start_values = cluster_dataframe[cluster_dataframe.Chromosome == chromosome].Start
                    if len(fragment_start_values) > 1:
                        for start1, start2 in zip(fragment_start_values[:-1], fragment_start_values[1:]):
                            if abs(start1 - start2) <= phasing_threshold:
                                phasing += 1

        # Store as tuple
        data = (cluster, fragments, chromosomes, phasing)

        # Apped to list
        cluster_data.append(data)

    return pd.DataFrame(data=cluster_data, columns=['Cluster', 'Fragments', 'Chromosomes', 'Phased_regions'])


def sort_dataframe(dataframe, info_dataframe, sort_type):
    """
    Sort and index dataframe based on type of sort
    """
    ascending = True
    if sort_type == 'Fragments' or sort_type == 'Phased regions':
        ascending = False

    clusters_index = {row.Cluster: nr for nr, (i, row) in enumerate(
        info_dataframe.sort_values(by=sort_type, ascending=ascending).iterrows())}

    dataframe['Index'] = pd.Series((np.zeros(len(dataframe)), 1))
    for i, df_row in dataframe.iterrows():
        dataframe.set_value(i, 'Index', clusters_index[df_row.Cluster])

    return dataframe.sort_values(by='Index'), clusters_index


def main(bed_file_name, genome_file_name, sort_cluster_type):
    #
    # Get initial information from bed file and genome file
    #

    print(strftime("%Y-%m-%d %H:%M:%S ", localtime()), "Reading genome file: ", genome_file_name)

    chrom_order, translate_key = genome_file_handling(genome_file_name)

    print(strftime("%Y-%m-%d %H:%M:%S ", localtime()), "Creating cluster dataframe from bed file: ", bed_file_name)

    cluster_df = get_cluster_df(bed_file_name, translate_key, chrom_order)

    print(strftime("%Y-%m-%d %H:%M:%S ", localtime()), "Collecting cluster info to dataframe")

    cluster_info = get_cluster_info(cluster_df, sort_cluster_type,
                                    include_chromosomes=(sort_cluster_type == 'Chromosomes'),
                                    include_fragments=(sort_cluster_type == 'Fragments'),
                                    include_phasing=(sort_cluster_type == 'Phased_regions'))

    print(strftime("%Y-%m-%d %H:%M:%S ", localtime()), "Sorting cluster dataframe based on: ", sort_cluster_type)

    cluster_df_sorted, cluster_index = sort_dataframe(cluster_df, cluster_info, sort_cluster_type)

    cluster_df_sorted.Cluster = cluster_df_sorted.Cluster.astype(str)

    print(strftime("%Y-%m-%d %H:%M:%S ", localtime()), "Setting up data source for plotting.")

    source = ColumnDataSource(cluster_df_sorted)

    #
    # Plotting
    #
    p = figure(plot_width=1500,
               plot_height=800,
               toolbar_location='right')

    cr = p.square(y='Index',
                    x='xmin',
                    line_width=4,
                    source=source,
                    fill_alpha=0.3,
                    line_alpha=0.3,
                    hover_line_color='red',
                    hover_fill_color='red',
                    hover_line_alpha=0.4,
                    hover_fill_alpha=0.4)

    p.add_tools(HoverTool(tooltips=[("Position", "(@Chromosome, @Start, @End)")],
                          renderers=[cr],
                          mode='hline'))

    # p.ygrid.grid_line_color = 'black'
    # p.ygrid.grid_line_width = 15
    # p.ygrid.grid_line_alpha = 0.1
    # p.ygrid = list(range(0, 20))

    p.xaxis.axis_label = None
    p.yaxis.axis_label = 'Cluster id'

    p.y_range = Range1d(0, 20)
    p.yaxis.ticker = list(range(0, len(cluster_index)))
    p.yaxis.major_label_overrides = {str(i): '#' + str(i) + ' = ' + str(key).rjust(8) for key, i in cluster_index.items() if i in range(0, len(cluster_index))}

    p.outline_line_color = None

    p.xaxis.ticker = [translate_key[chrom] for chrom in chrom_order]
    p.xaxis.major_label_overrides = {translate_key[chrom]: chrom for chrom in chrom_order}
    p.xaxis.major_label_orientation = 3.14/2

    p.xgrid.grid_line_color = None

    p.title.text = "Cluster overview " #+ str(current_pos) + ' - ' + str(current_pos + 20) + ' of total (' + str(len(cluster_index)) + ')'

    #show(p)

    # button_fwd = Button(label='Forward 20', button_type='success')
    # button_rev = Button(label='Reverse 20', button_type='success')
    #
    # def update_y_axis(change, cluster_index):
    #     global current_pos
    #     position = current_pos + change
    #     if 0 <= position:
    #         p.y_range.start += change
    #         p.y_range.end += change
    #         p.yaxis.ticker = list(range(position-1, position + 20))
    #         p.yaxis.major_label_overrides = {str(i): '#' + str(i) + ' = ' + str(key).rjust(8) for key, i in cluster_index.items() if
    #                                          i in range(position-1, position + 20)}
    #         current_pos = position
    #         p.title.text = "Cluster overview " + str(current_pos) + ' - ' + str(current_pos + 20) + ' of total (' + str(
    #             len(cluster_index)) + ')'
    #
    #     print(current_pos)
    #
    # def cb_fwd():
    #     print("FORWARD")
    #     update_y_axis(20, cluster_index)
    #
    # def cb_rev():
    #     print("REVERSE")
    #     update_y_axis(-20, cluster_index)
    #
    # button_fwd.on_click(cb_fwd)
    # button_rev.on_click(cb_rev)
    #
    # curdoc().add_root(layout([[p], [button_fwd, button_rev]]))
    index = TextInput(value="O", title="Range index:")

    def update_range(attr, old, new):
        p.y_range.start = int(index.value)
        p.y_range.end = int(index.value) + 20

    index.on_change('value', update_range)

    curdoc().add_root(layout([[p], [index]]))
    curdoc().title = "Test"

#
# Startup options
#

# Genome file
genome_file = 'data/hg38.genome'

if not sys.argv[1]:
    # List available data files
    data_files = os.listdir('data')
    bed_files = [file_name for file_name in data_files if '.bed' in file_name]

    print('AVAILABLE DATA FILES:')
    print('-'*30)
    for nr, file in enumerate(bed_files):
        print(nr, ':', file)
    print('-'*30)

    file_index = int(input('Give selected file index:'))
    bed_file = 'data/' + bed_files[file_index]

else:
    bed_file = sys.argv[1]

print('CURRENT FILE = ', bed_file)

print('AVAILABLE SORT OPTIONS:')
print('-'*30)
options = ['Cluster', 'Fragments', 'Chromosomes', 'Phased_regions']
for nr, option in enumerate(options):
    print(nr, ':', option)
print('-'*30)

option_index = int(input('Give selected sort option index:'))
sorting_type = options[option_index]

# Start main program
main(bed_file_name=bed_file, genome_file_name=genome_file, sort_cluster_type=sorting_type)
