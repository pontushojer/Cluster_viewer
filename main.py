import argparse
import pandas as pd
import numpy as np
from bokeh.io import show
from bokeh.models import ColumnDataSource, HoverTool, Button, Range1d
from bokeh.plotting import figure
from bokeh.layouts import layout
from bokeh.io import curdoc
import sys


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


def genome_file_handling(genome_file):
    genome_lens = {}
    with open(genome_file, 'r') as genome:
        for entry in genome.readlines():
            chrom_name, chrom_len = entry.strip().split('\t')

            if any(s in chrom_name for s in ['chrM', 'chrEBV', 'random', 'chrUn']):
                continue

            genome_lens[chrom_name] = int(chrom_len)

    chrom_order = list(genome_lens.keys())
    chrom_order = sorted(chrom_order, key=lambda x: chrom_number(x))

    translate_key = {}
    current_len = 0

    for chrom in chrom_order:
        translate_key[chrom] = current_len
        current_len += genome_lens[chrom]

    return chrom_order, translate_key


def get_cluster_df(bed_file, translate_key, chrom_order):
    cluster_regions = []

    with open(bed_file, 'r') as bed:
        for interval in bed.readlines():

            chrom, start, end, cluster_id, weight, direction = interval.strip().split('\t')
            if chrom in chrom_order:
                cluster_region = (int(cluster_id), str(chrom), int(start), int(end), int(start) + translate_key[chrom], int(end) + translate_key[chrom])
                cluster_regions.append(cluster_region)

    return pd.DataFrame(data=cluster_regions, columns=['Cluster', 'Chromosome', 'Start', 'End', 'xmin', 'xmax'])


def sort_dataframe(dataframe, sort_type):
    if sort_type == 'Cluster':
        clusters = dataframe.Cluster.unique()
        clusters = sorted(clusters, key=lambda x: int(x))
        clusters_index = {cluster: i for i, cluster in enumerate(clusters)}

        dataframe['Index'] = pd.Series((np.zeros(len(dataframe)), 1))
        for i, df_row in dataframe.iterrows():
            dataframe.set_value(i, 'Index', clusters_index[df_row.Cluster])

    return dataframe.sort_values(by='Index'), clusters_index


def main(bed_file, genome_file):

    chrom_order, translate_key = genome_file_handling(genome_file)

    cluster_df = get_cluster_df(bed_file, translate_key, chrom_order)

    cluster_df_sorted, cluster_index = sort_dataframe(cluster_df, 'Cluster')

    cluster_df_sorted.Cluster = cluster_df_sorted.Cluster.astype(str)

    selection = cluster_df_sorted.loc[cluster_df_sorted['Index'].isin(range(0, 20))]
    #source = ColumnDataSource(selection)

    source = ColumnDataSource(cluster_df_sorted)
    #
    # Plotting
    #

    #output_file('test.html')
    p = figure(#y_range=sorted(list(set(source.data['Cluster'])), key=lambda x: int(x), reverse=True),
               plot_width=1500,
               plot_height=800,
               toolbar_location='right')

    # cr = p.hbar(y='Cluster',
    #             left='xmin',
    #             right='xmax',
    #             height=0.3,
    #             line_width=4,
    #             source=source,
    #             fill_alpha=0.5,
    #             line_alpha=0.5,
    #             hover_line_color='red',
    #             hover_fill_color='red',
    #             hover_line_alpha=0.9,
    #             hover_fill_alpha=0.9)

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

    current_pos = 0
    global current_pos
    p.y_range = Range1d(current_pos - 1, current_pos + 21)
    p.yaxis.ticker = list(range(current_pos, current_pos + 21))
    p.yaxis.major_label_overrides = {str(i): '#' + str(i) + ' = ' + str(key) for key, i in cluster_index.items() if i in range(current_pos, current_pos + 20)}

    p.outline_line_color = None

    p.xaxis.ticker = [translate_key[chrom] for chrom in chrom_order]
    p.xaxis.major_label_overrides = {translate_key[chrom]: chrom for chrom in chrom_order}
    p.xaxis.major_label_orientation = 3.14/2

    p.xgrid.grid_line_color = None

    p.title.text = "Cluster overview " + str(current_pos) + ' - ' + str(current_pos + 20) + ' of total (' + str(len(cluster_index)) + ')'

    show(p)

    button_fwd = Button(label='Fwd', button_type='success')
    button_rev = Button(label='Rev', button_type='success')

    def update_axis(change, cluster_index):
        global current_pos
        position = current_pos + change
        if 0 <= position:
            p.y_range.start += change
            p.y_range.end += change
            p.yaxis.ticker = list(range(position-1, position + 20))
            p.yaxis.major_label_overrides = {str(i): '#' + str(i) + ' = ' + str(key) for key, i in cluster_index.items() if
                                             i in range(position-1, position + 21)}
            current_pos = position
            p.title.text = "Cluster overview " + str(current_pos) + ' - ' + str(current_pos + 20) + ' of total (' + str(
                len(cluster_index)) + ')'

        print(current_pos)

    def cb_fwd():
        print("FORWARD")
        update_axis(20, cluster_index)

    def cb_rev():
        print("REVERSE")
        update_axis(-20, cluster_index)

    button_fwd.on_click(cb_fwd)
    button_rev.on_click(cb_rev)

    curdoc().add_root(layout([[p], [button_fwd, button_rev]]))
    curdoc().title = "Test"
#load_data(join(dirname(__file__), 'test.bed'))
#load_data(join(dirname(__file__), 'hg38.genome'))

bed_file = 'data/test.bed'
genome_file = 'data/hg38.genome'
main(bed_file=bed_file, genome_file=genome_file)
