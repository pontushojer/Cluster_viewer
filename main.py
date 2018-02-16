import argparse
import pandas as pd
import numpy as np
from bokeh.models import ColumnDataSource, HoverTool, Button, Range1d, TextInput, Select, FixedTicker, Paragraph
from bokeh.models import DataTable, NumberFormatter, TableColumn, RangeSlider
from bokeh.plotting import figure
from bokeh.layouts import layout
from bokeh.models.callbacks import CustomJS
from bokeh.io import curdoc, show
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

    return chrom_order, translate_key, genome_lens


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
                                  int(end) + translate_key[chrom],
                                  int(end) - int(start),
                                  (int(end) + int(start))/2 + translate_key[chrom])

                cluster_regions.append(cluster_region)

    return pd.DataFrame(data=cluster_regions, columns=['Cluster', 'Chromosome', 'Start', 'End', 'xmin', 'xmax', 'width',
                                                       'xmid'])


def get_cluster_info(dataframe, sort_type, phasing_threshold=100000):
    """
    Get basic cluster information from dataframe
    """

    cluster_data = []
    # Loop over each cluster
    for cluster in dataframe.Cluster.unique():
        # Extract dataframe for cluster
        cluster_dataframe = dataframe[dataframe.Cluster == cluster]

        # Get number of fragments
        fragments = len(cluster_dataframe)

        # Get number of chromosomes
        chromosomes = len(cluster_dataframe.Chromosome.unique())

        # Get number of phased reads i.e pairs of neigbouring fragments within threshold
        phasing = 0
        if sort_type == 'Phased_regions':
            if fragments > 1:
                for chromosome in cluster_dataframe.Chromosome.unique():
                    fragment_start_values = cluster_dataframe[cluster_dataframe.Chromosome == chromosome].Start
                    if len(fragment_start_values) > 1:
                        for start1, start2 in zip(fragment_start_values[:-1], fragment_start_values[1:]):
                            if abs(start1 - start2) <= phasing_threshold:
                                phasing += 1

        # Store as tuple
        if sort_type == 'Phased_regions':
            data = (cluster, fragments, chromosomes, phasing)
        else:
            data = (cluster, fragments, chromosomes)

        # Apped to list
        cluster_data.append(data)

    if sort_type == 'Phased_regions':
        return pd.DataFrame(data=cluster_data, columns=['Cluster', 'Chromosomes', 'Fragments', 'Phased_regions'])
    else:
        return pd.DataFrame(data=cluster_data, columns=['Cluster', 'Chromosomes', 'Fragments'])


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
    dataframe['Index_list'] = pd.Series((np.zeros(len(dataframe)), 2))

    for i, df_row in dataframe.iterrows():
        dataframe.set_value(i, 'Index', clusters_index[df_row.Cluster])
        dataframe.set_value(i, 'Index_list', np.array((clusters_index[df_row.Cluster], clusters_index[df_row.Cluster])))

    return dataframe.sort_values(by='Index'), clusters_index


def main(bed_file_name, genome_file_name, sort_cluster_type):
    ########################################################
    # Get initial information from bed file and genome file
    ########################################################



    print(strftime("%Y-%m-%d %H:%M:%S ", localtime()), "Reading genome file: ", genome_file_name)

    chrom_order, translate_key, chrom_lens = genome_file_handling(genome_file_name)

    print(strftime("%Y-%m-%d %H:%M:%S ", localtime()), "Creating cluster dataframe from bed file: ", bed_file_name)

    cluster_df = get_cluster_df(bed_file_name, translate_key, chrom_order)

    print(strftime("%Y-%m-%d %H:%M:%S ", localtime()), "Collecting cluster info to dataframe")

    # Get informatio about cluster based on sort type. Will only gather phasing information relevant to sort.
    cluster_info = get_cluster_info(cluster_df, sort_cluster_type)

    print(strftime("%Y-%m-%d %H:%M:%S ", localtime()), "Sorting cluster dataframe based on: ", sort_cluster_type)

    cluster_df_sorted, cluster_index = sort_dataframe(cluster_df, cluster_info, sort_cluster_type)

    # Convert cluster type to string to enable factorial representation.
    cluster_df_sorted.Cluster = cluster_df_sorted.Cluster.astype(str)

    print(strftime("%Y-%m-%d %H:%M:%S ", localtime()), "Setting up data source for plotting.")

    source = ColumnDataSource(cluster_df_sorted)

    ###########
    # Plotting
    ###########
    min_visable_range = 5000  # To power requirements
    max_range = sum(chrom_lens.values())
    tools = "pan,wheel_zoom,box_zoom,reset,save"

    p = figure(x_range=(0, max_range),
               plot_width=1400,
               plot_height=750,
               toolbar_location='above',
               tools=tools)

    cluster_bar = p.rect(x=max_range/2,
                         width=max_range,
                         height=0.5,
                         y=list(cluster_index.values()),
                         fill_alpha=0.01,
                         line_alpha=0.01,
                         hover_line_color='red',
                         hover_fill_color='red',
                         hover_line_alpha=0.1,
                         hover_fill_alpha=0.1)

    glyph = p.rect(x='xmid',
                    y='Index',
                    width='width',
                    height=0.5,
                    source=source,
                    fill_alpha=0.3,
                    line_alpha=0.3,
                    line_width=1,
                    fill_color='black',
                    line_color='black',
                    hover_line_color='red',
                    hover_fill_color='red',
                    hover_line_alpha=0.5,
                    hover_fill_alpha=0.5)

    hovertool = HoverTool(mode='mouse',
                          tooltips=[("Position", "(@Chromosome: @Start{0,0} - @End{0,0})")],
                          renderers=[glyph],
                          line_policy='next')

    hovertool2 = HoverTool(mode='mouse', tooltips=None,
                           renderers=[cluster_bar])

    p.add_tools(hovertool, hovertool2)

    p.ygrid.grid_line_color = 'blue'
    p.ygrid.grid_line_width = 15
    p.ygrid.grid_line_alpha = 0.1
    p.ygrid.ticker = FixedTicker(ticks=np.arange(0, len(cluster_index), 1))
    p.xaxis.axis_label = None
    p.yaxis.axis_label = 'Cluster id'

    p.y_range = Range1d(-1, 20)
    p.yaxis.ticker = list(range(0, len(cluster_index)))
    p.yaxis.major_label_overrides = {str(i): '#' + str(key) + '(' + str(i).rjust(3) + ')'
                                     for key, i in cluster_index.items() if i in range(0, len(cluster_index))}

    p.outline_line_color = None

    p.xaxis.ticker = [translate_key[chrom] for chrom in chrom_order]
    p.xaxis.major_label_overrides = {translate_key[chrom]: chrom for chrom in chrom_order}
    p.xaxis.major_label_orientation = 3.14/4

    p.xgrid.grid_line_color = None

    #p.title.text = "Cluster overview " #+ str(current_pos) + ' - ' + str(current_pos + 20) + ' of total (' + str(len(cluster_index)) + ')'

    print(strftime("%Y-%m-%d %H:%M:%S ", localtime()), "Plotting figure.")

    #show(p)
    #sys.exit()

    ###########################
    # Interactivity section
    #########################

    def update_ticks(chrom, start, end):
        if chrom == 'All':
            # Ticks displayed
            p.xaxis.ticker = [translate_key[chrom] for chrom in chrom_order]
            p.xaxis.major_label_overrides = {translate_key[chrom]: chrom for chrom in chrom_order}
            p.xaxis.axis_label = None
        else:
            # Ticker setup
            window = end - start
            step = int(float(0.5 * 10 ** (len(str(window)) - 1)))
            ticks = list(range(translate_key[chrom],
                               translate_key[chrom] + chrom_lens[chrom], step))
            tick_labels = list(range(0, chrom_lens[chrom], step))

            # Ticks displayed
            p.xaxis.ticker = ticks
            p.xaxis.major_label_overrides = {tick: '{:,}'.format(label) for tick, label in zip(ticks, tick_labels)}
            p.xaxis.axis_label = chrom

        x_range_start.value = str(start)
        x_range_end.value = str(end)

    #
    # Select from which index to display clusters.
    #
    y_range_index = TextInput(value="0", title="Range index:")

    def update_y_range(attr, old, new):
        p.y_range.start = int(y_range_index.value) - 1
        p.y_range.end = int(y_range_index.value) + 20

    y_range_index.on_change('value', update_y_range)

    #
    # Select which chromosome to display, All refers to base state with all chroms visable
    #
    chr_options = ['All'] + chrom_order
    chr_select = Select(value='All', title='Select chromosome', options=chr_options)

    def update_chr(attr, old, new):
        if chr_select.value == 'All':
            # Range displayed
            p.x_range.start = 0
            p.x_range.end = max_range

        else:
            # Range displayed
            p.x_range.start = translate_key[chr_select.value]
            p.x_range.end = translate_key[chr_select.value] + chrom_lens[chr_select.value]

        update_ticks(chr_select.value, p.x_range.start, p.x_range.end)

    chr_select.on_change('value', update_chr)

    #
    # Select which regions to display based on input coordinates
    #
    x_range_start = TextInput(value='0', title='Start x:')
    x_range_end = TextInput(value=str(max_range), title='End x:')

    def update_x_range(attr, old, new):
        if int(x_range_end.value) - int(x_range_start.value) > min_visable_range:
            # Range displayed
            p.x_range.start = int(float(x_range_start.value))
            p.x_range.end = int(float(x_range_end.value))

            if chr_select.value != 'All':
                update_ticks(chr_select.value, p.x_range.start, p.x_range.end)
        else:
            print('Zoom in limit reached.')
    x_range_start.on_change('value', update_x_range)
    x_range_end.on_change('value', update_x_range)

    #
    # Zoom in button
    #
    zoom_in = Button(label='ZOOM IN', button_type='success')

    def update_zoom_in():
        diff = p.x_range.end - p.x_range.start
        if diff/2 > min_visable_range:
            p.x_range.start += int(float(diff/4))
            p.x_range.end -= int(float(diff / 4))
            update_ticks(chr_select.value, p.x_range.start, p.x_range.end)
        else:
            print('Zoom in limit reached.')

    zoom_in.on_click(update_zoom_in)

    #
    # Zoom out button
    #
    zoom_out = Button(label='ZOOM OUT', button_type='success')

    def update_zoom_out():
        diff = p.x_range.end - p.x_range.start
        chromosome = chr_select.value

        if p.x_range.start - int(float(diff / 4)) > 0:
            p.x_range.start -= int(float(diff / 4))
        else:
            if chromosome == 'All':
                p.x_range.start = 0
            else:
                p.x_range.start = translate_key[chromosome]

        if chromosome == 'All':
            if p.x_range.end + int(float(diff / 4)) < max_range:
                p.x_range.end += int(float(diff / 4))
            else:
                p.x_range.end = max_range

        else:
            if p.x_range.end + int(float(diff / 4)) < translate_key[chromosome] + chrom_lens[chromosome]:
                p.x_range.end += int(float(diff / 4))
            else:
                p.x_range.end = translate_key[chromosome] + chrom_lens[chromosome]

        update_ticks(chromosome, p.x_range.start, p.x_range.end)

    zoom_out.on_click(update_zoom_out)

    #
    # Cluster specific table + info
    #

    start_cluster = cluster_df_sorted[cluster_df_sorted.Index == 0].Cluster.values[0]

    text = Paragraph(text=str(start_cluster), width=250, height=10)

    cluster_data = cluster_df_sorted[cluster_df_sorted.Cluster == start_cluster]

    cluster_source = ColumnDataSource(cluster_data)

    columns = [TableColumn(field='Chromosome', title='Chrom'),
               TableColumn(field='Start', title='Start', formatter=NumberFormatter()),
               TableColumn(field='End', title='End', formatter=NumberFormatter())]

    data_table = DataTable(source=cluster_source, columns=columns, width=250, height=750, row_headers=False,
                           fit_columns=True, sortable=True, editable=False)

    cluster_select = TextInput(value=start_cluster, title='Cluster:', width=100)

    def update_cluster_info(attr, old, new):
        if int(cluster_select.value) in cluster_index.keys():
            text.text = 'Cluster: ' + cluster_select.value
            new_data = cluster_df_sorted[cluster_df_sorted.Cluster == cluster_select.value]
            cluster_source.data = {
                'Chromosome': new_data.Chromosome,
                'Start': new_data.Start,
                'End': new_data.End}
        else:
            text.text = 'Cluster: ' + cluster_select.value + ' not available!'
            print('Cluster ', cluster_select.value, 'not in data')

    cluster_select.on_change('value', update_cluster_info)

    #
    # HTML Document setup
    #
    curdoc().add_root(layout([[[p], [text, data_table]],
                              [zoom_in, zoom_out],
                              [y_range_index, chr_select, x_range_start, x_range_end, cluster_select]]))
    curdoc().title = "CLUSTER VIEWER"

#################
# Startup options
#################

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


# print('AVAILABLE SORT OPTIONS:')
# print('-'*30)
# options = ['Cluster', 'Fragments', 'Chromosomes', 'Phased_regions']
# for nr, option in enumerate(options):
#     print(nr, ':', option)
# print('-'*30)
#
# option_index = int(input('Give selected sort option index:'))
# sorting_type = options[option_index]
sorting_type = 'Cluster'

# Start main program
main(bed_file_name=bed_file, genome_file_name=genome_file, sort_cluster_type=sorting_type)
