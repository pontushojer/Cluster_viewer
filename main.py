import argparse
import pandas as pd
import numpy as np
from bokeh.models import ColumnDataSource, HoverTool, Button, Range1d, TextInput, Select, FixedTicker, Paragraph
from bokeh.models import DataTable, NumberFormatter, TableColumn
from bokeh.plotting import figure
from bokeh.layouts import layout
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


def neared_cluster_bin(dataframe, chroms, chroms_lenghts, bin_size=1000000):
    clusters = list(dataframe.Cluster.unique())
    translate_key = {}
    bins = 0

    for nr, chrom in enumerate(chroms):
        chrom_bins = int(chroms_lenghts[chrom] / bin_size + 1)
        translate_key[chrom] = bins
        bins += chrom_bins

    print(get_time(), 'Creating matrix for clustering')
    cluster_dict = {}
    bin_matrix = []
    count = 0
    for cluster in clusters:

        cluster_dataframe = dataframe[dataframe.Cluster == cluster]
        bin_array = [0] * bins
        for i, row in cluster_dataframe.iterrows():
            cluster_bin = translate_key[row.Chromosome] + int(row.Start / bin_size)
            bin_array[cluster_bin] += 1

        array_norm = np.array([bin_freq / sum(bin_array) for bin_freq in bin_array])
        bin_matrix.append(array_norm)
        cluster_dict[cluster] = count

        count += 1
        if count % 1000 == 0:
            print(get_time(), count, 'clusters processed.')

    bin_matrix = np.array(bin_matrix)

    import scipy.cluster.hierarchy as hcluster

    hclusters = hcluster.fclusterdata(bin_matrix, 0.1, criterion='distance', metric='euclidean')
    hclsts = [(i, hclst) for i, hclst in enumerate(hclusters)]
    clusters_index = {}

    for nr, (i, hclst) in enumerate(sorted(hclsts, key=lambda x: list(hclusters).count(x[1]), reverse=True)):
        clusters_index[clusters[i]] = nr



    # nbrs = NearestNeighbors(n_neighbors=5, algorithm='ball_tree').fit(bin_matrix)
    # dist, indices = nbrs.kneighbors(bin_matrix)

    # for i, index_list in enumerate(indices):
    #     if i != index_list[0]:

    # print(dist)
    # print(indices)

    #
    # clusters_index = {}
    # cluster = clusters.pop()
    # clusters_index[cluster] = 0
    # index = 0
    #
    # print(get_time(), 'Finding nearest neighbour.')
    # iter = 0
    # while len(clusters) > 1:
    #     iter += 1
    #     index += 1
    #     row = cluster_dict[cluster]
    #     dist = 1
    #     neighbour = None
    #
    #     for comp_cluster in clusters:
    #         comp_row = cluster_dict[comp_cluster]
    #         comp_dist = np.linalg.norm(row - comp_row)
    #
    #         if comp_dist == 0:
    #             neighbour = comp_cluster
    #             break
    #         elif comp_dist < dist:
    #             dist = comp_dist
    #             neighbour = comp_cluster
    #
    #     if neighbour:
    #         if dist < 0.01:
    #             clusters_index[neighbour] = index
    #
    #             cluster = neighbour
    #             clusters.remove(cluster)
    #         else:
    #             clusters.append(cluster)
    #             cluster = clusters.pop(randint(0, len(clusters) - 1))
    #             index -= 1
    #
    #     if iter % 1000 == 0:
    #         print(get_time(), iter, 'iteraction. Clusters remaining: ', len(clusters), 'Current cluster= ', cluster)
    #
    # clusters_index[clusters[0]] = index + 1

    dataframe['Index'] = pd.Series((np.zeros(len(dataframe)), 1))

    for i, df_row in dataframe.iterrows():
        # print(df_row.Cluster)
        try:
            dataframe.set_value(i, 'Index', clusters_index[df_row.Cluster])
        except KeyError:
            print('Keyerror', df_row.Cluster)
            print('In clusters', df_row.Cluster in list(dataframe.Cluster.unique()))
            print('In clusters index', df_row.Cluster in clusters_index.keys())
    return dataframe.sort_values(by='Index'), clusters_index


def neared_cluster_chrom(dataframe, chroms):
    clusters = list(dataframe.Cluster.unique())
    cluster_dict = {}

    for cluster in clusters:
        chrom_data = [0] * len(chroms)
        cluster_chroms = dataframe[dataframe.Cluster == cluster].Chromosome.values

        for chrom in cluster_chroms:
            chrom_data[chrom_number(chrom)-1] += 1

        chrom_data_norm = np.array([data / sum(chrom_data) for data in chrom_data])
        cluster_dict[cluster] = chrom_data_norm

    clusters_index = {}
    cluster = clusters.pop()
    clusters_index[cluster] = 0
    index = 0

    while len(clusters) > 1:
        index += 1
        row = cluster_dict[cluster]
        dist = 1
        neighbour = None

        for comp_cluster in clusters:
            comp_row = cluster_dict[comp_cluster]
            comp_dist = np.linalg.norm(row - comp_row)

            if comp_dist == 0:
                neighbour = comp_cluster
                break
            elif comp_dist < dist:
                dist = comp_dist
                neighbour = comp_cluster

        if neighbour:
            clusters_index[neighbour] = index

            cluster = neighbour
            clusters.remove(cluster)
        else:
            print(1, cluster)
    clusters_index[clusters[0]] = index + 1

    dataframe['Index'] = pd.Series((np.zeros(len(dataframe)), 1))

    for i, df_row in dataframe.iterrows():
        #print(df_row.Cluster)
        try:
            dataframe.set_value(i, 'Index', clusters_index[df_row.Cluster])
        except KeyError:
            print('Keyerror', df_row.Cluster)
            print('In clusters', df_row.Cluster in list(dataframe.Cluster.unique()))
            print('In clusters index', df_row.Cluster in clusters_index.keys())
    return dataframe.sort_values(by='Index'), clusters_index


def get_cluster_info(dataframe, sort_type, phasing_threshold=100000):
    """
    Get basic cluster information from dataframe
    """

    cluster_data = []
    clusters = dataframe.Cluster.unique()
    # Loop over each cluster
    for cluster in clusters:
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
            data = (cluster, chromosomes, fragments, phasing)
        else:
            data = (cluster, chromosomes, fragments)
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


def get_time():
    return strftime("%Y-%m-%d %H:%M:%S ", localtime())


def main(bed_file_name, genome_file_name, sort_cluster_type):
    ########################################################
    # Get initial information from bed file and genome file
    ########################################################

    print(get_time(), "Reading genome file: ", genome_file_name)

    chrom_order, translate_key, chrom_lens = genome_file_handling(genome_file_name)

    print(get_time(), "Creating cluster dataframe from bed file: ", bed_file_name)

    cluster_df = get_cluster_df(bed_file_name, translate_key, chrom_order)

    print(get_time(), "Collecting cluster info to dataframe")

    # Get informatio about cluster based on sort type. Will only gather phasing information relevant to sort.
    cluster_info = get_cluster_info(cluster_df, sort_cluster_type)

    print(get_time(), "Sorting cluster dataframe based on: ", sort_cluster_type)

    if sorting_type == 'Neighbour':
        cluster_df_sorted, cluster_index = neared_cluster_bin(cluster_df, chrom_order, chrom_lens, bin_size=1000000)
        #cluster_df_sorted, cluster_index = neared_cluster_chrom(cluster_df, chrom_order)
    else:
        cluster_df_sorted, cluster_index = sort_dataframe(cluster_df, cluster_info, sort_cluster_type)

    # Convert cluster type to string to enable factorial representation.
    cluster_df_sorted.Cluster = cluster_df_sorted.Cluster.astype(str)

    print(get_time(), "Setting up data source for plotting.")

    source = ColumnDataSource(cluster_df_sorted)

    ###########
    # Plotting
    ###########
    min_visable_range = 500  # To power requirements
    max_range = sum(chrom_lens.values())
    visable_clusters = 20

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
                    line_width=2,
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

    p.y_range = Range1d(-1, visable_clusters)
    p.yaxis.ticker = list(range(0, len(cluster_index)))
    p.yaxis.major_label_overrides = {str(i): '#' + str(key) + '(' + str(i).rjust(3) + ')'
                                     for key, i in cluster_index.items() if i in range(0, len(cluster_index))}

    p.outline_line_color = None

    p.xaxis.ticker = [translate_key[chrom] for chrom in chrom_order]
    p.xaxis.major_label_overrides = {translate_key[chrom]: chrom for chrom in chrom_order}
    p.xaxis.major_label_orientation = 3.14/4

    p.xgrid.grid_line_color = None

    p.title.text = "Cluster overview for data = " + str(bed_file_name)

    print(get_time(), "Plotting figure.")

    #show(p)
    #sys.exit()

    ###########################
    # Interactivity section
    #########################

    def update_ticks(chrom, start, end):
        print(get_time(), 'update_ticks')
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

            # Limit tick window to visable range + one more window.
            ticks__dict = {tick: '{:,}'.format(label) for tick, label in zip(ticks, tick_labels)
                           if tick + window > start and tick - window < end}

            # Ticks displayed
            p.xaxis.ticker = list(ticks__dict.keys())
            p.xaxis.major_label_overrides = ticks__dict
            p.xaxis.axis_label = chrom

            # if str(start - translate_key[chrom]) != x_range_start.value:
            #     x_range_start.value = str(start - translate_key[chrom])
            #
            # if str(end - translate_key[chrom]) != x_range_end.value:
            #     x_range_end.value = str(end - translate_key[chrom])

            if start - translate_key[chrom] != int(float(x_range_start.value.replace(',', ''))):
                x_range_start.value = '{:,}'.format(start - translate_key[chrom])

            if end - translate_key[chrom] != int(float(x_range_end.value.replace(',', ''))):
                x_range_end.value = '{:,}'.format(end - translate_key[chrom])

    #
    # Select from which index to display clusters.
    #
    y_range_index = TextInput(value="0", title="Range index:")

    def update_y_range(attr, old, new):
        print(get_time(), 'update_y_range')
        p.y_range.start = int(y_range_index.value) - 1
        p.y_range.end = int(y_range_index.value) + visable_clusters

    y_range_index.on_change('value', update_y_range)

    #
    # Select which chromosome to display, All refers to base state with all chroms visable
    #
    chr_options = ['All'] + chrom_order
    chr_select = Select(value='All', title='Select chromosome', options=chr_options)

    def update_chr(attr, old, new):
        print(get_time(), 'update_chr')
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
        print(get_time(), 'update_x_range')

        start = int(float(x_range_start.value.replace(',', '')))
        end = int(float(x_range_end.value.replace(',', '')))

        if end - start > min_visable_range:
            # Range displayed
            if chr_select.value != 'All':
                if p.x_range.start != start + translate_key[chr_select.value]:
                    p.x_range.start = start + translate_key[chr_select.value]

                if p.x_range.end != end + translate_key[chr_select.value]:
                    p.x_range.end = end + translate_key[chr_select.value]

                update_ticks(chr_select.value, p.x_range.start, p.x_range.end)
            else:
                p.x_range.start = start
                p.x_range.end = end

        else:
            print('Zoom in limit reached.')

        # if int(float(x_range_end.value)) - int(float(x_range_start.value)) > min_visable_range:
        #     # Range displayed
        #     if chr_select.value != 'All':
        #         if p.x_range.start != int(float(x_range_start.value)) + translate_key[chr_select.value]:
        #             p.x_range.start = int(float(x_range_start.value))
        #
        #         if p.x_range.end != int(float(x_range_end.value)) + translate_key[chr_select.value]:
        #             p.x_range.end = int(float(x_range_end.value))
        #
        #         update_ticks(chr_select.value, p.x_range.start, p.x_range.end)
        #     else:
        #         p.x_range.start = int(float(x_range_start.value))
        #         p.x_range.end = int(float(x_range_end.value))
        #
        # else:
        #     print('Zoom in limit reached.')

    x_range_start.on_change('value', update_x_range)
    x_range_end.on_change('value', update_x_range)

    #
    # Zoom in button
    #
    zoom_in = Button(label='ZOOM IN', button_type='success')

    def update_zoom_in():
        print(get_time(), 'update_zoom_in')
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
        print(get_time(), 'update_zoom_out')
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

    text = Paragraph(text='CLUSTER: ' + str(start_cluster) +
                          '. Frags = ' + str(cluster_info[cluster_info['Cluster'] == int(start_cluster)].Fragments.values[0]) +
                          ' Chroms = ' + str(cluster_info[cluster_info['Cluster'] == int(start_cluster)].Chromosomes.values[0]),
                     width=250, height=20)

    print(cluster_info.sort_values(by='Chromosomes').head(5))
    print(cluster_info.sort_values(by='Chromosomes').tail(5))

    cluster_data = cluster_df_sorted[cluster_df_sorted.Cluster == start_cluster]

    cluster_source = ColumnDataSource(cluster_data)

    columns = [TableColumn(field='Chromosome', title='Chrom'),
               TableColumn(field='Start', title='Start', formatter=NumberFormatter()),
               TableColumn(field='End', title='End', formatter=NumberFormatter())]

    data_table = DataTable(source=cluster_source, columns=columns, width=250, height=750, row_headers=False,
                           fit_columns=True, sortable=True, editable=False)

    cluster_select = TextInput(value=start_cluster, title='Cluster:', width=100)

    def update_cluster_info(attr, old, new):
        print(get_time(), 'update_cluster_info')
        print(cluster_info[cluster_info['Cluster'] == int(cluster_select.value)])
        if int(cluster_select.value) in cluster_index.keys():
            text.text = 'CLUSTER: ' + cluster_select.value + \
                        '. Frags = ' + str(cluster_info[cluster_info['Cluster'] == int(cluster_select.value)].Fragments.values[0]) + \
                        ' Chroms = ' + str(cluster_info[cluster_info['Cluster'] == int(cluster_select.value)].Chromosomes.values[0])
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
genome_file = '/Users/pontushojer/data_analysis/scripts/python_scripts/Cluster_viewer/data/hg38.genome'

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
options = ['Cluster', 'Fragments', 'Chromosomes', 'Phased_regions', 'Neighbour']
for nr, option in enumerate(options):
    print(nr, ':', option)
print('-'*30)

option_index = int(input('Give selected sort option index:'))
sorting_type = options[option_index]

# Start main program
main(bed_file_name=bed_file, genome_file_name=genome_file, sort_cluster_type=sorting_type)
