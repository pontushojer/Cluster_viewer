

import argparse
import pandas as pd
from sys import exit


def main(info_csv, bed_file_name, filter_strings, output_number, output_dir, lowest_values, thresholds, output_name,
         cluster_header):
    output_file = output_dir + '/' + output_name

    info_df = pd.DataFrame.from_csv(info_csv, index_col=None)

    print('Clusters at start:', len(info_df))
    if thresholds:
        for filt_s, thres in zip(filter_strings, thresholds):
            if '<' in thres:
                info_df = info_df[info_df[filt_s] < int(thres.replace('<', ''))]
            elif '>' in thres:
                info_df = info_df[info_df[filt_s] > int(thres.replace('>', ''))]
            elif '=' in thres:
                info_df = info_df[info_df[filt_s] == int(thres.replace('=', ''))]
            elif 'none' == thres:
                continue
            else:
                print('Wrong threshold given for: ', thres)
                exit()
    print('Clusters after filtering: ', len(info_df))

    if not lowest_values:
        lowest_values = [False]*len(filter_strings)

    info_df = info_df.sort_values(by=filter_strings, ascending=lowest_values)

    output_clusters = info_df[cluster_header].values
    output_clusters = output_clusters[:output_number]
    output_clusters = {str(output_cluster) for output_cluster in output_clusters}

    print('Clusters to output:', len(output_clusters))
    count = 0
    clusters = set()

    with open(bed_file_name, 'r') as bed_file:
        with open(output_file, 'a') as output:
            for line in bed_file.readlines():
                cluster = line.strip().split('\t')[3]
                clusters.add(cluster)
                if cluster in output_clusters:
                    count += 1
                    output.write(line)
    print('Clusters in bedfile: ', len(clusters))
    print('Clusters in bedfile to output:', len(output_clusters.intersection(clusters)))
    print('Fragment written to output: ', count)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Reads cluster information file and bed file containing reads to reduce size',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('cluster_info',
                        help='CSV file containing information relevant to clusters',
                        default=False)

    parser.add_argument('cluster_bed',
                        help='BED file containing fragments with cluster id on fourth column.',
                        default=False)

    parser.add_argument('-f', '--filter',
                        help='State which CSV header to filer on.',
                        default=False, required=True, nargs='+')

    parser.add_argument('-t', '--threshold',
                        help='Set threshold for equal e.g =3 less than e.g <3 or more than value e.g >3. If no threshold set to "none".',
                        nargs='+')

    parser.add_argument('-n', '--number',
                        help='Number of clusters to limit to',
                        default=100, type=int)

    parser.add_argument('-l', '--lowest',
                        help='Which values to output for the filter parameter.',
                        nargs='+')

    parser.add_argument('-o', '--output-name', help='Set name for output file', type=str, required=True)

    parser.add_argument('-c', '--cluster-header', help='Set cluster header for csv-file', type=str, default='Cluster id')

    parser.add_argument('--outdir',
                        help='Path to output dir',
                        default='/Users/pontushojer/data_analysis/scripts/python_scripts/Cluster_viewer/data/', type=str)

    args = parser.parse_args()

    print('INPUTS: ', args)

    main(args.cluster_info, args.cluster_bed, args.filter, args.number, args.outdir, args.lowest, args.threshold,
         args.output_name, args.cluster_header)
