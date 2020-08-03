#!/usr/bin/env python

import matplotlib
matplotlib.use('Cairo')
import re
import matplotlib.pyplot as plt
import argparse
import pandas as pd
import sys

def plot_hist(filename, fln, plot_perc, divider, log, hist_type, bins,
              col_array, x_label, y_label):

    if col_array:
        data = pd.read_csv(filename, header=None)
    else:
        data = pd.read_csv(filename)

    if plot_perc:
        if divider != None:
            ct = 0
            for i in data.columns.values:
                if ct != 0 and i != divider:
                    data[i] = data[i]/data[divider]
                ct += 1
            data = data.drop(columns=[divider])
        else:
            print("--plot-perc provided but --divider is not!")
            sys.exit(1)

    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111)
    ax.grid()

    if col_array:
        values = [i for i in data.columns.values]
        heights = [i for i in data.values[0]]
        plt.bar(values, heights, width=0.8, log=log)

    else:
        columns = data.columns.values
        if len(columns) == 1:
            try:
                values = [i[0] for i in data.values]
                bins = int(max(values))
                if bins == 1:
                    plt.hist(values, bins=100, width=0.005, log=log)
                else:
                    plt.hist(values, bins=bins-4, width=0.8, log=log)
            except:
                legend = ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
        else:
            first_col = data.columns.values[0]
            data = data.set_index(first_col)
            ax = data.plot.bar()
            legend = ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)

    if x_label:
        ax.set_xlabel(x_label)
    else:
        ax.set_xlabel('')
    if y_label:
        ax.set_ylabel(y_label)
    else:
        ax.set_ylabel('')

    plt.savefig(fln, format="png", bbox_inches='tight', dpi=300)
    plt.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Plot histogram. Sample names ' +
        'should be in a first column. One data column also can be used.')
    parser.add_argument('--input', dest='filename', type=str,
                        help='input csv file with data array')
    parser.add_argument('--output', dest='out_file', type=str,
                        help='Output file.')
    parser.add_argument('--plot-perc', dest='plot_perc', action='store_true',
                        help='Plot percentaged. Provide --divider. All columns will be devided.')
    parser.add_argument('--divider', dest='divider', type=str,
                        help='Column name of divider.')
    parser.add_argument('--log', dest='log_scale', action='store_true',
                        help='Use to plot in log scale.')
    parser.add_argument('--hist-type', dest='hist_type', type=str, default='bar',
                        help='Matplotlib hist(histtype=["bar", "barstacked", "step", "stepfilled"]).')
    parser.add_argument('--bins', dest='bins', type=int, default=100,
                        help='Matplotlib hist bins.')
    parser.add_argument('--use-columns-as-array', dest='col_array', action='store_true',
                        help='Use to plot columns as an array.')
    parser.add_argument('--x-label', dest='x_label', type=str,
                        help='Label for x axis.')
    parser.add_argument('--y-label', dest='y_label', type=str,
                        help='Label for y axis.')

    args = parser.parse_args()
    plot_hist(args.filename, args.out_file,
        args.plot_perc, args.divider, args.log_scale,
        args.hist_type, args.bins, args.col_array,
        args.x_label, args.y_label)
