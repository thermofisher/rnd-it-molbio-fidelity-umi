#!/usr/bin/env python

import matplotlib

matplotlib.use("Cairo")
import matplotlib.pyplot as plt
import argparse
import pandas as pd
import scipy.stats as ss
import numpy as np
from scipy.interpolate import make_interp_spline, BSpline
from scipy.signal import savgol_filter


def plot_line(filename, fln, x_name, y_name, z_name, sample_name, chromosome_name, x_label, y_label, colors, kwargs):

    if kwargs is None:
        kwargs = {}
    else:
        kwargs = dict(x.split('=') for x in kwargs)

    for kwarg in kwargs:
        try:
            kwargs[kwarg] = int(kwargs[kwarg])
        except:
            try:
                kwargs[kwarg] = float(kwargs[kwarg])
            except:
                pass

    if colors:
        colors = dict(x.split('=') for x in colors)

    data = pd.read_csv(filename)

    fig = plt.figure(figsize=(14, 10))
    total_chroms = len(data[chromosome_name].unique())

    data['frequency'] = data[z_name]/data[y_name]
    data['reference'] = data[x_name].str[0]
    data['mutation'] = data[x_name].str[-1]
    data['position'] =  data[x_name].str.replace(r"[\-a-zA-Z.*]", "").astype(int)
    csv_out = open(fln.replace('.png', '')+'_graph_data.csv', 'w')

    def ecdf(data):
        """ Compute ECDF """
        x = np.sort(data)
        n = x.size
        y = np.arange(1, n+1) / n
        return(x,y)

    dfout = pd.DataFrame({'Sample': [], 'X': [], 'Y': []})
    dfout.to_csv(csv_out)

    chr_counter = 0
    for chr, chr_group in data.groupby(chromosome_name):
        chr_counter += 1
        ax = fig.add_subplot(total_chroms, 1, chr_counter)
        ax.grid()
        ax.title.set_text(f'{chr} chromosome')

        # group tables by samples
        for name, group in chr_group.groupby(sample_name):
            grouped_by_position = group.groupby(['position']).agg({y_name: sum, z_name: np.mean})
            grouped_by_position['frequency'] = grouped_by_position[y_name] / grouped_by_position[z_name]
            grouped_by_position = grouped_by_position.reset_index()
            yser = grouped_by_position['frequency']
            xser, yser = ecdf(yser)

            if colors:
                for i in colors:
                    if i in name:
                        kwargs['color'] = colors[i]

            plt.plot(xser, yser, label=name, **kwargs)

            if x_label:
                ax.set_xlabel(x_label)
            else:
                ax.set_xlabel("")
            if y_label:
                ax.set_ylabel(y_label)
            else:
                ax.set_ylabel("")

            if chr_counter == 1:
                legend = ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0)

            dfout = pd.DataFrame({'Sample': name, 'X': xser, 'Y': yser})
            dfout.to_csv(csv_out, header=False)

    legend = ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0)
    if x_label:
        ax.set_xlabel(x_label)
    else:
        ax.set_xlabel("")
    if y_label:
        ax.set_ylabel(y_label)
    else:
        ax.set_ylabel("")

    plt.savefig(fln, format="png", bbox_inches="tight", dpi=300)
    plt.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Plot line. Sample names "
        + "should be in a first column. One data column also can be used."
    )
    parser.add_argument(
        "--input",
        dest="filename",
        type=str,
        required=True,
        help="input csv file with data array",
    )
    parser.add_argument(
        "--output", dest="out_file", type=str, required=True, help="Output file."
    )
    parser.add_argument(
        "--x-name",
        dest="x_name",
        type=str,
        required=True,
        help="Column name coresponding to the X data set.",
    )
    parser.add_argument(
        "--y-name",
        dest="y_name",
        type=str,
        help="Column name coresponding to the Y data set.",
    )
    parser.add_argument(
        "--z-name",
        dest="z_name",
        type=str,
        help="Column name coresponding to the Z data set.",
    )
    parser.add_argument(
        "--sample-name",
        dest="sample_name",
        type=str,
        required=True,
        help="Column name coresponding to the sample name column.",
    )
    parser.add_argument(
        "--chromosome-name",
        dest="chromosome_name",
        type=str,
        required=True,
        help="Column name coresponding to the chromosome name column.",
    )
    parser.add_argument(
        "--matplotlib-kwargs",
        dest="kwargs",
        type=str,
        nargs='+'
    )
    parser.add_argument(
        "--color-by-sample",
        dest="colors",
        type=str,
        nargs='+',
        help="Colors for samples: Sample=color.",
    )
    parser.add_argument("--x-label", dest="x_label", type=str, help="Label for x axis.")
    parser.add_argument("--y-label", dest="y_label", type=str, help="Label for y axis.")

    args = parser.parse_args()
    plot_line(
        args.filename,
        args.out_file,
        args.x_name,
        args.y_name,
        args.z_name,
        args.sample_name,
        args.chromosome_name,
        args.x_label,
        args.y_label,
        args.colors,
        args.kwargs
    )
