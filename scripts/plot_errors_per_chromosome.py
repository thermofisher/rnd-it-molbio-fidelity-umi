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

def unique(a1, a2):

    a3 = np.array([])
    a4 = np.array([])

    idx = 0
    for i in a1:
        if i not in a3:
            a3 = np.append(a3, i)
            a4 = np.append(a4, a2[idx])

        idx += 1
    return a3, a4

def plot(filename, fln, x_name, y_name, z_name, sample_name, repeat_name, chromosome_name, x_label, y_label, s, kwargs):

    data = pd.read_csv(filename)

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

    total_chroms = len(data[chromosome_name].unique())

    fig = plt.figure(figsize=(14, 10))
    data['frequency'] = data[z_name]/data[y_name]
    data['reference'] = data[x_name].str[0]
    data['mutation'] = data[x_name].str[-1]
    data['position'] =  data[x_name].str.replace(r"[\-a-zA-Z.*]", "").astype(int)
    csv_out = open(fln.replace('.png', '')+'_graph_data.csv', 'w')
    dfout = pd.DataFrame({'Sample': [], 'Chromosome': [], 'X': [], 'Y': [], 'STD': []})
    dfout.to_csv(csv_out)

    chr_counter = 0

    # group tables by chromosomes
    for chr, chr_group in data.groupby(chromosome_name):
        chr_counter += 1
        ax = fig.add_subplot(total_chroms, 1, chr_counter)
        ax.grid()
        ax.title.set_text(f'{chr} chromosome')

        # group tables by samples
        for name, group in chr_group.groupby(sample_name):
            # if repeats are provided in the table
            if repeat_name:
                # initiate "empty" numpy arrays. Actualy it has random first value. Later will be removed
                freqeuncies = {'position': np.empty(1, dtype=int), 'frequency': np.empty(1, dtype=float)}
                # group table by repeats in samples
                for r_name, repeat_group in group.groupby(repeat_name):
                    # agregate values by positions in chromosome since many mutations at the same position are posible
                    grouped_by_position = repeat_group.groupby(['position']).agg({y_name: sum, z_name: np.mean})
                    # calculate frequency of an error for each position
                    grouped_by_position['frequency'] = grouped_by_position[y_name] / grouped_by_position[z_name]
                    # reset index to find position  column by key
                    grouped_by_position = grouped_by_position.reset_index()

                    # collect information for further grouping
                    for i in grouped_by_position['position']:
                        freqeuncies['position'] = np.append(freqeuncies['position'], i)
                    for i in grouped_by_position['frequency']:
                        freqeuncies['frequency'] = np.append(freqeuncies['frequency'], i)

                # remove first numpy random value
                freqeuncies['frequency'] = np.delete(freqeuncies['frequency'], 0)
                freqeuncies['position'] = np.delete(freqeuncies['position'], 0)
                freqeuncies = pd.DataFrame(freqeuncies, dtype=float)

                # agregate data by position from all repeats and compute mean frequency and std for each
                ag_repeats = freqeuncies.groupby(['position']).agg(
                    frequency=pd.NamedAgg(column='frequency', aggfunc=np.mean),
                    std=pd.NamedAgg(column='frequency', aggfunc=np.std)
                )
                ag_repeats = ag_repeats.reset_index()
                xser = ag_repeats['position']
                yser = ag_repeats['frequency']
                errser = ag_repeats['std']
                dfout = pd.DataFrame({'Sample': name, 'Chromosome': chr, 'X': xser, 'Y': yser, 'STD': errser})
                dfout.to_csv(csv_out, header=False)
                plt.errorbar(xser, yser, yerr=errser, label=name, **kwargs)

            else:
                grouped_by_position = group.groupby(['position']).agg({y_name: sum, z_name: np.mean})
                grouped_by_position['frequency'] = grouped_by_position[y_name] / grouped_by_position[z_name]
                grouped_by_position = grouped_by_position.reset_index()
                xser = grouped_by_position['position']
                yser = grouped_by_position['frequency']
                dfout = pd.DataFrame({'Sample': name, 'Chromosome': chr, 'X': xser, 'Y': yser})
                dfout.to_csv(csv_out, header=False)

                plt.scatter(xser, yser, label=name, s=s, **kwargs)
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

    plt.savefig(fln, format="png", bbox_inches="tight", dpi=300)
    plt.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Plot errors per reference chromosome."
    )
    parser.add_argument(
        "--input",
        dest="filename",
        type=str,
        required=True,
        help="input csv file with Sample,Chromosome,ClustOfVariants,Count,Coverage,ConfidentCoverage columns",
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
        help="Column name coresponding to the Y data set which is errors found.",
    )
    parser.add_argument(
        "--z-name",
        dest="z_name",
        type=str,
        help="Column name coresponding to the Z data set which is confident bases.",
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
        help="Column name coresponding to the chromosome names column.",
    )
    parser.add_argument(
        "--repeat-name",
        dest="repeat_name",
        type=str,
        help="Column name coresponding to the repeat names column.",
    )
    parser.add_argument(
        "--scatter-size",
        dest="scatter_size",
        type=float,
        default=1
    )
    parser.add_argument(
        "--matplotlib-kwargs",
        dest="kwargs",
        type=str,
        nargs='+'
    )

    parser.add_argument("--x-label", dest="x_label", type=str, help="Label for x axis.")
    parser.add_argument("--y-label", dest="y_label", type=str, help="Label for y axis.")

    args = parser.parse_args()
    plot(
        args.filename,
        args.out_file,
        args.x_name,
        args.y_name,
        args.z_name,
        args.sample_name,
        args.repeat_name,
        args.chromosome_name,
        args.x_label,
        args.y_label,
        args.scatter_size,
        args.kwargs
    )
