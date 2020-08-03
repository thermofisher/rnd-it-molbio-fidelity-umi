if __name__ == '__main__':

    import argparse, os
    import pandas as pd
    import matplotlib
    matplotlib.use('Cairo')
    import matplotlib.pyplot as plt
    plt.style.use('bmh')

    parser = argparse.ArgumentParser(description='gets clean bases from bbbduk logs')
    parser.add_argument('--input', dest='input_files', type=str, nargs='+', help='input_file')
    parser.add_argument('--out', dest='out', type=str,  help='output_file for clean bases png')
    args = parser.parse_args()

    input_files = args.input_files
    out = args.out

    data = dict()
    info = ""
    for fl in input_files:
        smp = fl.split("/")[-1].replace("_trim.settings","")
        f = open(fl,'r')
        for l in f:
            if l.find("Result:") > -1:
                t_bases = float(l.split("\t")[2].split(" ")[-1].split(")")[0].replace("%","").replace("(",""))
                data[smp] = t_bases
        f.close()

    df=pd.DataFrame(data = list(data.items()), columns=['SAMPLE','CLEAN_BASES'])
    df=df.set_index('SAMPLE')

    width = 1
    plt.figure(figsize=(20, 20), dpi=300)
    font = {'family': 'DejaVu Sans',
            'size': 4}
    matplotlib.rc('font', **font)
    ax = df.plot(kind='bar')
    ax.xaxis.grid(False)
    ax.legend_.remove()
    ax.set_ylim(0, 100)

    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False)      # ticks along the bottom edge are off

    plt.tick_params(
        axis='y',          # changes apply to the y-axis
        which='both',      # both major and minor ticks are affected
        left=False)      # ticks along the left edge are off

    plt.title("Clean Bases")

    plt.savefig(out, format="png", bbox_inches='tight', dpi=300)
    csvName = out.replace(".png",".csv")
    df.to_csv(csvName, sep=",")
