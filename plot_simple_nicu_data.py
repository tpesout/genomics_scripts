#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

plt.style.use('ggplot')
text_fontsize = 8
# plt.rcParams['ytick.labelsize']=text_fontsize+4
plt.rcParams.update({'font.size': text_fontsize})
plt.rcParams['pdf.fonttype'] = 42
plt.switch_backend('agg')

# gsutil ls gs://ultra_rapid_analysis/small_variant_call/ | grep -v annotated | xargz gsutil cp {} .
# ls *vcf.gz | sed 's/.vcf.gz//' | xargz bash -c "rtg vcfstats {}.vcf.gz >{}.stats.txt &"
# ls *txt | sed 's/_pmd.stats.txt//' | xargz bash -c "echo \"'{}': [\$(cat {}_pmd.stats.txt | grep 'SNP Transitions/Transversions:' | sed 's#SNP Transitions/Transversions:.*(##' | sed 's#/#, #' | sed 's/)/],/' )\""
# ls *txt | sed 's/_pmd.stats.txt//' | xargz bash -c "echo \"'{}': [\$(cat {}_pmd.stats.txt | grep 'Total Het' | sed 's#Total Het/Hom ratio          :.*(##' | sed 's#/#, #' | sed 's/)/],/')\""

HET_HOM_DATA = {
    'Fast_001': [2658285, 1784920],
    'Fast_002': [2658001, 1784279],
    'Fast_003': [2764275, 1714075],
    'Fast_004': [2512953, 1954227],
    'Fast_005': [2783424, 1808957],
    'Fast_006': [2770999, 1729294],
    'Fast_007': [2574593, 1907721],
    'Fast_108': [2586607, 1917060],
    'Fast_109': [2840769, 1778498],
    'Fast_110': [2590270, 1773955],
    'Fast_111': [2428049, 1887499],
    'Fast_112': [3155839, 1614610],
    # 'Fast_113': [2745701, 1735191],
    # 'HG002_BC04': [2676530, 1784662],
    # 'HG002_No_BC': [2657802, 1780983]
}

TI_TV_DATA = {
    'Fast_001': [3703177, 1844994],
    'Fast_002': [3702538, 1844184],
    'Fast_003': [3692354, 1836180],
    'Fast_004': [3819353, 1904904],
    'Fast_005': [3802421, 1896421],
    'Fast_006': [3709590, 1846910],
    'Fast_007': [3814216, 1899969],
    'Fast_108': [3810580, 1900257],
    'Fast_109': [3811862, 1897808],
    'Fast_110': [3643634, 1805392],
    'Fast_111': [3683079, 1825657],
    'Fast_112': [3810437, 1880573],
    # 'Fast_113': [3702438, 1832496],
    # 'HG002_BC04': [3719435, 1855055],
    # 'HG002_No_BC': [3715178, 1851787]
}

QV_DATA = {
    'Fast_001': 12.5,
    'Fast_002': 12,
    'Fast_003': 11.5,
    'Fast_004': 11.2,
    'Fast_005': 11.5,
    'Fast_006': 11.1,
    'Fast_007': 10.7,
    'Fast_108': 11.9,
    'Fast_109': 11.8,
    'Fast_110': 11.5,
    'Fast_111': 11.4,
    'Fast_112': 11.2,
}



def parse_args(args = None):
    parser = argparse.ArgumentParser("Plots information from haplotagging_stats tsv")
    parser.add_argument('--figure_name', '-f', dest='figure_name', default="output", required=False, type=str,
                       help='Figure name (default to "output")')

    return parser.parse_args() if args is None else parser.parse_args(args)


def plot_simple_ratio(data, title, xlabel, ylabel, data_key, ylim=None, output_name=None):

    fig, ((ax1)) = plt.subplots(nrows=1,ncols=1)

    xs = sorted(list(data.keys()))
    ys = [data_key(data[x]) for x in xs]

    ax1.scatter(range(len(xs)), ys, marker = "D", c=cm.rainbow(np.linspace(1, 0, len(xs))))
    ax1.axhline(np.mean(ys), linestyle='dashed', color = 'black')
    print("Mean: {}".format(np.mean(ys)))

    plt.xticks(range(len(xs)), [x[6:] for x in xs])#, rotation = 45)

    plt.title(title)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    if ylim is None:
        ax1.set_ylim([0, 1.05 * max(ys)])
    else:
        ax1.set_ylim(ylim)

    plt.tight_layout()
    if output_name is not None:
        plt.savefig(output_name + ".pdf", dpi=300)
    plt.show()
    plt.close()


def plot_two_simple_ratio(data_l, data_r, title, xlabel, yllabel, yrlabel, yllim, yrlim, output_name=None):

    fig, ((ax1)) = plt.subplots(nrows=1,ncols=1)
    ax2 = ax1.twinx()

    xs = sorted(list(set(data_l.keys()).intersection(data_r.keys())))
    yls = [data_l[x][0] / data_l[x][1] for x in xs]
    yrs = [data_r[x][0] / data_r[x][1] for x in xs]

    l_line = ax1.scatter(range(len(xs)), yls, label=yllabel, color='blue')
    r_line = ax2.scatter(range(len(xs)), yrs, label=yrlabel, color='orange')

    plt.xticks(range(len(xs)), xs, rotation = 45)
    for tick in ax1.get_xticklabels():
        tick.set_rotation(45)

    plt.title(title)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(yllabel)
    ax2.set_ylabel(yrlabel)
    ax1.set_ylim(yllim)
    ax2.set_ylim(yrlim)

    lns = [l_line, r_line]
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=2)
    plt.tight_layout()
    plt.show()
    plt.close()

def main():
    args = parse_args()

    plot_simple_ratio(TI_TV_DATA, "Transition/Transversion Ratio", "Sample", "Ti/Tv Ratio", lambda x: x[0] / x[1], [1.9,2.1], "{}.titv".format(args.figure_name))
    plot_simple_ratio(HET_HOM_DATA, "Het/Hom Ratio", "Sample", "Het/Hom Ratio", lambda x: x[0] / x[1], [0.95, 2.05], "{}.hethom".format(args.figure_name))
    plot_simple_ratio(QV_DATA, "Avg Read QV", "Sample", "Avg Read QV", lambda x: x, [8, 13], "{}.readqv".format(args.figure_name))

if __name__ == "__main__":
    main()
