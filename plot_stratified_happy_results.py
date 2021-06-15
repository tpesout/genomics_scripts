import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

text_fontsize = 8
plt.style.use('ggplot')
# plt.rcParams['ytick.labelsize']=text_fontsize+4
plt.rcParams['pdf.fonttype'] = 42
plt.switch_backend('agg')

SNP_INDEX=2
INDEL_INDEX=0

REGION="Subset"
SAMPLE="Sample"
TYPE="Type"
RECALL="METRIC.Recall"
PRECISION="METRIC.Precision"
F1="METRIC.F1_Score"

# read inputs / CSV from hap.py
DATA = pd.read_csv("bcrows_nobcrows_nobcpara/HG002_Stratifications.csv")

SAMPLES = {"Barcoded": "HG002_BC04_rows", "Non-Barcoded": "HG002_No_BC_rows"}
SAMPLE_KEYS = ["Barcoded", "Non-Barcoded"]
OUTPUT_NAME="./HG002_barcode_nobc_stratifications.pdf"

# SAMPLES = {"Parabricks": "HG002_No_BC_parabricks", "Rows": "HG002_No_BC_rows"}
# SAMPLE_KEYS = ["Parabricks", "Rows"]
# OUTPUT_NAME="./HG002_parabricks_vs_rows_stratifications.pdf"


REGIONS = {
    "Homopolymer": "GRCh37_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz",
    "Non-Homopolymer": "GRCh37_notinAllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz",
    "Difficult Regions": "GRCh37_alldifficultregions.bed.gz",
    "Non-Difficult Regions": "GRCh37_notinalldifficultregions.bed.gz",
    "All Regions": "*",
    "Exons": "all_exons"
}
REGION_KEYS = ["All Regions", "Homopolymer", "Non-Homopolymer", "Difficult Regions", "Non-Difficult Regions", "Exons"]
TYPE_KEYS = ["SNP","INDEL"]

# x_lims = [0.9, 0.9, 0.9, 0.5]
# x_tick_frq = [0.02, 0.02, 0.02, 0.1]
x_lims_snp = [0.925 for _ in REGION_KEYS]
x_tick_frq_snp = [0.015 for _ in REGION_KEYS]
x_lims_indel = [0.0 for _ in REGION_KEYS]
x_tick_frq_indel = [0.2 for _ in REGION_KEYS]

fig, axes = plt.subplots(nrows=2, ncols=len(REGION_KEYS), figsize=(2.5*len(REGION_KEYS), 7))

for i in range(len(REGION_KEYS)):
    SNP_SCORES = [
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[0]]][DATA[TYPE]==TYPE_KEYS[0]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][RECALL].values[0],
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[0]]][DATA[TYPE]==TYPE_KEYS[0]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][PRECISION].values[0],
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[0]]][DATA[TYPE]==TYPE_KEYS[0]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][F1].values[0],
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[1]]][DATA[TYPE]==TYPE_KEYS[0]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][RECALL].values[0],
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[1]]][DATA[TYPE]==TYPE_KEYS[0]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][PRECISION].values[0],
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[1]]][DATA[TYPE]==TYPE_KEYS[0]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][F1].values[0],
    ]
    INDEL_SCORES = [
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[0]]][DATA[TYPE]==TYPE_KEYS[1]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][RECALL].values[0],
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[0]]][DATA[TYPE]==TYPE_KEYS[1]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][PRECISION].values[0],
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[0]]][DATA[TYPE]==TYPE_KEYS[1]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][F1].values[0],
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[1]]][DATA[TYPE]==TYPE_KEYS[1]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][RECALL].values[0],
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[1]]][DATA[TYPE]==TYPE_KEYS[1]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][PRECISION].values[0],
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[1]]][DATA[TYPE]==TYPE_KEYS[1]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][F1].values[0],
    ]

    # Set colors
    BARCODE_COLOR = [(238/256.0, 124/256.0, 8/256.0, 0.7), (238/256.0, 124/256.0, 8/256.0, 0.8), (238/256.0, 124/256.0, 8/256.0, 1.0)]
    NON_BARCODE_COLOR = [(10/256.0, 118/256.0, 148/256.0, 0.7), (10/256.0, 118/256.0, 148/256.0, 0.8), (10/256.0, 118/256.0, 148/256.0, 1.0)]

    SNP_COLORS = BARCODE_COLOR + NON_BARCODE_COLOR
    INDEL_COLORS = BARCODE_COLOR + NON_BARCODE_COLOR

    # Bar charts
    N = 4
    ind = np.arange(N)
    width = 0.98

    bar_positions = [0, 1, 2, 3.25, 4.25, 5.25]

    bar_positions.reverse()
    # no need to reverse the color

    x_lim_max = 1.00
    first_start = 0
    p1 = axes[0][i].barh(bar_positions, SNP_SCORES, width, color = SNP_COLORS)
    p2 = axes[1][i].barh(bar_positions, INDEL_SCORES, width, color = INDEL_COLORS)


    # Add text labels
    # plt.text(6, -950, "Assembler", fontsize = text_fontsize+3)
    # axes[0][i].set_xlabel('Value', fontsize = text_fontsize + 1, color='black')
    axes[1][i].set_xlabel('Value', fontsize = text_fontsize + 1, color='black')
    Category = ["R", "P", "F1"]

    x_lim_min_snp = x_lims_snp[i]
    value_offset_snp = (x_lim_max - x_lim_min_snp) / 10
    for ii, v in enumerate(SNP_SCORES):

        if v < 0.962: # and x_lim_min_snp == 0.9 or v < 0.6:
            loc = v
            ha_align = 'left'
            va_align = 'center'
            color = 'black'
        else:
            text_start = x_lim_min_snp + (x_lim_max - x_lim_min_snp) * 0.20
            value_start = x_lim_min_snp + (x_lim_max - x_lim_min_snp) * 0.39
            loc = v - value_offset_snp
            ha_align = 'right'
            va_align = 'center'
            color = 'white'

        axes[0][i].text(v, bar_positions[ii], str(v)[0:6], va=va_align, ha=ha_align, color=color, fontsize=text_fontsize)

    x_lim_min_indel = x_lims_indel[i]
    value_offset_indel = (x_lim_max - x_lim_min_indel) / 10
    for ii, v in enumerate(INDEL_SCORES):

        if v < 0.91 and x_lim_min_indel == 0.9 or v < 0.6:
            loc = v
            ha_align = 'left'
            va_align = 'center'
            color = 'black'
        else:
            text_start = x_lim_min_indel + (x_lim_max - x_lim_min_indel) * 0.20
            value_start = x_lim_min_indel + (x_lim_max - x_lim_min_indel) * 0.39
            loc = v - value_offset_indel
            ha_align = 'right'
            va_align = 'center'
            color = 'white'
        axes[1][i].text(v, bar_positions[ii], str(v)[0:6], va=va_align, ha=ha_align, color=color, fontsize=text_fontsize)

    if i == len(REGIONS) - 1:
        # caller_name_offser = x_lim_max + (0.08 * (1-x_lim_max)) # -0.13 for 0.0
        type_name_offset = x_lim_max + (0.10 * (1-x_lim_max)) # -0.07 for 0.0
        axes[0][i].text(type_name_offset, bar_positions[0] + ((bar_positions[-1] - bar_positions[0]) / 2), TYPE_KEYS[0], va='center', fontsize=text_fontsize, rotation=270)
        axes[1][i].text(type_name_offset, bar_positions[0] + ((bar_positions[-1] - bar_positions[0]) / 2), TYPE_KEYS[1], va='center', fontsize=text_fontsize, rotation=270)

    if i == 0:
        categories = SAMPLE_KEYS
        mid_points = [bar_positions[1], bar_positions[4],
                      bar_positions[1], bar_positions[4]]
        all_categories = categories + categories
        axes[0][i].set_yticks(mid_points)
        axes[0][i].set_yticklabels(all_categories, fontsize=text_fontsize+1, color='black', rotation=90, va='center')
        axes[1][i].set_yticks(mid_points)
        axes[1][i].set_yticklabels(all_categories, fontsize=text_fontsize+1, color='black', rotation=90, va='center')
    else:
        axes[0][i].set_yticks([])
        axes[0][i].set_yticklabels([])
        axes[1][i].set_yticks([])
        axes[1][i].set_yticklabels([])

    # Samples = ["ONT", "HG004"]
    # plt.text(bar_positions1[1], sample_name_offset, Samples[0], ha='center', color='black', fontsize=text_fontsize)
    # plt.text(bar_positions2[1], sample_name_offset, Samples[0], ha='center', color='black', fontsize=text_fontsize)

    # plt.text(bar_positions1[4], sample_name_offset, Samples[1], ha='center', color='black', fontsize=text_fontsize)
    # plt.text(bar_positions2[4], sample_name_offset, Samples[1], ha='center', color='black', fontsize=text_fontsize)

    categories = ["Recall", "Precision", "F1-Score"]
    if i == 0:
        # caller_name_offset = x_lim_min - (0.08 * (1-x_lim_max)) # -0.13 for 0.0

        for ii, v in enumerate(SNP_SCORES):
            axes[0][i].text(x_lim_min_snp + 0.001, bar_positions[ii], categories[ii % 3], color='white', va='center', ha='left', fontsize=text_fontsize)
        for ii, v in enumerate(INDEL_SCORES):
            axes[1][i].text(x_lim_min_indel + 0.001, bar_positions[ii], categories[ii % 3], color='white', va='center', ha='left', fontsize=text_fontsize)

    axes[0][i].tick_params(axis=u'both', which=u'both',length=0)
    axes[0][i].set_xlim(x_lim_min_snp, x_lim_max)
    axes[0][i].set_title(REGION_KEYS[i], fontsize=text_fontsize + 1)

    axes[1][i].tick_params(axis=u'both', which=u'both',length=0)
    axes[1][i].set_xlim(x_lim_min_indel, x_lim_max)
    # axes[1][i].set_title(REGION_KEYS[i], fontsize=text_fontsize + 1)

    x_ticks_snp = np.arange(x_lims_snp[i] + x_tick_frq_snp[i], 1.01, x_tick_frq_snp[i] * 2)
    x_lables_snp = ["%.2lf" % round(x,2) for x in x_ticks_snp]
    axes[0][i].set_xticks(x_ticks_snp)
    axes[0][i].set_xticklabels(x_lables_snp, fontsize=text_fontsize, color='black')

    x_ticks_indel = np.arange(x_lims_indel[i] + x_tick_frq_indel[i], 1.01, x_tick_frq_indel[i] * 2)
    x_lables_indel = ["%.2lf" % round(x,2) for x in x_ticks_indel]
    axes[1][i].set_xticks(x_ticks_indel)
    axes[1][i].set_xticklabels(x_lables_indel, fontsize=text_fontsize, color='black')


    # Save the figure
# plt.tight_layout()
# plt.margins(0, 0)
plt.show()
plt.savefig(OUTPUT_NAME, format='pdf', dpi=300)
