import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# plt.style.use('ggplot')
text_fontsize = 8
# # plt.rcParams['ytick.labelsize']=text_fontsize+4
# plt.rcParams['pdf.fonttype'] = 42
# plt.switch_backend('agg')

SNP_INDEX=2
INDEL_INDEX=0

REGION="Subset"
SAMPLE="Sample"
TYPE="Type"
RECALL="METRIC.Recall"
PRECISION="METRIC.Precision"
F1="METRIC.F1_Score"

# read inputs / CSV from hap.py
DATA = pd.read_csv("/home/tpesout/work/nicu/stratification_analysis/HG002_barcode_nobc_stratifications.csv")

SAMPLES = {"HG002 Barcoded": "HG002_BC04", "HG002 Non-Barcoded": "HG002_No_BC"}
SAMPLE_KEYS = ["HG002 Barcoded", "HG002 Non-Barcoded"]
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

x_lims = [0.9, 0.9, 0.9, 0.5]
x_tick_frq = [0.02, 0.02, 0.02, 0.1]
x_lims = [0.6 for _ in REGION_KEYS]
x_tick_frq = [0.1 for _ in REGION_KEYS]

fig, axes = plt.subplots(nrows=2, ncols=len(REGION_KEYS), figsize=(2*len(REGION_KEYS), 7))

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

    bar_positions2 = [0, 1, 2, 3.25, 4.25, 5.25]
    bar_positions1 = [7.25, 8.25, 9.25, 10.5, 11.5, 12.5]
    # bar_positions2 = [0, 1, 2, 3.25, 4.25, 5.25, 6.50, 7.50, 8.50]
    # bar_positions1 = [10.50, 11.50, 12.50, 13.75, 14.75, 15.75, 17, 18, 19]

    bar_positions1.reverse()
    bar_positions2.reverse()
    # HG003_SCORES.reverse()
    # HG004_SCORES.reverse()
    # HG003_COLOR.reverse()
    # HG004_COLOR.reverse()
    # no need to reverse the color

    x_lim_min = x_lims[i]
    x_lim_max = 1.00
    first_start = 0
    p1 = axes[i].barh(bar_positions1, SNP_SCORES, width, color = SNP_COLORS)
    axes[i].hlines((bar_positions1[-1] + bar_positions2[0]) / 2, x_lim_min, x_lim_max, colors='black', linestyles='solid')
    p2 = axes[i].barh(bar_positions2, INDEL_SCORES, width, color = INDEL_COLORS)



    # Add text labels
    # plt.text(6, -950, "Assembler", fontsize = text_fontsize+3)
    axes[i].set_xlabel('Value', fontsize = text_fontsize + 1, color='black')
    Category = ["R", "P", "F1"]


    value_offset = (x_lim_max - x_lim_min) / 10
    for ii, v in enumerate(SNP_SCORES):
        # axes[i].text(text_start, bar_positions1[ii], Category[ii % 3], va='center', color='black', fontsize=text_fontsize)
        if v < 0.91 and x_lim_min == 0.9 or v < 0.6:
            loc = v
            ha_align = 'left'
            va_align = 'center'
            color = 'black'
        else:
            text_start = x_lim_min + (x_lim_max - x_lim_min) * 0.20
            value_start = x_lim_min + (x_lim_max - x_lim_min) * 0.39
            loc = v - value_offset
            ha_align = 'right'
            va_align = 'center'
            color = 'white'

        axes[i].text(v, bar_positions1[ii], str(v)[0:5], va=va_align, ha=ha_align, color=color, fontsize=text_fontsize)
    for ii, v in enumerate(INDEL_SCORES):
        # axes[i].text(text_start, bar_positions2[ii], Category[ii % 3], va='center', color='black', fontsize=text_fontsize)
        if v < 0.91 and x_lim_min == 0.9 or v < 0.6:
            loc = v
            ha_align = 'left'
            va_align = 'center'
            color = 'black'
        else:
            text_start = x_lim_min + (x_lim_max - x_lim_min) * 0.20
            value_start = x_lim_min + (x_lim_max - x_lim_min) * 0.39
            loc = v - value_offset
            ha_align = 'right'
            va_align = 'center'
            color = 'white'
        axes[i].text(v, bar_positions2[ii], str(v)[0:5], va=va_align, ha=ha_align, color=color, fontsize=text_fontsize)

    if i == len(REGIONS) - 1:
        # caller_name_offser = x_lim_max + (0.08 * (1-x_lim_max)) # -0.13 for 0.0
        type_name_offset = x_lim_max + (0.10 * (1-x_lim_max)) # -0.07 for 0.0
        axes[i].text(type_name_offset, bar_positions1[0] + ((bar_positions1[-1] - bar_positions1[0]) / 2), TYPE_KEYS[0], va='center', fontsize=text_fontsize, rotation=270)
        axes[i].text(type_name_offset, bar_positions2[0] + ((bar_positions2[-1] - bar_positions2[0]) / 2), TYPE_KEYS[1], va='center', fontsize=text_fontsize, rotation=270)

    if i == 0:
        categories = SAMPLE_KEYS
        mid_points = [bar_positions1[1], bar_positions1[4],
                      bar_positions2[1], bar_positions2[4]]
        all_categories = categories + categories
        axes[i].set_yticks(mid_points)
        axes[i].set_yticklabels(all_categories, fontsize=text_fontsize+1, color='black')
    else:
        axes[i].set_yticks([])
        axes[i].set_yticklabels([])

    # Samples = ["ONT", "HG004"]
    # plt.text(bar_positions1[1], sample_name_offset, Samples[0], ha='center', color='black', fontsize=text_fontsize)
    # plt.text(bar_positions2[1], sample_name_offset, Samples[0], ha='center', color='black', fontsize=text_fontsize)

    # plt.text(bar_positions1[4], sample_name_offset, Samples[1], ha='center', color='black', fontsize=text_fontsize)
    # plt.text(bar_positions2[4], sample_name_offset, Samples[1], ha='center', color='black', fontsize=text_fontsize)

    categories = ["Recall", "Precision", "F1-Score"]
    if i == 0:
        caller_name_offset = x_lim_min - (0.08 * (1-x_lim_max)) # -0.13 for 0.0

        for ii, v in enumerate(SNP_SCORES):
            axes[i].text(x_lim_min + 0.001, bar_positions2[ii], categories[ii % 3], color='white', va='center', ha='left', fontsize=text_fontsize)
            axes[i].text(x_lim_min + 0.001, bar_positions1[ii], categories[ii % 3], color='white', va='center', ha='left', fontsize=text_fontsize)

    axes[i].tick_params(axis=u'both', which=u'both',length=0)
    axes[i].set_xlim(x_lim_min, x_lim_max)
    axes[i].set_title(REGION_KEYS[i], fontsize=text_fontsize + 1)

    x_ticks = np.arange(x_lims[i] + x_tick_frq[i], 1.01, x_tick_frq[i] * 2)
    x_lables = ["%.2lf" % round(x,2) for x in x_ticks]
    axes[i].set_xticks(x_ticks)
    axes[i].set_xticklabels(x_lables, fontsize=text_fontsize, color='black')
    # Save the figure
plt.tight_layout()
# plt.margins(0, 0)
plt.show()
# plt.savefig('./outputs/2b_Illumina_difficult.pdf', format='pdf', dpi=300)
