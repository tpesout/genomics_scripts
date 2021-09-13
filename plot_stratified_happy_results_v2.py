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

# SAMPLES = {"Barcoded": "HG002_BC04_rows", "Non-Barcoded": "HG002_No_BC_rows"}
# SAMPLE_KEYS = ["Barcoded", "Non-Barcoded"]
# OUTPUT_NAME="./HG002_barcode_nobc_stratifications_v2.pdf"

SAMPLES = {"Parabricks": "HG002_No_BC_parabricks", "Rows": "HG002_No_BC_rows"}
SAMPLE_KEYS = ["Parabricks", "Rows"]
OUTPUT_NAME="./HG002_parabricks_vs_rows_stratifications_v2.pdf"


REGIONS = {
    "Homopolymer": "GRCh37_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz",
    "Non-Homopolymer": "GRCh37_notinAllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz",
    "Difficult Regions": "GRCh37_alldifficultregions.bed.gz",
    "Non-Difficult Regions": "GRCh37_notinalldifficultregions.bed.gz",
    "All Regions": "*",
    "Exons": "all_exons"
}
REGION_KEYS = ["All Regions", "Exons", "Homopolymer", "Non-Homopolymer", "Difficult Regions", "Non-Difficult Regions"]
VARIANT_TYPE_KEYS = ["SNP", "INDEL"]
ANALYSIS_TYPE_KEYS = ["Precision", "Recall", "F1-Score"]

# x_lims = [0.9, 0.9, 0.9, 0.5]
# x_tick_frq = [0.02, 0.02, 0.02, 0.1]
# x_lims_prec = [0.6 for _ in REGION_KEYS]
# x_tick_frq_prec = [0.15 for _ in REGION_KEYS]
# x_lims_recall = [0.2 for _ in REGION_KEYS]
# x_tick_frq_recall = [0.3 for _ in REGION_KEYS]
# x_lims_f1 = [0.3 for _ in REGION_KEYS]
# x_tick_frq_f1 = [0.25 for _ in REGION_KEYS]

x_lims = [0.0 for _ in REGION_KEYS]
x_tick_frq = [0.25 for _ in REGION_KEYS]
# maintain separation in case we want to revert
x_lims_prec = x_lims
x_tick_frq_prec = x_tick_frq
x_lims_recall = x_lims
x_tick_frq_recall = x_tick_frq
x_lims_f1 = x_lims
x_tick_frq_f1 = x_tick_frq

on_or_off_bar_threshold = 0.7

fig, axes = plt.subplots(nrows=3, ncols=len(REGION_KEYS), figsize=(2.5*len(REGION_KEYS), 7))

# Set colors
# BARCODE_COLOR = [(238 / 256.0, 124 / 256.0, 8 / 256.0, 0.7), (238 / 256.0, 124 / 256.0, 8 / 256.0, 0.8),
#                  (238 / 256.0, 124 / 256.0, 8 / 256.0, 1.0)]
# NON_BARCODE_COLOR = [(10 / 256.0, 118 / 256.0, 148 / 256.0, 0.7), (10 / 256.0, 118 / 256.0, 148 / 256.0, 0.8),
#                      (10 / 256.0, 118 / 256.0, 148 / 256.0, 1.0)]
SAMPLE0_COLOR = (238 / 256.0, 124 / 256.0, 8 / 256.0, 0.8)
SAMPLE1_COLOR = (10 / 256.0, 118 / 256.0, 148 / 256.0, 0.8)

for i in range(len(REGION_KEYS)):
    # SNP_SCORES = [
    #     DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[0]]][DATA[TYPE]==TYPE_KEYS[0]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][RECALL].values[0],
    #     DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[0]]][DATA[TYPE]==TYPE_KEYS[0]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][PRECISION].values[0],
    #     DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[0]]][DATA[TYPE]==TYPE_KEYS[0]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][F1].values[0],
    #     DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[1]]][DATA[TYPE]==TYPE_KEYS[0]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][RECALL].values[0],
    #     DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[1]]][DATA[TYPE]==TYPE_KEYS[0]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][PRECISION].values[0],
    #     DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[1]]][DATA[TYPE]==TYPE_KEYS[0]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][F1].values[0],
    # ]
    # INDEL_SCORES = [
    #     DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[0]]][DATA[TYPE]==TYPE_KEYS[1]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][RECALL].values[0],
    #     DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[0]]][DATA[TYPE]==TYPE_KEYS[1]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][PRECISION].values[0],
    #     DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[0]]][DATA[TYPE]==TYPE_KEYS[1]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][F1].values[0],
    #     DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[1]]][DATA[TYPE]==TYPE_KEYS[1]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][RECALL].values[0],
    #     DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[1]]][DATA[TYPE]==TYPE_KEYS[1]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][PRECISION].values[0],
    #     DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[1]]][DATA[TYPE]==TYPE_KEYS[1]][DATA[REGION]==REGIONS[REGION_KEYS[i]]][F1].values[0],
    # ]

    PRECISION_SCORES = [
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[0]]][DATA[TYPE] == VARIANT_TYPE_KEYS[0]][DATA[REGION] == REGIONS[REGION_KEYS[i]]][PRECISION].values[0],
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[1]]][DATA[TYPE] == VARIANT_TYPE_KEYS[0]][DATA[REGION] == REGIONS[REGION_KEYS[i]]][PRECISION].values[0],
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[0]]][DATA[TYPE] == VARIANT_TYPE_KEYS[1]][DATA[REGION] == REGIONS[REGION_KEYS[i]]][PRECISION].values[0],
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[1]]][DATA[TYPE] == VARIANT_TYPE_KEYS[1]][DATA[REGION] == REGIONS[REGION_KEYS[i]]][PRECISION].values[0],
    ]
    RECALL_SCORES = [
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[0]]][DATA[TYPE] == VARIANT_TYPE_KEYS[0]][DATA[REGION] == REGIONS[REGION_KEYS[i]]][RECALL].values[0],
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[1]]][DATA[TYPE] == VARIANT_TYPE_KEYS[0]][DATA[REGION] == REGIONS[REGION_KEYS[i]]][RECALL].values[0],
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[0]]][DATA[TYPE] == VARIANT_TYPE_KEYS[1]][DATA[REGION] == REGIONS[REGION_KEYS[i]]][RECALL].values[0],
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[1]]][DATA[TYPE] == VARIANT_TYPE_KEYS[1]][DATA[REGION] == REGIONS[REGION_KEYS[i]]][RECALL].values[0],

    ]
    F1_SCORES = [
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[0]]][DATA[TYPE] == VARIANT_TYPE_KEYS[0]][DATA[REGION] == REGIONS[REGION_KEYS[i]]][F1].values[0],
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[1]]][DATA[TYPE] == VARIANT_TYPE_KEYS[0]][DATA[REGION] == REGIONS[REGION_KEYS[i]]][F1].values[0],
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[0]]][DATA[TYPE] == VARIANT_TYPE_KEYS[1]][DATA[REGION] == REGIONS[REGION_KEYS[i]]][F1].values[0],
        DATA[DATA[SAMPLE]==SAMPLES[SAMPLE_KEYS[1]]][DATA[TYPE] == VARIANT_TYPE_KEYS[1]][DATA[REGION] == REGIONS[REGION_KEYS[i]]][F1].values[0],

    ]


    # SNP_COLORS = BARCODE_COLOR + NON_BARCODE_COLOR
    # INDEL_COLORS = BARCODE_COLOR + NON_BARCODE_COLOR

    # Bar charts
    N = 4
    ind = np.arange(N)
    width = 0.98

    bar_positions = [0, 1, 3.25, 4.25]
    # bar_positions = [0, 1, 2, 3.25, 4.25, 5.25]

    # bar_positions.reverse()
    # no need to reverse the color

    x_lim_max = 1.00
    first_start = 0
    p1 = axes[0][i].bar(bar_positions, PRECISION_SCORES, width, color=[SAMPLE0_COLOR, SAMPLE1_COLOR, SAMPLE0_COLOR, SAMPLE1_COLOR])
    p2 = axes[1][i].bar(bar_positions, RECALL_SCORES, width, color=[SAMPLE0_COLOR, SAMPLE1_COLOR, SAMPLE0_COLOR, SAMPLE1_COLOR])
    p3 = axes[2][i].bar(bar_positions, F1_SCORES, width, color=[SAMPLE0_COLOR, SAMPLE1_COLOR, SAMPLE0_COLOR, SAMPLE1_COLOR])


    # Add text labels
    # plt.text(6, -950, "Assembler", fontsize = text_fontsize+3)
    # axes[0][i].set_xlabel('Value', fontsize = text_fontsize + 1, color='black')
    # axes[1][i].set_xlabel('Value', fontsize = text_fontsize + 1, color='black')
    # Category = ["R", "P", "F1"]

    x_lim_min_prec = x_lims_prec[i]
    # value_offset_snp = (x_lim_max - x_lim_min_prec) / 10
    for ii, v in enumerate(PRECISION_SCORES):

        if v < on_or_off_bar_threshold: #0.962: # and x_lim_min_snp == 0.9 or v < 0.6:
            loc = v + 0.025
            va_align = 'bottom'
            ha_align = 'center'
            color = 'black'
        else:
            text_start = x_lim_min_prec + (x_lim_max - x_lim_min_prec) * 0.20
            value_start = x_lim_min_prec + (x_lim_max - x_lim_min_prec) * 0.39
            loc = v - 0.025
            va_align = 'top'
            ha_align = 'center'
            color = 'white'

        axes[0][i].text(bar_positions[ii], loc, str(v)[0:6], va=va_align, ha=ha_align, rotation=90, color=color, fontsize=text_fontsize)

    x_lim_min_recall = x_lims_recall[i]
    # value_offset_indel = (x_lim_max - x_lim_min_recall) / 10
    for ii, v in enumerate(RECALL_SCORES):

        if v < on_or_off_bar_threshold: #0.91 and x_lim_min_recall == 0.9 or v < 0.6:
            loc = v + 0.025
            va_align = 'bottom'
            ha_align = 'center'
            color = 'black'
        else:
            text_start = x_lim_min_recall + (x_lim_max - x_lim_min_recall) * 0.20
            value_start = x_lim_min_recall + (x_lim_max - x_lim_min_recall) * 0.39
            loc = v - 0.025
            va_align = 'top'
            ha_align = 'center'
            color = 'white'
        axes[1][i].text(bar_positions[ii], loc, str(v)[0:6], va=va_align, ha=ha_align, rotation=90, color=color, fontsize=text_fontsize)

    x_lim_min_f1 = x_lims_f1[i]
    # value_offset_indel = (x_lim_max - x_lim_min_f1) / 10
    for ii, v in enumerate(F1_SCORES):

        if v < on_or_off_bar_threshold: #0.91 and x_lim_min_f1 == 0.9 or v < 0.6:
            loc = v + 0.025
            va_align = 'bottom'
            ha_align = 'center'
            color = 'black'
        else:
            text_start = x_lim_min_f1 + (x_lim_max - x_lim_min_f1) * 0.20
            value_start = x_lim_min_f1 + (x_lim_max - x_lim_min_f1) * 0.39
            loc = v - 0.025
            va_align = 'top'
            ha_align = 'center'
            color = 'white'
        axes[2][i].text(bar_positions[ii], loc, str(v)[0:6], va=va_align, ha=ha_align, rotation=90, color=color, fontsize=text_fontsize)

    # if i == len(REGIONS) - 1:
    #     # caller_name_offser = x_lim_max + (0.08 * (1-x_lim_max)) # -0.13 for 0.0
    #     type_name_offset = x_lim_max + (0.10 * (1-x_lim_max)) # -0.07 for 0.0
    #     axes[0][i].text(type_name_offset, bar_positions[0] + ((bar_positions[-1] - bar_positions[0]) / 2), ANALYSIS_TYPE_KEYS[0], va='center', fontsize=text_fontsize, rotation=270)
    #     axes[1][i].text(type_name_offset, bar_positions[0] + ((bar_positions[-1] - bar_positions[0]) / 2), ANALYSIS_TYPE_KEYS[1], va='center', fontsize=text_fontsize, rotation=270)
    #     axes[2][i].text(type_name_offset, bar_positions[0] + ((bar_positions[-1] - bar_positions[0]) / 2), ANALYSIS_TYPE_KEYS[2], va='center', fontsize=text_fontsize, rotation=270)

    # if i == 2:
    categories = VARIANT_TYPE_KEYS
    mid_points = [np.mean([bar_positions[0], bar_positions[1]]),
                  np.mean([bar_positions[2], bar_positions[3]])]
    all_categories = categories# + categories
    # axes[0][i].set_xticks(mid_points)
    # axes[0][i].set_xticklabels(all_categories, fontsize=text_fontsize+1, color='black', va='top')
    # axes[1][i].set_xticks(mid_points)
    # axes[1][i].set_xticklabels(all_categories, fontsize=text_fontsize+1, color='black', va='top')
    axes[2][i].set_xticks(mid_points)
    axes[2][i].set_xticklabels(all_categories, fontsize=text_fontsize+1, color='black', va='top')
    # else:
    axes[0][i].set_xticks([])
    axes[0][i].set_xticklabels([])
    axes[1][i].set_xticks([])
    axes[1][i].set_xticklabels([])
    # axes[2][i].set_xticks([])
    # axes[2][i].set_xticklabels([])

    # Samples = ["ONT", "HG004"]
    # plt.text(bar_positions1[1], sample_name_offset, Samples[0], ha='center', color='black', fontsize=text_fontsize)
    # plt.text(bar_positions2[1], sample_name_offset, Samples[0], ha='center', color='black', fontsize=text_fontsize)

    # plt.text(bar_positions1[4], sample_name_offset, Samples[1], ha='center', color='black', fontsize=text_fontsize)
    # plt.text(bar_positions2[4], sample_name_offset, Samples[1], ha='center', color='black', fontsize=text_fontsize)

    # categories = [" BC", " No BC"]
    # categories = ["Recall", "Precision", "F1-Score"]
    # if i == 0:
    # caller_name_offset = x_lim_min - (0.08 * (1-x_lim_max)) # -0.13 for 0.0

    # for ii, v in enumerate(PRECISION_SCORES):
    #     axes[0][i].text(bar_positions[ii], x_lim_min_prec + 0.001, categories[ii % 2], color='white', va='bottom', rotation=90, ha='center', fontsize=text_fontsize)
    # for ii, v in enumerate(RECALL_SCORES):
    #     axes[1][i].text(bar_positions[ii], x_lim_min_recall + 0.001, categories[ii % 2], color='white', va='bottom', rotation=90, ha='center', fontsize=text_fontsize)
    # for ii, v in enumerate(F1_SCORES):
    #     axes[2][i].text(bar_positions[ii], x_lim_min_f1 + 0.001, categories[ii % 2], color='white', va='bottom', rotation=90, ha='center', fontsize=text_fontsize)

    axes[0][i].tick_params(axis=u'both', which=u'both',length=0)
    axes[0][i].set_ylim(x_lim_min_prec, x_lim_max)
    # axes[0][i].set_title(REGION_KEYS[i], fontsize=text_fontsize + 1)

    axes[1][i].tick_params(axis=u'both', which=u'both',length=0)
    axes[1][i].set_ylim(x_lim_min_recall, x_lim_max)

    axes[2][i].tick_params(axis=u'both', which=u'both',length=0)
    axes[2][i].set_ylim(x_lim_min_f1, x_lim_max)
    axes[2][i].set_xlabel(REGION_KEYS[i])

    x_ticks_prec = np.arange(x_lims_prec[i] + x_tick_frq_prec[i], 1.01, x_tick_frq_prec[i] * 2)
    x_lables_prec = ["%.2lf" % round(x,2) for x in x_ticks_prec]
    axes[0][i].set_yticks(x_ticks_prec)

    x_ticks_recall = np.arange(x_lims_recall[i] + x_tick_frq_recall[i], 1.01, x_tick_frq_recall[i] * 2)
    x_lables_recall = ["%.2lf" % round(x,2) for x in x_ticks_recall]
    axes[1][i].set_yticks(x_ticks_recall)

    x_ticks_f1 = np.arange(x_lims_f1[i] + x_tick_frq_f1[i], 1.01, x_tick_frq_f1[i] * 2)
    x_lables_f1 = ["%.2lf" % round(x,2) for x in x_ticks_f1]
    axes[2][i].set_yticks(x_ticks_f1)

    if i == 0:
        axes[0][i].set_yticklabels(x_lables_prec, fontsize=text_fontsize, color='black')

        axes[1][i].set_yticklabels(x_lables_recall, fontsize=text_fontsize, color='black')

        axes[2][i].set_yticklabels(x_lables_f1, fontsize=text_fontsize, color='black')
    else:
        axes[0][i].set_yticklabels([])
        axes[1][i].set_yticklabels([])
        axes[2][i].set_yticklabels([])


    if i == 0:
        axes[0][i].set_ylabel(ANALYSIS_TYPE_KEYS[0])
        axes[1][i].set_ylabel(ANALYSIS_TYPE_KEYS[1])
        axes[2][i].set_ylabel(ANALYSIS_TYPE_KEYS[2])



# add legend
def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique), loc=8)
axes[2][len(REGION_KEYS) - 1].plot(0, 0, 0, color=SAMPLE0_COLOR, label=SAMPLE_KEYS[0])
axes[2][len(REGION_KEYS) - 1].plot(0, 0, 0, color=SAMPLE1_COLOR, label=SAMPLE_KEYS[1])
legend_without_duplicate_labels(axes[2][len(REGION_KEYS) - 1])

plt.show()
plt.savefig(OUTPUT_NAME, format='pdf', dpi=300)
