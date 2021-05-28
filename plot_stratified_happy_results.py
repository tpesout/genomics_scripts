import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# plt.style.use('ggplot')
# text_fontsize = 8
# # plt.rcParams['ytick.labelsize']=text_fontsize+4
# plt.rcParams['pdf.fonttype'] = 42
# plt.switch_backend('agg')

SNP_INDEX=2
INDEL_INDEX=0

# read inputs / CSV from hap.py

HG003_CCS_DV = pd.read_csv("./figure_inputs/Figure_2bc_difficult_regions/HG003_DV_v110_DV_margin_DV.extended.csv")
HG004_CCS_DV = pd.read_csv("./figure_inputs/Figure_2bc_difficult_regions/HG004_DV_v110_DV_margin_DV.extended.csv")
HG003_ILM_DV = pd.read_csv("./figure_inputs/Figure_2bc_difficult_regions/HG003_illumina_deepvariant_v1.extended.csv")
HG004_ILM_DV = pd.read_csv("./figure_inputs/Figure_2bc_difficult_regions/HG004_illumina_deepvariant_v1.extended.csv")
HG003_ONT_DV = pd.read_csv("./figure_inputs/Figure_2bc_difficult_regions/HG003_guppy422_90x_pmd_giab_421.extended.csv")
HG004_ONT_DV = pd.read_csv("./figure_inputs/Figure_2bc_difficult_regions/HG004_guppy422_90x_pmd_giab_421.extended.csv")

Regions = ["Major Histocompatibility Complex (MHC)", "Segmental Duplications", "Low mappability regions", "250bp+ non-unique regions"]
Region_index = [2752 - 2, 2854 - 2, 2826 - 2, 2830 - 2]
x_lims = [0.9, 0.9, 0.9, 0.5]
x_tick_frq = [0.02, 0.02, 0.02, 0.1]

fig, axes = plt.subplots(nrows=1, ncols=len(Regions), figsize=(8, 3.5))

for i in range(len(Regions)):
    print("PROCESSING: ", HG003_ONT_DV['Subset'][Region_index[i]])
    # HG003 RPF HG004 RPF
    HG003_SCORES = [HG003_ONT_DV['METRIC.Recall'][Region_index[i]],
                    HG003_ONT_DV['METRIC.Precision'][Region_index[i]],
                    HG003_ONT_DV['METRIC.F1_Score'][Region_index[i]],
                    HG003_ILM_DV['METRIC.Recall'][Region_index[i]],
                    HG003_ILM_DV['METRIC.Precision'][Region_index[i]],
                    HG003_ILM_DV['METRIC.F1_Score'][Region_index[i]],
                    HG003_CCS_DV['METRIC.Recall'][Region_index[i]],
                    HG003_CCS_DV['METRIC.Precision'][Region_index[i]],
                    HG003_CCS_DV['METRIC.F1_Score'][Region_index[i]]]

    HG004_SCORES = [HG004_ONT_DV['METRIC.Recall'][Region_index[i]],
                    HG004_ONT_DV['METRIC.Precision'][Region_index[i]],
                    HG004_ONT_DV['METRIC.F1_Score'][Region_index[i]],
                    HG004_ILM_DV['METRIC.Recall'][Region_index[i]],
                    HG004_ILM_DV['METRIC.Precision'][Region_index[i]],
                    HG004_ILM_DV['METRIC.F1_Score'][Region_index[i]],
                    HG004_CCS_DV['METRIC.Recall'][Region_index[i]],
                    HG004_CCS_DV['METRIC.Precision'][Region_index[i]],
                    HG004_CCS_DV['METRIC.F1_Score'][Region_index[i]]]

    Headers = ['Region', 'Sample', 'Sequencing technology', 'Total truth', 'True positives', 'False negatives', 'False positives', 'Recall', 'Precision', 'F1-score']
    print(','.join(Headers))
    print(','.join([Regions[i], 'HG003', 'ONT', str(HG003_ONT_DV['TRUTH.TOTAL'][Region_index[i]]), str(HG003_ONT_DV['TRUTH.TP'][Region_index[i]]), str(HG003_ONT_DV['TRUTH.FN'][Region_index[i]]), str(HG003_ONT_DV['QUERY.FP'][Region_index[i]]), str(HG003_ONT_DV['METRIC.Recall'][Region_index[i]]), str(HG003_ONT_DV['METRIC.Precision'][Region_index[i]]), str(HG003_ONT_DV['METRIC.F1_Score'][Region_index[i]])]))
    print(','.join([Regions[i], 'HG003', 'Illumina',        str(HG003_ILM_DV['TRUTH.TOTAL'][Region_index[i]]), str(HG003_ILM_DV['TRUTH.TP'][Region_index[i]]), str(HG003_ILM_DV['TRUTH.FN'][Region_index[i]]), str(HG003_ILM_DV['QUERY.FP'][Region_index[i]]), str(HG003_ILM_DV['METRIC.Recall'][Region_index[i]]), str(HG003_ILM_DV['METRIC.Precision'][Region_index[i]]), str(HG003_ILM_DV['METRIC.F1_Score'][Region_index[i]])]))
    print(','.join([Regions[i], 'HG003', 'PacBio',     str(HG003_CCS_DV['TRUTH.TOTAL'][Region_index[i]]), str(HG003_CCS_DV['TRUTH.TP'][Region_index[i]]), str(HG003_CCS_DV['TRUTH.FN'][Region_index[i]]), str(HG003_CCS_DV['QUERY.FP'][Region_index[i]]), str(HG003_CCS_DV['METRIC.Recall'][Region_index[i]]), str(HG003_CCS_DV['METRIC.Precision'][Region_index[i]]), str(HG003_CCS_DV['METRIC.F1_Score'][Region_index[i]])]))
    print(','.join([Regions[i], 'HG004', 'ONT', str(HG004_ONT_DV['TRUTH.TOTAL'][Region_index[i]]), str(HG004_ONT_DV['TRUTH.TP'][Region_index[i]]), str(HG004_ONT_DV['TRUTH.FN'][Region_index[i]]), str(HG004_ONT_DV['QUERY.FP'][Region_index[i]]), str(HG004_ONT_DV['METRIC.Recall'][Region_index[i]]), str(HG004_ONT_DV['METRIC.Precision'][Region_index[i]]), str(HG004_ONT_DV['METRIC.F1_Score'][Region_index[i]])]))
    print(','.join([Regions[i], 'HG004', 'Illumina',        str(HG004_ILM_DV['TRUTH.TOTAL'][Region_index[i]]), str(HG004_ILM_DV['TRUTH.TP'][Region_index[i]]), str(HG004_ILM_DV['TRUTH.FN'][Region_index[i]]), str(HG004_ILM_DV['QUERY.FP'][Region_index[i]]), str(HG004_ILM_DV['METRIC.Recall'][Region_index[i]]), str(HG004_ILM_DV['METRIC.Precision'][Region_index[i]]), str(HG004_ILM_DV['METRIC.F1_Score'][Region_index[i]])]))
    print(','.join([Regions[i], 'HG004', 'PacBio',     str(HG004_CCS_DV['TRUTH.TOTAL'][Region_index[i]]), str(HG004_CCS_DV['TRUTH.TP'][Region_index[i]]), str(HG004_CCS_DV['TRUTH.FN'][Region_index[i]]), str(HG004_CCS_DV['QUERY.FP'][Region_index[i]]), str(HG004_CCS_DV['METRIC.Recall'][Region_index[i]]), str(HG004_CCS_DV['METRIC.Precision'][Region_index[i]]), str(HG004_CCS_DV['METRIC.F1_Score'][Region_index[i]])]))

    # Set colors
    ONT_COLOR = [(238/256.0, 39/256.0, 8/256.0, 0.7), (238/256.0, 39/256.0, 8/256.0, 0.8), (238/256.0, 39/256.0, 8/256.0, 1.0)]
    ILM_COLOR = [(238/256.0, 124/256.0, 8/256.0, 0.7), (238/256.0, 124/256.0, 8/256.0, 0.8), (238/256.0, 124/256.0, 8/256.0, 1.0)]
    CCS_COLOR = [(10/256.0, 118/256.0, 148/256.0, 0.7), (10/256.0, 118/256.0, 148/256.0, 0.8), (10/256.0, 118/256.0, 148/256.0, 1.0)]

    HG003_COLOR = ONT_COLOR + ILM_COLOR + CCS_COLOR
    HG004_COLOR = ONT_COLOR + ILM_COLOR + CCS_COLOR

    # Bar charts
    N = 4
    ind = np.arange(N)
    width = 0.98
    bar_positions2 = [0, 1, 2, 3.25, 4.25, 5.25, 6.50, 7.50, 8.50]
    bar_positions1 = [10.50, 11.50, 12.50, 13.75, 14.75, 15.75, 17, 18, 19]

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
    p1 = axes[i].barh(bar_positions1, HG003_SCORES, width, color = HG003_COLOR)
    axes[i].hlines(9.50, x_lim_min, x_lim_max, colors='black', linestyles='solid')
    p2 = axes[i].barh(bar_positions2, HG004_SCORES, width, color = HG004_COLOR)



    # Add text labels
    # plt.text(6, -950, "Assembler", fontsize = text_fontsize+3)
    axes[i].set_xlabel('Value', fontsize = text_fontsize + 1, color='black')
    Category = ["R", "P", "F1"]


    value_offset = (x_lim_max - x_lim_min) / 10
    for ii, v in enumerate(HG003_SCORES):
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
    for ii, v in enumerate(HG004_SCORES):
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

    if i == len(Regions) - 1:
        caller_name_offser = x_lim_max + (0.08 * (1-x_lim_max)) # -0.13 for 0.0
        sample_name_offset = x_lim_max + (0.10 * (1-x_lim_max)) # -0.07 for 0.0
        samples = ["HG003 (SNP)", "HG004 (SNP)"]
        axes[i].text(sample_name_offset, bar_positions1[0] + ((bar_positions1[-1] - bar_positions1[0]) / 2), samples[0], va='center', fontsize=text_fontsize, rotation=270)
        axes[i].text(sample_name_offset, bar_positions2[0] + ((bar_positions2[-1] - bar_positions2[0]) / 2), samples[1], va='center', fontsize=text_fontsize, rotation=270)

    if i == 0:
        categories = ["Nanopore", "Illumina", "PacBio\nHiFi"]
        mid_points = [bar_positions1[1], bar_positions1[4], bar_positions1[7],
                      bar_positions2[1], bar_positions2[4], bar_positions2[7]]
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

        for ii, v in enumerate(HG003_SCORES):
            axes[i].text(x_lim_min + 0.001, bar_positions2[ii], categories[ii % 3], color='white', va='center', ha='left', fontsize=text_fontsize)
            axes[i].text(x_lim_min + 0.001, bar_positions1[ii], categories[ii % 3], color='white', va='center', ha='left', fontsize=text_fontsize)

    axes[i].tick_params(axis=u'both', which=u'both',length=0)
    axes[i].set_xlim(x_lim_min, x_lim_max)
    axes[i].set_title(Regions[i], fontsize=text_fontsize + 1)

    x_ticks = np.arange(x_lims[i] + x_tick_frq[i], 1.01, x_tick_frq[i] * 2)
    x_lables = ["%.2lf" % round(x,2) for x in x_ticks]
    axes[i].set_xticks(x_ticks)
    axes[i].set_xticklabels(x_lables, fontsize=text_fontsize, color='black')
    # Save the figure
plt.tight_layout()
# plt.margins(0, 0)
plt.show()
# plt.savefig('./outputs/2b_Illumina_difficult.pdf', format='pdf', dpi=300)
