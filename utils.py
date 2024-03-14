import matplotlib
import matplotlib.pyplot as plt
import collections
from matplotlib import colors
import pandas as pd
import numpy as np
import seaborn as sns
# from statannotations.Annotator import Annotator
from scipy import stats
import statistics
from matplotlib import rc  # used to change font


def PiePlot_df1_byPatient(data1, title, order, angle=220, save=False, nolegend=False, by="MVB", casetype=None, subtype=None):

    if not casetype is None:
        data1 = data1.loc[data1['CaseType'] == casetype, :]

    if not subtype is None:
        data1 = data1.loc[data1['subtype'] == subtype, :]

    fc_p1 = data1.groupby([by]).apply(lambda x: x['phenotype_combined'].value_counts(
        sort=False, normalize=True).rename_axis('Type').reset_index(name='Percent'))

    df_fc1 = fc_p1.reset_index()

    df_percent1 = df_fc1.set_index([by, 'Type'])['Percent'].unstack(
    ).reset_index()  # pivot dataframe unstack used to reshape dataframe

    df_percent1 = df_percent1.set_index(by)

    percents_ave = df_percent1.mean()

    s = percents_ave[order]
    new_percents = np.array(s.values)
    new_labels = ["{:.2%}".format(ele) for ele in s]

    plt.pie(new_percents, startangle=angle, labels=new_labels, textprops={'fontsize': 11},
            colors=[get_marker_color(m)
                    for m in order],  # 'tab:blue', 'tab:pink',
            wedgeprops={"edgecolor": "black",
                        'linewidth': 0.1,
                        'antialiased': True})
    plt.ylabel('')

    if not nolegend:

        legend_labels = [get_marker_label(
            marker, ret=False) for marker in order]

        # legend_labels = ['B cells(CD20)', '$\mathregular{CD4^{+}}$ T cells',
        #               '$\mathregular{CD8^{+}}$ T cells','Treg cells (FOXP3)',
        #               'Macrophages (CD68)', 'Dendritic cells (CD11c)',
        #               'NK cells (CD56)', 'Pan CK cells', 'Unidentified cells']

        plt.legend(bbox_to_anchor=(1.9, 1), loc='upper right',
                   prop={'size': 11},
                   # labels=["{0}, {1:.2%}".format(l, s) for l, s in zip(new_labels, new_percents)])
                   labels=legend_labels,
                   frameon=False,
                   labelspacing=1)

    plt.title(title)

    if save:
        # plt.tight_layout()
        # plt.subplots_adjust(left=0.3, right=1, bottom=0.1, top=0.2)
        plt.savefig(
            f"../../Paper/Figure/Composition/PiePlot {title}.pdf", format='pdf', bbox_inches="tight", dpi=500)

    plt.show()


def PiePlot_combined_byPatient(data1, data2, title, order, angle=220, save=False, nolegend=False, by="MVB", casetype=None, subtype=None):

    # change_font()

    if not casetype is None:
        data1 = data1.loc[data1['CaseType'] == casetype, :]
        data2 = data2.loc[data2['CaseType'] == casetype, :]

    if not subtype is None:
        data1 = data1.loc[data1['subtype'] == subtype, :]
        data2 = data2.loc[data2['subtype'] == subtype, :]

    fc_p1 = data1.groupby([by]).apply(lambda x: x['phenotype_combined'].value_counts(
        sort=False, normalize=True).rename_axis('Type').reset_index(name='Percent'))
    fc_p2 = data2.groupby([by]).apply(lambda x: x['phenotype_combined'].value_counts(
        sort=False, normalize=True).rename_axis('Type').reset_index(name='Percent'))

    df_fc1 = fc_p1.reset_index()
    df_fc2 = fc_p2.reset_index()

    df_percent1 = df_fc1.set_index([by, 'Type'])['Percent'].unstack(
    ).reset_index()  # pivot dataframe unstack used to reshape dataframe
    df_percent2 = df_fc2.set_index([by, 'Type'])['Percent'].unstack(
    ).reset_index()  # colnames ['MVB', 'CD11c', 'CD56', 'Unidentified']
    df_percent2.drop(['other'], axis=1, inplace=True)  # delete column other
    # combine df_percent and df_percent2, deduct CD56 and CD11c percents from other in df_percent

    df_percent = pd.merge(df_percent1, df_percent2, on=by)

    df_percent['Unidentified'] = df_percent['Unidentified'] - \
        df_percent['CD56'] - df_percent['CD11c']
    # remove MVB that Unidentified < 0
    df_percent = df_percent.loc[df_percent['Unidentified'] >= 0, :]

    df_percent = df_percent.set_index(by)
    
    # df_percent.to_csv("../../Paper/Figure/NEW_ORGANIZATION/Phenotype/Composition/phenotype_percent_MVB.csv")

    percents_ave = df_percent.mean()

    s = percents_ave[order]
    new_percents = np.array(s.values)
    new_labels = ["{:.2%}".format(ele) for ele in s]

    plt.pie(new_percents, startangle=angle, labels=new_labels, textprops={'fontsize': 11},
            colors=[get_marker_color(m) for m in order],
            wedgeprops={"edgecolor": "black",
                        'linewidth': 0,
                        'antialiased': True})
    plt.ylabel('')

    if not nolegend:

        legend_labels = ['B cells(CD20)', '$\mathregular{CD4^{+}}$ T cells',
                         '$\mathregular{CD8^{+}}$ T cells', 'Treg cells (FOXP3)',
                         'Macrophages (CD68)', 'Dendritic cells (CD11c)',
                         'NK cells (CD56)', 'Pan CK cells', 'Unidentified cells']

        plt.legend(bbox_to_anchor=(1.9, 1), loc='upper right',
                   prop={'size': 11},
                   # fontname='Arial',
                   # labels=["{0}, {1:.2%}".format(l, s) for l, s in zip(new_labels, new_percents)])
                   labels=legend_labels,
                   frameon=False,
                   labelspacing=1)
    csfont = {'fontname': 'Arial'}
    plt.title(title, **csfont)

    if save:
        # plt.tight_layout()
        # plt.subplots_adjust(left=0.3, right=1, bottom=0.1, top=0.2)
        plt.savefig(
            f"../../Paper/Figure/Composition/PiePlot {title}.pdf", format='pdf', bbox_inches="tight", dpi=500)
        # plt.savefig(f"../../Paper/Figure/Composition/PiePlot {title}.png", format='png', bbox_inches = "tight", dpi=500)

    plt.show()


def PiePlot_combined(types1, types2, title, order, angle=220, save=True):
    type_count1 = pd.Series(types1).value_counts(normalize=True)
    type_count2 = pd.Series(types2).value_counts(normalize=True)

    newPercentOther = type_count1['Unidentified'] - \
        type_count2['CD11c'] - type_count2['CD56']

    type_count = pd.concat([type_count1, type_count2[['CD11c', 'CD56']]])
    # type_count.drop(['Unidentified'], inplace=True) #drop the old one
    type_count['Unidentified'] = newPercentOther

    s = type_count[order]
    # s.sort_values(inplace=True)
    new_percents = np.array(s.values)
    # new_labels = np.array(s.index)
    new_labels = ["{:.2%}".format(ele) for ele in s]
    legend_labels = ['B cells(CD20)', '$\mathregular{CD4^{+}}$ T cells',
                     '$\mathregular{CD8^{+}}$ T cells', 'Treg cells (FOXP3)',
                     'Macrophages (CD68)', 'Dendritic cells (CD11c)',
                     'NK cells (CD56)', 'Pan CK cells', 'Unidentified cells']

    plt.pie(new_percents, startangle=angle, labels=new_labels, textprops={'fontsize': 11},
            colors=['tab:brown', 'tab:olive', 'tab:purple', 'tab:red',
                    'tab:green', 'tab:blue', 'tab:pink', 'tab:orange', 'tab:gray'],
            wedgeprops={"edgecolor": "black",
                        'linewidth': 0.1,
                        'antialiased': True})
    plt.ylabel('')

    plt.legend(bbox_to_anchor=(1.9, 1), loc='upper right',
               prop={'size': 11},
               # labels=["{0}, {1:.2%}".format(l, s) for l, s in zip(new_labels, new_percents)])
               labels=legend_labels,
               frameon=False,
               labelspacing=1)

    plt.title(title)
    if save:
        # plt.tight_layout()
        # plt.subplots_adjust(left=0.3, right=1, bottom=0.1, top=0.2)
        plt.savefig(
            f"../figure/composite/PiePlot {title}.pdf", format='pdf', bbox_inches="tight", dpi=500)


def PiePlot(types, title, angle=0, save=True):
    type_count = collections.Counter(types)
    labels = type_count.keys()
    sizes = type_count.values()
    percents = np.array(list(sizes)) * 100 / np.array(list(sizes)).sum()

    # plt.pie(sizes,autopct='%1.1f%%', startangle=0)

    # plt.legend(labels, bbox_to_anchor=(0.8,0.7),  fontsize=10,
    #         bbox_transform=plt.gcf().transFigure)

    plt.pie(percents, startangle=angle,
            labels=labels, textprops={'fontsize': 7})
    plt.ylabel('')
    #--------------------------------------------------------------------

    plt.legend(bbox_to_anchor=(1.7, 0.9), loc='upper right',
               prop={'size': 8},
               labels=['%s, %1.1f %%' % (l, s) for l, s in zip(labels, percents)])
    plt.title(title)

    if save:
        # plt.tight_layout()
        # plt.subplots_adjust(left=0.3, right=1, bottom=0.1, top=0.2)
        plt.savefig(f"../figure/PiePlot {title}.pdf",
                    format='pdf', bbox_inches="tight")
    plt.title(title)

    if save:
        # plt.tight_layout()
        # plt.subplots_adjust(left=0.3, right=1, bottom=0.1, top=0.2)
        plt.savefig(f"../figure/PiePlot {title}.pdf",
                    dformat='pdf', bbox_inches="tight")


def PiePlot2(types, title, save=True):
    type_count = collections.Counter(types)
    labels = type_count.keys()
    sizes = type_count.values()
    plt.pie(sizes, startangle=0, labels=labels)
    plt.ylabel('')
    percents = np.array(list(sizes)) * 100 / np.array(list(sizes)).sum()
    plt.legend(bbox_to_anchor=(1.35, 1.1), loc='upper right',
               labels=['%s, %1.1f %%' % (l, s) for l, s in zip(labels, percents)])


def GetTypeCompForCol(data_merged, col, phenotype):
    fc = data_merged.groupby(col).apply(lambda x: x[phenotype].value_counts(
        sort=False, normalize=False).rename_axis('Type').reset_index(name='Counts'))
    fc_p = data_merged.groupby(col).apply(lambda x: x[phenotype].value_counts(
        sort=False, normalize=True).rename_axis('Type').reset_index(name='Percent'))
    fc['Percent'] = fc_p['Percent']
    df_fc = fc.reset_index()
    df_melt = df_fc.drop(['level_1'], axis=1)
# those above lines can be simplified to the following, will do the same, not need to reset_index twice
    # fc = data_merged.groupby(col).apply(lambda x: x['Cluster'].value_counts(sort=False, normalize=False))
    # fc_p = data_merged.groupby(col).apply(lambda x: x['Cluster'].value_counts(sort=False, normalize=True))
    # fc.columns = [[col, 'Type', 'Counts']]
    # fc_p.columns = [[col, 'Type', 'Percent']]
    # fc['Percent'] = fc_p['Percent']
    # df_fc = fc.reset_index()

    # df_melt.to_csv("./output/Type_composition_melted_"+col+".csv", index=False)

    df_count = df_melt.set_index([col, 'Type'])[
        'Counts'].unstack().reset_index()
    df_count.index = df_count[col]
    df_count.drop([col], axis=1, inplace=True)

    # df_count.to_csv("./output/Type_composition_"+col.lower()+"_counts.csv", index=False)

    df_percent = df_melt.set_index([col, 'Type'])[
        'Percent'].unstack().reset_index()
    df_percent.index = df_percent[col]
    df_percent.drop([col], axis=1, inplace=True)

    # df_percent.to_csv("./output/Type_composition_"+col.lower()+"_percent.csv", index=False)

    return df_fc, df_count, df_percent


def get_marker_color(marker):
    marker_color = {"CD8": "#E41A1C",
                    "CK": "#377EB8",
                    "CD20": "#4DAF4A",
                    "CD4": "#984EA3",
                    "FOXP3": "#FF7F00",
                    "CD68": "#FFFF33",
                    "CD11c": "#A65628",
                     "CD56": "#F781BF",
                     "Unidentified": "#999999"
                    }
    return(marker_color[marker])


def get_marker_label(marker, ret=True, short=True):

    marker_label_ret_short = {"CD4": "$\mathregular{CD4^{+}}$\nT cells",
                        "CD8": "$\mathregular{CD8^{+}}$\nT cells",
                        "CD20": "B cells\n(CD20)",
                        "FOXP3": "Treg cells\n(FOXP3)",
                        "CD68": "MΦ\n(CD68)",
                        "CK": "Pan CK\ncells",
                        "CD11c": "Dendritic\n(CD11c)",
                        "CD56": "NK cells\n(CD56)",
                        "Unidentified": "Unidentified\ncells"}
    
    marker_label_ret = {"CD4": "$\mathregular{CD4^{+}}$\nT cells",
                        "CD8": "$\mathregular{CD8^{+}}$\nT cells",
                        "CD20": "B cells\n(CD20)",
                        "FOXP3": "Treg cells\n(FOXP3)",
                        "CD68": "Macrophages\n(CD68)",
                        "CK": "Pan CK\ncells",
                        "CD11c": "Dendritic\n(CD11c)",
                        "CD56": "NK cells\n(CD56)",
                        "Unidentified": "Unidentified\ncells"}

    marker_label = {"CD4": "$\mathregular{CD4^{+}}$ T cells",
                    "CD8": "$\mathregular{CD8^{+}}$ T cells",
                    "CD20": "B cells (CD20)",
                    "FOXP3": "Treg cells (FOXP3)",
                    "CD68": "Macrophages (CD68)",
                    "CK": "Pan CK cells",
                    "CD11c": "Dendritic (CD11c)",
                        "CD56": "NK cells (CD56)",
                        "Unidentified": "Unidentified cells"}
    if not marker in marker_label.keys():
        return(marker)
    else:
        if ret:
           if short:
               return(marker_label_ret_short[marker])
           else:
               return(marker_label_ret[marker])

        else:
            return(marker_label[marker])


def BarPlot(df_plot, title, angle=270, save=True, changeColor=False):

    from matplotlib import cm

    if changeColor:
        # , rot=270 ) add this do the same as rotation xticks
        df_plot.plot.bar(stacked=True, cmap=cm.get_cmap('tabS'))
    else:
        df_plot.plot.bar(stacked=True)
    plt.xticks(rotation=angle, horizontalalignment="center")
    plt.legend(bbox_to_anchor=(1.0, 1.0))
    plt.title(title)
    plt.xlabel("")

    if save:
        plt.savefig(f"../figure/BarPlot {title}.pdf",
                    format='pdf', bbox_inches="tight")


palette_mapper = {
    'CD4': colors.to_hex('seagreen'),
    'CD8': colors.to_hex('salmon'),
    'CD20': colors.to_hex('lightskyblue'),
    'CD68': colors.to_hex('red'),
    'FOXP3': colors.to_hex('magenta'),
    'CK': colors.to_hex('orange'),
    'Unidentified': colors.to_hex('blueviolet'),
}


def make_palette(column, groups):
    return [colors.to_rgb(column.color_picker(g, palette_mapper.get(g.lower(), colors.to_hex('black')))) for g in groups]


def label(x, color, label):
    ax = plt.gca()
    ax.text(0, .2, label, color='black', fontsize=35,
            ha="left", va="center", transform=ax.transAxes)
    ax.tick_params(axis='x', labelsize=30)


def ridge_plot(df_plt: pd.DataFrame, xlabel: str, group: str = 'phenotype_CD8'):
    groups_ = np.sort(np.unique(df_plt[group].values))
    groups_ = np.array([i.lower() for i in groups_])

    # ## dirty sorting
    # groups = []
    # for g in severity_order:
    #     if g in groups_:
    #         groups.append(g)

    groups = groups_

    palette = sns.color_palette('pastel')

    g = sns.FacetGrid(df_plt,
                      palette=palette,
                      row=group, hue=group, aspect=9, height=1.2)

    g.map_dataframe(sns.kdeplot, x=xlabel, fill=True, alpha=1)
    g.map_dataframe(sns.kdeplot, x=xlabel, color='black')

    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, color='black', fontsize=13,
                ha="left", va="center", transform=ax.transAxes)

    g.map(label, group)
    g.fig.subplots_adjust(hspace=-.15)
    g.set_titles("")
    g.set(yticks=[], xlabel=xlabel)
    g.despine(left=True)
    g.set(yticks=[], ylabel="")
    # g.despine(bottom=True, left=True)

    # plt.xlabel(xlabel, fontsize=10, labelpad=5)
    # plt.gcf().subplots_adjust(bottom=0.25)
    # plt.tight_layout()


# sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0), 'axes.linewidth':2})
# palette = sns.color_palette("Set2", 12)
# g = sns.FacetGrid(df_filtered, palette=palette, row="Language", hue="Language", aspect=9, height=1.2)
# g.map_dataframe(sns.kdeplot, x="IMDB Score", fill=True, alpha=1)
# g.map_dataframe(sns.kdeplot, x="IMDB Score", color='black')
# def label(x, color, label):
#     ax = plt.gca()
#     ax.text(0, .2, label, color='black', fontsize=13,
#             ha="left", va="center", transform=ax.transAxes)

# g.map(label, "Language")
# g.fig.subplots_adjust(hspace=-.5)
# g.set_titles("")
# g.set(yticks=[], xlabel="IMDB Score")
# g.despine( left=True)
# plt.suptitle('Netflix Originals - IMDB Scores by Language', y=0.98)

    return g


def violin_plot(title, data_group, x, y, hue):
    fig, ax = plt.subplots(figsize=(25, 8))

    # if x=="Celltype":
    #     sns.violinplot(x=x, y=y, data=data_group, ax=ax, order=[ "LT-HSC HLF","HSC HIST1H2AC","HSC WNT11","HSC CACNB2","HSC MYADM-CD97","ST-HSC PBX1","LMPP CDK6-FLT3","MPP SPINK2-CD99","MPP Ribo-high","pre-MEP","MEP-MKP","ERP","BMCP","ML-Gran","MultiLin-ATAC","pre-Gran CP","LMPP PRSSI","LMPP LSAMP","LMPP Naive T-cell"])
    # else:
    #     sns.violinplot(x=x, y=y, data=data_group, ax=ax, order=["ND251_34","ND251_HS","ND251_LMPP","ND251_MPP"])

    sns.violinplot(x=x, y=y, hue=hue, data=data_group, ax=ax)

    ax.set_title("{} ".format(title),  fontsize=26)
    ax.set_ylabel("density", fontsize=20)
    ax.set_xlabel("")
    plt.xticks(rotation=45, fontsize=16)
    plt.setp(ax.xaxis.get_majorticklabels(), ha='right')
    return (fig)


def violin_plot_plus(title, data, x, y, hue, markers, hue_order, colors=None, bbox_to_anchor=None, save=True, nolegend=False):
     sns.set_theme(style='white')

     colors = [get_color_code(it) for it in hue_order]
     ax = sns.violinplot(
        data=data,
        x=x, y=y,
        hue=hue,
        order=markers,
        palette=colors,
        hue_order=hue_order
     )

    # Get the legend from just the bar chart
     handles, labels = ax.get_legend_handles_labels()

     ax.set_title(title,  fontsize=30)
     # ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
     # ax.set_ylabel("log(density)")
     if y == "density":
        ax.set_ylabel(
            "Number of cells per $\mathregular{mm^{2}}$", fontsize=30)
     if y == "percent":
        ax.set_ylabel("Number of cells per total cells", fontsize=30)
     ax.set(xlabel=None)
    # ax.set_ylim([0,10])
     ax.tick_params(axis='both', which='major', labelsize=25)
     ax.legend_.remove()

     if not nolegend:
        ax.legend(
            handles,
            labels,
            loc=(1.02, 0.8),
            fontsize="30",
            bbox_to_anchor=bbox_to_anchor,
            frameon=False

        )

     # labels = [item.get_text() for item in ax.get_xticklabels()]
     labels = [get_marker_label(marker, ret=True) for marker in markers]
     # labels = ['B cells\n(CD20)', '$\mathregular{CD4^{+}}$\nT cells',
     #               '$\mathregular{CD8^{+}}$\nT cells','Treg cells\n(FOXP3)',
     #               'Macrophages\n(CD68)', 'Dendritic\n(CD11c)',
     #               'NK cells\n(CD56)', 'Pan CK\ncells', 'Unidentified\ncells']

     ax.set_xticklabels(labels)

     # ax.spines[['right', 'top']].set_visible(False)
     ax.spines['top'].set_visible(False)
     ax.spines['right'].set_visible(False)
     plt.show()

     fig = ax.get_figure()
     # fig = plt.gcf()
     fig.set_size_inches(23, 10.5)

     # fig = ax.get_figure()
     # fig(figsize=(8, 6), dpi=80)
     if save:
         if ('Tumor' in title) or ('Stroma' in title):
             fig.savefig(
                 f"../../Paper/Figure/BarPlot/Violin_Tumor&Stroma/{title}.pdf", format='pdf',  bbox_inches="tight",  dpi=500)
         else:
             fig.savefig(
                 f"../../Paper/Figure/BarPlot/Violin/{title}.pdf", format='pdf',  bbox_inches="tight",  dpi=500)


def box_plot(outpath, title, data, x, y, hue, markers, hue_order, colors=None, bbox_to_anchor=None, save=True, nolegend=False, size=(30, 10.5)):
    sns.set_theme(style='white')

    colors = [get_color_code(it) for it in hue_order]
    ax = sns.boxplot(
        data=data,
        x=x, y=y,
        hue=hue,
        order=markers,
        palette=colors,
        hue_order=hue_order
    )

    # Get the legend from just the bar chart
    handles, labels = ax.get_legend_handles_labels()

    ax.set_title(title,  fontsize=30)

    if y == "density":
         ax.set_ylabel(
             "Number of cells per $\mathregular{mm^{2}}$", fontsize=30)
    if y == "percent":
         ax.set_ylabel("Number of cells per total cells", fontsize=30)
    ax.set(xlabel=None)
     # ax.set_ylim([0,10])
    ax.tick_params(axis='both', which='major', labelsize=25)
    ax.legend_.remove()

    if not nolegend:
         ax.legend(
             handles,
             labels,
             loc=(1.02, 0.8),
             fontsize="30",
             bbox_to_anchor=bbox_to_anchor,
             frameon=False

         )

    labels = [get_marker_label(marker, ret=True) for marker in markers]
    ax.set_xticklabels(labels)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.show()

    fig = ax.get_figure()

    fig.set_size_inches(size)

    if save:
        fig.savefig(f"{outpath}/{title}.pdf", format='pdf', bbox_inches="tight",  dpi=500)
                        



def box_plot_plus(outpath, title, data, x, y, hue, markers, hue_order, colors=None, bbox_to_anchor=None, save=False, nolegend=False, size=(30, 10.5), short=True):
    sns.set_theme(style='white')

    colors = [get_color_code(it) for it in hue_order]
    ax = sns.boxplot(
        data=data,
        x=x, y=y,
        hue=hue,
        order=markers,
        palette=colors,
        hue_order=hue_order
    )

    # Get the legend from just the bar chart
    handles, labels = ax.get_legend_handles_labels()

    # Draw the stripplot
    sns.stripplot(
        data=data,
        x=x, y=y,
        hue=hue,
        dodge=True,
        jitter=True,
        ax=ax,
        palette=colors,
        # color="white",
         edgecolor="#484f4d",  # "#808A87",   #"white",
         marker='o', size=2, linewidth=1,  # marker='o', size=7, linewidth=0.5, #,
         alpha=0.5,
         order=markers,
         hue_order=hue_order
        )

    ax.set_title(title,  fontsize=20)
     # ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
     # ax.set_ylabel("log(density)")
    # if y == "density":
    #     ax.set_ylabel(
    #         "Number of cells per $\mathregular{mm^{2}}$", fontsize=20)
    # if y == "percent":
    #     ax.set_ylabel("Number of cells / total cells", fontsize=20)
    ax.set(xlabel=None)
    # ax.set_ylim([0,10])
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.legend_.remove()

    if not nolegend:
        ax.legend(
            handles,
            labels,
            loc=(1.02, 0.8),
            fontsize="30",
            bbox_to_anchor=bbox_to_anchor,
            frameon=False

        )

    # labels = [item.get_text() for item in ax.get_xticklabels()]
    labels = [get_marker_label(marker, ret=True, short=short) for marker in markers]
    # labels = ['B cells\n(CD20)', '$\mathregular{CD4^{+}}$\nT cells',
    #               '$\mathregular{CD8^{+}}$\nT cells','Treg cells\n(FOXP3)',
    #               'Macrophages\n(CD68)', 'Dendritic\n(CD11c)',
    #               'NK cells\n(CD56)', 'Pan CK\ncells', 'Unidentified\ncells']

    ax.set_xticklabels(labels)

    # ax.spines[['right', 'top']].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.show()

    fig = ax.get_figure()
    # fig = plt.gcf()
    fig.set_size_inches(size)

    # fig = ax.get_figure()
    # fig(figsize=(8, 6), dpi=80)
    if save:
        fig.savefig(f"{outpath}/{title}.pdf", format='pdf',
                        bbox_inches="tight",  dpi=500)
                


def box_plot_plus_horizontal(outpath, title, data, x, y, hue, markers, hue_order, colors=None, bbox_to_anchor=None, save=True, nolegend=True, size=(30, 10.5), short=False):
    sns.set_theme(style='white')

    colors = [get_color_code(it) for it in hue_order]
    ax = sns.boxplot(
        data=data,
        x=x, y=y,
        hue=hue,
        order=markers,
        palette=colors,
        hue_order=hue_order
    )

    # Get the legend from just the bar chart
    handles, labels = ax.get_legend_handles_labels()

    # Draw the stripplot
    sns.stripplot(
        data=data,
        x=x, y=y,
        hue=hue,
        dodge=True,
        jitter=True,
        ax=ax,
        palette=colors,
        # color="white",
         edgecolor="#484f4d",  # "#808A87",   #"white",
         marker='o', size=7, linewidth=1,  # marker='o', size=7, linewidth=0.5, #,
         alpha=0.5,
         order=markers,
         hue_order=hue_order
        )

    ax.set_title(title,  fontsize=20)
     # ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
     # ax.set_ylabel("log(density)")
    if x == "density":
         ax.set_xlabel(
             "Number of cells per $\mathregular{mm^{2}}$", fontsize=10)
    if x == "percent":
         ax.set_xlabel("Number of cells / total cells", fontsize=10)
    ax.set(ylabel=None)
     # ax.set_ylim([0,10])
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.legend_.remove()

    if not nolegend:
         ax.legend(
             handles,
             labels,
             loc=(1.02, 0.8),
             fontsize="30",
             bbox_to_anchor=bbox_to_anchor,
             frameon=False

         )

     # labels = [item.get_text() for item in ax.get_xticklabels()]
    labels = [get_marker_label(marker, ret=True, short=short) for marker in markers]
     # labels = ['B cells\n(CD20)', '$\mathregular{CD4^{+}}$\nT cells',
     #               '$\mathregular{CD8^{+}}$\nT cells','Treg cells\n(FOXP3)',
     #               'Macrophages\n(CD68)', 'Dendritic\n(CD11c)',
     #               'NK cells\n(CD56)', 'Pan CK\ncells', 'Unidentified\ncells']

    ax.set_yticklabels(labels)

     # ax.spines[['right', 'top']].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.show()

    fig = ax.get_figure()
     # fig = plt.gcf()
    fig.set_size_inches(size)

     # fig = ax.get_figure()
     # fig(figsize=(8, 6), dpi=80)
    if save:
         if ('Tumor' in title) or ('Stroma' in title):
             fig.savefig(
                 f"../../Paper/Figure/BarPlot/BoxPlot_Tumor&Stroma/{title}.pdf", format='pdf',  bbox_inches="tight",  dpi=500)
         else:
             fig.savefig(f"{outpath}/{title}.pdf", format='pdf',
                         bbox_inches="tight",  dpi=500)
             fig.savefig(f"{outpath}/{title}.png", format='png',
                         bbox_inches="tight",  dpi=500)

def get_pvalue_pairwise(x, y, z, non_parametric=True):
    if non_parametric == False:
        if np.nanmedian(x) >= np.nanmedian(y):
            alt = "greater"
        else:
            alt = "less"
        xy = stats.ttest_ind(x, y, alternative=alt, equal_var=False)[1]
        
        if np.nanmedian(x) >= np.nanmedian(z):
            alt = "greater"
        else:
            alt = "less"
        xz = stats.ttest_ind(x, z, alternative=alt, equal_var=False)[1]
        
        if np.nanmedian(y) >= np.nanmedian(z):
            alt = "greater"
        else:
            alt = "less"
        yz = stats.ttest_ind(y, z, alternative=alt, equal_var=False)[1]
     
    else:
        
        if np.nanmedian(x) >= np.nanmedian(y):
            alt = "greater"
        else:
            alt = "less"
        xy = stats.ranksums(x, y, alternative=alt)[1]
        
        if np.nanmedian(x) >= np.nanmedian(z):
            alt = "greater"
        else:
            alt = "less"
        xz = stats.ranksums(x, z, alternative=alt)[1]
        
        if np.nanmedian(y) >= np.nanmedian(z):
            alt = "greater"
        else:
            alt = "less"
        yz = stats.ranksums(y, z, alternative=alt)[1]
    
    ret = [xy, xz, yz]
    ret = ["{:.2E}".format(pvalue) for pvalue in ret]
    return (ret)
        

def get_pvalue(x, y, z=None, one_sided=True, non_parametric=True):

    if non_parametric == False:
        if z is None:
            if one_sided:
                if np.nanmedian(x) >= np.nanmedian(y):
                    alt = "greater"
                else:
                    alt = "less"
                pvalue = stats.ttest_ind(x, y, alternative=alt, equal_var=False)[1]
            else:
                pvalue = stats.ttest_ind(x, y, equal_var=False)[1]
        else:
            pvalue = stats.f_oneway(x, y, z)[1]  # , nan_policy='omit'
    else:
        if z is None:
            if one_sided:
                if np.nanmedian(x) >= np.nanmedian(y):
                    alt = "greater"
                else:
                    alt = "less"
                pvalue = stats.ranksums(x, y, alternative=alt)[1]
            else:
                pvalue = stats.ranksums(x, y)[1]
        else:
            pvalue = stats.kruskal(x, y, z)[1]
    return("{:.2E}".format(pvalue))


# def get_pvalue_oneSided(x,y, z=None, non_parametric=True):
#     if non_parametric == False:
#         if z is None:
#             pvalue_greater = stats.ttest_ind(x, y, alternative="greater")[1]
#             pvalue_less = stats.ttest_ind(x, y, alternative="less")[1]

#         else:
#             pvalue = stats.f_oneway(x,y,z)[1]  #, nan_policy='omit'
#             pvalue_greater = pvalue
#             pvalue_less = pvalue
#     else:
#         if z is None:

#             pvalue_greater = stats.ranksums(x, y, alternative="greater")[1]
#             pvalue_less = stats.ranksums(x, y, alternative="less")[1]

#         else:
#             pvalue = stats.kruskal(x,y,z)[1]
#             pvalue_greater = pvalue
#             pvalue_less = pvalue

#     pvalue_greater = ("{:.2E}".format(pvalue_greater))
#     pvalue_less = ("{:.2E}".format(pvalue_less))

#     return([pvalue_greater, pvalue_less])

  # casetype="subtype"
  # groups=hue_order

def min_max_norm(s):
    return((s - s.min()) / (s.max() - s.min()))


def generate_pvalues(df_plot, markers, casetype, groups, grouptype, outpath, title, non_parametric=True, save=False):
    pvalues = dict()
    pvalues_pairwise = dict()
    for marker in markers:
        x = df_plot.loc[(df_plot[casetype] == groups[0]) & (
            df_plot['phenotype'] == marker), grouptype].dropna().values
        y = df_plot.loc[(df_plot[casetype] == groups[1]) & (
            df_plot['phenotype'] == marker), grouptype].dropna().values
        if len(groups) == 3:
            z = df_plot.loc[(df_plot[casetype] == groups[2]) & (
                df_plot['phenotype'] == marker), grouptype].dropna().values
                 
            pvalues[marker] = get_pvalue(x, y, z, non_parametric=non_parametric) 
            pvalues_pairwise[marker] = get_pvalue_pairwise( x, y, z, non_parametric=non_parametric)
            
        elif len(groups) == 2:
            pvalues[marker] = get_pvalue(x, y, non_parametric=non_parametric)

    df_pvalue = pd.DataFrame(pvalues, index=[0])
    
    
    if save:
        if non_parametric:
            df_pvalue.to_csv(
                f"{outpath}/{title} pvalue non_parametric.csv", index=True)
        else:
            df_pvalue.to_csv(f"{outpath}/{title} pvalue.csv", index=True)


    if pvalues_pairwise != {}:
        df_pvalue_pairwise = pd.DataFrame(pvalues_pairwise, index=[groups[0] + " ~ " + groups[1], groups[0] + " ~ " + groups[2], groups[1] + " ~ " + groups[2]])
        if save:
            if non_parametric:
                df_pvalue_pairwise.to_csv(
                   f"{outpath}/{title} pvalue non_parametric pairwise.csv", index=True)
            else:
                df_pvalue_pairwise.to_csv(
                   f"{outpath}/{title} pvalue pairwise.csv", index=True)
                
    return(df_pvalue)



# this function is not used anymore
def generate_barplot_and_pvalues(df_feature, outpath, title, casetype, markers, hue_order, grouptype, non_parametric=True, size=(30, 10.5)):
    if casetype is None:
        df_plot = df_feature
        # generate_pvalues(df_plot, markers, "casetype", hue_order,
        #                  grouptype, outpath, title, non_parametric=non_parametric)
        # bar_plot_plus(title,df_plot, 'phenotype', grouptype, 'casetype', markers, hue_order=hue_order)
        box_plot_plus(outpath, title, df_plot, 'phenotype', grouptype,
                 'casetype', markers, hue_order=hue_order)
        # violin_plot_plus(title,df_plot, 'phenotype', grouptype, 'casetype', markers, hue_order=hue_order)

    else:
        df_plot = df_feature.loc[df_feature['casetype'] == casetype, :]
        df_plot = df_plot.loc[df_plot['subtype'].isin(hue_order), :]
        # generate_pvalues(df_plot, markers, "subtype", hue_order,
        #                  grouptype, outpath, title, non_parametric=non_parametric)
        # bar_plot_plus(title,df_plot, 'phenotype', grouptype, 'subtype', markers, hue_order=hue_order)
        box_plot_plus(outpath, title, df_plot, 'phenotype', grouptype,
                 'subtype', markers, hue_order=hue_order, size=size)
        # violin_plot_plus(title,df_plot, 'phenotype', grouptype, 'subtype', markers, hue_order=hue_order)


def get_color_code(item):
    colors = {"epithelioid": "#738fc0",
              "biphasic": "#7db988",
              "sarcomatoid": "#e59d76",
              "non epithelioid": "tab:red",
              "pleural": "#8f6798",
              "peritoneal": "#fef17c"
              }
    return(colors[item.lower()])


def bar_plot_plus(title, data, x, y, hue, markers, hue_order, colors=None, bbox_to_anchor=None, save=True, nolegend=False):

    sns.set_theme(style='white')

    colors = [get_color_code(it) for it in hue_order]
    # colors = ['tab:blue','tab:red'] # ['#00E5EE','#BF3EFF'] #Pleural epithelioid vs non epithelioid
    # colors = ['tab:blue','tab:orange'] #["#00E5EE","#FF3E96"] #Peritoneal epithlial vs biphasic
    # colors = ['darkcyan', 'darkkhaki']
    #['darkcyan', 'cornsilk'] #['#00E5EE','#BF3EFF'] #['C2', 'C3']##[ "#00E5EE","#FF3E96","#F0F8FF"] # #['darkcyan', 'cornsilk']#['gold', 'seagreen']# ['gold', 'teal'] #
    # Draw the bar chart
    ax = sns.barplot(
        data=data,
        x=x, y=y,
        hue=hue,
        alpha=1,
        ci=None,
        # errorbar=None,
        order=markers,
        palette=colors,
        edgecolor="black",
        hue_order=hue_order,
        # width=0.7

    )

    # Get the legend from just the bar chart
    handles, labels = ax.get_legend_handles_labels()

    # Draw the stripplot
    sns.stripplot(
        data=data,
        x=x, y=y,
        hue=hue,
        dodge=True,
        jitter=True,
        ax=ax,
        palette=colors,
        # color="white",
        edgecolor="#484f4d",  # "#808A87",   #"white",
        marker='o', size=7, linewidth=1,  # marker='o', size=7, linewidth=0.5, #,
        alpha=0.5,
        order=markers,
        hue_order=hue_order
    )

    ax.set_title(title,  fontsize=30)
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    # ax.set_ylabel("log(density)")
    if y == "density":
        ax.set_ylabel(
            "Number of cells per $\mathregular{mm^{2}}$", fontsize=20)
    if y == "percent":
        ax.set_ylabel("Number of cells per total cells", fontsize=20)
    ax.set(xlabel=None)
    # ax.set_ylim([0,10])
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.legend_.remove()

    if not nolegend:
        ax.legend(
            handles,
            labels,
            loc=(1.02, 0.8),
            fontsize="30",
            bbox_to_anchor=bbox_to_anchor,
            frameon=False

        )

    # labels = [item.get_text() for item in ax.get_xticklabels()]
    labels = [get_marker_label(marker, ret=True) for marker in markers]
    # labels = ['B cells\n(CD20)', '$\mathregular{CD4^{+}}$\nT cells',
    #               '$\mathregular{CD8^{+}}$\nT cells','Treg cells\n(FOXP3)',
    #               'Macrophages\n(CD68)', 'Dendritic\n(CD11c)',
    #               'NK cells\n(CD56)', 'Pan CK\ncells', 'Unidentified\ncells']

    ax.set_xticklabels(labels)

    # ax.spines[['right', 'top']].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.show()

    fig = ax.get_figure()
    # fig = plt.gcf()
    fig.set_size_inches(23, 10.5)

    # fig = ax.get_figure()
    # fig(figsize=(8, 6), dpi=80)
    if save:
        if ('Tumor' in title) or ('Stroma' in title):
            fig.savefig(
                f"../../Paper/Figure/BarPlot_Tumor&Stroma/{title}.pdf", format='pdf',  bbox_inches="tight",  dpi=500)
        else:
            fig.savefig(
                f"../../Paper/Figure/BarPlot/{title}.pdf", format='pdf',  bbox_inches="tight",  dpi=500)


def bar_plot_plus_withAnno(title, data, x, y, hue, markers=None, hue_order=None, colors=None, pvalues=None, bbox_to_anchor=None, save=True, nolegend=False):
    sns.set_theme(style='white')
    colors = ['#00E5EE', '#BF3EFF']  # Pleural epithelioid vs non epithelioid
    colors = ["#00E5EE", "#FF3E96"]  # Peritoneal epithlial vs biphasic
    colors = ['darkcyan', 'darkkhaki']
    #['darkcyan', 'cornsilk'] #['#00E5EE','#BF3EFF'] #['C2', 'C3']##[ "#00E5EE","#FF3E96","#F0F8FF"] # #['darkcyan', 'cornsilk']#['gold', 'seagreen']# ['gold', 'teal'] #
    # Draw the bar chart
    ax = sns.barplot(
        data=data,
        x=x, y=y,
        hue=hue,
        alpha=1,
        errorbar=None,
        order=markers,
        palette=colors,
        edgecolor="black",
        hue_order=hue_order,
        width=0.7

    )

    # Get the legend from just the bar chart
    handles, labels = ax.get_legend_handles_labels()

    # Draw the stripplot
    sns.stripplot(
        data=data,
        x=x, y=y,
        hue=hue,
        dodge=True,
        jitter=True,
        ax=ax,
        palette=colors,
        # color="white",
        edgecolor="#484f4d",  # "#808A87",   #"white",
        marker='o', size=7, linewidth=1,  # marker='o', size=7, linewidth=0.5, #,
        alpha=0.5,
        order=markers,
        hue_order=hue_order
    )

    ax.set_title(title,  fontsize=30)
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    # ax.set_ylabel("log(density)")
    ax.set_ylabel("Number of cells per $\mathregular{mm^{2}}$", fontsize=30)
    ax.set(xlabel=None)
    ax.set_ylim([0, 10])
    ax.tick_params(axis='both', which='major', labelsize=25)
    ax.legend_.remove()

    if not nolegend:
        ax.legend(
            handles,
            labels,
            loc=(1.05, 0.8),
            fontsize="30",
            # (1.6, .8) #bbox_to_anchor=(1.3, .8) #, # ,
            bbox_to_anchor=bbox_to_anchor

        )

    labels = [item.get_text() for item in ax.get_xticklabels()]

    labels = [get_marker_label(marker, ret=True) for marker in labels]

    # ['B cells\n(CD20)', '$\mathregular{CD4^{+}}$\nT cells',
    #               '$\mathregular{CD8^{+}}$\nT cells','Treg cells\n(FOXP3)',
    #               'Macrophages\n(CD68)', 'Dendritic\n(CD11c)',
    #               'NK cells\n(CD56)', 'Pan CK\ncells', 'Unidentified\ncells']

    ax.set_xticklabels(labels)
    ax.spines[['right', 'top']].set_visible(False)
    plt.show()

    # add annotation
    box_pairs = [

        (("CD8", "Pleural"), ("CD8", "Peritoneal")),
        (("CK", "Pleural"), ("CK", "Peritoneal")),
        (("CD11c", "Pleural"), ("CD11c", "Peritoneal")),
        (("CD56", "Pleural"), ("CD56", "Peritoneal")),
        (("other", "Pleural"), ("other", "Peritoneal")),
    ]

    annot = Annotator(ax, box_pairs, data=data, x=x,
                      y=y, hue=hue, order=markers)

    annot.new_plot(ax, pairs=box_pairs, plot='barplot',
                   data=data, x=x, y=y, hue=hue, hue_order=hue_order, seed=2021)
    # annot.configure(loc='outside')
    pvalues = ['%.1e' % x for x in pvalues]

    annot.set_custom_annotations(pvalues)
    annot.annotate()

    if save:
        plt.savefig(
            f"../figure/DensityPlot {title}_anno.pdf", format='pdf', bbox_inches="tight")

    plt.show()


def bar_plot(title, data, x, y, hue, xorder=None, save=True, nolegend=False):
    sns.set_theme(style='white')
    # colors = [ "#00E5EE","#FF3E96"] #[ "#00E5EE","#FF3E96","#F0F8FF"]  #['darkcyan', 'cornsilk']#['gold', 'seagreen']# ['gold', 'teal'] #
    # Draw the bar chart
    ax = sns.barplot(
        data=data,
        x=x, y=y,
        hue=hue,
        alpha=1,
        ci=None,
        order=xorder,
        # palette=colors,
        # edgecolor = "black",
        # hue_order=hue_order
    )
    ax.legend(
        loc=7,
        bbox_to_anchor=(1.45, .8),  # bbox_to_anchor=(1.3, .5),

    )
    ax.set_title(title,  fontsize=14)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    # ax.set_ylabel("log(density)")
    ax.set_ylabel("Contact score")
    # ax.set(ylim=(0, 0.5))
    # ax.set_ylim([0,10])
    if save:
        plt.savefig(f"../figure/BarPlot {title}.pdf",
                    format='pdf', bbox_inches="tight")


def heatmap_plot(data, outpath, title, vmin=0, vmax=0, mask=None, contain0=False, save=True):
    sns.set_style("white")

    platte = sns.light_palette("seagreen", as_cmap=True)

    if contain0:
        palette = sns.diverging_palette(220, 15, as_cmap=True)
    else:
        palette = sns.cubehelix_palette(
            start=2, rot=0, dark=0, light=.95,  as_cmap=True)  # reverse=True,
        # palette  = sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True)
        palette = palette.reversed()
        # palette = sns.color_palette("rocket", as_cmap=True)#sns.color_palette("viridis", as_cmap=True)
    sns.set(font_scale=1)
    if (vmin == 0 and vmax == 0):
        ax = sns.heatmap(data, cmap=palette, annot=True, fmt='.2f',
                         mask=mask, center=0, square=True, linewidths=.5, cbar_kws={"shrink": .5})  # Force cells to be square

    else:
        if contain0:
            ax = sns.heatmap(data, cmap=palette, annot=True, vmin=vmin, vmax=vmax,
                             fmt='.2f', mask=mask, center=0)  # annot_kws={"size": 10}
        else:
            ax = sns.heatmap(data, cmap=palette, annot=True, vmin=vmin,
                             vmax=vmax, fmt='.2f', mask=mask)  # annot_kws={"size": 10}

    plt.title(title, fontsize=14)

    labels = [get_marker_label(marker, ret=False) for marker in data.columns]
    # labels = ['B cells (CD20)', '$\mathregular{CD4^{+}}$ T cells',
    #               '$\mathregular{CD8^{+}}$ T cells','Treg cells (FOXP3)',
    #               'Macrophages (CD68)', 'Dendritic (CD11c)',
    #               'NK cells (CD56)', 'Pan CK cells']
    # plt.set_xticklabels(labels)
    # ax.set_xticklabels(labels=labels, rotation=90)
    # ax.set_yticklabels(labels=labels, rotation=0)
    plt.xlabel("")
    plt.ylabel("")
    # plt.yticks(rotation=0)

    if save:
        # plt.tight_layout()
        # plt.subplots_adjust(left=0.3, right=1, bottom=0.1, top=0.2)
        # plt.savefig(f"../../Paper/Figure/Correlation/heatmap {title}.pdf", format='pdf', bbox_inches = "tight", dpi=500)

        plt.savefig(f"{outpath}heatmap {title}.pdf",
                    format='pdf', bbox_inches="tight", dpi=500)


def bar_plot_series(data, title):
    data.plot.bar(color=['C0', 'C1', 'C2', 'C3'])
    plt.ylabel("Percentage")
    plt.title(title)
    plt.xticks(rotation=45)


def scatter_plot(data, x, y, grouptype, hue,  subs, title, corr, outpath, save=True):
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Arial'

    data = data.loc[data[hue].isin(subs), :]

    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)
    # sns.regplot(x=x+grouptype, y=y+grouptype, data=data);
    sns.lmplot(x=x+grouptype, y=y+grouptype, data=data,
               hue=hue, fit_reg=False, legend=False)

    plt.legend(bbox_to_anchor=(1.6, 0.96), loc='upper right',
               prop={'size': 16},

               # frameon=True,
               # labelspacing=1
               )
    plt.title(title, size=20)
    plt.xlabel(get_marker_label(x) + " ave " +
               grouptype + " of patient", size=16)
    plt.ylabel(y + " ave " + grouptype + " of patient", size=16)
    _, xright = plt.xlim()
    _, ytop = plt.ylim()
    if corr > 0:
        plt.text(xright*0.1, ytop*0.8, "corr = {:.2f}".format(corr), size=16)
    else:
        plt.text(xright*0.6, ytop*0.8, "corr = {:.2f}".format(corr), size=16)

    if save:
        # plt.savefig(f"../figure/correlation/{title}.png",dpi=100, bbox_inches = "tight")
        plt.savefig(f"{outpath}/{title}.pdf", format="pdf",
                    bbox_inches="tight", dpi=500)


def waterfall_plot(data1, data2, casetype, subtype=None, by="MVB"):

    fc_p1 = data1.groupby([by]).apply(lambda x: x['phenotype_combined'].value_counts(
        sort=False, normalize=True).rename_axis('Type').reset_index(name='Percent'))
    fc_p2 = data2.groupby([by]).apply(lambda x: x['phenotype_combined'].value_counts(
        sort=False, normalize=True).rename_axis('Type').reset_index(name='Percent'))

    df_fc1 = fc_p1.reset_index()
    df_fc2 = fc_p2.reset_index()

    df_percent1 = df_fc1.set_index([by, 'Type'])['Percent'].unstack(
    ).reset_index()  # pivot dataframe unstack used to reshape dataframe
    df_percent2 = df_fc2.set_index([by, 'Type'])['Percent'].unstack(
    ).reset_index()  # colnames ['MVB', 'CD11c', 'CD56', 'Unidentified']
    df_percent2.drop(['other'], axis=1, inplace=True)  # delete column other
    # combine df_percent and df_percent2, deduct CD56 and CD11c percents from other in df_percent

    df_percent = pd.merge(df_percent1, df_percent2, on=by)

    df_percent['Unidentified'] = df_percent['Unidentified'] - \
        df_percent['CD56'] - df_percent['CD11c']
    # remove MVB that Unidentified < 0
    df_percent = df_percent.loc[df_percent['Unidentified'] >= 0, :]

    df_plot = df_percent
    df_plot.index = df_plot[by]
    df_plot.drop([by], axis=1, inplace=True)
    df_plot.index.name = None
    df_plot.columns.name = None

    df_plot = df_plot.fillna(0)
    sortby = "CK"
    df_plot2 = df_plot.sort_values(by=[sortby])

    if casetype is None:
        title = "Cell type composition for all patients by " + by + " sortby " + sortby
    else:
        title = f"Cell type composition for {casetype} by " + \
            by + " sortby " + sortby
        if not subtype is None:
            title = f"{title}-{subtype}"

    df_plot2.to_csv(f"./output/figure/waterfall {title}.csv")

    ############################################################################
    # Stacked bar plot
    # barplot with matplotlib
    markers = ["Unidentified", "CD4", "CD8", "CD68",
               "FOXP3", "CD11c", "CD20", "CD56", "CK"]

    x = list(df_plot2.index)
    y1 = df_plot2['Unidentified']
    y2 = df_plot2['CD4']
    y3 = df_plot2['CD8']
    y4 = df_plot2['CD68']
    y5 = df_plot2['FOXP3']
    y6 = df_plot2['CD11c']
    y7 = df_plot2['CD20']
    y8 = df_plot2['CD56']
    y9 = df_plot2['CK']

    fig, ax = plt.subplots(figsize=(0.1*df_plot2.shape[0], 6))

    ax.set_title(title)
    ax.axis('off')
    ax.set_ylim((0, 1.05))

    # ax.legend(loc="upper right")
    # ax.xticks(rotation='vertical')
    ax.set_xticks([])
    ax.bar(x, y1, color='tab:gray', label='Unidentified\ncells')
    ax.bar(x, y2, bottom=y1, color='tab:olive',
           label='$\mathregular{CD4^{+}}$\nT cells')  # CD4
    ax.bar(x, y3, bottom=y1+y2, color='tab:purple',
           label='$\mathregular{CD8^{+}}$\nT cells',)  # CD8
    ax.bar(x, y4, bottom=y1+y2+y3, color='tab:green',
           label='$\mathregular{CD8^{+}}$\nT cells')  # CD68
    ax.bar(x, y5, bottom=y1+y2+y3+y4, color='tab:red',
           label='Treg cells\n(FOXP3)')  # FOXP3
    ax.bar(x, y6, bottom=y1+y2+y3+y4+y5, color='tab:blue',
           label='Dendritic\n(CD11c)')  # CD11c
    ax.bar(x, y7, bottom=y1+y2+y3+y4+y5+y6,
           color='tab:brown', label='B cells\n(CD20)')  # CD20
    ax.bar(x, y8, bottom=y1+y2+y3+y4+y5+y6+y7,
           color='tab:pink', label='NK cells\n(CD56)')  # CD56
    ax.bar(x, y9, bottom=y1+y2+y3+y4+y5+y6+y7+y8,
           color='tab:orange', label='Pan CK\ncells')  # CK

    legend_labels = ['Pan CK\ncells', 'NK cells\n(CD56)', 'B cells\n(CD20)', 'Dendritic\n(CD11c)',
                     'Treg cells\n(FOXP3)', 'Macrophages\n(CD68)', '$\mathregular{CD8^{+}}$\nT cells',
                     '$\mathregular{CD4^{+}}$\nT cells',  'Unidentified\ncells'
                     ]

    handles, labels = fig.gca().get_legend_handles_labels()
    # specify order
    order = [8, 7, 6, 5, 4, 3, 2,1,0]

    # pass handle & labels lists along with order as below
    plt.legend(handles=[handles[i] for i in order],
               labels=legend_labels,
               scatterpoints=1,
               loc="lower center",  # "upper center" puts it below the line
               ncol=9,
               prop={"size": 10},
               bbox_to_anchor=(0.5, 0.95),
               bbox_transform=fig.transFigure
               )

    # ax.legend(
    #     labels= legend_labels,
    #     scatterpoints=1,
    #     loc="lower center", # "upper center" puts it below the line
    #     ncol=9,
    #     prop = { "size": 10 },
    #     bbox_to_anchor=(0.5, 0.95),
    #     bbox_transform=fig.transFigure
    # )

    fig.savefig(
        f"../../Paper/Figure/Composition/waterfall {title}.pdf", bbox_inches='tight', format='pdf', dpi=500)
    # fig.savefig(f"./output/figure/waterfall {title}.eps", bbox_inches='tight', format='eps', dpi=1000)


def plot_with_anno(marker, df_comp, x, y, pairs, pvalues, ylim=[0, 1], yloc=None):
    ax = sns.boxplot(data=df_comp, x=x, y=y)
    ax.xaxis.set_tick_params(rotation=90)
    ax.set_ylabel("Composition")
    ax.set(xlabel=None)

    ax.set_title(marker,  fontsize=14, y=yloc)
    annot = Annotator(ax, pairs, data=df_comp, x=x, y=y)
    annot.new_plot(ax, pairs=pairs, plot='boxplot',
                   data=df_comp, x=x, y=y)
    # annot.configure(loc='outside')
    pvalues = ['%.1e' % x for x in pvalues]

    annot.set_custom_annotations(pvalues)
    annot.annotate()

    ax.set_ylim([0, 1])
    ax.spines[['right', 'top']].set_visible(False)
    plt.show()


def change_font(font="Arial", font_family="sans-serif"):

    font = {'size': 16}
    # matplotlib.rc('font', **font)
    #
    # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    #rc('font',**{'family':'serif','serif':['Palatino']})
    #rc('text', usetex=False

    # change font
    # matplotlib.rcParams['font.family'] = font_family
    # matplotlib.rcParams['font.sans-serif'] = font

# def generate_pvalues_oneSided(df_plot, markers, casetype, groups,grouptype, title, non_parametric=False):
#     pvalues = dict()
#     for marker in markers:
#         x = df_plot.loc[(df_plot[casetype]==groups[0]) & (df_plot['phenotype']==marker),grouptype].dropna().values
#         y = df_plot.loc[ (df_plot[casetype]==groups[1]) & (df_plot['phenotype']==marker),grouptype].dropna().values
#         if len(groups) == 3:
#             z = df_plot.loc[ (df_plot[casetype]==groups[2]) & (df_plot['phenotype']==marker),grouptype].dropna().values
#             pvalues[marker] = get_pvalue_oneSided(x,y,z, non_parametric=non_parametric)
#         elif len(groups) == 2:
#             pvalues[marker] = get_pvalue_oneSided(x,y, non_parametric=non_parametric)

#     df_pvalue = pd.DataFrame(pvalues,index=["greater", "less"])
#     if non_parametric:
#         df_pvalue.to_csv(f"../../Paper/Figure/BarPlot/{title} pvalue one-sided_non_parametric.csv",index=True)
#     else:
#         df_pvalue.to_csv(f"../../Paper/Figure/BarPlot/{title} pvalue one-sided.csv",index=True)
