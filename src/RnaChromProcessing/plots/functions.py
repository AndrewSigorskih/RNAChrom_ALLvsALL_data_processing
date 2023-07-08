import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np 
import pandas as pd
import seaborn as sns

def set_style_white() -> None:
    sns.set_style('white')
    sns.set_palette('husl')
    plt.rc('font',  
           serif = 'Ubuntu',
           monospace = 'Ubuntu Mono',
           size = 10)
    plt.rc('axes', 
           labelsize = 16,
           labelweight = 'bold',
           labelpad = 10,
           titlesize = 22,
           titlepad = 10,
           titleweight = 'bold')
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14
    plt.rcParams['legend.fontsize'] = 14


def rna_strand_barplot(wins: pd.DataFrame,
                       total_genes: int,
                       out_dir: str,
                       prefix: str) -> None:
    # init figure
    fig, ax = plt.subplots()
    fig.set_size_inches(11.7, 8.27)
    ax.set_ylim(-total_genes-2, total_genes+2)
    # variables
    labels = wins.index
    x = np.arange(len(labels))  # the label locations
    width = 0.5  # the width of the bars
    groups = pd.Series([x.split('_')[0] for x in labels])
    # groups coloring
    colors =[]
    patches=[]
    palette = sns.color_palette("husl", len(groups.unique()))
    for i, name in enumerate(groups.unique()):
        colors += [palette[i]] * groups[groups==name].shape[0]
        patches.append(
            mpatches.Patch(
                color=palette[i], 
                label=name)
            )
    # plot bars
    rects = []
    negrects = []
    for i in x:
        rects += ax.bar(i, wins.iat[i, 0], width, color=colors[i], alpha=.5)
        negrects += ax.bar(i, -wins.iat[i, 1], width, color=colors[i], alpha=.5)
    # labels and legend
    ax.set_ylabel('Numbers of wins', fontsize=20)
    ax.set_title(f'Numbers of wins and losses\nout of {total_genes} genes', fontsize=16)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=90)
    plt.legend(handles=patches, fontsize=14)
    
    ax.axes.get_yaxis().set_ticks([])
    
    # exact numbers over rects
    def autolabel(rects,neg=1):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in rects:
            height = round(rect.get_height(),2)
            ax.annotate(f'{np.abs(height)}',
                    xy=(rect.get_x()+width/2, height),
                    xytext=(0, 3*neg),  # 3 points vertical offset
                    textcoords='offset points',
                    ha='center', va='bottom')
    autolabel(rects)
    autolabel(negrects,-3)
    # save
    ax.axhline(color='grey')
    fig.tight_layout()
    plt.savefig(f'{out_dir}/{prefix}_wins.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{out_dir}/{prefix}_wins.svg', dpi=300, bbox_inches='tight', format='svg')
