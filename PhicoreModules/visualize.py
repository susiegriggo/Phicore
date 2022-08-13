"""
Plotting


"""

__author__ = 'Vijini'


def visualise(df, length, labels, style="o"):

    """
    df: dataframe
    length: length of genomes
    lables: lables of stats
    """

    axes = df.plot.line(subplots=True, figsize=(16, 8), legend=False, style=style)

    for i in range(len(axes)):
        axes[i].set_xlim(0, length)
        axes[i].yaxis.grid()
        axes[i].set_axisbelow(True)
        axes[i].set(ylabel=labels[i])


