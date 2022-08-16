import pandas as pd
from typing import List

__author__ = 'Przemyslaw Decewicz'

def write_df_to_artemis(df : pd.DataFrame, filename : str, colors : List[str]):
    """
    Write and clean dataframe to Artemis graph file
    :param df: Pandas dataframe
    :param filename:
    :return:
    """
    df.to_csv(filename, sep=' ', index=False)
    with open(filename, 'r') as f:
        lines = f.readlines()
        # BASE TAA TAG TGA 
        # colour 255:0:0 0:255:0 0:0:250

    colors = ['255 0 0', '0 255 0', '0 0 255', '255 255 0', '0 255 255', '255 0 255', '255 255 255']
    with open(filename, 'w') as f:
        header = lines[0].split()[1:]
        f.write('# BASE {}\n'.format(' '.join(header)))
        f.write('# COLOUR {}\n'.format (' '.join(colors[:len(header)])))                                                                      
        for line in lines[1:]:
            if len(line.split()) == 1:
                continue
            f.write(line)
