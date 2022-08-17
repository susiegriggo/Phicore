import pandas as pd
from typing import List, Tuple

__author__ = 'Przemyslaw Decewicz'

def colours() -> List[Tuple[str, str]]:
    """
    Get list of colours for Artemis graph file
    :return:
    """
    
    colours = [
        ('0 0 0','black'),
        ('230 25 75','red'),
        ('60 180 75','green'),
        ('0 130 200','blue'),
        ('255 225 25','yellow'),
        ('145 30 180','purple'),
        ('70 240 240','cyan'),
        ('245 130 48','orange'),
        ('240 50 230','pink'),
        ('210 245 60','lime'),
        ('250 190 212','peach'),
        ('0 128 128','teal'),
        ('220 190 255','lavender'),
        ('170 110 40','brown'),
        ('255 250 200','beige'),
        ('128 0 0','maroon'),
        ('170 255 195','mint'),
        ('128 128 0','olive'),
        ('255 215 180','coral'),
        ('0 0 128','navy'),
        ('128 128 128','gray'),
        ('255 255 255','white'),
    ]
    return colours

def write_df_to_artemis(df : pd.DataFrame,  filename : str,  colours : List[str] = colours()):
    """
    Write and clean dataframe to Artemis graph file
    :param df: Pandas dataframe with base position in first column and values in other columns
    :param filename: name of the file to write to
    :return:
    """
    df.to_csv(filename, sep=' ', index=False)
    with open(filename, 'r') as f:
        lines = f.readlines()

    with open(filename, 'w') as f:
        header = lines[0].split()[1:]
        f.write('# BASE {}\n'.format(' '.join(header)))

        colours_rgb = []
        colours_names = []
        for rgb in colours[:len(header)]:
            colours_rgb.append(':'.join(rgb[0].split()))
            colours_names.append(rgb[1])
        f.write('# colour {}\n'.format (' '.join(colours_rgb)))
        f.write('# label {}\n'.format(' '.join(header))) # the label must go after colour
        f.write('# name {}\n'.format (' '.join(colours_names)))
        for line in lines[1:]:
            if len(line.split()) == 1:
                continue
            f.write(line)
