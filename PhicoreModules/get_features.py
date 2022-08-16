import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from typing import List, Union

__author__ = 'Przemyslaw Decewicz'

def get_features_of_type(seqiorec: SeqRecord, ftype: str) -> List[SeqFeature]:
    """
    Get features of a given type from SeqRecord
    :param seqiorec: a SeqRecord object
    :param ftype: type of a feature
    :return:
    """

    flist = []
    for feature in seqiorec.features:
        if feature.type == ftype:
            flist.append(feature)
    
    return flist

def get_gc_content(seq: Union[str, Seq, SeqRecord]) -> float:
    """
    Calculate GC content of a nucleotide sequence
    :param seq: a nucleotide sequence
    :return:
    """

    gc = 0
    for i in seq:
        if i == 'G' or i == 'C':
            gc += 1
    return gc / len(seq)

def get_features_lengths(seqiorec: SeqRecord, ftype: str) -> List[float]:
    """
    Get average length of SeqFeatures of a given type
    :param seqiorec: a SeqRecord object
    :param ftype: type of a feature
    :return:
    """

    lengths = 0
    for feature in seqiorec.features:
        if feature.type == ftype:
            lengths += float(len(feature.location.extract(seqiorec).seq))

    if ftype == 'CDS':
        return [x / 3 for x in lengths]
    else:
        return lengths

def get_coding_density(seqiorec: SeqRecord, ftypes: List[str] = ['CDS', 'tRNA', 'rRNA']) -> float:
    """
    Get coding density for a SeqRecord considering given features types
    :param seqiorec: SeqRecord object
    :param ftypes: a list of feature types
    :return:
    """

    cdcov = np.zeros(len(seqiorec.seq))
    for feature in seqiorec.features:
        if feature.type in ftypes:
            start, stop = map(int, sorted([feature.location.start, feature.location.end]))
            cdcov[start:stop] += 1
    return sum([1 if x > 0 else 0 for x in cdcov]) / len(seqiorec.seq)

def get_distribution_of_stops(seqiorec: SeqRecord, window: int = 210, step: int = 1) -> pd.DataFrame:
    """
    Get distribution of STOP codons in a sequence
    :param seqiorec: SeqRecord object
    :param window: window size
    :param step: step size
    :return:
    """

    stops = ['TAA', 'TAG', 'TGA']

    stops_distr = {
        'x': range(1, len(seqiorec.seq) + 1),
        'TAA': [np.NAN]*int(window/2),
        'TAG': [np.NAN]*int(window/2),
        'TGA': [np.NAN]*int(window/2)
    }
    
    i = 0
    while i + window/2 + 1 <= len(seqiorec.seq) - window/2:
        window_seq = seqiorec.seq[i : i + window]
        taa = window_seq.count('TAA')
        tag = window_seq.count('TAG')
        tga = window_seq.count('TGA')
        stops_distr['TAA'].extend([taa]*(step))
        stops_distr['TAG'].extend([tag]*(step))
        stops_distr['TGA'].extend([tga]*(step))
        i += step
        
    i -= step
    left = len(seqiorec.seq) - len(stops_distr['TAA'])
    if left > 0:   
        stops_distr['TAA'].extend([np.NAN]*left)
        stops_distr['TAG'].extend([np.NAN]*left)
        stops_distr['TGA'].extend([np.NAN]*left)

    return pd.DataFrame(stops_distr)
