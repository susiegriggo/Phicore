"""
the modules
"""

from .read_genbank import parse_genbank
from .get_features import get_features_of_type, get_gc_content, get_features_lengths, get_coding_density, get_distribution_of_stops, get_mean_cds_length_rec_window, get_rolling_gc, get_rolling_mean_cds
from .stats import mean, median, mode, stdev
from .visualize import visualise
from .write import write_df_to_artemis


__all__ = [
    'parse_genbank', 'mean', 'median', 'mode', 'stdev',
    'get_features_of_type', 'get_gc_content', 'get_features_lengths', 'get_coding_density', 'get_distribution_of_stops', 'get_mean_cds_length_rec_window', 'get_rolling_gc', 'get_rolling_mean_cds'
    'visualise',
    'write_df_to_artemis'
]

