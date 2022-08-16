"""
the modules
"""

from .read_genbank import parse_genbank
from .get_features import get_features_of_type, get_gc_content, get_features_lengths, get_coding_density, get_distribution_of_stops
from .stats import mean, median, mode, stdev
from .visualize import visualise
from .write import write_df_to_artemis




__all__ = [
    'parse_genbank', 'mean', 'median', 'mode', 'stdev',
    'get_features_of_type', 'get_gc_content', 'get_features_lengths', 'get_coding_density', 'get_distribution_of_stops',
    'visualise',
    'write_df_to_artemis'
]

