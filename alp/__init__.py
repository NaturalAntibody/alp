"""Gumbel (Karlin-Altschul) parameters for local alignment via NCBI ALP."""

from ._core import (
    ALPError,
    GumbelParams,
    BLOSUM62_FREQS,
    ROBINSON_FREQS,
    compute_gumbel_params,
    match_mismatch_matrix,
)
from ._matrices import AA_ALPHABET, BLOSUM

__all__ = [
    "AA_ALPHABET",
    "ALPError",
    "BLOSUM",
    "GumbelParams",
    "BLOSUM62_FREQS",
    "ROBINSON_FREQS",
    "compute_gumbel_params",
    "match_mismatch_matrix",
]
__version__ = "0.1.0"
