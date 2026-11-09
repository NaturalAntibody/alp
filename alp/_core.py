"""ctypes wrapper over ``libalp.so`` (NCBI ALP library)."""

from __future__ import annotations

import ctypes
import functools
import math
import os
import pathlib
from dataclasses import dataclass, field
from typing import Dict, Mapping, Optional, Sequence, Tuple, Union

from ._matrices import AA_ALPHABET, BLOSUM

# ---------------------------------------------------------------------------
# Library loading
# ---------------------------------------------------------------------------

_LIB_NAME = "libalp.so"


def _load_library() -> ctypes.CDLL:
    here = pathlib.Path(__file__).resolve().parent
    candidates = [here / _LIB_NAME, here.parent / "build" / _LIB_NAME]
    env = os.environ.get("ALP_LIB")
    if env:
        candidates.insert(0, pathlib.Path(env))
    for p in candidates:
        if p.exists():
            return ctypes.CDLL(str(p))
    raise OSError(
        f"Cannot find {_LIB_NAME}. Build it with `pip install .` (or "
        "`cmake -B build && cmake --build build`) or set ALP_LIB."
    )


_lib = _load_library()

_PARAM_NAMES = (
    "lambda_", "K", "a_I", "b_I", "a_J", "b_J",
    "alpha_I", "beta_I", "alpha_J", "beta_J", "sigma", "tau",
)


class _CParams(ctypes.Structure):
    _fields_ = [
        (f, ctypes.c_double)
        for name in _PARAM_NAMES
        for f in (name, name.rstrip("_") + "_error")
    ]


_ERRBUF_LEN = 1024

_lib.gumbel_compute.argtypes = [
    ctypes.POINTER(ctypes.c_long), ctypes.c_int,
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
    ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.c_long,
    ctypes.c_int,
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_long,                                   # seed
    ctypes.c_long, ctypes.c_long, ctypes.POINTER(ctypes.c_long),  # fixed-plan samples
    ctypes.POINTER(_CParams), ctypes.c_char_p, ctypes.c_int,
]
_lib.gumbel_compute.restype = ctypes.c_int

_lib.gumbel_evalue.argtypes = [
    ctypes.POINTER(_CParams), ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double), ctypes.c_char_p, ctypes.c_int,
]
_lib.gumbel_evalue.restype = ctypes.c_int


class ALPError(RuntimeError):
    """Raised when ALP fails (e.g. invalid scoring system or no convergence)."""


# ---------------------------------------------------------------------------
# Background frequencies
# ---------------------------------------------------------------------------

def _normalised(d: Dict[str, float]) -> Dict[str, float]:
    total = sum(d.values())
    return {k: v / total for k, v in d.items()}


#: Robinson & Robinson (1991) amino-acid frequencies, exactly as BLAST uses
#: them for its statistics (``Robinson_prob`` in ``blast_stat.c``).
ROBINSON_FREQS: Dict[str, float] = _normalised({
    "A": 78.05, "R": 51.29, "N": 44.87, "D": 53.64,
    "C": 19.25, "Q": 42.64, "E": 62.95, "G": 73.77,
    "H": 21.99, "I": 51.42, "L": 90.19, "K": 57.44,
    "M": 22.43, "F": 38.56, "P": 52.03, "S": 71.20,
    "T": 58.41, "W": 13.30, "Y": 32.16, "V": 64.41,
})

#: BLOSUM62 background frequencies (Henikoff & Henikoff), as listed in
#: MMseqs2's ``blosum62.out`` and used by MMseqs2 for its E-values.
BLOSUM62_FREQS: Dict[str, float] = {
    "A": 0.07422, "R": 0.05161, "N": 0.04465, "D": 0.05363,
    "C": 0.02469, "Q": 0.03426, "E": 0.05431, "G": 0.07415,
    "H": 0.02621, "I": 0.06792, "L": 0.09891, "K": 0.05815,
    "M": 0.02499, "F": 0.04742, "P": 0.03854, "S": 0.05723,
    "T": 0.05089, "W": 0.01303, "Y": 0.03228, "V": 0.07292,
}

_FREQ_PRESETS = {
    "robinson": ROBINSON_FREQS, "blast": ROBINSON_FREQS,
    "blosum62": BLOSUM62_FREQS, "mmseqs": BLOSUM62_FREQS,
}

#: Default fixed-plan sample size (see ``compute_gumbel_params``).
DEFAULT_N_SAMPLES = 20000

FreqSpec = Union[None, str, Mapping[str, float], Sequence[float]]
MatrixSpec = Union[str, Mapping[str, Mapping[str, int]], Sequence[Sequence[int]]]


# ---------------------------------------------------------------------------
# Result object
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class GumbelParams:
    """Gumbel parameters of the local alignment score distribution.

    ``I`` refers to sequence #1 (query), ``J`` to sequence #2 (subject).
    ``errors`` holds ALP's estimated standard error of each parameter.
    """

    lambda_: float
    K: float
    a_I: float
    b_I: float
    a_J: float
    b_J: float
    alpha_I: float
    beta_I: float
    alpha_J: float
    beta_J: float
    sigma: float
    tau: float
    gapped: bool = True
    errors: Dict[str, float] = field(default_factory=dict, compare=False, repr=False)
    #: Realizations ALP used for its (ALP, killing) stages; (0, 0) if gapless
    #: or not computed here.
    samples_used: Tuple[int, int] = field(default=(0, 0), compare=False, repr=False)

    # -- conversions -------------------------------------------------------

    def to_dict(self) -> Dict[str, float]:
        """The 12 parameters as a plain dict."""
        return {name: getattr(self, name) for name in _PARAM_NAMES}

    def as_tuple(self) -> Tuple[float, ...]:
        """The 12 parameters in ALP order (lambda, K, a_I, b_I, a_J, ...)."""
        return tuple(getattr(self, name) for name in _PARAM_NAMES)

    # -- statistics --------------------------------------------------------

    def bit_score(self, raw_score: float) -> float:
        """``(lambda * S - ln K) / ln 2``."""
        return (self.lambda_ * raw_score - math.log(self.K)) / math.log(2.0)

    def evalue(self, raw_score: float, query_length: float, subject_length: float) -> float:
        """E-value with ALP's finite-size (area) correction.

        For a database search pass the total database length (residues) as
        ``subject_length``.
        """
        return self._calc(raw_score, query_length, subject_length)[0]

    def pvalue(self, raw_score: float, query_length: float, subject_length: float) -> float:
        """P-value ``1 - exp(-E)``."""
        return self._calc(raw_score, query_length, subject_length)[1]

    def area(self, raw_score: float, query_length: float, subject_length: float) -> float:
        """Finite-size corrected search space (~ ``query_length * subject_length``)."""
        return self._calc(raw_score, query_length, subject_length)[2]

    def _calc(self, score: float, len1: float, len2: float) -> Tuple[float, float, float]:
        c = _CParams(**{n: getattr(self, n) for n in _PARAM_NAMES})
        e, p, a = ctypes.c_double(), ctypes.c_double(), ctypes.c_double()
        err = ctypes.create_string_buffer(_ERRBUF_LEN)
        rc = _lib.gumbel_evalue(
            ctypes.byref(c), float(score), float(len1), float(len2),
            ctypes.byref(e), ctypes.byref(p), ctypes.byref(a), err, _ERRBUF_LEN,
        )
        if rc != 0:
            raise ALPError(err.value.decode(errors="replace") or f"ALP error code {rc}")
        return e.value, p.value, a.value


# ---------------------------------------------------------------------------
# Input normalisation
# ---------------------------------------------------------------------------

def _as_int(value, what: str) -> int:
    try:
        f = float(value)
    except (TypeError, ValueError):
        raise TypeError(f"{what} must be a number, got {value!r}") from None
    if not math.isfinite(f) or f != int(f):
        raise ValueError(f"{what} must be an integer (ALP works on integer scores), got {value!r}")
    return int(f)


def _resolve_matrix(
    matrix: MatrixSpec, alphabet: Optional[Sequence[str]]
) -> Tuple[Tuple[Tuple[int, ...], ...], Optional[Tuple[str, ...]]]:
    """Return (square int matrix, alphabet or None)."""
    if isinstance(matrix, str):
        key = matrix.lower().replace("-", "").replace("_", "")
        if key not in BLOSUM:
            raise ValueError(f"unknown matrix {matrix!r}; built-ins: {sorted(BLOSUM)}")
        if alphabet is not None and tuple(alphabet) != tuple(AA_ALPHABET):
            raise ValueError("custom alphabet is not supported with a named built-in matrix")
        return BLOSUM[key], tuple(AA_ALPHABET)

    if isinstance(matrix, Mapping):
        if alphabet is None:
            keys = list(matrix.keys())
            # e.g. blosum.BLOSUM(62) also contains B, Z, X, J, U, O, *.  Those
            # are ambiguity/stop codes that would distort the statistics, so
            # restrict to the standard amino acids when they are all present.
            alphabet = AA_ALPHABET if set(AA_ALPHABET) <= set(keys) else keys
        alphabet = tuple(alphabet)
        try:
            rows = tuple(
                tuple(_as_int(matrix[r][c], f"matrix[{r!r}][{c!r}]") for c in alphabet)
                for r in alphabet
            )
        except KeyError as e:
            raise ValueError(f"matrix has no entry for symbol {e.args[0]!r}") from None
        return rows, alphabet

    # 2-D array-like (list of lists, numpy array, ...)
    if hasattr(matrix, "tolist"):
        matrix = matrix.tolist()
    rows = tuple(
        tuple(_as_int(v, f"matrix[{i}][{j}]") for j, v in enumerate(row))
        for i, row in enumerate(matrix)
    )
    n = len(rows)
    if n == 0 or any(len(r) != n for r in rows):
        raise ValueError("matrix must be a non-empty square 2-D array")
    if alphabet is not None:
        alphabet = tuple(alphabet)
        if len(alphabet) != n:
            raise ValueError(f"alphabet has {len(alphabet)} symbols but matrix is {n}x{n}")
    return rows, alphabet


def _resolve_freqs(
    spec: FreqSpec, n: int, alphabet: Optional[Tuple[str, ...]]
) -> Tuple[float, ...]:
    if spec is None:
        spec = "robinson" if alphabet is not None and set(alphabet) == set(AA_ALPHABET) else "uniform"
    if isinstance(spec, str):
        name = spec.lower()
        if name == "uniform":
            return tuple([1.0 / n] * n)
        if name in _FREQ_PRESETS:
            spec = _FREQ_PRESETS[name]
            if alphabet is None or not set(alphabet) <= set(spec):
                raise ValueError(f"{name!r} frequencies need an amino-acid alphabet")
        else:
            raise ValueError(f"unknown background frequency preset {spec!r}")
    if isinstance(spec, Mapping):
        if alphabet is None:
            raise ValueError("bg_freqs given as a mapping requires a known alphabet")
        missing = [a for a in alphabet if a not in spec]
        if missing:
            raise ValueError(f"bg_freqs missing symbols: {missing}")
        values = [float(spec[a]) for a in alphabet]
    else:
        if hasattr(spec, "tolist"):
            spec = spec.tolist()
        values = [float(v) for v in spec]
        if len(values) != n:
            raise ValueError(f"bg_freqs has {len(values)} entries, expected {n}")
    if any(v < 0 or not math.isfinite(v) for v in values) or sum(values) <= 0:
        raise ValueError("bg_freqs must be non-negative and sum to a positive value")
    total = sum(values)
    return tuple(v / total for v in values)


def _check_scoring_system(mat, f1, f2) -> None:
    """Karlin-Altschul requirements: some positive score, negative expectation."""
    n = len(mat)
    if not any(mat[i][j] > 0 and f1[i] > 0 and f2[j] > 0 for i in range(n) for j in range(n)):
        raise ValueError("scoring system has no positive score with non-zero probability")
    expected = sum(f1[i] * f2[j] * mat[i][j] for i in range(n) for j in range(n))
    if expected >= 0:
        raise ValueError(
            f"expected score per aligned pair is {expected:.4g} >= 0; local alignment "
            "statistics require a negative expected score"
        )


# ---------------------------------------------------------------------------
# Core computation (cached)
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=128)
def _compute(
    mat: Tuple[Tuple[int, ...], ...],
    f1: Tuple[float, ...],
    f2: Tuple[float, ...],
    gapped: bool,
    gap_open: int,
    gap_extend: int,
    insertions_after_deletions: bool,
    eps_lambda: float,
    eps_k: float,
    max_seconds: float,
    max_mem_mb: float,
    seed: int,
    n_samples: int,
    n_prelim: int,
) -> GumbelParams:
    n = len(mat)
    flat = (ctypes.c_long * (n * n))(*(v for row in mat for v in row))
    c_f1 = (ctypes.c_double * n)(*f1)
    c_f2 = (ctypes.c_double * n)(*f2)
    out = _CParams()
    used = (ctypes.c_long * 2)()
    err = ctypes.create_string_buffer(_ERRBUF_LEN)
    rc = _lib.gumbel_compute(
        flat, n, c_f1, c_f2, int(gapped),
        gap_open, gap_extend, gap_open, gap_extend,
        int(insertions_after_deletions),
        eps_lambda, eps_k, max_seconds, max_mem_mb, seed,
        n_samples, n_prelim, used,
        ctypes.byref(out), err, _ERRBUF_LEN,
    )
    if rc != 0:
        msg = err.value.decode(errors="replace").strip() or f"error code {rc}"
        raise ALPError(f"ALP failed: {msg}")
    values = {name: getattr(out, name) for name in _PARAM_NAMES}
    errors = {name: getattr(out, name.rstrip("_") + "_error") for name in _PARAM_NAMES}
    return GumbelParams(**values, gapped=gapped, errors=errors,
                        samples_used=(int(used[0]), int(used[1])))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def match_mismatch_matrix(match: int, mismatch: int, size: int = 4) -> Tuple[Tuple[int, ...], ...]:
    """Build a ``size x size`` matrix with ``match`` on the diagonal and
    ``-abs(mismatch)`` elsewhere."""
    match = _as_int(match, "match")
    mismatch = -abs(_as_int(mismatch, "mismatch"))
    if match <= 0:
        raise ValueError("match score must be positive")
    if size < 2:
        raise ValueError("alphabet size must be >= 2")
    return tuple(tuple(match if i == j else mismatch for j in range(size)) for i in range(size))


def compute_gumbel_params(
    matrix: Optional[MatrixSpec] = None,
    *,
    match: Optional[int] = None,
    mismatch: Optional[int] = None,
    gap_open: Optional[int] = None,
    gap_extend: Optional[int] = None,
    alphabet: Optional[Sequence[str]] = None,
    bg_freqs: FreqSpec = None,
    bg_freqs2: FreqSpec = None,
    gap_convention: str = "ncbi",
    insertions_after_deletions: bool = False,
    lambda_tolerance: float = 0.01,
    k_tolerance: float = 0.05,
    max_seconds: Optional[float] = None,
    max_mem_mb: float = 500.0,
    seed: int = 42,
    n_samples: Optional[int] = DEFAULT_N_SAMPLES,
    n_preliminary_samples: int = 200,
) -> GumbelParams:
    """Estimate Gumbel parameters (lambda, K, finite-size corrections) with ALP.

    Scoring system — give exactly one of:

    * ``matrix``: a built-in name (``"blosum45"``, ``"blosum50"``,
      ``"blosum62"``, ``"blosum80"``, ``"blosum90"``), a dict-of-dicts
      ``{row: {col: score}}`` (e.g. ``blosum.BLOSUM(62)``), or a square 2-D
      array-like of integer scores.
    * ``match`` / ``mismatch``: a uniform scoring scheme over ``alphabet``
      (default ``"ACGT"``).  ``mismatch`` may be given as a cost (``3``) or a
      score (``-3``); it is always applied as ``-abs(mismatch)``.

    Parameters
    ----------
    gap_open, gap_extend : int, optional
        Positive gap penalties.  When both are omitted the (analytic, fast)
        gapless parameters are computed; otherwise the gapped parameters are
        estimated by importance-sampling simulation (seconds to a minute).
    alphabet : sequence of str, optional
        Symbol order.  Needed when ``bg_freqs`` is a mapping and ``matrix`` is
        an array.  For a dict-of-dicts matrix it selects which symbols take
        part; by default the 20 standard amino acids if all are present
        (ambiguity codes such as B/Z/X/* are dropped), otherwise every key.
    bg_freqs : None, preset name, mapping or sequence
        Letter background frequencies (normalised to sum to 1).  Presets:
        ``"uniform"``, ``"robinson"``/``"blast"`` (Robinson & Robinson, as in
        BLAST) and ``"blosum62"``/``"mmseqs"`` (BLOSUM62 background, as in
        MMseqs2).  Default: Robinson for the 20-amino-acid alphabet,
        otherwise uniform.
    bg_freqs2 : same as ``bg_freqs``, optional
        Frequencies for sequence #2 if they differ from sequence #1.
    gap_convention : {"ncbi", "open_includes_extend"}
        ``"ncbi"`` (ALP/BLAST/MMseqs2): a gap of length k costs
        ``gap_open + k*gap_extend``.  ``"open_includes_extend"`` (SSW,
        parasail, EMBOSS-style): it costs ``gap_open + (k-1)*gap_extend``;
        the values are converted before calling ALP.
    insertions_after_deletions : bool
        Allow an insertion to directly follow a deletion.
    lambda_tolerance, k_tolerance : float
        Target relative errors of lambda and K (adaptive mode only).
    max_seconds : float or None, default None
        Time limit.  ``None`` (default): no limit in fixed-plan and gapless
        mode, ALP's 60 s budget in adaptive mode.  In fixed-plan and gapless
        mode a limit is a hard timeout that raises ``ALPError`` and never
        changes the result.  In adaptive mode it is the budget ALP plans its
        sampling around.
    max_mem_mb : float
        Memory limit for the simulation.
    seed : int
        Random seed.
    n_samples : int or None, default 20000
        Number of importance-sampling realizations for the gapped simulation
        (fixed-plan mode).  The result is a pure function of the inputs and
        ``seed``: bit-identical across runs, machine speed and load.
        ``lambda_tolerance``/``k_tolerance`` are ignored and ``max_seconds``
        only acts as a hard timeout (``ALPError``).  Larger values reduce the
        statistical error roughly as 1/sqrt(n) and cost linear time.
        ``None`` selects ALP's adaptive mode: sample sizes are chosen from the
        tolerances *and* from measured run time, so results can vary between
        machines or under load whenever ``max_seconds`` comes into play.
    n_preliminary_samples : int
        Realizations per preliminary stage in fixed-plan mode.

    Returns
    -------
    GumbelParams

    Raises
    ------
    ValueError
        Invalid inputs, or a scoring system that has no valid local
        alignment statistics (no positive score / non-negative expectation).
    ALPError
        ALP itself failed (e.g. did not reach the requested accuracy in
        ``max_seconds``).

    Examples
    --------
    >>> p = compute_gumbel_params("blosum62", gap_open=11, gap_extend=1)
    >>> round(p.lambda_, 3), round(p.K, 3)
    (0.274, 0.045)
    >>> p = compute_gumbel_params(match=1, mismatch=2, gap_open=5, gap_extend=2)
    """
    if (matrix is None) == (match is None and mismatch is None):
        raise ValueError("give either `matrix` or `match`/`mismatch`, not both / neither")

    if matrix is None:
        if match is None or mismatch is None:
            raise ValueError("both `match` and `mismatch` are required")
        alphabet = tuple(alphabet) if alphabet is not None else tuple("ACGT")
        mat = match_mismatch_matrix(match, mismatch, len(alphabet))
    else:
        mat, alphabet = _resolve_matrix(matrix, alphabet)

    n = len(mat)
    f1 = _resolve_freqs(bg_freqs, n, alphabet)
    f2 = f1 if bg_freqs2 is None else _resolve_freqs(bg_freqs2, n, alphabet)
    _check_scoring_system(mat, f1, f2)

    if (gap_open is None) != (gap_extend is None):
        raise ValueError("give both `gap_open` and `gap_extend` (or neither for gapless)")
    gapped = gap_open is not None
    go = ge = 0
    if gapped:
        go = abs(_as_int(gap_open, "gap_open"))
        ge = abs(_as_int(gap_extend, "gap_extend"))
        if gap_convention == "open_includes_extend":
            go -= ge
            if go < 0:
                raise ValueError("gap_open must be >= gap_extend with 'open_includes_extend'")
        elif gap_convention != "ncbi":
            raise ValueError(f"unknown gap_convention {gap_convention!r}")
        if ge <= 0:
            raise ValueError("gap_extend must be positive for gapped statistics")
    if n_samples is not None and (int(n_samples) <= 0 or int(n_preliminary_samples) <= 0):
        raise ValueError("n_samples and n_preliminary_samples must be positive")

    return _compute(
        mat, f1, f2, gapped, go, ge, bool(insertions_after_deletions),
        float(lambda_tolerance), float(k_tolerance),
        float(max_seconds) if max_seconds is not None else 0.0, float(max_mem_mb), int(seed),
        int(n_samples or 0), int(n_preliminary_samples),
    )
