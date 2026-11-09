# alp

Standalone Python wrapper over the NCBI **ALP** library (Sheetlin & Frith,
public domain) for computing the Gumbel / Karlin–Altschul parameters
(λ, K and the finite-size corrections a, b, α, β, σ, τ) of local alignment
scores. Use these parameters to turn raw Smith–Waterman scores into bit-scores
and E-values.

The directory is self-contained: it ships its own copy of ALP (`lib/alp`), a
small C API (`src/alp_api.{h,cpp}`), a CMake build and a ctypes package
with **no runtime Python dependencies**.

## Install

```bash
pip install .            # needs a C++17 compiler; cmake is pulled in automatically
# or, for development:
pip install cmake && python scripts/build_ext.py && pip install -e '.[test]'
pytest
```

The shared library is looked up next to the package, then in `./build/`;
`ALP_LIB=/path/to/libalp.so` overrides both.

## Usage

```python
from alp import compute_gumbel_params

# Built-in BLOSUM (45/50/62/80/90) with Robinson & Robinson background freqs
p = compute_gumbel_params("blosum62", gap_open=11, gap_extend=1)
p.lambda_, p.K                    # (0.2741, 0.0494)

# Any dict-of-dicts matrix, e.g. from the `blosum` package
import blosum
p = compute_gumbel_params(blosum.BLOSUM(62), gap_open=11, gap_extend=1)

# Square integer array (list of lists / numpy) + optional alphabet
p = compute_gumbel_params(my_matrix, alphabet="ACGT", bg_freqs=[.3, .2, .2, .3],
                          gap_open=5, gap_extend=2)

# Match / mismatch scoring (default alphabet ACGT, uniform frequencies)
p = compute_gumbel_params(match=2, mismatch=3, gap_open=5, gap_extend=2)

# Omit gap penalties for (exact, analytic) gapless statistics
p = compute_gumbel_params("blosum62")

# Statistics
p.bit_score(55)
p.evalue(55, query_length=120, subject_length=1_000_000)  # ALP finite-size corrected
p.pvalue(55, 120, 1_000_000)
p.to_dict(), p.as_tuple(), p.errors                        # errors = ALP std. errors
```

`GumbelParams(*stored_tuple)` rebuilds a parameter set from saved values, so
the (sub-second to ~1 min) simulation only needs to run once per scoring
system. Within a process, results are also memoised.

### Options

| argument | default | meaning |
|---|---|---|
| `alphabet` | auto | symbol order; for dict matrices the 20 standard amino acids are used when present (B/Z/X/\* are dropped) |
| `bg_freqs`, `bg_freqs2` | Robinson for 20 AA, else uniform | `"uniform"`, `"robinson"`/`"blast"` (BLAST's Robinson & Robinson), `"blosum62"`/`"mmseqs"` (BLOSUM62 background, as in MMseqs2), a mapping, or a sequence; `bg_freqs2` for sequence #2 |
| `gap_convention` | `"ncbi"` | `"ncbi"`: gap of length k costs `open + k·extend` (ALP, BLAST, MMseqs2). `"open_includes_extend"`: `open + (k−1)·extend` (SSW, parasail); converted before calling ALP |
| `n_samples` | 20000 | fixed number of simulation samples, so results are bit-identical for given inputs and seed. `None` = ALP's adaptive mode (see below) |
| `n_preliminary_samples` | 200 | samples per preliminary stage (fixed mode) |
| `lambda_tolerance`, `k_tolerance` | 0.01, 0.05 | target relative errors (adaptive mode only) |
| `max_seconds` | `None` | no time limit by default. If set, it is a hard timeout (`ALPError`) in fixed and gapless mode; in adaptive mode it is ALP's planning budget (default 60 s) |
| `max_mem_mb` | 500 | memory limit |
| `seed` | 42 | random seed |

`ValueError` is raised for invalid input, including scoring systems with no
valid local statistics (no positive score or non-negative expected score).
`ALPError` is raised when ALP itself fails.

## Notes

* **Reproducibility.** By default ALP chooses its sample sizes adaptively,
  from the tolerances *and* from measured run time. With a binding time
  limit, results then depend on machine speed and load: BLOSUM62 9/1 on a
  loaded machine gave K and α values 3.5% and 10% different. This package
  therefore defaults to a fixed plan (`n_samples=20000`). Results are
  bit-identical across runs, load and `max_seconds`. They depend only on
  the inputs, the seed and the build, because floating-point differences
  between compilers or math libraries could still change them.
* **Accuracy vs `n_samples`.** Statistical error is roughly ∝ 1/√n and
  runtime is linear in n. At 20000, λ is good to ≤0.5% and K to ≤2%. The
  finite-size α/σ converge more slowly and are biased low at small n for
  low gap costs. Raise `n_samples` (e.g. 100000) when E-values for short
  sequences matter. Runtime at 20000 is ~2.5 s for BLOSUM62 11/1, ~10 s for
  9/1 and ~36 s for BLOSUM45 16/1.
* `n_samples=None` (seed 42, tolerances 0.01/0.05, 60 s) reproduces
  MMseqs2's ALP settings.
* `lib/alp` is vendored NCBI ALP with one patch (`sls_alignment_evaluer.cpp`,
  marked "alp-python patch"). In fixed-plan mode, the gapless pre-step
  originally had a hard-coded 120 s cap; it now honours `max_seconds`, and has
  no limit when that is unset.
* Gapless parameters are exact. They reproduce BLAST's ungapped rows with
  Robinson frequencies and MMseqs2's gapless BLOSUM62 entry with
  `bg_freqs="mmseqs"`.
* C API: see `src/alp_api.h` (`gumbel_compute`, `gumbel_evalue`).

## Comparison with BLAST and MMseqs2 tables

`scripts/compare_blast.py` recomputes every precomputed row in BLAST's
`blast_stat.c` (BLOSUM45–90, all blastn reward/penalty tables). It also adds
MMseqs2's hard-coded entries from `EvalueComputation.h` where they exist:

```bash
python scripts/compare_blast.py path/to/blast_stat.c \
    --mmseqs path/to/EvalueComputation.h --csv cmp.csv --md cmp.md \
    [--n-samples 100000]   # or --n-samples 0 for ALP adaptive mode
```
