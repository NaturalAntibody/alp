import math

import pytest

from alp import (
    AA_ALPHABET,
    ALPError,
    BLOSUM,
    GumbelParams,
    ROBINSON_FREQS,
    compute_gumbel_params,
    match_mismatch_matrix,
)

# Reference values from MMseqs2's EvalueComputation.h precomputed table.
MMSEQS_BLOSUM62_GAPLESS = (0.3207378152604042354, 0.13904657125294345166, 0.76221128839920349041)
MMSEQS_NUCLEOTIDE_5_2 = (0.62092274139392822363, 0.35177597988201619872, 0.74528059208662511548)


def blosum62_dict():
    m = BLOSUM["blosum62"]
    d = {r: {c: m[i][j] for j, c in enumerate(AA_ALPHABET)} for i, r in enumerate(AA_ALPHABET)}
    # extra ambiguity symbols, as in the `blosum` package, must be ignored
    for r in list(d):
        d[r]["X"] = -1
    d["X"] = {c: -1 for c in list(AA_ALPHABET) + ["X"]}
    return d


class TestMatrixInput:
    def test_gapless_blosum62_matches_mmseqs_with_its_freqs(self):
        p = compute_gumbel_params("blosum62", bg_freqs="mmseqs")
        assert not p.gapped
        lam, k, a = MMSEQS_BLOSUM62_GAPLESS
        assert p.lambda_ == pytest.approx(lam, rel=1e-9)
        assert p.K == pytest.approx(k, rel=1e-9)
        assert p.a_I == pytest.approx(a, rel=1e-9)

    def test_gapless_blosum62_matches_blast_by_default(self):
        # blast_stat.c blosum62 ungapped row: lambda .3176, K .134, alpha .7916, alpha_v 4.96466
        p = compute_gumbel_params("blosum62")
        assert p.lambda_ == pytest.approx(0.3176, abs=5e-5)
        assert p.K == pytest.approx(0.134, abs=5e-4)
        assert p.a_I == pytest.approx(0.7916, abs=5e-5)
        assert p.alpha_I == pytest.approx(4.96466, abs=5e-6)

    def test_gapped_blosum62_close_to_blast(self):
        # blast_stat.c BLOSUM62 11/1: lambda=0.267, K=0.041 (tabulated to 2-3
        # significant digits; ALP at default tolerances is within ~1% / ~10%).
        p = compute_gumbel_params("blosum62", gap_open=11, gap_extend=1)
        assert p.gapped
        assert p.lambda_ == pytest.approx(0.267, rel=0.01)
        assert p.K == pytest.approx(0.041, rel=0.12)
        assert set(p.errors) == set(p.to_dict())

    def test_all_input_forms_agree(self):
        named = compute_gumbel_params("BLOSUM62", gap_open=11, gap_extend=1)
        as_dict = compute_gumbel_params(blosum62_dict(), gap_open=11, gap_extend=1)
        as_array = compute_gumbel_params(
            [list(r) for r in BLOSUM["blosum62"]], alphabet=AA_ALPHABET,
            bg_freqs="robinson", gap_open=11, gap_extend=1,
        )
        assert named == as_dict == as_array

    def test_array_without_alphabet_defaults_to_uniform(self):
        uniform = compute_gumbel_params(BLOSUM["blosum62"])
        named_uniform = compute_gumbel_params("blosum62", bg_freqs="uniform")
        robinson = compute_gumbel_params("blosum62")
        assert uniform == named_uniform
        assert uniform.lambda_ != robinson.lambda_

    def test_numpy_input(self):
        np = pytest.importorskip("numpy")
        a = compute_gumbel_params(np.array(BLOSUM["blosum62"]), bg_freqs="uniform")
        assert a == compute_gumbel_params("blosum62", bg_freqs="uniform")

    def test_fixed_plan_is_default_and_reproducible(self):
        from alp import _core
        a = compute_gumbel_params("blosum62", gap_open=11, gap_extend=1, n_samples=3000)
        _core._compute.cache_clear()
        b = compute_gumbel_params("blosum62", gap_open=11, gap_extend=1, n_samples=3000)
        assert a.as_tuple() == b.as_tuple()
        assert a.samples_used == (3000, 3000)
        assert compute_gumbel_params("blosum62", gap_open=11, gap_extend=1).samples_used == (20000, 20000)

    def test_adaptive_mode(self):
        p = compute_gumbel_params("blosum62", gap_open=11, gap_extend=1, n_samples=None)
        assert p.lambda_ == pytest.approx(0.267, rel=0.01)
        assert p.samples_used[0] > 0

    def test_deterministic_for_seed(self):
        kw = dict(gap_open=10, gap_extend=1, bg_freqs="uniform")
        a = compute_gumbel_params("blosum80", seed=1, **kw)
        from alp import _core
        _core._compute.cache_clear()
        b = compute_gumbel_params("blosum80", seed=1, **kw)
        assert a.as_tuple() == b.as_tuple()

    @pytest.mark.parametrize("name", sorted(BLOSUM))
    def test_all_builtins_symmetric_and_gapless_ok(self, name):
        m = BLOSUM[name]
        assert all(m[i][j] == m[j][i] for i in range(20) for j in range(20))
        assert compute_gumbel_params(name).lambda_ > 0


class TestMatchMismatch:
    def test_reproduces_mmseqs_nucleotide_5_2(self):
        # MMseqs2 ran ALP in adaptive mode (tolerances 0.01/0.05, 60 s, seed 42)
        p = compute_gumbel_params(match=2, mismatch=3, gap_open=5, gap_extend=2, n_samples=None)
        lam, k, a = MMSEQS_NUCLEOTIDE_5_2
        assert p.lambda_ == pytest.approx(lam, rel=1e-9)
        assert p.K == pytest.approx(k, rel=1e-9)
        assert p.a_I == pytest.approx(a, rel=1e-9)

    def test_mismatch_sign_is_irrelevant(self):
        a = compute_gumbel_params(match=1, mismatch=-2, gap_open=5, gap_extend=2)
        b = compute_gumbel_params(match=1, mismatch=2, gap_open=5, gap_extend=2)
        assert a == b

    def test_gapless_analytic_lambda(self):
        # Gapless, uniform ACGT, +1/-1: lambda solves 1/4 e^l + 3/4 e^-l = 1  ->  e^l = 3
        p = compute_gumbel_params(match=1, mismatch=1)
        assert p.lambda_ == pytest.approx(math.log(3), rel=1e-6)

    def test_alphabet_size(self):
        p4 = compute_gumbel_params(match=1, mismatch=1)
        p20 = compute_gumbel_params(match=1, mismatch=1, alphabet=AA_ALPHABET, bg_freqs="uniform")
        assert p4.lambda_ != p20.lambda_

    def test_matrix_helper(self):
        assert match_mismatch_matrix(2, 3, 3) == ((2, -3, -3), (-3, 2, -3), (-3, -3, 2))


class TestGapConvention:
    def test_open_includes_extend_is_shifted_ncbi(self):
        a = compute_gumbel_params("blosum62", gap_open=12, gap_extend=1,
                                  gap_convention="open_includes_extend")
        b = compute_gumbel_params("blosum62", gap_open=11, gap_extend=1)
        assert a == b


class TestStatistics:
    def test_bit_score(self):
        p = compute_gumbel_params("blosum62", gap_open=11, gap_extend=1)
        s = 50
        assert p.bit_score(s) == pytest.approx((p.lambda_ * s - math.log(p.K)) / math.log(2))

    def test_evalue_decreases_with_score_and_scales_with_db(self):
        p = compute_gumbel_params("blosum62", gap_open=11, gap_extend=1)
        assert p.evalue(60, 200, 1e6) < p.evalue(40, 200, 1e6)
        e1, e2 = p.evalue(60, 200, 1e6), p.evalue(60, 200, 1e7)
        assert e2 / e1 == pytest.approx(10, rel=0.05)

    def test_evalue_close_to_naive_for_large_space(self):
        p = compute_gumbel_params("blosum62", gap_open=11, gap_extend=1)
        s, m, n = 80, 5000, 1e8
        naive = p.K * m * n * math.exp(-p.lambda_ * s)
        assert p.evalue(s, m, n) == pytest.approx(naive, rel=0.3)
        assert p.area(s, m, n) < m * n  # finite-size correction shrinks the area

    def test_pvalue(self):
        p = compute_gumbel_params("blosum62", gap_open=11, gap_extend=1)
        e = p.evalue(30, 100, 100)
        assert p.pvalue(30, 100, 100) == pytest.approx(1 - math.exp(-e))

    def test_roundtrip_from_stored_values(self):
        p = compute_gumbel_params("blosum62", gap_open=11, gap_extend=1)
        q = GumbelParams(*p.as_tuple())
        assert q.evalue(45, 300, 1e6) == p.evalue(45, 300, 1e6)


class TestValidation:
    def test_needs_exactly_one_scoring_system(self):
        with pytest.raises(ValueError):
            compute_gumbel_params()
        with pytest.raises(ValueError):
            compute_gumbel_params("blosum62", match=1, mismatch=1)
        with pytest.raises(ValueError):
            compute_gumbel_params(match=1)

    def test_gap_params_together(self):
        with pytest.raises(ValueError):
            compute_gumbel_params("blosum62", gap_open=11)

    def test_positive_expected_score_rejected(self):
        with pytest.raises(ValueError, match="expected score"):
            compute_gumbel_params(match=5, mismatch=1)

    def test_non_integer_scores_rejected(self):
        with pytest.raises(ValueError, match="integer"):
            compute_gumbel_params([[1.5, -1], [-1, 1]])

    def test_non_square_rejected(self):
        with pytest.raises(ValueError, match="square"):
            compute_gumbel_params([[1, -1, -1], [-1, 1, -1]])

    def test_unknown_name(self):
        with pytest.raises(ValueError, match="unknown matrix"):
            compute_gumbel_params("pam250")

    def test_freq_length_mismatch(self):
        with pytest.raises(ValueError, match="entries"):
            compute_gumbel_params(match=1, mismatch=1, bg_freqs=[0.5, 0.5])

    def test_robinson_freqs_sum(self):
        assert sum(ROBINSON_FREQS.values()) == pytest.approx(1.0, abs=1e-4)

    def test_time_limit_is_hard_timeout_in_fixed_mode(self):
        with pytest.raises(ALPError):
            compute_gumbel_params("blosum62", gap_open=11, gap_extend=1, n_samples=4000, max_seconds=1e-6)
        # a generous limit does not change the result
        a = compute_gumbel_params("blosum62", gap_open=11, gap_extend=1, n_samples=4000)
        b = compute_gumbel_params("blosum62", gap_open=11, gap_extend=1, n_samples=4000, max_seconds=600)
        assert a.as_tuple() == b.as_tuple()

    def test_alp_error_is_reported(self):
        # Absurdly low time budget for a gapped simulation should fail inside ALP.
        with pytest.raises(ALPError):
            compute_gumbel_params("blosum45", gap_open=1, gap_extend=1, max_seconds=1e-9)
