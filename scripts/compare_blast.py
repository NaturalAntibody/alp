"""3-way comparison: our ALP estimates vs BLAST vs MMseqs2 precomputed tables.

Usage:
    python scripts/compare_blast.py BLAST_STAT_C [--mmseqs EvalueComputation.h]
        [--csv out.csv] [--md out.md] [--n-samples 20000 | --n-samples 0 --lambda-tol 0.01 --k-tol 0.05]
        [--max-seconds S] [--jobs N]

Every row of BLAST's precomputed tables (``blosumNN_values`` and
``blastn_values_R_P`` in ``algo/blast/core/blast_stat.c``) is recomputed with
ALP using the same scores and background frequencies (BLAST's
``Robinson_prob`` for proteins, uniform ACGT for nucleotides).  Where MMseqs2
has a hard-coded entry (``defaultParameter`` in ``EvalueComputation.h``) its
values are added, and ALP is also run with MMseqs2's own background
frequencies so that each reference is compared on its own inputs.  MMseqs2
entries without a BLAST counterpart are appended at the end.

Column mapping (BLAST -> ALP):  lambda, K, alpha (col 7) -> a,
alpha_v (col 10) -> alpha, sigma (col 11) -> sigma.  BLAST's H and the old
Altschul-Gish beta have no ALP counterpart; its b/beta/tau are derived, not
tabulated.  Both use the NCBI gap convention (gap of length k costs
open + k*extend).

Special BLAST nucleotide rows:
  * ``{0, 0}`` rows are megablast's non-affine (greedy) statistics: linear gaps
    with extension cost reward/2 + penalty.  For odd rewards the scores are
    doubled so ALP gets integers, and lambda is rescaled.
  * Tables for even rewards (2/-7, 2/-5, 2/-3) and the megablast rows assume
    scores lie on a lattice of span 2 in ALP's units.  For a lattice of span
    d, K scales by d(1-e^-λ)/(1-e^-dλ), so ``alp_K_span2`` = K·2/(1+e^-λ) is the
    like-for-like value for those rows.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # run from a checkout
from alp import ALPError, BLOSUM62_FREQS, compute_gumbel_params  # noqa: E402

INF = 32767  # INT2_MAX
EVEN_SCORE_TABLES = {(2, 7), (2, 5), (2, 3)}  # "can only be applied to even scores"
# MMseqs2 matrix file -> scoring system name used here
MMSEQS_SYSTEMS = {"blosum62.out": "blosum62", "nucleotide.out": "+2/-3"}
PARAMS = ("lambda_", "K", "a", "alpha", "sigma")


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"//[^\n]*", "", text)


def _num(tok: str) -> float:
    return float(INF) if "INT2_MAX" in tok else float(tok)


def parse_blast(path: Path):
    text = _strip_comments(path.read_text())
    tables = {}
    for m in re.finditer(r"array_of_8\s+(blosum\d+|blastn_values_\d+_\d+)(?:_values)?\s*\[[^\]]*\]\s*=\s*\{(.*?)\};", text, re.S):
        rows = [[_num(t) for t in r.split(",") if t.strip()] for r in re.findall(r"\{([^{}]*)\}", m.group(2))]
        tables[m.group(1)] = rows
    rob = re.search(r"Robinson_prob\[\]\s*=\s*\{(.*?)\};", text, re.S).group(1)
    freqs = {c: float(p) for c, p in re.findall(r"\{\s*'(\w)'\s*,\s*([\d.]+)\s*\}", rob)}
    total = sum(freqs.values())
    return tables, {k: v / total for k, v in freqs.items()}


def parse_mmseqs(path: Path):
    """{(system, gap): {lambda_, K, a, alpha, sigma}} from defaultParameter[]."""
    text = _strip_comments(path.read_text())
    body = re.search(r"defaultParameter\[\]\s*=\s*\{(.*)\};", text, re.S).group(1)
    out = {}
    for m in re.finditer(r'\{\s*"([^"]+)"\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(true|false)\s*,\s*\{([^{}]*)\}\s*\}', body):
        name, go, ge, gapped, vals = m.groups()
        v = [float(x) for x in vals.split(",")]
        # AlignmentEvaluerParameters order: lambda, K, a1, b1, a2, b2,
        # alpha1, beta1, alpha2, beta2, sigma, tau (I/J symmetric here)
        gap = f"{go}/{ge}" if gapped == "true" else "ungapped"
        out[(MMSEQS_SYSTEMS.get(name, name), gap)] = dict(lambda_=v[0], K=v[1], a=v[2], alpha=v[6], sigma=v[10])
    return out


# ---------------------------------------------------------------------------
# ALP jobs
# ---------------------------------------------------------------------------

def _alp_job(job):
    key, kw, scale = job
    t0 = time.time()
    try:
        p = compute_gumbel_params(**kw)
    except (ALPError, ValueError) as e:
        return key, None, time.time() - t0, str(e).splitlines()[0][:80]
    # ALP ran on scores multiplied by `scale`; convert back to original units.
    res = dict(
        lambda_=p.lambda_ * scale, K=p.K, a=p.a_I * scale, alpha=p.alpha_I * scale,
        sigma=p.sigma * scale, lambda_err=p.errors["lambda_"] * scale, K_err=p.errors["K"],
        K_span2=p.K * 2 / (1 + math.exp(-p.lambda_)),  # lambda in ALP's units
    )
    return key, res, time.time() - t0, ""


def pct(a, b):
    if a is None or b is None or b == 0:
        return None
    return 100.0 * (a - b) / b


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("blast_stat_c", type=Path)
    ap.add_argument("--mmseqs", type=Path, help="MMseqs2 src/alignment/EvalueComputation.h")
    ap.add_argument("--csv", type=Path)
    ap.add_argument("--md", type=Path)
    ap.add_argument("--n-samples", type=int, default=20000,
                    help="fixed-plan sample size (reproducible); 0 = ALP adaptive mode")
    ap.add_argument("--lambda-tol", type=float, default=0.01, help="adaptive mode only")
    ap.add_argument("--k-tol", type=float, default=0.05, help="adaptive mode only")
    ap.add_argument("--max-seconds", type=float, default=None, help="default: no limit (fixed plan)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    args = ap.parse_args()

    tables, robinson = parse_blast(args.blast_stat_c)
    mmseqs = parse_mmseqs(args.mmseqs) if args.mmseqs else {}
    common = dict(lambda_tolerance=args.lambda_tol, k_tolerance=args.k_tol,
                  max_seconds=args.max_seconds, seed=args.seed,
                  n_samples=args.n_samples or None)

    rows, jobs = [], []

    def gaps_kw(go, ge):
        return {} if go >= INF else {"gap_open": int(go), "gap_extend": int(ge)}

    def add(system, gap, kind, blast, scoring, gkw, scale=1, note="", mm_freqs=None):
        key = len(rows)
        mm = mmseqs.pop((system, gap), None)
        rows.append(dict(system=system, gap=gap, kind=kind, note=note, blast=blast, mmseqs=mm))
        jobs.append(((key, "alp"), {**scoring, **gkw, **common}, scale))
        if mm is not None and mm_freqs is not None:
            jobs.append(((key, "alp_mm"), {**scoring, **gkw, **common, "bg_freqs": mm_freqs}, scale))

    for name in ("blosum45", "blosum50", "blosum62", "blosum80", "blosum90"):
        for r in tables[name]:
            blast = dict(zip(("lambda_", "K", "H", "a", "beta", "C", "alpha", "sigma"), r[3:11]))
            ungapped = r[0] >= INF
            add(name, "ungapped" if ungapped else f"{int(r[0])}/{int(r[1])}",
                "ungapped" if ungapped else "affine", blast,
                {"matrix": name, "bg_freqs": robinson}, gaps_kw(r[0], r[1]), mm_freqs=BLOSUM62_FREQS)

    for name, trows in tables.items():
        if not name.startswith("blastn_values_"):
            continue
        reward, penalty = map(int, name.split("_")[2:])
        system = f"+{reward}/-{penalty}"
        even = (reward, penalty) in EVEN_SCORE_TABLES
        for r in trows:
            blast = {"lambda_": r[2], "K": r[3]}
            if r[0] == 0 and r[1] == 0:
                # megablast non-affine: linear gap, extension = reward/2 + penalty
                s = 1 if reward % 2 == 0 else 2
                add(system, f"lin {reward / 2 + penalty:g}", "megablast", blast,
                    {"match": reward * s, "mismatch": penalty * s, "bg_freqs": "uniform"},
                    {"gap_open": 0, "gap_extend": (reward * s) // 2 + penalty * s}, scale=s,
                    note="megablast non-affine row; K in () is span-2 corrected")
            else:
                add(system, f"{int(r[0])}/{int(r[1])}", "affine", blast,
                    {"match": reward, "mismatch": penalty, "bg_freqs": "uniform"},
                    gaps_kw(r[0], r[1]), mm_freqs="uniform",
                    note="even-score table; K in () is span-2 corrected" if even else "")

    # MMseqs2 entries with no BLAST counterpart
    for (system, gap) in list(mmseqs):
        go, ge = (INF, INF) if gap == "ungapped" else map(int, gap.split("/"))
        if system == "+2/-3":
            scoring, mmf = {"match": 2, "mismatch": 3, "bg_freqs": "uniform"}, "uniform"
        elif system == "blosum62":
            scoring, mmf = {"matrix": "blosum62", "bg_freqs": robinson}, BLOSUM62_FREQS
        else:
            continue
        add(system, gap, "ungapped" if go >= INF else "affine", None, scoring, gaps_kw(go, ge),
            mm_freqs=mmf, note="MMseqs2 only (no BLAST row)")

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        for (key, which), res, dt, err in ex.map(_alp_job, jobs):
            rows[key][which] = res
            rows[key][f"{which}_seconds"] = round(dt, 2)
            if err:
                rows[key][f"{which}_error"] = err
    wall = time.time() - t0

    # ---- flatten -------------------------------------------------------------
    flat = []
    for r in rows:
        b, m, a, am = r["blast"] or {}, r["mmseqs"] or {}, r.get("alp") or {}, r.get("alp_mm") or {}
        span = "span-2" in r["note"]
        f = dict(system=r["system"], gaps=r["gap"], kind=r["kind"], note=r["note"])
        for src, d in (("blast", b), ("mmseqs", m), ("alp", a)):
            for k in PARAMS:
                f[f"{src}_{k.rstrip('_')}"] = d.get(k)
        f.update(alp_lambda_err=a.get("lambda_err"), alp_K_err=a.get("K_err"),
                 alp_K_span2=a.get("K_span2") if span else None)
        for k in PARAMS:
            f[f"alp_mmfreqs_{k.rstrip('_')}"] = am.get(k)
        f.update(
            d_lambda_alp_vs_blast=pct(a.get("lambda_"), b.get("lambda_")),
            d_K_alp_vs_blast=pct(a.get("K_span2") if span else a.get("K"), b.get("K")),
            d_a_alp_vs_blast=pct(a.get("a"), b.get("a")),
            d_alpha_alp_vs_blast=pct(a.get("alpha"), b.get("alpha")),
            d_sigma_alp_vs_blast=pct(a.get("sigma"), b.get("sigma")),
            d_lambda_alp_vs_mmseqs=pct(am.get("lambda_"), m.get("lambda_")),
            d_K_alp_vs_mmseqs=pct(am.get("K"), m.get("K")),
            d_lambda_mmseqs_vs_blast=pct(m.get("lambda_"), b.get("lambda_")),
            d_K_mmseqs_vs_blast=pct(m.get("K"), b.get("K")),
            alp_seconds=r.get("alp_seconds"),
            error=r.get("alp_error", "") or r.get("alp_mm_error", ""),
        )
        flat.append(f)

    # ---- report --------------------------------------------------------------
    def g(x, n=4):
        return "" if x is None else f"{x:.{n}g}"

    def p(x):
        return "" if x is None else f"{x:+.1f}%"

    cols = ["system", "gaps", "λ BLAST", "λ MMseqs2", "λ ALP", "Δλ ALP–BLAST",
            "K BLAST", "K MMseqs2", "K ALP", "ΔK ALP–BLAST",
            "a BLAST / ALP", "α BLAST / ALP", "σ BLAST / ALP",
            "ALP(MMseqs2 freqs) λ / K", "Δλ / ΔK ALP–MMseqs2", "note"]
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for f in flat:
        kcol = g(f["alp_K"], 3) + (f" ({g(f['alp_K_span2'], 3)})" if f["alp_K_span2"] is not None else "")

        def pair(k):
            return "" if f[f"blast_{k}"] is None else f"{g(f[f'blast_{k}'], 3)} / {g(f[f'alp_{k}'], 3)}"

        has_mm = f["mmseqs_lambda"] is not None
        lines.append("| " + " | ".join([
            f["system"], f["gaps"], g(f["blast_lambda"]), g(f["mmseqs_lambda"]), g(f["alp_lambda"]),
            p(f["d_lambda_alp_vs_blast"]), g(f["blast_K"], 3), g(f["mmseqs_K"], 3), kcol,
            p(f["d_K_alp_vs_blast"]), pair("a"), pair("alpha"), pair("sigma"),
            f"{g(f['alp_mmfreqs_lambda'])} / {g(f['alp_mmfreqs_K'], 3)}" if has_mm else "",
            f"{p(f['d_lambda_alp_vs_mmseqs'])} / {p(f['d_K_alp_vs_mmseqs'])}" if has_mm else "",
            f["note"] + (f" ERROR: {f['error']}" if f["error"] else ""),
        ]) + " |")

    def summary(sel, key):
        v = sorted(abs(f[key]) for f in flat if sel(f) and f[key] is not None)
        if not v:
            return None
        return f"median {v[len(v) // 2]:.1f}%, 90th pct {v[int(0.9 * (len(v) - 1))]:.1f}%, max {v[-1]:.1f}% (n={len(v)})"

    groups = {
        "protein ungapped": lambda f: f["system"].startswith("blosum") and f["kind"] == "ungapped",
        "protein gapped": lambda f: f["system"].startswith("blosum") and f["kind"] == "affine",
        "nucleotide affine": lambda f: f["system"].startswith("+") and f["kind"] == "affine",
        "nucleotide megablast": lambda f: f["kind"] == "megablast",
    }
    mode = (f"fixed plan n_samples={args.n_samples}" if args.n_samples
            else f"adaptive lambda_tol={args.lambda_tol}, k_tol={args.k_tol}")
    summ = [f"ALP settings: {mode}, seed={args.seed}; "
            f"{len(jobs)} ALP runs, {wall:.0f}s wall, {sum(1 for f in flat if f['error'])} errors\n"]
    for name, sel in groups.items():
        summ.append(f"**{name}**")
        for key, label in (("d_lambda_alp_vs_blast", "|Δλ|"), ("d_K_alp_vs_blast", "|ΔK|"),
                           ("d_a_alp_vs_blast", "|Δa|"), ("d_alpha_alp_vs_blast", "|Δα|"),
                           ("d_sigma_alp_vs_blast", "|Δσ|")):
            s = summary(sel, key)
            if s:
                summ.append(f"- {label} ALP vs BLAST: {s}")
        summ.append("")
    report = "\n".join(summ) + "\n" + "\n".join(lines) + "\n"
    print(report)

    if args.md:
        args.md.write_text("# ALP vs BLAST vs MMseqs2 precomputed Gumbel parameters\n\n" + report)
    if args.csv:
        with args.csv.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(flat[0]))
            w.writeheader()
            w.writerows(flat)


if __name__ == "__main__":
    main()
