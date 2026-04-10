"""
analysis.py — Core analysis pipeline for diffpeak.
Requires Python >= 3.8 (fixes multiprocessing 2GB pickle limit in 3.7).
"""

import os
import sys
import math
from multiprocessing import Pool

from .cli      import parse_args
from .io       import parse_peak_bed, parse_bgnorm, validate_files
from .peaks    import merge_peaks, compute_coverage_for_peak
from .stats    import compute_pseudocount, log2fc, log2_values, welch_t_test, bh_correction, compute_credibility
from .output   import write_tsv, write_bed, write_summary, Logger


def run():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    log = Logger(args.outdir, args.name)
    log.info("=" * 52)
    log.info("  diffpeak — Differential Peak Analyzer")
    log.info("=" * 52)

    try:
        validate_files([args.treat, args.ctrl, args.peaks])
    except FileNotFoundError as e:
        log.error(str(e))
        sys.exit(1)

    has_rep  = len(args.treat) > 1 or len(args.ctrl) > 1
    use_pval = has_rep and args.pval is not None
    use_qval = has_rep and args.qval is not None

    log.info("Mode      : {}".format("replicate (t-test + BH)" if has_rep else "single sample (log2FC only)"))
    log.info("Treat     : {} bgnorm file(s)".format(len(args.treat)))
    log.info("Control   : {} bgnorm file(s)".format(len(args.ctrl)))
    log.info("Threads   : {}".format(args.threads))
    pval_str = ", pval <= {}".format(args.pval) if use_pval else ""
    qval_str = ", qval <= {}".format(args.qval) if use_qval else ""
    log.info("Filters   : log2FC >= {}{}{}{}".format(
        args.log2fc, pval_str, qval_str,
        ", direction = {}".format(args.direction)))

    # ── Step 1: Load and merge all peak files ─────────────────
    log.info("Step 1/4: Loading and merging peak files...")
    all_peaks = []
    for fpath in args.peaks:
        pk = parse_peak_bed(fpath)
        all_peaks.extend(pk)
        log.info("  {}: {} peaks".format(os.path.basename(fpath), len(pk)))

    if not all_peaks:
        log.error("No peaks found in any peak file. Exiting.")
        sys.exit(1)

    ref_peaks = merge_peaks(all_peaks, gap=args.merge_gap)
    log.info("  Merged reference: {} peaks (gap={} bp)".format(
        len(ref_peaks), args.merge_gap))

    # ── Step 2: Load bgnorm signal files ──────────────────────
    log.info("Step 2/4: Loading normalized signal files...")
    treat_signals, ctrl_signals = [], []

    for fpath in args.treat:
        sig = parse_bgnorm(fpath)
        treat_signals.append(sig)
        n_iv = sum(len(v) for v in sig.values())
        log.info("  treat {}: {} intervals".format(os.path.basename(fpath), n_iv))

    for fpath in args.ctrl:
        sig = parse_bgnorm(fpath)
        ctrl_signals.append(sig)
        n_iv = sum(len(v) for v in sig.values())
        log.info("  ctrl  {}: {} intervals".format(os.path.basename(fpath), n_iv))

    pseudocount = compute_pseudocount(treat_signals + ctrl_signals)
    log.info("  Pseudocount (min non-zero signal / 2): {:.2e}".format(pseudocount))

    # ── Step 3: Parallel coverage computation ─────────────────
    log.info("Step 3/4: Computing signal coverage ({} peaks, {} threads)...".format(
        len(ref_peaks), args.threads))

    tasks = [
        (chr_, start, end, treat_signals, ctrl_signals, args.min_signal)
        for chr_, start, end in ref_peaks
    ]

    raw_results = []
    if args.threads > 1:
        chunk = max(1, len(tasks) // (args.threads * 4))
        with Pool(processes=args.threads) as pool:
            for i, res in enumerate(pool.imap(
                    compute_coverage_for_peak, tasks, chunksize=chunk)):
                if res is not None:
                    raw_results.append(res)
                if (i + 1) % 10000 == 0:
                    log.info("  {}/{} peaks processed...".format(i + 1, len(tasks)))
    else:
        for i, task in enumerate(tasks):
            res = compute_coverage_for_peak(task)
            if res is not None:
                raw_results.append(res)
            if (i + 1) % 10000 == 0:
                log.info("  {}/{} peaks processed...".format(i + 1, len(tasks)))

    n_skipped = len(ref_peaks) - len(raw_results)
    log.info("  Evaluated: {} peaks (skipped low-signal: {})".format(
        len(raw_results), n_skipped))

    # ── Compute log2FC ─────────────────────────────────────────
    for r in raw_results:
        r["log2fc"] = log2fc(r["treat_mean"], r["ctrl_mean"], pseudocount)
        r["pval"]   = 1.0
        r["qval"]   = 1.0

    # ── Step 4: Statistical testing ───────────────────────────
    log.info("Step 4/4: Statistical testing and filtering...")

    if has_rep:
        p_vals = [
            welch_t_test(
                log2_values(r["treat_vals"], pseudocount),
                log2_values(r["ctrl_vals"],  pseudocount)
            )[1]
            for r in raw_results
        ]
        q_vals = bh_correction(p_vals)
        for i, r in enumerate(raw_results):
            r["pval"] = p_vals[i]
            r["qval"] = q_vals[i]
        log.info("  Welch t-test + BH correction applied.")
        if not use_pval and not use_qval:
            log.info("  Note: p/q-value computed but not used as filters.")
        elif use_pval and not use_qval:
            log.info("  Filtering by p-value only (pval <= {}).".format(args.pval))
        elif use_qval and not use_pval:
            log.info("  Filtering by q-value only (qval <= {}).".format(args.qval))
        else:
            log.info("  Filtering by p-value <= {} and q-value <= {}.".format(
                args.pval, args.qval))
    else:
        log.info("  Single sample: p/q-value not computed (log2FC filter only).")

    # ── Compute credibility on full raw_results pool ───────────
    compute_credibility(raw_results)
    log.info("  Credibility scores computed from {} peaks (full distribution).".format(
        len(raw_results)))

    # ── Filter ────────────────────────────────────────────────
    def passes_filter(r):
        if abs(r["log2fc"]) < args.log2fc:
            return False
        if args.direction == "up"   and r["log2fc"] <= 0:
            return False
        if args.direction == "down" and r["log2fc"] >= 0:
            return False
        if use_pval and r["pval"] > args.pval:
            return False
        if use_qval and r["qval"] > args.qval:
            return False
        peak_len = r["end"] - r["start"]
        if args.min_len is not None and peak_len < args.min_len:
            return False
        if args.max_len is not None and peak_len > args.max_len:
            return False
        if args.cred is not None and r.get("credibility_p", 1.0) > args.cred:
            return False
        return True

    results_filtered = [r for r in raw_results if passes_filter(r)]
    up   = sum(1 for r in results_filtered if r["log2fc"] > 0)
    down = sum(1 for r in results_filtered if r["log2fc"] < 0)
    len_info = ""
    if args.min_len is not None or args.max_len is not None:
        lo = "{}bp".format(args.min_len) if args.min_len is not None else "any"
        hi = "{}bp".format(args.max_len) if args.max_len is not None else "any"
        len_info = "  (length filter: {} ~ {})".format(lo, hi)
    log.info("  Peaks passing filter: {}  (up {}, down {}){}".format(
        len(results_filtered), up, down, len_info))

    # ── Write outputs ─────────────────────────────────────────
    tsv_path = write_tsv(results_filtered, args.outdir, args.name, has_rep)
    bed_path = write_bed(results_filtered, args.outdir, args.name)
    elapsed  = log.elapsed()

    sum_path = write_summary(
        args=args,
        pseudocount=pseudocount,
        total_ref_peaks=len(ref_peaks),
        total_evaluated=len(raw_results),
        results_all=raw_results,
        results_filtered=results_filtered,
        has_rep=has_rep,
        outdir=args.outdir,
        name=args.name,
        elapsed_sec=elapsed,
    )

    log.info("Output TSV     : {}".format(tsv_path))
    log.info("Output BED     : {}".format(bed_path))
    log.info("Summary report : {}".format(sum_path))
    log.info("Log file       : {}".format(log.log_path))
    log.info("Finished in {:.1f} s".format(elapsed))
    log.close()
