"""
output.py — Write analysis results to TSV, BED, summary report, and log file.
"""

import os
import datetime


# ──────────────────────────────────────────────────────────────
# 1. TSV (full statistics table)
# ──────────────────────────────────────────────────────────────

def write_tsv(results, outdir, name, has_rep):
    """
    Write all significant peaks with statistics to a TSV file.
    Columns: chr  start  end  treat_mean  ctrl_mean  log2FC  [pval  qval]
             direction  dominant_signal  credibility_p
    """
    path = os.path.join(outdir, "{}_diff_peaks.tsv".format(name))
    headers = ["chr", "start", "end", "treat_mean", "ctrl_mean", "log2FC"]
    if has_rep:
        headers += ["pval", "qval"]
    headers += ["direction", "dominant_signal", "credibility_p"]

    with open(path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for r in results:
            is_increase  = r["log2fc"] > 0
            direction    = "increase" if is_increase else "decrease"
            dom_signal   = r["treat_mean"] if is_increase else r["ctrl_mean"]
            credibility_p = r.get("credibility_p", 1.0)
            row = [
                r["chr"], str(r["start"]), str(r["end"]),
                "{:.6f}".format(r["treat_mean"]),
                "{:.6f}".format(r["ctrl_mean"]),
                "{:.6f}".format(r["log2fc"]),
            ]
            if has_rep:
                row += ["{:.4e}".format(r["pval"]), "{:.4e}".format(r["qval"])]
            row += [direction, "{:.6f}".format(dom_signal), "{:.6f}".format(credibility_p)]
            f.write("\t".join(row) + "\n")
    return path


# ──────────────────────────────────────────────────────────────
# 2. BED (IGV-compatible, with log2FC in score/extra column)
# ──────────────────────────────────────────────────────────────

def write_bed(results, outdir, name):
    """
    Write significant peaks as an 8-column BED file:
    chr  start  end  name  igv_score  strand(+/-)  log2FC  credibility_p

    score column (col 5):
      igv_score = int(percentile_rank * 1000), 0-1000
      Derived from the signal strength of the dominant group relative to
      all peaks. Higher score = stronger signal = darker color in IGV.

    strand column (col 6):
      + = increase (treat > ctrl)  → IGV colors blue
      - = decrease (ctrl > treat)  → IGV colors red

    credibility_p (col 8):
      1 - percentile_rank of the dominant signal across all peaks.
      Smaller value = stronger signal = more likely a real peak.
      Treat like a p-value: credibility_p < 0.05 means top 5%% signal.
    """
    path = os.path.join(outdir, "{}_diff_peaks.bed".format(name))
    with open(path, "w") as f:
        f.write('track name="{}" description="diffpeak differential peaks" colorByStrand="0,0,255 255,0,0"\n'.format(name))
        for i, r in enumerate(results, 1):
            is_increase   = r["log2fc"] > 0
            direction     = "increase" if is_increase else "decrease"
            strand        = "+" if is_increase else "-"
            peak_name     = "{}_{}_{}".format(name, i, direction)
            igv_score     = r.get("igv_score", 0)
            credibility_p = r.get("credibility_p", 1.0)
            f.write(
                "{}\t{}\t{}\t{}\t{}\t{}\t{:.4f}\t{:.6f}\n".format(
                    r["chr"], r["start"], r["end"],
                    peak_name, igv_score, strand,
                    r["log2fc"], credibility_p
                )
            )
    return path


# ──────────────────────────────────────────────────────────────
# 3. Summary report (plain text)
# ──────────────────────────────────────────────────────────────

def write_summary(
    args,
    pseudocount,
    total_ref_peaks,
    total_evaluated,
    results_all,
    results_filtered,
    has_rep,
    outdir,
    name,
    elapsed_sec,
):
    """
    Write a human-readable summary report (.txt).
    """
    path = os.path.join(outdir, f"{name}_summary.txt")
    now  = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    up   = sum(1 for r in results_filtered if r["log2fc"] > 0)
    down = sum(1 for r in results_filtered if r["log2fc"] < 0)

    # log2FC distribution of filtered peaks
    lfc_vals = [r["log2fc"] for r in results_filtered]
    mean_lfc = sum(lfc_vals) / len(lfc_vals) if lfc_vals else 0.0
    max_lfc  = max(lfc_vals) if lfc_vals else 0.0
    min_lfc  = min(lfc_vals) if lfc_vals else 0.0

    with open(path, "w") as f:
        f.write("=" * 60 + "\n")
        f.write("  diffpeak — Differential Peak Analysis Summary\n")
        f.write("=" * 60 + "\n")
        f.write(f"  Run date   : {now}\n")
        f.write(f"  Run time   : {elapsed_sec:.1f} s\n")
        f.write(f"  Name       : {name}\n")
        f.write(f"  Output dir : {outdir}\n")
        f.write("\n")

        f.write("── Input files ──────────────────────────────────────────\n")
        f.write(f"  Treat bgnorm    ({len(args.treat)} file(s)):\n")
        for fn in args.treat:
            f.write(f"    {fn}\n")
        f.write(f"  Control bgnorm  ({len(args.ctrl)} file(s)):\n")
        for fn in args.ctrl:
            f.write(f"    {fn}\n")
        f.write("  Peak files      ({} file(s)):\n".format(len(args.peaks)))
        for fn in args.peaks:
            f.write("    {}\n".format(fn))
        f.write("\n")

        f.write("── Parameters ───────────────────────────────────────────\n")
        f.write("  Mode           : {}\n".format("replicate (t-test + BH)" if has_rep else "single sample (log2FC only)"))
        f.write("  Min |log2FC|   : {}\n".format(args.log2fc))
        if has_rep:
            f.write("  Max p-value    : {}\n".format(args.pval if args.pval is not None else "not filtered"))
            f.write("  Max q-value    : {}\n".format(args.qval if args.qval is not None else "not filtered"))
        f.write("  Direction      : {}\n".format(args.direction))
        f.write("  Peak merge gap : {} bp\n".format(args.merge_gap))
        f.write("  Min signal     : {}\n".format(args.min_signal))
        f.write("  Pseudocount    : {:.2e} (min non-zero signal / 2)\n".format(pseudocount))
        f.write("  Min peak len   : {}\n".format("{} bp".format(args.min_len) if args.min_len is not None else "not filtered"))
        f.write("  Max peak len   : {}\n".format("{} bp".format(args.max_len) if args.max_len is not None else "not filtered"))
        f.write("  Max credibility_p : {}\n".format(args.cred if args.cred is not None else "not filtered"))
        f.write(f"  Threads        : {args.threads}\n")
        f.write("\n")

        f.write("── Results ──────────────────────────────────────────────\n")
        f.write(f"  Merged ref peaks  : {total_ref_peaks}\n")
        f.write(f"  Peaks evaluated   : {total_evaluated}\n")
        f.write(f"  Peaks passing filter : {len(results_filtered)}\n")
        f.write(f"    ↑ Increase      : {up}\n")
        f.write(f"    ↓ Decrease      : {down}\n")
        f.write("\n")
        if lfc_vals:
            f.write("── log2FC distribution (filtered peaks) ─────────────────\n")
            f.write(f"  Mean log2FC : {mean_lfc:.4f}\n")
            f.write(f"  Max  log2FC : {max_lfc:.4f}\n")
            f.write(f"  Min  log2FC : {min_lfc:.4f}\n")
            f.write("\n")

        f.write("── Output files ─────────────────────────────────────────\n")
        f.write(f"  {name}_diff_peaks.tsv\n")
        f.write(f"  {name}_diff_peaks.bed\n")
        f.write(f"  {name}_summary.txt\n")
        f.write(f"  {name}.log\n")
        f.write("=" * 60 + "\n")

    return path


# ──────────────────────────────────────────────────────────────
# 4. Logger
# ──────────────────────────────────────────────────────────────

class Logger:
    """
    Simple dual logger — writes to stdout and a .log file simultaneously.
    """

    def __init__(self, outdir, name):
        os.makedirs(outdir, exist_ok=True)
        self.log_path = os.path.join(outdir, f"{name}.log")
        self._fh = open(self.log_path, "w")
        self._start = datetime.datetime.now()
        self._write_header(name)

    def _write_header(self, name):
        ts = self._start.strftime("%Y-%m-%d %H:%M:%S")
        self.info(f"diffpeak started at {ts}")
        self.info(f"Log file: {self.log_path}")

    def _ts(self):
        return datetime.datetime.now().strftime("%H:%M:%S")

    def info(self, msg):
        line = f"[{self._ts()}] {msg}"
        print(line)
        self._fh.write(line + "\n")
        self._fh.flush()

    def warn(self, msg):
        line = f"[{self._ts()}] WARNING: {msg}"
        print(line)
        self._fh.write(line + "\n")
        self._fh.flush()

    def error(self, msg):
        line = f"[{self._ts()}] ERROR: {msg}"
        print(line)
        self._fh.write(line + "\n")
        self._fh.flush()

    def elapsed(self):
        delta = datetime.datetime.now() - self._start
        return delta.total_seconds()

    def close(self):
        self._fh.close()
