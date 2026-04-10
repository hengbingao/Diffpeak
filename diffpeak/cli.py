import argparse
from . import __version__


def get_parser():
    p = argparse.ArgumentParser(
        prog="diffpeak",
        description=(
            "diffpeak v{} -- Differential peak analysis for CUT&Tag / ChIP-seq\n"
            "Computes log2FC (and p/q-value with replicates) from normalized\n"
            "bedgraph signal files against a merged peak reference."
        ).format(__version__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single sample: log2FC filter only
  diffpeak -t treat.bgnorm -c ctrl.bgnorm
           --peaks treat.bed ctrl.bed
           --log2fc 1 -o results/ -n H3K27ac_BEP396

  # Multiple peak files (treat + ctrl mixed freely)
  diffpeak -t t1.bgnorm t2.bgnorm
           -c c1.bgnorm c2.bgnorm
           --peaks t1.bed t2.bed c1.bed c2.bed
           --log2fc 1 --pval 0.05
           --threads 8 -o results/ -n H3K27ac_BEP396

  # Broad histone mark: merge nearby peaks before analysis
  diffpeak -t treat.bgnorm -c ctrl.bgnorm
           --peaks treat.bed ctrl.bed
           --merge-gap 500
           --log2fc 1 -o results/ -n H3K27me3_BEP396
        """,
    )

    # ── Required: signal files ────────────────────────────────
    sig = p.add_argument_group("signal files (required)")
    sig.add_argument(
        "-t", "--treat", nargs="+", required=True, metavar="FILE",
        help="Treat normalized bedgraph / bgnorm file(s). Multiple files = replicates.",
    )
    sig.add_argument(
        "-c", "--ctrl", nargs="+", required=True, metavar="FILE",
        help="Control normalized bedgraph / bgnorm file(s). Multiple files = replicates.",
    )

    # ── Required: peak files ──────────────────────────────────
    pk = p.add_argument_group("peak files (required)")
    pk.add_argument(
        "--peaks", nargs="+", required=True, metavar="FILE",
        help=(
            "Peak BED file(s) — any number, treat and ctrl mixed freely. "
            "All files are pooled and merged into a single non-overlapping "
            "reference peak set used for signal quantification."
        ),
    )

    # ── Peak merging ──────────────────────────────────────────
    merge = p.add_argument_group("peak merging  (applies to --peaks only)")
    merge.add_argument(
        "--merge-gap", type=int, default=0, metavar="INT",
        help=(
            "Merge peaks separated by <= INT bp into one region "
            "(equivalent to bedtools merge -d). "
            "Default: 0 (only merge overlapping peaks). "
            "Recommended: 0-100 for sharp marks (H3K4me3), "
            "200-500 for broad marks (H3K27me3, H3K9me3)."
        ),
    )

    # ── Output ────────────────────────────────────────────────
    out = p.add_argument_group("output")
    out.add_argument(
        "-o", "--outdir", default=".", metavar="DIR",
        help="Output directory (default: current directory).",
    )
    out.add_argument(
        "-n", "--name", default="diffpeak", metavar="NAME",
        help="Output file prefix (default: diffpeak).",
    )

    # ── Filter thresholds ─────────────────────────────────────
    filt = p.add_argument_group("filter thresholds")
    filt.add_argument(
        "--log2fc", type=float, default=1.0, metavar="FLOAT",
        help="Minimum |log2FC| to report a peak (default: 1.0).",
    )
    filt.add_argument(
        "--pval", type=float, default=None, metavar="FLOAT",
        help=(
            "Maximum p-value (e.g. 0.05). Only available with replicates. "
            "Independent of --qval — either or both can be specified. "
            "If omitted, p-value is not used as a filter."
        ),
    )
    filt.add_argument(
        "--qval", type=float, default=None, metavar="FLOAT",
        help=(
            "Maximum q-value / BH-FDR (e.g. 0.1). Only available with replicates. "
            "Independent of --pval — either or both can be specified. "
            "If omitted, q-value is not used as a filter."
        ),
    )
    filt.add_argument(
        "--direction", choices=["both", "up", "down"], default="both",
        help="Direction to report: both / up (treat>ctrl) / down (ctrl>treat). Default: both.",
    )
    filt.add_argument(
        "--min-signal", type=float, default=0.0, metavar="FLOAT",
        help="Minimum mean signal in either group to keep a peak (default: 0).",
    )
    filt.add_argument(
        "--cred", type=float, default=None, metavar="FLOAT",
        help=(
            "Maximum credibility_p to report a peak (e.g. 0.2). "
            "credibility_p reflects the signal strength of the dominant group "
            "(treat for increase, ctrl for decrease) relative to all peaks. "
            "Smaller value = stronger signal = more likely a real peak. "
            "Similar to a p-value: --cred 0.05 keeps only the top 5%% signal peaks. "
            "If omitted, no credibility filter is applied."
        ),
    )
    filt.add_argument(
        "--min-len", type=int, default=None, metavar="INT",
        help=(
            "Minimum peak length (bp) to report (e.g. 200). "
            "Peaks shorter than this are discarded. If omitted, no lower length filter."
        ),
    )
    filt.add_argument(
        "--max-len", type=int, default=None, metavar="INT",
        help=(
            "Maximum peak length (bp) to report (e.g. 10000). "
            "Peaks longer than this are discarded. If omitted, no upper length filter."
        ),
    )

    # ── Miscellaneous ─────────────────────────────────────────
    misc = p.add_argument_group("miscellaneous")
    misc.add_argument(
        "--threads", type=int, default=4, metavar="INT",
        help="Parallel threads for signal coverage computation (default: 4).",
    )
    misc.add_argument(
        "--version", action="version", version="%(prog)s {}".format(__version__),
    )

    return p


def parse_args():
    return get_parser().parse_args()
