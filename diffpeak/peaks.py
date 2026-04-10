"""
peaks.py — Peak merging (bedtools merge equivalent) and
           signal coverage computation from bedgraph intervals.
"""

from collections import defaultdict


# ──────────────────────────────────────────────────────────────
# 1. Peak merging
# ──────────────────────────────────────────────────────────────

def merge_peaks(all_peaks, gap=0):
    """
    Merge a list of (chr, start, end) tuples into a non-overlapping
    reference set, combining peaks separated by <= gap bp.

    Equivalent to: bedtools merge -i <peaks> -d <gap>

    Args:
        all_peaks : list of (chr, start, end)
        gap       : int, maximum allowed gap between peaks to merge

    Returns:
        list of (chr, start, end), sorted by chr then start
    """
    by_chr = defaultdict(list)
    for chr_, start, end in all_peaks:
        by_chr[chr_].append((start, end))

    merged = []
    for chr_ in sorted(by_chr.keys()):
        sorted_ivs = sorted(by_chr[chr_], key=lambda x: x[0])
        cur_start, cur_end = sorted_ivs[0]
        for start, end in sorted_ivs[1:]:
            if start <= cur_end + gap:
                cur_end = max(cur_end, end)
            else:
                merged.append((chr_, cur_start, cur_end))
                cur_start, cur_end = start, end
        merged.append((chr_, cur_start, cur_end))

    return merged


# ──────────────────────────────────────────────────────────────
# 2. Coverage computation
# ──────────────────────────────────────────────────────────────

def _binary_search_first(ivs, peak_start):
    """
    Binary search: find the index of the first interval whose
    end > peak_start (i.e. the first one that could overlap the peak).
    """
    lo, hi, first = 0, len(ivs) - 1, len(ivs)
    while lo <= hi:
        mid = (lo + hi) // 2
        if ivs[mid][1] > peak_start:
            first = mid
            hi = mid - 1
        else:
            lo = mid + 1
    return first


def weighted_mean_in_peak(intervals_by_chr, chr_, peak_start, peak_end):
    """
    Compute the length-weighted mean normalized signal value
    within [peak_start, peak_end) on chromosome chr_.

    Uses binary search to locate overlapping intervals efficiently.

    Returns:
        float — weighted mean signal (0.0 if no overlap found)
    """
    if chr_ not in intervals_by_chr:
        return 0.0

    ivs = intervals_by_chr[chr_]
    first = _binary_search_first(ivs, peak_start)

    total_bp = 0
    weighted_sum = 0.0

    for i in range(first, len(ivs)):
        iv_start, iv_end, iv_val = ivs[i]
        if iv_start >= peak_end:
            break
        ov_start = max(iv_start, peak_start)
        ov_end   = min(iv_end,   peak_end)
        if ov_end > ov_start:
            bp = ov_end - ov_start
            total_bp     += bp
            weighted_sum += iv_val * bp

    return weighted_sum / total_bp if total_bp > 0 else 0.0


def compute_coverage_for_peak(args_tuple):
    """
    Worker function for multiprocessing.Pool.map().
    args_tuple: (chr_, start, end, treat_signals, ctrl_signals)
    Returns a dict with raw values, or None if both means below min_signal.
    """
    chr_, start, end, treat_signals, ctrl_signals, min_signal = args_tuple

    treat_vals = [weighted_mean_in_peak(sig, chr_, start, end) for sig in treat_signals]
    ctrl_vals  = [weighted_mean_in_peak(sig, chr_, start, end) for sig in ctrl_signals]

    treat_mean = sum(treat_vals) / len(treat_vals)
    ctrl_mean  = sum(ctrl_vals)  / len(ctrl_vals)

    if treat_mean < min_signal and ctrl_mean < min_signal:
        return None

    return {
        "chr": chr_, "start": start, "end": end,
        "treat_mean": treat_mean, "ctrl_mean": ctrl_mean,
        "treat_vals": treat_vals, "ctrl_vals": ctrl_vals,
    }
