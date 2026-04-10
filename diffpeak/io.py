"""
io.py — File parsing for peak BED files and normalized bedgraph (bgnorm) files.
"""

import os
from collections import defaultdict


def parse_peak_bed(filepath):
    """
    Parse a BED / narrowPeak / broadPeak file.
    Only the first three columns (chr, start, end) are used.
    Returns: list of (chr, start, end)
    """
    peaks = []
    with open(filepath) as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            try:
                chr_, start, end = cols[0], int(cols[1]), int(cols[2])
                if end > start:
                    peaks.append((chr_, start, end))
            except ValueError:
                pass  # skip malformed lines silently
    return peaks


def parse_bgnorm(filepath):
    """
    Parse a 4-column bedgraph-style normalized signal file:
        chr  start  end  value
    Returns: dict { chr: sorted list of (start, end, value) }
    Intervals are sorted by start position for efficient binary-search lookup.
    """
    intervals = defaultdict(list)
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                continue
            try:
                chr_ = cols[0]
                start = int(cols[1])
                end = int(cols[2])
                val = float(cols[3])
                if end > start:
                    intervals[chr_].append((start, end, val))
            except ValueError:
                pass

    # Sort by start position for fast binary search during coverage queries
    for chr_ in intervals:
        intervals[chr_].sort(key=lambda x: x[0])

    return dict(intervals)


def validate_files(file_lists):
    """
    Check that all input files exist.
    Raises FileNotFoundError listing all missing files.
    """
    missing = [f for files in file_lists for f in files if not os.path.isfile(f)]
    if missing:
        msg = "Input file(s) not found:\n" + "\n".join(f"  {f}" for f in missing)
        raise FileNotFoundError(msg)
