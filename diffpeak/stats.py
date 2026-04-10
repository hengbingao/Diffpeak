"""
stats.py — Statistical routines for diffpeak:
           log2FC, Welch's t-test, Benjamini-Hochberg FDR correction.

Pure Python standard library — no numpy/scipy dependency.
"""

import math


# ──────────────────────────────────────────────────────────────
# 1. Quantile pseudocount + log2 fold-change
# ──────────────────────────────────────────────────────────────

def compute_pseudocount(all_signals):
    """
    Compute a data-driven pseudocount from all loaded bgnorm signal files.

    Strategy: collect all non-zero signal values across every file,
    then use half of the global minimum non-zero value as the pseudocount.
    This ensures the pseudocount is small relative to real signals
    and does not compress fold-change differences.

    Args:
        all_signals : list of dicts  {chr: [(start, end, value), ...]}
                      (treat_signals + ctrl_signals combined)

    Returns:
        float -- pseudocount value (fallback: 1e-6 if no non-zero values found)
    """
    min_nonzero = None
    for sig_dict in all_signals:
        for intervals in sig_dict.values():
            for _, _, val in intervals:
                if val > 0:
                    if min_nonzero is None or val < min_nonzero:
                        min_nonzero = val
    if min_nonzero is None:
        return 1e-6
    return min_nonzero / 2.0


def log2fc(treat_mean, ctrl_mean, pseudocount=1e-6):
    """
    log2FC = log2((treat_mean + pc) / (ctrl_mean + pc))

    Uses a data-driven pseudocount (min non-zero signal / 2) instead
    of fixed +1, so fold-change differences are not compressed.
    """
    return math.log2((treat_mean + pseudocount) / (ctrl_mean + pseudocount))


def log2_values(vals, pseudocount=1e-6):
    """
    Convert signal values to log2 scale with pseudocount.
    Used for t-test in log2 space (more appropriate for continuous signals).
    """
    return [math.log2(v + pseudocount) for v in vals]


# ──────────────────────────────────────────────────────────────
# 2. Special functions (Lanczos / incomplete beta)
# ──────────────────────────────────────────────────────────────

_LANCZOS_C = [
    76.18009172947146, -86.50532032941677,  24.01409824083091,
    -1.231739572450155,  1.208650973866179e-3, -5.395239384953e-6,
]


def _lgamma(z):
    """Log-gamma via Lanczos approximation."""
    x = y = z
    tmp = x + 5.5
    tmp -= (x + 0.5) * math.log(tmp)
    ser = 1.000000000190015
    for c in _LANCZOS_C:
        y += 1
        ser += c / y
    return -tmp + math.log(2.5066282746310005 * ser / x)


def _incomplete_beta(x, a, b):
    """Regularized incomplete beta function I_x(a, b) via series expansion."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    result = 0.0
    term = 1.0
    for i in range(200):
        term *= x * (a + i) / (a + b + i)
        result += term / (a + i + 1)
        if abs(term) < 1e-12:
            break
    log_b = _lgamma(a) + _lgamma(b) - _lgamma(a + b)
    coeff = math.exp(a * math.log(x) + b * math.log(1.0 - x) - log_b)
    return coeff * (1.0 / a + result)


def _t_cdf(t_abs, df):
    """CDF of the t-distribution evaluated at |t|."""
    x = df / (df + t_abs * t_abs)
    return 1.0 - 0.5 * _incomplete_beta(x, df / 2.0, 0.5)


# ──────────────────────────────────────────────────────────────
# 3. Welch's t-test
# ──────────────────────────────────────────────────────────────

def welch_t_test(treat_vals, ctrl_vals):
    """
    Two-sample Welch's t-test (unequal variance, two-sided).

    Args:
        treat_vals : list of floats (treat replicate means at a single peak)
        ctrl_vals  : list of floats (ctrl  replicate means at a single peak)

    Returns:
        (t_statistic, p_value)
        Returns (0.0, 1.0) when sample size is insufficient or SE == 0.
    """
    na, nb = len(treat_vals), len(ctrl_vals)
    if na < 2 or nb < 2:
        return 0.0, 1.0

    mean_a = sum(treat_vals) / na
    mean_b = sum(ctrl_vals)  / nb
    var_a  = sum((v - mean_a) ** 2 for v in treat_vals) / (na - 1)
    var_b  = sum((v - mean_b) ** 2 for v in ctrl_vals)  / (nb - 1)

    se = math.sqrt(var_a / na + var_b / nb)
    if se == 0.0:
        return 0.0, 1.0

    t = (mean_a - mean_b) / se

    # Welch–Satterthwaite degrees of freedom
    num   = (var_a / na + var_b / nb) ** 2
    denom = (var_a / na) ** 2 / (na - 1) + (var_b / nb) ** 2 / (nb - 1)
    df    = num / denom if denom > 0.0 else 1.0

    p = 2.0 * (1.0 - _t_cdf(abs(t), df))
    return t, max(0.0, min(1.0, p))


# ──────────────────────────────────────────────────────────────
# 4. Benjamini-Hochberg FDR correction
# ──────────────────────────────────────────────────────────────

def bh_correction(p_values):
    """
    Benjamini-Hochberg procedure for controlling the False Discovery Rate.

    Args:
        p_values : list of raw p-values (in original order)

    Returns:
        list of q-values in the same order as the input
    """
    n = len(p_values)
    if n == 0:
        return []

    # Sort by p-value, keeping original indices
    ranked = sorted(enumerate(p_values), key=lambda x: x[1])
    q_values = [1.0] * n
    min_q = 1.0

    # Traverse from largest to smallest rank (step-down)
    for step, (orig_idx, pval) in enumerate(reversed(ranked)):
        rank = n - step          # rank is 1-based
        q = min(1.0, pval * n / rank)
        min_q = min(min_q, q)
        q_values[orig_idx] = min_q

    return q_values


# ──────────────────────────────────────────────────────────────
# 5. Peak credibility score (Z-score based)
# ──────────────────────────────────────────────────────────────

def _norm_sf(z):
    """
    Survival function of the standard normal: P(Z > z).
    Approximation using the complementary error function.
    Returns a value in (0, 1].
    """
    # P(Z > z) = 0.5 * erfc(z / sqrt(2))
    # erfc approximation (Abramowitz & Stegun 7.1.26)
    if z < 0:
        return 1.0 - _norm_sf(-z)
    t = 1.0 / (1.0 + 0.3275911 * z)
    poly = t * (0.254829592
                + t * (-0.284496736
                       + t * (1.421413741
                              + t * (-1.453152027
                                     + t * 1.061405429))))
    return poly * math.exp(-z * z / 2.0)


def compute_credibility(results):
    """
    Compute a credibility p-value for each differential peak based on
    the Z-score of the dominant group signal relative to the global
    background distribution of that group across all peaks.

    Rationale:
      A high-confidence differential peak should have a strong signal
      in the dominant group, not just a large fold-change.
      Example: treat=2.0, ctrl=0.1  →  treat is the dominant group.
               Even though ctrl is low (as expected for a real diff peak),
               the credibility is driven entirely by how strong treat is,
               not by how low ctrl is.

    Method:
      - increase peak  →  dominant signal = mean of treat_vals per replicate
      - decrease peak  →  dominant signal = mean of ctrl_vals  per replicate

      For single sample:
        z = (dominant_mean - global_mean) / global_std
        credibility_p = P(Z > z)  [one-sided]

      For replicates:
        Compute z-score for each replicate independently,
        then average the z-scores before converting to p-value.
        This accounts for replicate consistency: a peak that is
        consistently high across all replicates gets a lower p-value
        than one that is high in only one replicate.

    Global distribution:
      Computed separately for treat and ctrl across ALL evaluated peaks
      (raw_results), so the background reflects the full signal landscape,
      not just the filtered subset.

    Args:
        results : list of peak dicts with treat_vals, ctrl_vals, log2fc

    Adds in-place:
        credibility_p : float, one-sided p-value (small = more credible)
        igv_score     : int 0-1000 for BED score column (larger = more credible)
    """
    if not results:
        return results

    # ── Build global distribution for treat and ctrl ──────────
    # Flatten all per-replicate values across all peaks
    all_treat_flat = []
    all_ctrl_flat  = []
    for r in results:
        all_treat_flat.extend(r["treat_vals"])
        all_ctrl_flat.extend(r["ctrl_vals"])

    def _mean_std(vals):
        if len(vals) < 2:
            m = vals[0] if vals else 0.0
            return m, 1.0
        m = sum(vals) / len(vals)
        var = sum((v - m) ** 2 for v in vals) / len(vals)
        s = math.sqrt(var) if var > 0 else 1.0
        return m, s

    treat_mean_global, treat_std_global = _mean_std(all_treat_flat)
    ctrl_mean_global,  ctrl_std_global  = _mean_std(all_ctrl_flat)

    # ── Compute credibility for each peak ─────────────────────
    for r in results:
        is_increase = r["log2fc"] > 0

        if is_increase:
            # Dominant group = treat; compute z-score per replicate
            dom_vals    = r["treat_vals"]
            glob_mean   = treat_mean_global
            glob_std    = treat_std_global
        else:
            # Dominant group = ctrl
            dom_vals    = r["ctrl_vals"]
            glob_mean   = ctrl_mean_global
            glob_std    = ctrl_std_global

        # Z-score per replicate, then average
        z_scores = [(v - glob_mean) / glob_std for v in dom_vals]
        z_mean   = sum(z_scores) / len(z_scores)

        # One-sided p-value: P(Z > z_mean)
        cred_p = _norm_sf(z_mean)

        r["credibility_p"] = round(max(0.0, min(1.0, cred_p)), 6)
        # IGV score: higher z → lower p → higher score → darker in IGV
        r["igv_score"]     = min(1000, max(0, int((1.0 - cred_p) * 1000)))

    return results
