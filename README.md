# diffpeak

Differential peak analysis for **CUT&Tag / ChIP-seq** from normalized bedgraph (bgnorm) signal files.

- Pure Python standard library — **no numpy / scipy / bedtools required**
- Automatic **replicate detection**: single sample → log2FC only; multiple samples → Welch t-test + BH FDR correction
- **Quantile pseudocount**: data-driven pseudocount (min non-zero signal / 2) avoids fold-change compression
- **t-test in log2 space**: more appropriate for continuous normalized signals
- **Credibility score**: Z-score based signal strength evaluation of the dominant group
- **Parallel coverage computation** via Python multiprocessing
- Outputs: TSV (full statistics) + BED (IGV-ready, +/- strand coloring) + summary report + log file

---

## Installation

```bash
unzip diffpeak.zip
cd diffpeak
pip install .

# Verify
diffpeak --version
```

### Upgrade from a previous version

```bash
pip uninstall diffpeak -y
cd diffpeak
pip install .
```

---

## Usage

```
diffpeak -t TREAT [TREAT ...] -c CTRL [CTRL ...]
         --peaks FILE [FILE ...]
         [options]
```

---

## Arguments

### Signal files (required)

| Argument | Description |
|----------|-------------|
| `-t / --treat` | Treat normalized bedgraph / bgnorm file(s). Multiple files = replicates. |
| `-c / --ctrl` | Control normalized bedgraph / bgnorm file(s). Multiple files = replicates. |

### Peak files (required)

| Argument | Description |
|----------|-------------|
| `--peaks` | Peak BED file(s) — treat and ctrl mixed freely, any number. All files are pooled and merged into a single non-overlapping reference peak set. |

### Peak merging

| Argument | Default | Description |
|----------|---------|-------------|
| `--merge-gap` | `0` | Merge peaks separated by <= N bp (equivalent to `bedtools merge -d`). Recommended: 0–100 for sharp marks (H3K4me3), 200–500 for broad marks (H3K27me3). |

### Output

| Argument | Default | Description |
|----------|---------|-------------|
| `-o / --outdir` | `.` | Output directory |
| `-n / --name` | `diffpeak` | Output file prefix |

### Filter thresholds

| Argument | Default | Description |
|----------|---------|-------------|
| `--log2fc` | `1.0` | Minimum \|log2FC\| |
| `--pval` | *(off)* | Maximum p-value (replicates only, independent of --qval) |
| `--qval` | *(off)* | Maximum q-value / BH-FDR (replicates only, independent of --pval) |
| `--direction` | `both` | `both` / `up` (treat > ctrl) / `down` (ctrl > treat) |
| `--min-signal` | `0` | Minimum mean signal to include a peak (filters background noise) |
| `--cred` | *(off)* | Maximum credibility_p (e.g. 0.05). Filters by signal strength of the dominant group. See below. |
| `--min-len` | *(off)* | Minimum peak length (bp) |
| `--max-len` | *(off)* | Maximum peak length (bp) |

### Miscellaneous

| Argument | Default | Description |
|----------|---------|-------------|
| `--threads` | `4` | Parallel threads for coverage computation |

---

## Examples

### Single sample (log2FC filter only)

```bash
diffpeak \
  -t treat.bgnorm \
  -c ctrl.bgnorm \
  --peaks treat.broadPeak ctrl.broadPeak \
  --log2fc 1 --min-len 200 \
  -o results/ -n H3K27ac_BEP396
```

### With credibility filter

```bash
diffpeak \
  -t treat.bgnorm \
  -c ctrl.bgnorm \
  --peaks treat.broadPeak ctrl.broadPeak \
  --log2fc 0.5 --min-len 200 --cred 0.05 \
  -o results/ -n H3K27ac_BEP396
```

### With biological replicates (t-test + BH correction)

```bash
diffpeak \
  -t t1.bgnorm t2.bgnorm t3.bgnorm \
  -c c1.bgnorm c2.bgnorm c3.bgnorm \
  --peaks t1.broadPeak t2.broadPeak t3.broadPeak \
          c1.broadPeak c2.broadPeak c3.broadPeak \
  --log2fc 1 --pval 0.05 --cred 0.1 \
  --threads 8 \
  -o results/ -n H3K27ac_BEP396
```

### Broad histone mark

```bash
diffpeak \
  -t treat.bgnorm -c ctrl.bgnorm \
  --peaks treat.broadPeak ctrl.broadPeak \
  --merge-gap 500 --log2fc 1 --min-len 200 \
  -o results/ -n H3K27me3_BEP396
```

### In a SLURM loop

```bash
ALL_PEAKS=$(ls ./peak_calling/*.broadPeak 2>/dev/null | tr '\n' ' ')

diffpeak \
    -t ${BEP_fragments} \
    -c ${GFP_fragments} \
    --peaks ${ALL_PEAKS} \
    --log2fc 0.5 --min-len 200 --cred 0.1 \
    -o ${OUTDIR} \
    -n ${histone}_${BEP}
```

---

## Input file formats

### Peak BED file (`--peaks`)
Standard BED / broadPeak / narrowPeak — only columns 1–3 (chr, start, end) are used:
```
chr1    10090   10465   peak_name   27   .   4.39   5.59   2.71
```

### Normalized signal file (`-t` / `-c`)
4-column bedgraph (chr, start, end, normalized_value):
```
chr1    10087   10115   0.077254
chr1    10153   10169   0.154508
chr1    10169   10176   0.309016
```

---

## Output files

| File | Description |
|------|-------------|
| `<n>_diff_peaks.tsv` | Full statistics table (see columns below) |
| `<n>_diff_peaks.bed` | BED format, IGV-ready with strand-based coloring |
| `<n>_summary.txt` | Human-readable run summary |
| `<n>.log` | Detailed run log |

### TSV columns

| Column | Description |
|--------|-------------|
| chr, start, end | Peak coordinates |
| treat_mean | Mean normalized signal across treat replicates |
| ctrl_mean | Mean normalized signal across ctrl replicates |
| log2FC | log2((treat+pc) / (ctrl+pc)), pc = quantile pseudocount |
| pval | Welch t-test p-value in log2 space *(replicates only)* |
| qval | BH-corrected FDR *(replicates only)* |
| direction | increase / decrease |
| dominant_signal | Signal of the dominant group (treat if increase, ctrl if decrease) |
| credibility_p | Z-score based credibility p-value (see below) |

### BED columns

```
chr  start  end  name  igv_score  strand  log2FC  credibility_p
```

- **strand**: `+` = increase (IGV blue), `-` = decrease (IGV red)
- **igv_score**: 0–1000, derived from credibility — higher = stronger signal = darker color in IGV

---

## Statistical methods

### log2FC

Uses a **data-driven quantile pseudocount** instead of the conventional +1:

```
pseudocount = min(all non-zero signal values across all files) / 2
log2FC = log2((treat_mean + pc) / (ctrl_mean + pc))
```

This prevents fold-change compression when signal values are small decimals
(e.g. normalized bgnorm values like 0.077), where +1 would shrink a true
5-fold difference to an apparent 0.45 log2FC.

### t-test (replicates only)

Welch's t-test is performed **in log2 space** on per-replicate signal values,
which is more appropriate for continuous normalized signals than testing
raw values directly. BH correction is applied across all peaks.

### Credibility score

Answers the question: *"Is this peak's signal genuinely strong, not just
relatively higher than the other group?"*

```
increase peak → dominant group = treat
decrease peak → dominant group = ctrl
```

For each peak, the Z-score of the dominant group signal is computed
relative to the **global background distribution** of that group across
all evaluated peaks:

```
z = (dominant_signal - global_mean) / global_std
credibility_p = P(Z > z)   [one-sided normal]
```

With replicates, z-scores are computed per replicate and averaged,
so a peak that is consistently strong across all replicates gets a
lower credibility_p than one that is high in only one replicate.

**All four scenarios are handled correctly:**

| Scenario | Dominant | credibility_p | Kept? |
|----------|----------|---------------|-------|
| treat high, ctrl low (increase) | treat | small | ✅ |
| treat high, ctrl also high (increase) | treat | medium | filtered |
| treat low, ctrl high (decrease) | ctrl | small | ✅ |
| treat low, ctrl also low (decrease) | ctrl | large | filtered |

**`--cred` threshold guide:**

| `--cred` | Meaning |
|----------|---------|
| `0.05` | Top ~5% signal strength (≥1.65 SD above global mean) — strict |
| `0.10` | Top ~10% (≥1.28 SD above mean) |
| `0.20` | Top ~20% (≥0.84 SD above mean) — lenient, good starting point |
| not set | No credibility filter applied |

> **Tip**: Run once without `--cred`, inspect the `credibility_p` column
> distribution in the TSV, then choose a threshold that fits your data.

---

## Recommended workflow

```bash
# Step 1: run without credibility filter to see full results
diffpeak -t treat.bgnorm -c ctrl.bgnorm \
         --peaks treat.broadPeak ctrl.broadPeak \
         --log2fc 0.5 --min-len 200 \
         -o results/ -n test_run

# Step 2: inspect credibility_p distribution in TSV
# e.g. in R: hist(read.table("test_run_diff_peaks.tsv", header=T)$credibility_p)

# Step 3: apply appropriate --cred threshold
diffpeak -t treat.bgnorm -c ctrl.bgnorm \
         --peaks treat.broadPeak ctrl.broadPeak \
         --log2fc 0.5 --min-len 200 --cred 0.1 \
         -o results/ -n final_run
```
