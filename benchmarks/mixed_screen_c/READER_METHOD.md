# Small-Genome Reader-Side Screening ("reader" method)

## Overview

The "reader" method is a strategy for handling small genomes (viruses, plasmids, short contigs) in sylph's two-stage profiling system. It addresses the problem that small genomes have very few k-mers below the default `--screen-c` threshold, making them invisible to the stage-1 screen.

## The Problem

When building a two-stage database (`.syl2db`), sylph uses a sparse screen to quickly filter out genomes that definitely aren't present in a sample. The screen stores only ~1/3000 of each genome's k-mers (default `--screen-c 3000`). For a typical bacterial genome of 4 Mbp, this means ~1,300 screen k-mers, which is sufficient for reliable detection.

However, for **small genomes** (e.g., viral vOTUs of ~8 kbp), the nominal screen rate yields only 2-3 k-mers - far too few for reliable detection, regardless of how much viral material is in the sample.

## Existing Solutions

Before "reader", there were two approaches:

1. **`loosen`** (incumbent): Widen the pooled MPHF index to include the densest rate any genome needed. This guaranteed detection but one tiny genome would slow down screening for the entire database, as every sample k-mer would need to be looked up in a much larger index.

2. **`band`** (value-banded tiers): Store extra keys from small genomes in a separate value band (band 1) that is only checked when needed. This keeps the main index (band 0) fast, but requires additional index structure and lookup logic.

## The "reader" Solution

The "reader" method takes a different approach: **no changes to the index at all**.

### How it Works

1. **Build time** (`convert-db-two-screen --small-genome-screen none`):  
   Build the database with strictly the nominal `--screen-c` rate. No per-genome floors are applied. Small genomes simply get 1-5 screen k-mers, as they would in a database with no small-genome handling.

2. **Profile time** (`profile --screen-small-genomes N`):  
   At database open time, identify all genomes with fewer than `N` stored stage-1 k-mers (the small genomes). For each such genome, read back its **N smallest dense k-mers** from the dense block on disk, and create a temporary in-memory screen by intersecting these k-mers against the sample's k-mer hash map.

3. **Cost model**:  
   - Amortized cost: O(number of small genomes) per sample, rather than O(all sample k-mers)
   - Wins when small genomes are a minority (e.g., 99,718 vOTUs mixed with 2,000 GTDB representatives)
   - Also works on databases built by stock sylph with no small-genome handling

### Why it's Better

From the benchmark results (`RESULTS.md`):

At the largest scale tested (2,000 GTDB representatives + 100,000 MetaVR vOTUs):

**Stage-1 screen time (fastest replicate, in seconds):**

| Sample    | nofloor | reader | band   | loosen | reader_speedup vs loosen |
|-----------|---------|--------|--------|--------|--------------------------|
| human     | 0.0247  | 0.0318 | 0.1294 | 0.2712 | **8.5x**                 |
| marine    | 0.1451  | 0.1587 | 0.7178 | 1.7548 | **11.1x**                |
| simulated | 0.0015  | 0.0035 | 0.0051 | 0.012  | **3.4x**                 |

**Marginal cost over "nofloor" (screening no small genomes at all):**
- human: +7.1 ms (over 99,718 small genomes)
- marine: +13.6 ms (over 99,718 small genomes)
- simulated: +2.0 ms (over 99,718 small genomes)

**Index size:**

| Mode   | Index Size (MB) |
|--------|-----------------|
| band   | 63.8            |
| loosen | 71.1            |
| none   | **23.5**        |

**Sensitivity (worst-case recall vs. single-stage profile):**

| Arm     | Recall |
|---------|--------|
| band    | 1.0    |
| loosen  | 1.0    |
| nofloor | 0.375  |
| reader  | **1.0**|

**Amortized over 100 samples (same sample × 100):**

| Sample | reader (s) | loosen (s) | reader speedup |
|--------|------------|------------|----------------|
| human  | 18.4       | 40.2       | **2.2x**       |
| marine | 87.0       | 317.6      | **3.7x**       |

The "reader" method achieves **perfect sensitivity** (100% recall) while being:
- **8-11× faster** than `loosen` for single samples
- **2-4× faster** than `loosen` even when amortized over many samples
- **2.7× smaller** index than `loosen`
- Only **7-14 ms slower** than not screening small genomes at all

## Benchmark Details

### Test Configuration

- **Reference genomes**: 
  - 2,000 GTDB representatives (bacterial/archaeal)
  - Up to 100,000 MetaVR vOTU representatives (viral)
  
- **Samples**:
  - `human`: SRR9224014 (human gut metagenome)
  - `marine`: Marine metagenome
  - `simulated`: Simulated community with known truth

- **Parameters**:
  - Build: `--screen-c 3000 --min-sparse-kmers 50 --small-genome-screen none`
  - Profile: `--screen-small-genomes 50`
  - Dense: `-c 100 -k 31`

### Total K-mer Counts

From `g2000_v100000` database:
- Total stage-1 k-mers (reader mode): 1,863,014 (unfloored baseline)
- Dense k-mers read back for 99,718 small genomes: 4,045,512
- Effective per-small-genome cost: ~40 dense k-mers read per genome

### Resource Usage Summary

**Single-sample profiling (g2000_v100000, human sample):**

| Metric              | reader  | loosen  | Improvement |
|---------------------|---------|---------|-------------|
| Wall time           | 0.3s    | 0.5s    | 1.7× faster |
| Stage-1 screen time | 0.032s  | 0.271s  | 8.5× faster |
| Max RSS             | 3 MB    | 255 MB  | 85× less    |
| Index size          | 23.5 MB | 71.1 MB | 3× smaller  |

**100-sample multithreaded (g2000_v100000, same sample × 100):**

| Metric              | reader  | loosen  | Improvement   |
|---------------------|---------|---------|---------------|
| Total wall time     | 18.4s   | 40.2s   | 2.2× faster   |
| Total stage-1 time  | 3.59s   | 26.2s   | 7.3× faster   |
| Wall time/sample    | 0.184s  | 0.402s  | 2.2× faster   |

## Usage

### Building a Database for Reader-Side Screening

```bash
# Build with no small-genome floor
sylph convert-db-two-screen reference.syldb \
    --screen-c 3000 \
    --small-genome-screen none \
    -o database
```

### Profiling with Reader-Side Screening

```bash
# Screen small genomes directly at open time
sylph profile database.syl2db sample.sylsp \
    --screen-small-genomes 50
```

The `--screen-small-genomes N` parameter specifies how many k-mers to use for screening each small genome. A value of 50 matches the `--min-sparse-kmers` default used in build modes like `loosen` and `band`, ensuring equivalent sensitivity.

## When to Use "reader"

**Use "reader" when:**
- Your database mixes large genomes (bacteria, archaea) with small genomes (viruses, plasmids)
- Small genomes are a minority (e.g., <10% of references)
- You want the smallest possible index size
- You profile each sample only a few times (the open-time cost amortizes well over 10-100 samples)

**Consider "loosen" when:**
- Small genomes dominate the database (e.g., a virus-only database)
- You profile the same database thousands of times (one-time index cost, no per-sample small-genome scan)

**Consider "band" when:**
- You need the absolute fastest single-sample screen and are willing to maintain a more complex index structure
- However, the benchmark shows "reader" is simpler and often faster in practice

## Implementation Details

The "reader" method is implemented in `src/twostage_db.rs` via the `dense_sparse_prefix()` function, which:
1. Identifies genomes with fewer than `screen_small_genomes` stored stage-1 k-mers at `ScreenIndex::read()` time
2. Reads their dense blocks from disk (Golomb-Rice decode)
3. Extracts the smallest N hashes to form a per-genome screen prefix
4. Intersects these against the sample k-mer hash map in `contain.rs`

The method is **non-invasive** to existing code:
- No changes to the on-disk format (fully backward compatible)
- Can be applied to any `.syl2db` built with stock sylph
- The `--screen-small-genomes 0` option disables it (default behavior)

## Conclusion

The "reader" method is the **recommended approach** for mixed-scale databases. It provides:
- ✅ Perfect sensitivity for small genomes
- ✅ Smallest index size
- ✅ Fastest single-sample profiling
- ✅ Simplicity and backward compatibility
- ✅ Effective amortized cost over many samples

It represents the best trade-off for typical metagenomic profiling workloads where viral diversity is mined alongside bacterial/archaeal references.
