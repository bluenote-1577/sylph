//! Two-stage seekable genome database (`.syl2db`), automatically used by
//! `sylph query`/`sylph profile` whenever given as input (no flag needed).
//!
//! A standard sylph database (`.syldb`) is a single bincoded `Vec<GenomeSketch>`
//! that must be loaded in full to profile a sample. For two-stage profiling we
//! only ever need the *dense* k-mers of the handful of genomes a sample actually
//! contains; loading every genome's dense k-mers up front is wasteful for a large
//! reference. This module re-packs a `.syldb` into a two-stage seekable layout,
//! mirroring (in spirit) the `.sylref` format of wwood/sylph#2 but deliberately
//! simpler:
//!
//!   * **No k-mer dereplication.** Each genome keeps its *own complete* k-mer set;
//!     a k-mer shared by several genomes is stored in each of them.
//!   * **No shared/pooled hashes.** There is no conserved-k-mer pool.
//!
//! ## Layout
//!
//! ```text
//! [4]  magic  "SY2D"
//! [1]  version
//! [8]  index_offset  (u64 LE)
//! [8]  footer_offset (u64 LE)
//! ---- body ----
//! dense block 0, dense block 1, ...      (each: Golomb-Rice compressed)
//! ---- index ----
//! ScreenIndex                            (pooled-MPHF stage-1 screen index)
//! ---- footer ----
//! bincode(Footer)                        (per-genome metadata)
//! ```
//!
//! **Stage 1 (pooled MPHF, loaded fully).** The `ScreenIndex` is one minimal
//! perfect hash over the *distinct* sparse (`screen_c`) k-mers of all genomes,
//! with a multi-owner CSR mapping each sparse k-mer to the list of genomes that
//! carry it. A sample is screened by a single pass over its k-mers (work ∝
//! sample, not reference): each sample k-mer is looked up once and its coverage
//! pushed to every owning genome. The contained genomes are exactly those a
//! per-genome `get_stats` screen would find -- this is the "Path B" organisation
//! (see `experiments/7_mphf_screen_again`), but multi-owner because `.syl2db`
//! keeps shared k-mers.
//!
//! ### Small genomes and mixed screen rates
//!
//! A genome of a few kb has almost no k-mers below the nominal `screen_c`
//! threshold (a mean-length viral vOTU at `-c 100 --screen-c 3000` gets ~2), so
//! the stage-1 screen misses it however much of it the sample contains. The fix
//! is a per-genome floor (`--min-sparse-kmers`): a genome whose nominal subsample
//! falls short instead uses its `min_sparse_kmers` *smallest* dense hashes, i.e.
//! a denser, genome-specific screen rate. `--small-genome-screen` chooses how
//! that mixture of rates is stored and screened. Every mode screens the same
//! per-genome k-mer set and therefore yields the same survivors; they differ only
//! in cost:
//!
//!   * `loosen` -- one pooled MPHF, its key set widened to the densest
//!     genome-specific threshold any genome needed. Simple, but one 8 kb genome
//!     drags the *whole* database's effective screen rate toward `c`: the MPHF
//!     grows toward the dense k-mer count and every sample k-mer below
//!     `MAX/c` (i.e. all of them) has to be looked up in it.
//!   * `band` -- value-banded tiers. Keys are partitioned by hash *value*, not by
//!     genome: band 0 (`h < MAX/screen_c`) is the pooled MPHF at exactly the
//!     nominal rate, so it costs what an unloosened build costs; band 1
//!     (`MAX/screen_c <= h < MAX/dense_c`) holds only the extra keys densified
//!     genomes needed, as a sorted `Vec<u64>` behind a bloom pre-filter -- no
//!     MPHF, no fingerprints. A sample k-mer picks its band with two
//!     comparisons. Band 1 is appended after the band-0 index block, which
//!     `ScreenIndex::read` consumes sequentially, so an older reader simply does
//!     not see it (and gets the `none` behaviour).
//!   * `none` -- no floor at all; the file is what a build with no small-genome
//!     handling would write, so it stays readable by stock sylph with no loss.
//!     Small-genome sensitivity then comes from the *reader* side:
//!     `--screen-small-genomes N` reconstructs the short genomes' sparse prefix
//!     from their dense blocks at open time (`dense_sparse_prefix`) and
//!     intersects each against the sample hashmap directly. That costs
//!     O(number of short references) per sample instead of O(sample), so it wins
//!     only while short genomes are a small minority -- but it also works on
//!     databases built by stock sylph.
//!
//! **Stage 2 (dense, Golomb-Rice, loaded on demand).** Each genome's *full*
//! `genome_kmers` and `pseudotax_tracked_nonused_kmers` are an independently
//! Golomb-Rice-coded block at a known offset. Only the genomes that pass the
//! stage-1 screen are decoded. The profiling path intentionally does not retain
//! decoded blocks across samples: its permissive screen can admit many genomes
//! that fail the dense pass, and caching that growing union would trade away the
//! format's bounded-RSS advantage. `load_dense` remains available to callers
//! whose workload has a known-small, repeatedly accessed working set.

use crate::cmdline::DbConvertArgs;
use crate::constants::*;
use crate::types::*;
use boomphf::Mphf;
use fxhash::FxHashMap;
use log::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::os::unix::fs::FileExt;
use std::path::Path;
use std::sync::{Arc, Mutex};

const MAGIC: &[u8; 4] = b"SY2D";
/// Format version (dense Golomb-Rice blocks + pooled-MPHF stage-1 screen index).
const VERSION: u8 = 2;
/// magic (4) + version (1) + index offset (8) + footer offset (8)
const HEADER_LEN: u64 = 21;
/// boomphf construction gamma (space/speed trade-off), matching the ref-delta
/// sparse index.
const MPHF_GAMMA: f64 = 2.0;

// --- primitive integer / bit coding -----------------------------------------

fn write_uvarint(w: &mut Vec<u8>, mut x: u64) {
    loop {
        let mut byte = (x & 0x7f) as u8;
        x >>= 7;
        if x != 0 {
            byte |= 0x80;
        }
        w.push(byte);
        if x == 0 {
            break;
        }
    }
}

fn read_uvarint<R: Read>(r: &mut R) -> io::Result<u64> {
    let mut result: u64 = 0;
    let mut shift: u32 = 0;
    loop {
        let mut b = [0u8; 1];
        r.read_exact(&mut b)?;
        result |= ((b[0] & 0x7f) as u64) << shift;
        if b[0] & 0x80 == 0 {
            break;
        }
        shift += 7;
        if shift >= 64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "uvarint overflow",
            ));
        }
    }
    Ok(result)
}

/// LSB-first bit writer accumulating into a byte buffer.
struct BitWriter {
    buf: Vec<u8>,
    cur: u8,
    nbits: u8,
}

impl BitWriter {
    fn new() -> Self {
        BitWriter {
            buf: Vec::new(),
            cur: 0,
            nbits: 0,
        }
    }
    #[inline]
    fn write_bit(&mut self, b: u32) {
        if b != 0 {
            self.cur |= 1 << self.nbits;
        }
        self.nbits += 1;
        if self.nbits == 8 {
            self.buf.push(self.cur);
            self.cur = 0;
            self.nbits = 0;
        }
    }
    #[inline]
    fn write_bits(&mut self, val: u64, n: u32) {
        for i in 0..n {
            self.write_bit(((val >> i) & 1) as u32);
        }
    }
    #[inline]
    fn write_unary(&mut self, q: u64) {
        for _ in 0..q {
            self.write_bit(1);
        }
        self.write_bit(0);
    }
    fn finish(mut self) -> Vec<u8> {
        if self.nbits > 0 {
            self.buf.push(self.cur);
        }
        self.buf
    }
}

/// LSB-first bit reader that decodes a word at a time: bytes are buffered into a
/// 64-bit accumulator so `read_bits`/`read_unary` extract many bits per
/// instruction (shift/mask, trailing_ones) instead of one bit per call.
struct BitReader<'a> {
    buf: &'a [u8],
    pos: usize,
    acc: u64,
    nbits: u32,
}

impl<'a> BitReader<'a> {
    fn new(buf: &'a [u8]) -> Self {
        BitReader {
            buf,
            pos: 0,
            acc: 0,
            nbits: 0,
        }
    }
    /// Pull bytes into the accumulator until it holds >= 56 bits (so it never
    /// exceeds 63, keeping shifts in range) or the input is exhausted.
    #[inline]
    fn refill(&mut self) {
        while self.nbits < 56 && self.pos < self.buf.len() {
            self.acc |= (self.buf[self.pos] as u64) << self.nbits;
            self.pos += 1;
            self.nbits += 8;
        }
    }
    #[inline]
    fn read_bits(&mut self, n: u32) -> io::Result<u64> {
        if n == 0 {
            return Ok(0);
        }
        if n > 32 {
            // Split so each half fits the (>=56 bit) accumulator comfortably.
            let lo = self.read_bits(32)?;
            let hi = self.read_bits(n - 32)?;
            return Ok(lo | (hi << 32));
        }
        if self.nbits < n {
            self.refill();
            if self.nbits < n {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "bitstream truncated",
                ));
            }
        }
        let v = self.acc & ((1u64 << n) - 1);
        self.acc >>= n;
        self.nbits -= n;
        Ok(v)
    }
    /// Unary = run of 1s terminated by a 0 (matches `BitWriter::write_unary`).
    #[inline]
    fn read_unary(&mut self) -> io::Result<u64> {
        let mut q = 0u64;
        loop {
            if self.nbits == 0 {
                self.refill();
                if self.nbits == 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "bitstream truncated",
                    ));
                }
            }
            let ones = (self.acc | (1u64 << self.nbits)).trailing_ones(); // <= nbits
            if ones >= self.nbits {
                // all buffered bits are 1s; consume them and continue
                q += self.nbits as u64;
                self.acc = 0;
                self.nbits = 0;
            } else {
                // `ones` 1-bits then the terminating 0
                q += ones as u64;
                self.acc >>= ones + 1;
                self.nbits -= ones + 1;
                return Ok(q);
            }
        }
    }
}

/// Sort + delta + Golomb-Rice encode a set of hashes onto `out`. Order is not
/// preserved (hash sets are order-independent); duplicates become zero gaps and
/// are preserved. The Rice parameter is chosen from the mean gap and written
/// inline, so the block is self-delimiting given the leading count.
fn write_hashes(out: &mut Vec<u8>, hashes: &[u64]) {
    let mut sorted = hashes.to_vec();
    sorted.sort_unstable();
    write_uvarint(out, sorted.len() as u64);
    if sorted.is_empty() {
        return;
    }
    let mut deltas = Vec::with_capacity(sorted.len());
    let mut prev = 0u64;
    for &h in &sorted {
        deltas.push(h - prev);
        prev = h;
    }
    // Rice parameter k ~ log2(mean gap): near-optimal for the geometric gap
    // distribution of uniformly random hashes.
    let sum: u128 = deltas.iter().map(|&d| d as u128).sum();
    let mean = (sum / deltas.len() as u128).max(1);
    let mut k = 0u32;
    while k < 63 && (1u128 << (k + 1)) <= mean {
        k += 1;
    }
    out.push(k as u8);
    let mut bw = BitWriter::new();
    for &d in &deltas {
        bw.write_unary(d >> k);
        if k > 0 {
            bw.write_bits(d & ((1u64 << k) - 1), k);
        }
    }
    let bits = bw.finish();
    write_uvarint(out, bits.len() as u64);
    out.extend_from_slice(&bits);
}

fn read_hashes<R: Read>(r: &mut R) -> io::Result<Vec<u64>> {
    let n = read_uvarint(r)? as usize;
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut kb = [0u8; 1];
    r.read_exact(&mut kb)?;
    let k = kb[0] as u32;
    let blen = read_uvarint(r)? as usize;
    let mut bits = vec![0u8; blen];
    r.read_exact(&mut bits)?;
    let mut br = BitReader::new(&bits);
    let mut out = Vec::with_capacity(n);
    let mut prev = 0u64;
    for _ in 0..n {
        let q = br.read_unary()?;
        let low = if k > 0 { br.read_bits(k)? } else { 0 };
        let d = (q << k) | low;
        prev = prev.wrapping_add(d);
        out.push(prev);
    }
    Ok(out)
}

// --- footer (stage-1 sparse index + metadata) -------------------------------

/// Per-genome metadata. Everything here is loaded into memory when the database
/// is opened; the dense block at `dense_offset` is decoded lazily. The stage-1
/// sparse k-mers live in the pooled `ScreenIndex`, not here.
#[derive(serde::Serialize, serde::Deserialize, Clone, Debug, PartialEq)]
pub struct GenomeMeta {
    pub file_name: String,
    pub first_contig_name: String,
    pub gn_size: usize,
    pub min_spacing: usize,
    pub has_pseudotax: bool,
    pub dense_offset: u64,
}

#[derive(serde::Serialize, serde::Deserialize, Clone, Debug, PartialEq)]
pub struct Footer {
    /// Dense rate: every k-mer is kept in the dense blocks at this `-c`.
    pub c: usize,
    pub k: usize,
    /// Sparse stage-1 screen rate (`screen_c >= c`).
    pub screen_c: usize,
    pub genomes: Vec<GenomeMeta>,
}

// --- stage-1 pooled-MPHF screen index ---------------------------------------

#[inline]
fn sparse_fingerprint(h: u64) -> u32 {
    (h ^ (h >> 32)) as u32
}

/// FracMinHash threshold for the stage-1 screen: a k-mer is "sparse" iff its
/// hash is `< u64::MAX / screen_c` (the same rule as `subsample_view` and the
/// `.syl2db` build).
#[inline]
pub(crate) fn screen_threshold(screen_c: usize) -> u64 {
    u64::MAX / screen_c.max(1) as u64
}

/// How a build stores the mixture of screen rates that small genomes force (see
/// the module docs). All modes screen the same per-genome k-mer set -- a genome
/// short of `min_sparse_kmers` at the nominal rate uses its `min_sparse_kmers`
/// smallest dense hashes -- except `None`, which applies no floor at all.
#[derive(clap::ArgEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum SmallGenomeMode {
    /// Widen the single pooled MPHF to the densest rate any genome needed.
    Loosen,
    /// Keep the pooled MPHF at the nominal rate; put the extra keys of
    /// densified genomes in a separate value band.
    Band,
    /// No floor: strictly the nominal `screen_c` (see
    /// `--screen-small-genomes` for the reader-side alternative).
    None,
}

impl std::fmt::Display for SmallGenomeMode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let s = match self {
            SmallGenomeMode::Loosen => "loosen",
            SmallGenomeMode::Band => "band",
            SmallGenomeMode::None => "none",
        };
        f.write_str(s)
    }
}

// --- band-1 index (value-banded screen tiers) -------------------------------

const BAND_MAGIC: &[u8; 4] = b"BND1";

/// Two independent bloom bit positions from one hash (splitmix-style mixing;
/// the input is already a hash, but its low bits are what the FracMinHash
/// threshold constrains, so it must be re-mixed).
#[inline]
fn bloom_bits(h: u64) -> (u64, u64) {
    let mut a = h.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    a ^= a >> 29;
    let mut b = h.wrapping_mul(0xBF58_476D_1CE4_E5B9);
    b ^= b >> 31;
    (a, b)
}

/// Stage-1 screen keys in the value band `[lo, hi)` above the nominal
/// `screen_c` threshold: the extra keys that densified (small) genomes need,
/// and nothing else. Deliberately not an MPHF -- there are few enough of these
/// in the intended (mostly-large-genomes) case that a sorted `Vec<u64>` behind
/// a bloom pre-filter is cache-resident and cheaper than a perfect hash plus
/// fingerprint array.
pub struct BandIndex {
    /// Inclusive lower bound: the nominal `screen_threshold(screen_c)`.
    pub lo: u64,
    /// Exclusive upper bound: one past the largest key any genome needed.
    pub hi: u64,
    /// Sorted, distinct keys.
    keys: Vec<u64>,
    /// CSR row offsets, length `keys.len() + 1`.
    owner_offsets: Vec<u32>,
    /// CSR owner genome ids (a multiset, as in `ScreenIndex`).
    owners: Vec<u32>,
    /// Bloom pre-filter over `keys`, `BAND_BLOOM_BITS_PER_KEY` bits per key.
    bloom: Vec<u64>,
    bloom_word_mask: u64,
    /// Per genome: number of band-1 keys it owns (absent = none).
    counts: FxHashMap<u32, u32>,
}

impl BandIndex {
    /// Build from `(key, genome)` pairs; `pairs` may be unsorted and may repeat
    /// a `(key, genome)` pair (a genome carrying the same k-mer twice), which is
    /// preserved so counts match a per-genome intersection.
    fn build(mut pairs: Vec<(u64, u32)>, lo: u64, hi: u64) -> BandIndex {
        pairs.sort_unstable();
        let mut counts: FxHashMap<u32, u32> = FxHashMap::default();
        for &(_, g) in &pairs {
            *counts.entry(g).or_insert(0) += 1;
        }
        let mut keys: Vec<u64> = Vec::new();
        let mut owners: Vec<u32> = Vec::with_capacity(pairs.len());
        let mut owner_offsets: Vec<u32> = vec![0];
        let mut i = 0;
        while i < pairs.len() {
            let h = pairs[i].0;
            keys.push(h);
            while i < pairs.len() && pairs[i].0 == h {
                owners.push(pairs[i].1);
                i += 1;
            }
            owner_offsets.push(owners.len() as u32);
        }
        let (bloom, bloom_word_mask) = Self::build_bloom(&keys);
        BandIndex {
            lo,
            hi,
            keys,
            owner_offsets,
            owners,
            bloom,
            bloom_word_mask,
            counts,
        }
    }

    fn build_bloom(keys: &[u64]) -> (Vec<u64>, u64) {
        let want_bits = (keys.len() * BAND_BLOOM_BITS_PER_KEY).max(64);
        let n_words = (want_bits / 64).next_power_of_two();
        let mut bloom = vec![0u64; n_words];
        let mask = (n_words - 1) as u64;
        for &h in keys {
            let (a, b) = bloom_bits(h);
            bloom[((a >> 6) & mask) as usize] |= 1u64 << (a & 63);
            bloom[((b >> 6) & mask) as usize] |= 1u64 << (b & 63);
        }
        (bloom, mask)
    }

    #[inline]
    fn bloom_maybe(&self, h: u64) -> bool {
        let (a, b) = bloom_bits(h);
        let m = self.bloom_word_mask;
        self.bloom[((a >> 6) & m) as usize] & (1u64 << (a & 63)) != 0
            && self.bloom[((b >> 6) & m) as usize] & (1u64 << (b & 63)) != 0
    }

    /// Genomes owning `h`, or empty if `h` is outside the band or not a key.
    /// Two comparisons and (almost always) one bloom word for a non-key.
    #[inline]
    fn owners_of(&self, h: u64) -> &[u32] {
        if h < self.lo || h >= self.hi || !self.bloom_maybe(h) {
            return &[];
        }
        match self.keys.binary_search(&h) {
            Ok(i) => {
                &self.owners[self.owner_offsets[i] as usize..self.owner_offsets[i + 1] as usize]
            }
            Err(_) => &[],
        }
    }

    pub fn num_keys(&self) -> usize {
        self.keys.len()
    }
    pub fn num_owners(&self) -> usize {
        self.owners.len()
    }
    pub fn num_genomes(&self) -> usize {
        self.counts.len()
    }
    #[inline]
    fn count(&self, g: u32) -> u32 {
        self.counts.get(&g).copied().unwrap_or(0)
    }

    /// Append the band block. Keys are delta+uvarint coded (they are sorted and
    /// dense within a narrow value range) and the bloom is rebuilt on read, so
    /// only keys/owners cost file space.
    fn write_to_vec(&self, out: &mut Vec<u8>) {
        out.extend_from_slice(BAND_MAGIC);
        out.extend_from_slice(&self.lo.to_le_bytes());
        out.extend_from_slice(&self.hi.to_le_bytes());
        write_uvarint(out, self.keys.len() as u64);
        write_uvarint(out, self.owners.len() as u64);
        write_uvarint(out, self.counts.len() as u64);
        let mut prev = self.lo;
        for &h in &self.keys {
            write_uvarint(out, h - prev);
            prev = h;
        }
        for i in 0..self.keys.len() {
            write_uvarint(
                out,
                (self.owner_offsets[i + 1] - self.owner_offsets[i]) as u64,
            );
        }
        for &o in &self.owners {
            out.extend_from_slice(&o.to_le_bytes());
        }
        let mut entries: Vec<(u32, u32)> = self.counts.iter().map(|(&g, &c)| (g, c)).collect();
        entries.sort_unstable();
        for (g, c) in entries {
            write_uvarint(out, g as u64);
            write_uvarint(out, c as u64);
        }
    }

    /// Parse a band block whose magic has already been confirmed present.
    fn read(r: &mut &[u8]) -> io::Result<BandIndex> {
        let mut magic = [0u8; 4];
        r.read_exact(&mut magic)?;
        debug_assert_eq!(&magic, BAND_MAGIC);
        let mut buf = [0u8; 8];
        r.read_exact(&mut buf)?;
        let lo = u64::from_le_bytes(buf);
        r.read_exact(&mut buf)?;
        let hi = u64::from_le_bytes(buf);
        let n_keys = read_uvarint(r)? as usize;
        let n_owners = read_uvarint(r)? as usize;
        let n_counts = read_uvarint(r)? as usize;
        let mut keys = Vec::with_capacity(n_keys);
        let mut prev = lo;
        for _ in 0..n_keys {
            prev += read_uvarint(r)?;
            keys.push(prev);
        }
        let mut owner_offsets = Vec::with_capacity(n_keys + 1);
        owner_offsets.push(0u32);
        let mut acc = 0u32;
        for _ in 0..n_keys {
            acc += read_uvarint(r)? as u32;
            owner_offsets.push(acc);
        }
        if acc as usize != n_owners {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "band-1 owner counts do not sum to the owner total",
            ));
        }
        let mut owners = Vec::with_capacity(n_owners);
        for _ in 0..n_owners {
            let mut b = [0u8; 4];
            r.read_exact(&mut b)?;
            owners.push(u32::from_le_bytes(b));
        }
        let mut counts: FxHashMap<u32, u32> =
            FxHashMap::with_capacity_and_hasher(n_counts, Default::default());
        for _ in 0..n_counts {
            let g = read_uvarint(r)? as u32;
            let c = read_uvarint(r)? as u32;
            counts.insert(g, c);
        }
        let (bloom, bloom_word_mask) = Self::build_bloom(&keys);
        Ok(BandIndex {
            lo,
            hi,
            keys,
            owner_offsets,
            owners,
            bloom,
            bloom_word_mask,
            counts,
        })
    }
}

/// Pooled stage-1 screen index: one MPHF over the distinct sparse k-mers of all
/// genomes, plus a multi-owner CSR (a k-mer may belong to several genomes).
/// Owners are a *multiset* -- a genome appears once per occurrence of the k-mer
/// in its sparse set -- so duplicate k-mers count exactly as the per-genome
/// `get_stats` loop would.
pub struct ScreenIndex {
    pub screen_c: usize,
    pub k: usize,
    mphf: Mphf<u64>,
    /// Per slot: fingerprint of the k-mer, to reject foreign (non-indexed) hashes
    /// that the MPHF would otherwise map to an arbitrary slot.
    fingerprints: Vec<u32>,
    /// CSR row offsets, length `n_slots + 1`.
    owner_offsets: Vec<u32>,
    /// CSR owner genome ids (flat); slot `s` owns `owners[off[s]..off[s+1]]`.
    owners: Vec<u32>,
    /// Per genome: number of sparse k-mers in *this* (band-0) structure. Use
    /// `screen_kmers` for the ANI denominator, which also counts band 1.
    pub sparse_count: Vec<u32>,
    /// Optional band-1 tier (`--small-genome-screen band`), holding the keys
    /// above the nominal threshold that densified genomes needed.
    pub band: Option<BandIndex>,
}

impl ScreenIndex {
    /// Build from each genome's sparse (`screen_c`) k-mers. Owners are kept as a
    /// multiset so the screen reproduces the per-genome `get_stats` counts
    /// exactly (including any duplicate k-mers within a genome).
    pub fn build(sparse_per_genome: &[Vec<u64>], screen_c: usize, k: usize) -> ScreenIndex {
        let sparse_count: Vec<u32> = sparse_per_genome.iter().map(|v| v.len() as u32).collect();
        let total: usize = sparse_per_genome.iter().map(|v| v.len()).sum();

        // (k-mer, genome) pairs, sorted so equal k-mers form contiguous runs.
        let mut pairs: Vec<(u64, u32)> = Vec::with_capacity(total);
        for (g, v) in sparse_per_genome.iter().enumerate() {
            for &h in v {
                pairs.push((h, g as u32));
            }
        }
        pairs.sort_unstable();

        // Distinct keys for the MPHF.
        let mut keys: Vec<u64> = Vec::new();
        for &(h, _) in &pairs {
            if keys.last() != Some(&h) {
                keys.push(h);
            }
        }
        let mphf = Mphf::new_parallel(MPHF_GAMMA, &keys, Some(0));
        let n_slots = keys.len();

        // Per-slot owner counts -> CSR offsets.
        let mut fingerprints = vec![0u32; n_slots];
        let mut owner_offsets = vec![0u32; n_slots + 1];
        let mut i = 0;
        while i < pairs.len() {
            let h = pairs[i].0;
            let mut j = i;
            while j < pairs.len() && pairs[j].0 == h {
                j += 1;
            }
            let slot = mphf.hash(&h) as usize;
            fingerprints[slot] = sparse_fingerprint(h);
            owner_offsets[slot + 1] = (j - i) as u32; // count, prefix-summed below
            i = j;
        }
        for s in 0..n_slots {
            owner_offsets[s + 1] += owner_offsets[s];
        }

        // Fill owners using a per-slot write cursor.
        let mut owners = vec![0u32; total];
        let mut cursor: Vec<u32> = owner_offsets[..n_slots].to_vec();
        let mut i = 0;
        while i < pairs.len() {
            let h = pairs[i].0;
            let slot = mphf.hash(&h) as usize;
            let mut j = i;
            while j < pairs.len() && pairs[j].0 == h {
                owners[cursor[slot] as usize] = pairs[j].1;
                cursor[slot] += 1;
                j += 1;
            }
            i = j;
        }

        ScreenIndex {
            screen_c,
            k,
            mphf,
            fingerprints,
            owner_offsets,
            owners,
            sparse_count,
            band: None,
        }
    }

    pub fn num_genomes(&self) -> usize {
        self.sparse_count.len()
    }

    /// Total number of stage-1 screen k-mers for genome `g` across all bands --
    /// the `n_kmers` ANI denominator, and what a per-genome `get_stats` screen
    /// would use.
    #[inline]
    pub fn screen_kmers(&self, g: u32) -> u32 {
        self.sparse_count[g as usize] + self.band.as_ref().map_or(0, |b| b.count(g))
    }

    /// Single inverted pass over the sample: for each sample k-mer below the
    /// screen threshold with non-zero count, look it up and push its coverage to
    /// every owning genome. Returns `genome -> matched coverage counts`, exactly
    /// the per-genome `covs` a `get_stats(.., None, ..)` screen would collect.
    /// With a band-1 tier present, each sample k-mer picks its band with two
    /// comparisons: below the nominal threshold it goes to the MPHF, otherwise
    /// (and below the band ceiling) to the small band-1 structure.
    pub fn gather_hits(&self, sample: &SequencesSketch) -> FxHashMap<u32, Vec<u32>> {
        let thresh = screen_threshold(self.screen_c);
        let band_hi = self.band.as_ref().map_or(0, |b| b.hi);
        let mut hits: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
        for (&h, &cnt) in sample.kmer_counts.iter() {
            if cnt == 0 {
                continue;
            }
            if h >= thresh {
                if h < band_hi {
                    for &g in self.band.as_ref().unwrap().owners_of(h) {
                        hits.entry(g).or_default().push(cnt);
                    }
                }
                continue;
            }
            if let Some(slot) = self.mphf.try_hash(&h) {
                let slot = slot as usize;
                if slot < self.fingerprints.len()
                    && self.fingerprints[slot] == sparse_fingerprint(h)
                {
                    let lo = self.owner_offsets[slot] as usize;
                    let hi = self.owner_offsets[slot + 1] as usize;
                    for &g in &self.owners[lo..hi] {
                        hits.entry(g).or_default().push(cnt);
                    }
                }
            }
        }
        hits
    }

    /// Serialize the index into `out` (raw little-endian blocks).
    fn write_to_vec(&self, out: &mut Vec<u8>) -> io::Result<()> {
        let mphf_bytes = bincode::serialize(&self.mphf).map_err(io::Error::other)?;
        write_uvarint(out, mphf_bytes.len() as u64);
        out.extend_from_slice(&mphf_bytes);
        write_uvarint(out, self.fingerprints.len() as u64); // n_slots
        write_uvarint(out, self.owners.len() as u64);
        write_uvarint(out, self.sparse_count.len() as u64); // n_genomes
        for &fp in &self.fingerprints {
            out.extend_from_slice(&fp.to_le_bytes());
        }
        for &o in &self.owner_offsets {
            out.extend_from_slice(&o.to_le_bytes());
        }
        for &o in &self.owners {
            out.extend_from_slice(&o.to_le_bytes());
        }
        for &c in &self.sparse_count {
            out.extend_from_slice(&c.to_le_bytes());
        }
        // Band 1 is *appended*: `read` consumes everything above sequentially and
        // stops, so a reader that predates banding sees a plain nominal-rate
        // index (and screens exactly band 0) rather than failing.
        if let Some(band) = &self.band {
            band.write_to_vec(out);
        }
        Ok(())
    }

    /// Parse an index block produced by `write_to_vec`. `screen_c`/`k` come from
    /// the footer (not duplicated in the block).
    fn read(mut r: &[u8], screen_c: usize, k: usize) -> io::Result<ScreenIndex> {
        let mphf_len = read_uvarint(&mut r)? as usize;
        let mut mphf_bytes = vec![0u8; mphf_len];
        r.read_exact(&mut mphf_bytes)?;
        let mphf: Mphf<u64> = bincode::deserialize(&mphf_bytes)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let n_slots = read_uvarint(&mut r)? as usize;
        let n_owners = read_uvarint(&mut r)? as usize;
        let n_genomes = read_uvarint(&mut r)? as usize;

        let read_u32_vec = |r: &mut &[u8], n: usize| -> io::Result<Vec<u32>> {
            let mut v = Vec::with_capacity(n);
            for _ in 0..n {
                let mut buf = [0u8; 4];
                r.read_exact(&mut buf)?;
                v.push(u32::from_le_bytes(buf));
            }
            Ok(v)
        };
        let fingerprints = read_u32_vec(&mut r, n_slots)?;
        let owner_offsets = read_u32_vec(&mut r, n_slots + 1)?;
        let owners = read_u32_vec(&mut r, n_owners)?;
        let sparse_count = read_u32_vec(&mut r, n_genomes)?;
        // Optional trailing band-1 block (absent in a `loosen`/`none` build).
        let band = if r.len() >= BAND_MAGIC.len() && &r[..BAND_MAGIC.len()] == BAND_MAGIC {
            Some(BandIndex::read(&mut r)?)
        } else {
            None
        };

        Ok(ScreenIndex {
            screen_c,
            k,
            mphf,
            fingerprints,
            owner_offsets,
            owners,
            sparse_count,
            band,
        })
    }
}

// --- writing ----------------------------------------------------------------

/// If even the adaptive floor (see `write_two_stage_db`'s `min_sparse_kmers`
/// parameter) can't be reached (i.e. the genome has fewer than this many
/// k-mers in total, dense included), warn loudly: this genome will also drag
/// the whole database's effective stage-1 screen rate down toward its own
/// dense rate (see `write_two_stage_db`), and its detection reliability at
/// this size is inherently poor.
const SPARSE_WARN_THRESHOLD: usize = 20;

/// Number of offending genomes named in an aggregated build diagnostic. The
/// counts matter; the names are only there to make the problem findable.
const DIAGNOSTIC_EXAMPLES: usize = 3;

/// Re-pack genome sketches into the two-stage seekable layout and write to `w`.
/// `screen_c` is the (coarser) stage-1 subsampling rate; it must be `>= c`.
/// Dense blocks are Golomb-Rice coded. Each genome's sparse stage-1 set is
/// selected at `screen_c` if that clears `min_sparse_kmers`, otherwise (unless
/// `mode` is `None`) it is that genome's `min_sparse_kmers` smallest dense
/// hashes -- a denser, genome-specific rate. `mode` decides how that mixture of
/// rates is stored (see `SmallGenomeMode` and the module docs):
///
///   * `Loosen` widens the single pooled MPHF, so the database-wide effective
///     screen rate in the footer is derived from whichever genome ended up
///     densest and the pooled index's early-exit filter stays correct for every
///     genome.
///   * `Band` keeps the pooled MPHF at exactly `screen_c` and puts the extra
///     keys in an appended band-1 block, leaving the footer's `screen_c`
///     nominal.
///   * `None` applies no floor, so nothing is densified at all.
///
/// `min_sparse_kmers` must be `>= 1` -- callers must enforce this (see
/// `run_db_convert`'s validation), since 0 would let a genome's sparse set end
/// up empty, making it silently invisible to the stage-1 screen forever.
pub fn write_two_stage_db<W: Write>(
    mut w: W,
    sketches: &[GenomeSketch],
    screen_c: usize,
    min_sparse_kmers: usize,
    min_contain: usize,
    mode: SmallGenomeMode,
) -> io::Result<()> {
    let c = sketches.first().map(|s| s.c).unwrap_or(0);
    let k = sketches.first().map(|s| s.k).unwrap_or(0);
    // FracMinHash threshold for the nominal sparse subsample (same rule as
    // subsample_view), and for this database's own dense rate -- every dense
    // k-mer already satisfies `h < dense_thresh` by construction, so it is a
    // safe (and tight) upper bound for a genome forced to use all of them.
    let nominal_thresh = screen_threshold(screen_c);
    let dense_thresh = screen_threshold(c);

    let mut body: Vec<u8> = Vec::new();
    let mut genomes: Vec<GenomeMeta> = Vec::with_capacity(sketches.len());
    let mut sparse_per_genome: Vec<Vec<u64>> = Vec::with_capacity(sketches.len());
    // Densest per-genome threshold used by any genome; the database-wide
    // effective screen threshold must be at least this permissive so the
    // pooled index's early-exit never prunes a real stored key.
    let mut thresh_needed: u64 = 0;
    // `Band` mode only: (key, genome) pairs above the nominal threshold.
    let mut band_pairs: Vec<(u64, u32)> = Vec::new();
    // Per-genome diagnostics are counted and reported once at the end: a viral
    // database densifies millions of genomes, and one `warn!` each would bury
    // everything else (and dwarf the database itself in log bytes).
    let mut n_densified = 0usize;
    let mut n_too_few = 0usize;
    let mut too_few_examples: Vec<String> = Vec::new();
    let mut skipped_examples: Vec<String> = Vec::new();
    let mut n_skipped = 0usize;

    for gs in sketches {
        let dense_total = gs.genome_kmers.len();
        if dense_total < min_contain {
            n_skipped += 1;
            if skipped_examples.len() < DIAGNOSTIC_EXAMPLES {
                skipped_examples.push(format!(
                    "'{}' (file {}, {} dense k-mers)",
                    gs.first_contig_name, gs.file_name, dense_total
                ));
            }
            trace!(
                "genome '{}' (file {}) has only {} dense k-mers (< --min-contain={}); excluded",
                gs.first_contig_name,
                gs.file_name,
                dense_total,
                min_contain
            );
            continue;
        }

        let dense_offset = HEADER_LEN + body.len() as u64;
        write_hashes(&mut body, &gs.genome_kmers);
        match &gs.pseudotax_tracked_nonused_kmers {
            Some(p) => {
                body.push(1);
                write_hashes(&mut body, p);
            }
            None => body.push(0),
        }

        // The nominal (`screen_c`) subsample, which is all any mode ever puts in
        // band 0 / the pooled MPHF.
        let nominal: Vec<u64> = gs
            .genome_kmers
            .iter()
            .copied()
            .filter(|&h| h < nominal_thresh)
            .collect();
        // Extra keys above the nominal threshold needed to reach the floor: the
        // difference between the `min_sparse_kmers` smallest dense hashes and
        // `nominal`. Empty unless this genome is too small (or `mode` is `None`).
        let mut extra: Vec<u64> = Vec::new();
        let mut genome_thresh = nominal_thresh;
        if mode != SmallGenomeMode::None && nominal.len() < min_sparse_kmers {
            n_densified += 1;
            // O(n) partial selection of the `want` smallest hashes; `want` is the
            // whole dense set when the genome has fewer than the floor to begin
            // with.
            let want = min_sparse_kmers.min(dense_total);
            let mut v = gs.genome_kmers.clone();
            if want < v.len() {
                v.select_nth_unstable(want - 1);
                v.truncate(want);
            }
            let kth_val = v.iter().copied().max().unwrap_or(0);
            debug_assert!(kth_val < dense_thresh || dense_total == 0);
            genome_thresh = kth_val.saturating_add(1);
            extra = v.into_iter().filter(|&h| h >= nominal_thresh).collect();
            trace!(
                "genome '{}' (file {}) has only {} sparse k-mers at --screen-c {}; densified to \
                 {} using a genome-specific rate",
                gs.first_contig_name,
                gs.file_name,
                nominal.len(),
                screen_c,
                want
            );
        }

        // "Didn't reach your own configured target" -- stays at the fixed
        // SPARSE_WARN_THRESHOLD for any target >= that (the common case),
        // shrinks only if min_sparse_kmers itself was set below it, so it never
        // misfires on a genome that hit its own (deliberately small) target.
        let effective_warn_threshold = SPARSE_WARN_THRESHOLD.min(min_sparse_kmers);
        if nominal.len() + extra.len() < effective_warn_threshold {
            n_too_few += 1;
            if too_few_examples.len() < DIAGNOSTIC_EXAMPLES {
                too_few_examples.push(format!(
                    "'{}' (file {}, {} screen k-mers)",
                    gs.first_contig_name,
                    gs.file_name,
                    nominal.len() + extra.len()
                ));
            }
        }

        thresh_needed = thresh_needed.max(genome_thresh);
        let g_id = genomes.len() as u32;
        match mode {
            // One pooled structure: the extra keys join the nominal ones and the
            // database-wide screen rate is loosened to cover them.
            SmallGenomeMode::Loosen => {
                let mut sparse = nominal;
                sparse.extend_from_slice(&extra);
                sparse_per_genome.push(sparse);
            }
            // Partition by hash value: nominal keys to band 0, extra keys to
            // band 1, which leaves band 0 exactly as an unfloored build.
            SmallGenomeMode::Band => {
                band_pairs.extend(extra.iter().map(|&h| (h, g_id)));
                sparse_per_genome.push(nominal);
            }
            SmallGenomeMode::None => sparse_per_genome.push(nominal),
        }
        genomes.push(GenomeMeta {
            file_name: gs.file_name.clone(),
            first_contig_name: gs.first_contig_name.clone(),
            gn_size: gs.gn_size,
            min_spacing: gs.min_spacing,
            has_pseudotax: gs.pseudotax_tracked_nonused_kmers.is_some(),
            dense_offset,
        });
    }

    if genomes.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "all {} input genome(s) have fewer than --min-contain={} dense k-mers; nothing to write",
                sketches.len(),
                min_contain
            ),
        ));
    }

    if n_skipped > 0 {
        warn!(
            "{} of {} genome(s) had fewer than --min-contain={} dense k-mers and were excluded \
             from this two-stage database (they could never pass `profile`/`query`'s hit \
             threshold at that setting anyway, or are likely erroneous/fragmentary). Examples: {}",
            n_skipped,
            sketches.len(),
            min_contain,
            skipped_examples.join(", ")
        );
    }
    if n_densified > 0 {
        let extra_note = match mode {
            SmallGenomeMode::Loosen => {
                "they force the WHOLE database's stage-1 screen rate down toward the dense -c, \
                 which slows the screen for every sample (see --small-genome-screen band)"
            }
            SmallGenomeMode::Band => {
                "their extra keys live in the band-1 tier, leaving the band-0 screen rate nominal"
            }
            SmallGenomeMode::None => unreachable!("no genome is densified in `none` mode"),
        };
        info!(
            "{} of {} genome(s) had fewer than --min-sparse-kmers={} k-mers at --screen-c {} and \
             use a denser genome-specific stage-1 screen rate; {}",
            n_densified,
            genomes.len(),
            min_sparse_kmers,
            screen_c,
            extra_note
        );
    }
    if n_too_few > 0 {
        warn!(
            "{} of {} genome(s) have fewer than {} stage-1 screen k-mers even after \
             densification: their whole dense sketch is that small, so detection reliability at \
             this size is inherently poor. Consider excluding tiny genomes/contigs/plasmids, or \
             check that they were sketched as intended. Examples: {}",
            n_too_few,
            genomes.len(),
            SPARSE_WARN_THRESHOLD.min(min_sparse_kmers),
            too_few_examples.join(", ")
        );
    }

    // Database-wide effective screen rate: safe by the number-theory identity
    // floor(a / floor(a/b)) >= b for positive integers a, b -- so
    // screen_threshold(effective_screen_c) >= thresh_needed, and every
    // genome's actually-stored k-mers (each < their own genome_thresh <=
    // thresh_needed) pass the pooled index's early-exit filter. Only `loosen`
    // pays this: the other modes keep every band-0 key below the nominal
    // threshold, so their footer rate stays nominal (and an older reader that
    // ignores band 1 still filters correctly).
    let effective_screen_c = if thresh_needed == 0 || mode != SmallGenomeMode::Loosen {
        screen_c
    } else {
        (u64::MAX / thresh_needed).max(1) as usize
    };
    if effective_screen_c != screen_c {
        info!(
            "effective stage-1 screen -c adjusted from {} to {} due to small genome(s), \
             see warnings above",
            screen_c, effective_screen_c
        );
    }

    // Pooled stage-1 screen index (band 0), then the optional band-1 tier.
    let mut screen_index = ScreenIndex::build(&sparse_per_genome, effective_screen_c, k);
    drop(sparse_per_genome);
    if mode == SmallGenomeMode::Band && !band_pairs.is_empty() {
        let n_pairs = band_pairs.len();
        let band = BandIndex::build(band_pairs, nominal_thresh, thresh_needed);
        info!(
            "band-1 tier: {} distinct keys ({} owners) over {} genome(s), value band \
             [MAX/{}, MAX/{})",
            band.num_keys(),
            n_pairs,
            band.num_genomes(),
            screen_c,
            (u64::MAX / thresh_needed.max(1)).max(1)
        );
        screen_index.band = Some(band);
    }
    let mut index_block: Vec<u8> = Vec::new();
    screen_index.write_to_vec(&mut index_block)?;

    let footer = Footer {
        c,
        k,
        screen_c: effective_screen_c,
        genomes,
    };
    let footer_bytes = bincode::serialize(&footer).map_err(io::Error::other)?;
    let index_offset = HEADER_LEN + body.len() as u64;
    let footer_offset = index_offset + index_block.len() as u64;

    // One machine-readable line summarising what this build cost, so a
    // benchmark harness can compare modes without re-reading the database.
    // Debug-level: the format is for parsers, not for people.
    debug!(
        "bench_build mode={} screen_c={} effective_screen_c={} min_sparse_kmers={} genomes={} \
         densified={} band_keys={} band_owners={} dense_bytes={} index_bytes={} footer_bytes={}",
        mode,
        screen_c,
        effective_screen_c,
        min_sparse_kmers,
        footer.genomes.len(),
        n_densified,
        screen_index.band.as_ref().map_or(0, |b| b.num_keys()),
        screen_index.band.as_ref().map_or(0, |b| b.num_owners()),
        body.len(),
        index_block.len(),
        footer_bytes.len(),
    );

    w.write_all(MAGIC)?;
    w.write_all(&[VERSION])?;
    w.write_all(&index_offset.to_le_bytes())?;
    w.write_all(&footer_offset.to_le_bytes())?;
    w.write_all(&body)?;
    w.write_all(&index_block)?;
    w.write_all(&footer_bytes)?;
    Ok(())
}

// --- opened database --------------------------------------------------------

/// Backing store for the dense blocks.
///   * `File`  - positional `read_at` (pread) of just the requested block bytes.
///     No shared cursor, so concurrent reads from any number of threads need no
///     lock; only the touched block bytes (plus reclaimable OS page cache) cost
///     memory, so RSS stays low. This is the path used for `.syl2db` files.
///   * `Owned` - whole file in memory; for in-memory readers / tests.
enum DenseData {
    File(File),
    Owned(Vec<u8>),
}

impl DenseData {
    /// Whole-file bytes; only valid for the in-memory backing (used to parse the
    /// header/footer). The `File` backing is read positionally via `with_block`.
    #[inline]
    fn bytes(&self) -> &[u8] {
        match self {
            DenseData::Owned(v) => &v[..],
            DenseData::File(_) => unreachable!("File-backed db is read via with_block"),
        }
    }
}

/// An opened two-stage database. Construction loads only the stage-1 sparse
/// index (the bincoded footer); per-genome dense blocks are decoded on demand,
/// in parallel (each decode positionally reads its own block, no shared cursor
/// or lock).
pub struct TwoStageDb {
    pub c: usize,
    pub k: usize,
    pub screen_c: usize,
    /// Path this database was opened from (empty for in-memory readers); used
    /// only for diagnostics.
    pub path: String,
    /// File offset where the dense-block region ends (start of the screen index);
    /// used to bound the last genome's block for positional reads.
    index_offset: u64,
    genomes: Vec<GenomeMeta>,
    /// Pooled stage-1 screen index (Path B). Replaces the per-genome sparse
    /// sketches; querying a sample against it yields the contained genomes.
    pub screen_index: ScreenIndex,
    /// Reader-side small-genome screen (`--screen-small-genomes`); see
    /// `SmallGenomeScreen`.
    pub small_screen: Option<SmallGenomeScreen>,
    data: DenseData,
    cache: Mutex<FxHashMap<u32, Arc<GenomeSketch>>>,
}

/// Reader-side small-genome screen: the stage-1 screen k-mers of genomes the
/// *stored* index left short are reconstructed from their dense blocks when the
/// database is opened, and intersected against the sample directly instead of
/// through the pooled index.
///
/// This needs no index change at all -- a database written with
/// `--small-genome-screen none` is what a build with no small-genome handling
/// writes, so this also gives small-genome sensitivity on stock-sylph databases.
/// The trade-off is where the work goes: cost is O(retained reference k-mers)
/// per sample rather than O(sample), which is a win only while the short genomes
/// are a small minority of the reference.
/// Retained k-mers are stored flat (CSR), not as a `Vec` per genome: a viral
/// reference can put millions of genomes in here, and a `Vec` each would add an
/// allocation header per genome and a pointer chase per genome per sample.
pub struct SmallGenomeScreen {
    /// Genomes the stored index left below `floor`, ascending.
    genome_ids: Vec<u32>,
    /// Row offsets into `hashes`, length `genome_ids.len() + 1`.
    offsets: Vec<u32>,
    /// Each genome's `floor` smallest dense hashes -- exactly the set a
    /// densified build stores for it.
    hashes: Vec<u64>,
    /// Screen k-mer counts (the ANI denominator) for those genomes.
    counts: FxHashMap<u32, u32>,
    /// One past the largest retained hash, i.e. the densest rate this screen
    /// reaches (the reader's analogue of a `loosen` build's widened threshold).
    densest_thresh: u64,
}

impl SmallGenomeScreen {
    pub fn num_genomes(&self) -> usize {
        self.genome_ids.len()
    }
    pub fn num_kmers(&self) -> usize {
        self.hashes.len()
    }
}

/// Parse the magic + version header; return `(index_offset, footer_offset)`.
fn parse_header(hdr: &[u8]) -> io::Result<(u64, u64)> {
    if hdr.len() < HEADER_LEN as usize || &hdr[0..4] != MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "not a sylph two-stage database",
        ));
    }
    if hdr[4] != VERSION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unsupported two-stage database version",
        ));
    }
    let index_offset = u64::from_le_bytes(hdr[5..13].try_into().unwrap());
    let footer_offset = u64::from_le_bytes(hdr[13..21].try_into().unwrap());
    Ok((index_offset, footer_offset))
}

/// Assemble a `TwoStageDb` from its parsed footer + screen index + backing store.
fn build_db(
    footer: Footer,
    index_offset: u64,
    screen_index: ScreenIndex,
    data: DenseData,
    path: String,
) -> TwoStageDb {
    TwoStageDb {
        c: footer.c,
        k: footer.k,
        screen_c: footer.screen_c,
        path,
        index_offset,
        genomes: footer.genomes,
        screen_index,
        small_screen: None,
        data,
        cache: Mutex::new(FxHashMap::default()),
    }
}

/// Parse the header + index + footer of a `.syl2db` already resident in `data`.
fn from_bytes(data: DenseData) -> io::Result<TwoStageDb> {
    let bytes = data.bytes();
    let (index_offset, footer_offset) = parse_header(bytes)?;
    if index_offset > footer_offset || footer_offset as usize > bytes.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "two-stage database offsets out of range",
        ));
    }
    let footer: Footer = bincode::deserialize(&bytes[footer_offset as usize..])
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let screen_index = ScreenIndex::read(
        &bytes[index_offset as usize..footer_offset as usize],
        footer.screen_c,
        footer.k,
    )?;
    Ok(build_db(
        footer,
        index_offset,
        screen_index,
        data,
        String::new(),
    ))
}

/// Open a `.syl2db` from an in-memory reader (reads it all into memory).
pub fn open<R: Read>(mut r: R) -> io::Result<TwoStageDb> {
    let mut v = Vec::new();
    r.read_to_end(&mut v)?;
    from_bytes(DenseData::Owned(v))
}

impl TwoStageDb {
    pub fn len(&self) -> usize {
        self.genomes.len()
    }
    pub fn is_empty(&self) -> bool {
        self.genomes.is_empty()
    }

    /// Source file name of genome `g` (for `--screen-dump` and diagnostics).
    pub fn genome_file_name(&self, g: u32) -> &str {
        &self.genomes[g as usize].file_name
    }

    /// First contig name of genome `g` (for duplicate-genome detection across
    /// combined sources -- see `contain()`).
    pub fn genome_first_contig_name(&self, g: u32) -> &str {
        &self.genomes[g as usize].first_contig_name
    }

    /// Full metadata for genome `g` (for `sylph inspect`) -- cheap, no dense
    /// decode needed, `GenomeMeta` is fully loaded when the database is opened.
    pub fn genome_meta(&self, g: u32) -> &GenomeMeta {
        &self.genomes[g as usize]
    }

    /// End offset of genome `g`'s region (start of the next genome's block, or
    /// the screen index for the last genome). Genomes are stored in ascending
    /// offset, and the index block immediately follows the last dense block.
    fn block_end(&self, g: u32) -> u64 {
        let gi = g as usize;
        if gi + 1 < self.genomes.len() {
            self.genomes[gi + 1].dense_offset
        } else {
            self.index_offset
        }
    }

    /// Run `f` on genome `g`'s whole on-disk region (`genome_kmers` block, the
    /// pseudotax flag, and the optional pseudotax block). For the file backing
    /// this is one positional read of just that region; for the owned backing it
    /// is a zero-copy slice. No shared cursor, so it is safe to call concurrently
    /// from many threads.
    fn with_block<T>(&self, g: u32, f: impl FnOnce(&[u8]) -> io::Result<T>) -> io::Result<T> {
        let start = self.genomes[g as usize].dense_offset as usize;
        match &self.data {
            DenseData::Owned(v) => f(&v[start..]),
            DenseData::File(file) => {
                let end = self.block_end(g) as usize;
                let mut buf = vec![0u8; end - start];
                file.read_exact_at(&mut buf, start as u64)?;
                f(&buf)
            }
        }
    }

    /// The `n` smallest hashes of genome `g`'s dense block -- exactly the sparse
    /// set a densified (`--min-sparse-kmers n`) build would have stored for it,
    /// since `read_hashes` returns the block in ascending order. Reads only the
    /// `genome_kmers` block, not the (longer) pseudotax one.
    pub fn dense_sparse_prefix(&self, g: u32, n: usize) -> io::Result<Vec<u64>> {
        let mut hashes = self.with_block(g, |bytes| {
            let mut cur = bytes;
            read_hashes(&mut cur)
        })?;
        hashes.truncate(n);
        Ok(hashes)
    }

    /// Build the reader-side small-genome screen (see `SmallGenomeScreen`):
    /// every genome with fewer than `floor` *stored* screen k-mers gets its
    /// `floor` smallest dense hashes read back and retained in memory.
    pub fn enable_small_genome_screen(&mut self, floor: usize) -> io::Result<()> {
        let genome_ids: Vec<u32> = (0..self.len() as u32)
            .filter(|&g| (self.screen_index.screen_kmers(g) as usize) < floor)
            .collect();
        let per_genome: Vec<Vec<u64>> = genome_ids
            .par_iter()
            .map(|&g| self.dense_sparse_prefix(g, floor))
            .collect::<io::Result<Vec<_>>>()?;

        // Row offsets are u32, as everywhere else in this format. 12.7 M viral
        // genomes at the default floor need 0.6 G of them, but a large enough
        // `floor` on a large enough reference would wrap silently, so refuse.
        let total_kmers: usize = per_genome.iter().map(Vec::len).sum();
        if total_kmers > u32::MAX as usize {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "--screen-small-genomes {} would retain {} k-mers across {} genomes, more than \
                     this index can address; use a smaller value",
                    floor,
                    total_kmers,
                    genome_ids.len()
                ),
            ));
        }
        let mut offsets: Vec<u32> = Vec::with_capacity(genome_ids.len() + 1);
        offsets.push(0);
        let mut hashes: Vec<u64> = Vec::with_capacity(total_kmers);
        let mut counts: FxHashMap<u32, u32> =
            FxHashMap::with_capacity_and_hasher(genome_ids.len(), Default::default());
        let mut densest_thresh = 0u64;
        for (&g, v) in genome_ids.iter().zip(&per_genome) {
            counts.insert(g, v.len() as u32);
            // `dense_sparse_prefix` returns the block in ascending order, so the
            // last hash is this genome's own (densest) screen threshold.
            if let Some(&last) = v.last() {
                densest_thresh = densest_thresh.max(last.saturating_add(1));
            }
            hashes.extend_from_slice(v);
            offsets.push(hashes.len() as u32);
        }
        self.small_screen = Some(SmallGenomeScreen {
            genome_ids,
            offsets,
            hashes,
            counts,
            densest_thresh,
        });
        Ok(())
    }

    /// Stage-1 screen k-mer count for genome `g` (the ANI denominator), taking
    /// whichever screen actually covers it.
    #[inline]
    pub fn screen_kmers(&self, g: u32) -> u32 {
        if let Some(s) = &self.small_screen {
            if let Some(&c) = s.counts.get(&g) {
                return c;
            }
        }
        self.screen_index.screen_kmers(g)
    }

    /// The densest stage-1 screen rate anything in this database is screened at.
    /// Coarse-rate assumptions (e.g. scaling `-M` from the dense rate to the
    /// screen rate) must use this, not the footer's nominal `screen_c`, or they
    /// over-reject the genomes that needed a denser rate.
    pub fn effective_screen_c(&self) -> usize {
        let mut thresh = screen_threshold(self.screen_c);
        if let Some(b) = &self.screen_index.band {
            thresh = thresh.max(b.hi);
        }
        if let Some(s) = &self.small_screen {
            thresh = thresh.max(s.densest_thresh);
        }
        (u64::MAX / thresh.max(1)).max(1) as usize
    }

    /// Stage 1: matched coverage counts per genome, from the pooled index (plus
    /// its band-1 tier) and, if enabled, the reader-side small-genome screen.
    pub fn screen(&self, sample: &SequencesSketch) -> FxHashMap<u32, Vec<u32>> {
        let mut hits = self.screen_index.gather_hits(sample);
        if let Some(s) = &self.small_screen {
            // The retained prefix is a superset of whatever the stored index
            // holds for these genomes, so it replaces (never adds to) their
            // entry -- otherwise shared k-mers would be counted twice.
            let extra: Vec<(u32, Vec<u32>)> = (0..s.genome_ids.len())
                .into_par_iter()
                .filter_map(|i| {
                    let (lo, hi) = (s.offsets[i] as usize, s.offsets[i + 1] as usize);
                    let mut covs: Vec<u32> = Vec::new();
                    for h in &s.hashes[lo..hi] {
                        if let Some(&cnt) = sample.kmer_counts.get(h) {
                            if cnt != 0 {
                                covs.push(cnt);
                            }
                        }
                    }
                    if covs.is_empty() {
                        None
                    } else {
                        Some((s.genome_ids[i], covs))
                    }
                })
                .collect();
            for (g, covs) in extra {
                hits.insert(g, covs);
            }
        }
        hits
    }

    /// Decode genome `g`'s full dense `GenomeSketch` without touching the cache.
    /// Concurrent calls from different threads do not contend on any shared
    /// cursor (each does its own positional read). The two-stage pass-1 uses this
    /// to decode each survivor into a short-lived sketch, probe it, and drop it
    /// unless it passes -- so the discarded majority is never cached.
    pub fn decode_dense(&self, g: u32) -> io::Result<GenomeSketch> {
        let meta = &self.genomes[g as usize];
        let (genome_kmers, pseudotax) = self.with_block(g, |bytes| {
            let mut cur = bytes;
            let gk = read_hashes(&mut cur)?;
            let mut flag = [0u8; 1];
            cur.read_exact(&mut flag)?;
            let pt = if flag[0] != 0 {
                Some(read_hashes(&mut cur)?)
            } else {
                None
            };
            Ok((gk, pt))
        })?;
        Ok(GenomeSketch {
            genome_kmers,
            pseudotax_tracked_nonused_kmers: pseudotax,
            file_name: meta.file_name.clone(),
            first_contig_name: meta.first_contig_name.clone(),
            c: self.c,
            k: self.k,
            gn_size: meta.gn_size,
            min_spacing: meta.min_spacing,
        })
    }

    /// Decode genome `g`'s full dense `GenomeSketch`, caching it across calls.
    ///
    /// This cache is unbounded. The profiling path deliberately uses
    /// `decode_dense` instead, because caching every permissive screen survivor
    /// across a large or diverse sample cohort can eventually retain a large
    /// fraction of the dense database in memory.
    pub fn load_dense(&self, g: u32) -> io::Result<Arc<GenomeSketch>> {
        if let Some(a) = self.cache.lock().unwrap().get(&g) {
            return Ok(a.clone());
        }
        let sketch = Arc::new(self.decode_dense(g)?);
        self.cache.lock().unwrap().insert(g, sketch.clone());
        Ok(sketch)
    }
}

/// Open a `.syl2db` file from a path. Only the header + footer (the stage-1
/// sparse index) are read up front; dense blocks are read positionally on demand
/// during profiling, so opening is cheap and RSS stays low.
pub fn open_file(path: &str) -> io::Result<TwoStageDb> {
    let file = File::open(path)?;
    let mut hdr = [0u8; HEADER_LEN as usize];
    file.read_exact_at(&mut hdr, 0)?;
    let (index_offset, footer_offset) = parse_header(&hdr)?;
    let flen = file.metadata()?.len();
    if index_offset > footer_offset || footer_offset > flen {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "two-stage database offsets out of range",
        ));
    }
    let mut fbytes = vec![0u8; (flen - footer_offset) as usize];
    file.read_exact_at(&mut fbytes, footer_offset)?;
    let footer: Footer =
        bincode::deserialize(&fbytes).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let mut ibytes = vec![0u8; (footer_offset - index_offset) as usize];
    file.read_exact_at(&mut ibytes, index_offset)?;
    let screen_index = ScreenIndex::read(&ibytes, footer.screen_c, footer.k)?;
    Ok(build_db(
        footer,
        index_offset,
        screen_index,
        DenseData::File(file),
        path.to_string(),
    ))
}

// --- CLI handler ------------------------------------------------------------

fn load_genome_sketches(path: &str) -> Vec<GenomeSketch> {
    let file = File::open(path).unwrap_or_else(|_| panic!("Could not open {}", path));
    let reader = BufReader::with_capacity(10_000_000, file);
    bincode::deserialize_from(reader)
        .unwrap_or_else(|_| panic!("{} is not a valid database sketch (.syldb)", path))
}

pub fn run_db_convert(args: DbConvertArgs) {
    let level = if args.trace {
        log::LevelFilter::Trace
    } else if args.debug {
        log::LevelFilter::Debug
    } else {
        log::LevelFilter::Info
    };
    simple_logger::SimpleLogger::new()
        .with_level(level)
        .init()
        .unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .ok();

    if args.files.is_empty() {
        error!("No genome database sketches (*.syldb) supplied; exiting");
        std::process::exit(1);
    }

    let mut sketches: Vec<GenomeSketch> = Vec::new();
    for f in &args.files {
        info!("Loading genome sketches from {}", f);
        sketches.extend(load_genome_sketches(f));
    }
    if sketches.is_empty() {
        error!("No genome sketches found in input; exiting");
        std::process::exit(1);
    }

    let c = sketches[0].c;
    let k = sketches[0].k;
    for s in &sketches {
        if s.c != c || s.k != k {
            error!("Input sketches have inconsistent -c/-k; exiting");
            std::process::exit(1);
        }
    }
    if sketches
        .iter()
        .any(|s| s.pseudotax_tracked_nonused_kmers.is_none())
    {
        error!(
            "Some input genomes were sketched with --disable-profiling (no profiling k-mers), \
             which a two-stage database always needs (it is used by both `query` and \
             `profile`); re-sketch without --disable-profiling. Exiting"
        );
        std::process::exit(1);
    }
    if args.screen_c < c {
        error!(
            "--screen-c ({}) must be >= the database -c ({}); the screen can only be made sparser, never denser. Exiting",
            args.screen_c, c
        );
        std::process::exit(1);
    }
    if args.min_sparse_kmers == 0 {
        error!(
            "--min-sparse-kmers must be >= 1 (0 can silently make small genomes invisible to the stage-1 screen). Exiting"
        );
        std::process::exit(1);
    }

    let out = if args.output.ends_with(TWO_STAGE_DB_SUFFIX) {
        args.output.clone()
    } else {
        format!("{}{}", args.output, TWO_STAGE_DB_SUFFIX)
    };
    if let Some(parent) = Path::new(&out).parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).ok();
        }
    }
    info!(
        "Converting {} genomes (dense -c {}, stage-1 screen -c {}, small-genome screen {}) -> {}",
        sketches.len(),
        c,
        args.screen_c,
        args.small_genome_screen,
        out
    );
    let w =
        BufWriter::new(File::create(&out).unwrap_or_else(|_| panic!("Could not create {}", out)));
    write_two_stage_db(
        w,
        &sketches,
        args.screen_c,
        args.min_sparse_kmers,
        args.min_contain,
        args.small_genome_screen,
    )
    .unwrap_or_else(|e| panic!("Failed to write {}: {}", out, e));
    info!("Wrote two-stage database to {}", out);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip_hashes(input: &[u64]) {
        let mut expected = input.to_vec();
        expected.sort_unstable();
        let mut buf = Vec::new();
        write_hashes(&mut buf, input);
        let mut r = &buf[..];
        assert_eq!(read_hashes(&mut r).unwrap(), expected);
        assert!(r.is_empty(), "read_hashes left trailing bytes");
    }

    #[test]
    fn hashes_roundtrip_various() {
        roundtrip_hashes(&[]);
        roundtrip_hashes(&[0]);
        roundtrip_hashes(&[42]);
        roundtrip_hashes(&[5, 5, 5]); // duplicates -> zero gaps
        roundtrip_hashes(&[u64::MAX, 0, 1, u64::MAX / 2]);
        roundtrip_hashes(&[10, 9, 8, 7, 6, 5, 4, 3, 2, 1]);
        // many uniformly-spread hashes (the realistic FracMinHash case)
        let mut v = Vec::new();
        let mut x = 0xdead_beef_u64;
        for _ in 0..5000 {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            v.push(x >> 8); // keep them in a sub-range like a fracminhash threshold
        }
        roundtrip_hashes(&v);
    }

    fn gsketch(name: &str, kmers: Vec<u64>, pt: Option<Vec<u64>>) -> GenomeSketch {
        GenomeSketch {
            genome_kmers: kmers,
            pseudotax_tracked_nonused_kmers: pt,
            file_name: name.to_string(),
            first_contig_name: format!("{}_c1", name),
            c: 50,
            k: 31,
            gn_size: 12345,
            min_spacing: 30,
        }
    }

    #[test]
    fn db_write_open_load_roundtrip() {
        // screen_c = 200 (coarser than c = 50): the sparse subset keeps hashes
        // below u64::MAX/200. Each genome has > SPARSE_TARGET_MIN_DEFAULT dense
        // k-mers below thresh, so both take the nominal (non-adaptive) branch
        // and the stored screen_c should round-trip to the requested nominal rate.
        let thresh = u64::MAX / 200;
        let mut g0_kmers: Vec<u64> = (1..=60u64).collect();
        g0_kmers.extend([thresh - 1, thresh + 10, thresh * 3, 9_000_000_000]);
        let mut g1_kmers: Vec<u64> = (100..=161u64).collect();
        g1_kmers.extend([thresh + 1, thresh * 2, 123_456_789_000, 7]);
        let sketches = vec![
            gsketch("g0.fa", g0_kmers.clone(), Some(vec![100, 200, 300])),
            gsketch("g1.fa", g1_kmers.clone(), Some(vec![])),
        ];

        let mut buf = Vec::new();
        write_two_stage_db(
            &mut buf,
            &sketches,
            200,
            SPARSE_TARGET_MIN_DEFAULT,
            0,
            SmallGenomeMode::Loosen,
        )
        .unwrap();
        let db = open(std::io::Cursor::new(buf)).unwrap();

        assert_eq!(db.c, 50);
        assert_eq!(db.k, 31);
        // Both genomes take the nominal branch (>= SPARSE_TARGET_MIN_DEFAULT
        // below thresh), so the effective screen_c round-trips from the
        // nominal one via the same formula the implementation uses.
        let expected_screen_c = (u64::MAX / screen_threshold(200)).max(1) as usize;
        assert_eq!(db.screen_c, expected_screen_c);
        assert_eq!(db.len(), 2);

        // stage-1 index: per-genome sparse count = fracminhash subset size at screen_c
        let expect_sparse = |ks: &[u64]| -> Vec<u64> {
            let mut v: Vec<u64> = ks.iter().copied().filter(|&h| h < thresh).collect();
            v.sort_unstable();
            v
        };
        assert_eq!(
            db.screen_index.sparse_count[0] as usize,
            expect_sparse(&g0_kmers).len()
        );
        assert_eq!(
            db.screen_index.sparse_count[1] as usize,
            expect_sparse(&g1_kmers).len()
        );

        // stage-2 dense block reconstructs the exact genome k-mers + pseudotax
        let d0 = db.load_dense(0).unwrap();
        let mut got = d0.genome_kmers.clone();
        got.sort_unstable();
        let mut exp = g0_kmers.clone();
        exp.sort_unstable();
        assert_eq!(got, exp);
        assert_eq!(d0.c, 50);
        assert_eq!(d0.k, 31);
        assert_eq!(d0.gn_size, 12345);
        assert_eq!(d0.file_name, "g0.fa");
        assert_eq!(
            d0.pseudotax_tracked_nonused_kmers,
            Some(vec![100, 200, 300])
        );

        let d1 = db.load_dense(1).unwrap();
        let mut got1 = d1.genome_kmers.clone();
        got1.sort_unstable();
        let mut exp1 = g1_kmers.clone();
        exp1.sort_unstable();
        assert_eq!(got1, exp1);
        assert_eq!(d1.pseudotax_tracked_nonused_kmers, Some(vec![]));

        // second load hits the cache and returns the same data
        let d0b = db.load_dense(0).unwrap();
        assert_eq!(d0b.genome_kmers, d0.genome_kmers);
    }

    fn sample_from(counts: &[(u64, u32)]) -> SequencesSketch {
        let mut s = SequencesSketch::new(String::new(), 50, 31, false, None, 0.0);
        for &(h, c) in counts {
            s.kmer_counts.insert(h, c);
        }
        s
    }

    /// `gather_hits` must reproduce, per genome, the exact coverage multiset that
    /// the per-genome `get_stats` loop collects: intersection of the genome's
    /// sparse k-mers with the sample, with duplicate k-mers counted per
    /// occurrence and shared k-mers credited to every owner.
    #[test]
    fn screen_index_matches_per_genome_intersection() {
        let screen_c = 100usize;
        let thresh = screen_threshold(screen_c);
        // All below threshold so every k-mer is "sparse". g0 has a duplicate (5);
        // 5 and 9 are shared across genomes.
        let sparse = vec![
            vec![5u64, 5, 9, 11],  // g0: duplicate 5
            vec![9u64, 11, 20],    // g1
            vec![5u64, 30, 40, 9], // g2
        ];
        for v in &sparse {
            assert!(v.iter().all(|&h| h < thresh));
        }
        let idx = ScreenIndex::build(&sparse, screen_c, 31);
        assert_eq!(idx.sparse_count, vec![4u32, 3, 4]);

        // Sample: some matching k-mers (with counts), a zero-count k-mer (ignored),
        // a k-mer above threshold (ignored), and a foreign k-mer (no owner).
        let sample = sample_from(&[
            (5, 7),
            (9, 3),
            (11, 2),
            (40, 9),
            (99, 0),         // zero count -> ignored
            (thresh + 1, 5), // above screen threshold -> ignored
            (123456, 4),     // foreign -> no owner
        ]);

        let hits = idx.gather_hits(&sample);

        // Brute-force per-genome reference (mirrors get_stats winner_map=None).
        let mut expected: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
        for (g, v) in sparse.iter().enumerate() {
            for &h in v {
                if h < thresh {
                    if let Some(&c) = sample.kmer_counts.get(&h) {
                        if c != 0 {
                            expected.entry(g as u32).or_default().push(c);
                        }
                    }
                }
            }
        }
        let norm = |m: &FxHashMap<u32, Vec<u32>>| -> Vec<(u32, Vec<u32>)> {
            let mut out: Vec<(u32, Vec<u32>)> = m
                .iter()
                .map(|(&g, v)| {
                    let mut v = v.clone();
                    v.sort_unstable();
                    (g, v)
                })
                .collect();
            out.sort();
            out
        };
        assert_eq!(norm(&hits), norm(&expected));
        // g0 sees 5 twice (duplicate) + 9 + 11 -> covs {7,7,3,2}
        let mut g0 = hits[&0].clone();
        g0.sort_unstable();
        assert_eq!(g0, vec![2, 3, 7, 7]);
    }

    /// The index survives a serialize/parse round-trip through the file format.
    #[test]
    fn screen_index_roundtrips_through_db() {
        let thresh = u64::MAX / 200;
        // Filler below-thresh k-mers push each genome past SPARSE_TARGET_MIN_DEFAULT
        // so both take the nominal (non-adaptive) branch, matching the probe values
        // used against `sample` below (which assume plain screen_c=200 filtering).
        let mut g0 = vec![5u64, 5, thresh - 1, thresh + 9]; // last is above screen thresh
        g0.extend(1000..1062u64);
        let mut g1 = vec![5u64, 7, thresh - 2];
        g1.extend(2000..2062u64);
        let sketches = vec![
            gsketch("g0.fa", g0.clone(), Some(vec![1])),
            gsketch("g1.fa", g1.clone(), Some(vec![2])),
        ];
        let mut buf = Vec::new();
        write_two_stage_db(
            &mut buf,
            &sketches,
            200,
            SPARSE_TARGET_MIN_DEFAULT,
            0,
            SmallGenomeMode::Loosen,
        )
        .unwrap();
        let db = open(std::io::Cursor::new(buf)).unwrap();

        let sample = sample_from(&[(5, 4), (7, 6), (thresh - 1, 1), (thresh - 2, 9)]);
        let hits = db.screen_index.gather_hits(&sample);
        // g0: 5 twice + (thresh-1) once -> {4,4,1}; g1: 5 + 7 + (thresh-2) -> {4,6,9}
        let mut g0h = hits[&0].clone();
        g0h.sort_unstable();
        assert_eq!(g0h, vec![1, 4, 4]);
        let mut g1h = hits[&1].clone();
        g1h.sort_unstable();
        assert_eq!(g1h, vec![4, 6, 9]);
    }

    /// Genomes with <= SPARSE_TARGET_MIN_DEFAULT total dense k-mers use all of
    /// them as their sparse set (regardless of the nominal screen_c), and loosen
    /// the database-wide effective screen_c to cover them.
    #[test]
    fn write_two_stage_db_small_genomes_use_all_kmers() {
        let small_a: Vec<u64> = (1..=15u64).collect(); // 15 total: < SPARSE_WARN_THRESHOLD
        let small_b: Vec<u64> = (1000..=1034u64).collect(); // 35 total: >= warn, < target
        let sketches = vec![
            gsketch("a.fa", small_a.clone(), None),
            gsketch("b.fa", small_b.clone(), None),
        ];
        let mut buf = Vec::new();
        write_two_stage_db(
            &mut buf,
            &sketches,
            3000,
            SPARSE_TARGET_MIN_DEFAULT,
            0,
            SmallGenomeMode::Loosen,
        )
        .unwrap();
        let db = open(std::io::Cursor::new(buf)).unwrap();

        assert_eq!(db.screen_index.sparse_count[0] as usize, small_a.len());
        assert_eq!(db.screen_index.sparse_count[1] as usize, small_b.len());

        // Both genomes are far below SPARSE_TARGET_MIN_DEFAULT, so both use every
        // dense k-mer they have and the database-wide effective screen_c is
        // loosened to whatever covers the largest of those (here b's 1034) --
        // never coarser, or the pooled index's early-exit would prune real keys.
        let largest_stored = *small_b.iter().max().unwrap();
        let expected_screen_c = (u64::MAX / (largest_stored + 1)).max(1) as usize;
        assert_eq!(db.screen_c, expected_screen_c);
        assert!(
            screen_threshold(db.screen_c) > largest_stored,
            "the effective screen threshold must cover every stored k-mer"
        );
        // ... and the screen really does find them.
        let sample = sample_from(&[(small_a[0], 2), (largest_stored, 5)]);
        let hits = db.screen_index.gather_hits(&sample);
        assert_eq!(hits[&0], vec![2]);
        assert_eq!(hits[&1], vec![5]);
    }

    /// A genome whose nominal (screen_c) sparse subset falls short of
    /// SPARSE_TARGET_MIN_DEFAULT, but which has enough total dense k-mers to
    /// clear the target, adaptively selects exactly the SPARSE_TARGET_MIN_DEFAULT
    /// smallest hashes -- not just however many happened to be below the
    /// nominal threshold.
    #[test]
    fn write_two_stage_db_adaptive_floor_selects_smallest() {
        let screen_c = 3000usize;
        let nominal_thresh = screen_threshold(screen_c);
        let dense_thresh = screen_threshold(50); // c=50, matching `gsketch`
        // 200 k-mers, all strictly between nominal_thresh and dense_thresh, so
        // the nominal filter keeps none of them (forcing the adaptive branch),
        // evenly spaced so the 50 smallest are exactly known.
        let step = (dense_thresh - nominal_thresh) / 250;
        let kmers: Vec<u64> = (0..200u64).map(|i| nominal_thresh + 1 + i * step).collect();
        assert!(kmers.iter().all(|&h| h > nominal_thresh && h < dense_thresh));

        let sketches = vec![gsketch("c.fa", kmers.clone(), None)];
        let mut buf = Vec::new();
        write_two_stage_db(
            &mut buf,
            &sketches,
            screen_c,
            SPARSE_TARGET_MIN_DEFAULT,
            0,
            SmallGenomeMode::Loosen,
        )
        .unwrap();
        let db = open(std::io::Cursor::new(buf)).unwrap();

        assert_eq!(db.screen_index.sparse_count[0] as usize, SPARSE_TARGET_MIN_DEFAULT);

        // The effective screen_c is derived exactly from the 50th-smallest
        // (index 49) k-mer's threshold -- denser (smaller) than the nominal
        // screen_c, since this genome needed a finer rate to reach the floor.
        let kth_val = kmers[49]; // kmers is already sorted ascending by construction
        let expected_screen_c = (u64::MAX / (kth_val + 1)).max(1) as usize;
        assert_eq!(db.screen_c, expected_screen_c);
        assert!(db.screen_c < screen_c, "adaptive floor should force a denser effective screen_c");

        // Confirm it's specifically the 50 *smallest* that were kept: a sample
        // containing the 50th-smallest k-mer (index 49) must hit; one
        // containing only the 51st-smallest (index 50, just past the floor)
        // must not.
        let sample_in = sample_from(&[(kmers[49], 3)]);
        let hits_in = db.screen_index.gather_hits(&sample_in);
        assert_eq!(hits_in.get(&0).map(|v| v.len()), Some(1));

        let sample_out = sample_from(&[(kmers[50], 3)]);
        let hits_out = db.screen_index.gather_hits(&sample_out);
        assert!(hits_out.get(&0).is_none(), "51st-smallest k-mer should not have been kept");
    }

    // --- small-genome screen modes ------------------------------------------

    const MODE_SCREEN_C: usize = 3000;

    /// One genome big enough to have plenty of k-mers at the nominal screen rate,
    /// and one so small that *all* of its k-mers sit above the nominal threshold
    /// (the viral/plasmid case that motivates the floor).
    fn big_and_tiny() -> (Vec<GenomeSketch>, Vec<u64>, Vec<u64>) {
        let nominal_thresh = screen_threshold(MODE_SCREEN_C);
        let dense_thresh = screen_threshold(50); // c=50, matching `gsketch`
        let big: Vec<u64> = (1..=100u64).map(|i| i * (nominal_thresh / 200)).collect();
        let step = (dense_thresh - nominal_thresh) / 100;
        let tiny: Vec<u64> = (0..40u64).map(|i| nominal_thresh + 1 + i * step).collect();
        assert!(big.iter().all(|&h| h < nominal_thresh));
        assert!(tiny.iter().all(|&h| h > nominal_thresh && h < dense_thresh));
        let sketches = vec![
            gsketch("big.fa", big.clone(), None),
            gsketch("tiny.fa", tiny.clone(), None),
        ];
        (sketches, big, tiny)
    }

    fn build_with_mode(sketches: &[GenomeSketch], mode: SmallGenomeMode) -> Vec<u8> {
        let mut buf = Vec::new();
        write_two_stage_db(
            &mut buf,
            sketches,
            MODE_SCREEN_C,
            SPARSE_TARGET_MIN_DEFAULT,
            0,
            mode,
        )
        .unwrap();
        buf
    }

    fn sorted_hits(hits: &FxHashMap<u32, Vec<u32>>) -> Vec<(u32, Vec<u32>)> {
        let mut out: Vec<(u32, Vec<u32>)> = hits
            .iter()
            .map(|(&g, v)| {
                let mut v = v.clone();
                v.sort_unstable();
                (g, v)
            })
            .collect();
        out.sort();
        out
    }

    /// The three modes are meant to be different *storage* for the same screen,
    /// not different screens: `loosen`, `band` and `none` + the reader-side
    /// small-genome screen must all give the same per-genome k-mer counts and the
    /// same matched coverages. `none` on its own is the sensitivity loss they all
    /// exist to avoid, so it must (and does) miss the tiny genome entirely.
    #[test]
    fn small_genome_modes_screen_identically() {
        let (sketches, big, tiny) = big_and_tiny();
        let sample = sample_from(&[(big[0], 4), (tiny[0], 7), (tiny[39], 9)]);

        let loosen = open(std::io::Cursor::new(build_with_mode(
            &sketches,
            SmallGenomeMode::Loosen,
        )))
        .unwrap();
        let band = open(std::io::Cursor::new(build_with_mode(
            &sketches,
            SmallGenomeMode::Band,
        )))
        .unwrap();
        let mut none = open(std::io::Cursor::new(build_with_mode(
            &sketches,
            SmallGenomeMode::None,
        )))
        .unwrap();

        // Without a floor the tiny genome has no screen k-mers at all and is
        // invisible, however much of it the sample contains.
        assert_eq!(none.screen_kmers(1), 0);
        assert!(none.screen(&sample).get(&1).is_none());

        none.enable_small_genome_screen(SPARSE_TARGET_MIN_DEFAULT)
            .unwrap();
        let reader = none;

        for db in [&loosen, &band, &reader] {
            assert_eq!(db.screen_kmers(0) as usize, big.len());
            assert_eq!(db.screen_kmers(1) as usize, tiny.len());
        }
        let expected = sorted_hits(&loosen.screen(&sample));
        assert_eq!(expected, vec![(0, vec![4]), (1, vec![7, 9])]);
        assert_eq!(sorted_hits(&band.screen(&sample)), expected);
        assert_eq!(sorted_hits(&reader.screen(&sample)), expected);

        // Only `loosen` pays for the tiny genome database-wide; the other two
        // keep the nominal rate for everything the tiny genome does not own, but
        // all three report the same densest (effective) rate, so downstream
        // k-mer-count scaling behaves identically.
        assert!(loosen.screen_c < MODE_SCREEN_C);
        assert_eq!(band.screen_c, MODE_SCREEN_C);
        assert_eq!(reader.screen_c, MODE_SCREEN_C);
        assert_eq!(band.effective_screen_c(), loosen.effective_screen_c());
        assert_eq!(reader.effective_screen_c(), loosen.effective_screen_c());
    }

    /// Band 1 is *appended* to the band-0 index block, so a `band` database's
    /// index is byte-for-byte a `none` database's index plus a trailing block: a
    /// reader that predates banding consumes band 0 and stops, seeing a valid
    /// nominal-rate index rather than an error.
    #[test]
    fn band_index_extends_unfloored_index() {
        let (sketches, _, tiny) = big_and_tiny();
        let band_bytes = build_with_mode(&sketches, SmallGenomeMode::Band);
        let none_bytes = build_with_mode(&sketches, SmallGenomeMode::None);

        let index_block = |b: &[u8]| -> Vec<u8> {
            let (i, f) = parse_header(b).unwrap();
            b[i as usize..f as usize].to_vec()
        };
        let (band_index, none_index) = (index_block(&band_bytes), index_block(&none_bytes));
        assert!(
            band_index.starts_with(&none_index),
            "band-0 keys/fingerprints/counts must be identical to an unfloored build"
        );
        assert_eq!(
            &band_index[none_index.len()..none_index.len() + 4],
            BAND_MAGIC
        );

        // The dense blocks are untouched by banding, so only the index differs
        // (and, in the header, the footer offset the longer index pushes along).
        let dense_body = |b: &[u8]| {
            let (i, _) = parse_header(b).unwrap();
            b[HEADER_LEN as usize..i as usize].to_vec()
        };
        assert_eq!(dense_body(&band_bytes), dense_body(&none_bytes));

        // Band 1 holds exactly the tiny genome's keys, and nothing else.
        let db = open(std::io::Cursor::new(band_bytes)).unwrap();
        let band = db.screen_index.band.as_ref().unwrap();
        assert_eq!(band.num_keys(), tiny.len());
        assert_eq!(band.num_genomes(), 1);
        assert_eq!(band.count(1) as usize, tiny.len());
    }

    /// `none` means no floor at all, so the file it writes cannot depend on
    /// `--min-sparse-kmers` -- that is what makes it byte-identical to a build
    /// with no small-genome handling, and therefore still readable by stock
    /// sylph.
    #[test]
    fn none_mode_ignores_the_sparse_floor() {
        let (sketches, _, _) = big_and_tiny();
        let mut with_1 = Vec::new();
        write_two_stage_db(
            &mut with_1,
            &sketches,
            MODE_SCREEN_C,
            1,
            0,
            SmallGenomeMode::None,
        )
        .unwrap();
        let mut with_500 = Vec::new();
        write_two_stage_db(
            &mut with_500,
            &sketches,
            MODE_SCREEN_C,
            500,
            0,
            SmallGenomeMode::None,
        )
        .unwrap();
        assert_eq!(with_1, with_500);
    }

    /// The reader-side screen reconstructs a genome's sparse prefix from its
    /// dense block: the `n` smallest hashes, which is exactly what a densified
    /// build stores. Enabling it on a database that already densified must
    /// therefore change nothing -- it either has nothing to do, or (for a genome
    /// whose whole dense sketch is smaller than the floor, so no build could
    /// reach it) re-derives the same set and replaces it with itself.
    #[test]
    fn reader_side_screen_matches_stored_prefix() {
        let (sketches, big, tiny) = big_and_tiny();
        let db = open(std::io::Cursor::new(build_with_mode(
            &sketches,
            SmallGenomeMode::None,
        )))
        .unwrap();

        let mut expected: Vec<u64> = big.clone();
        expected.sort_unstable();
        expected.truncate(10);
        assert_eq!(db.dense_sparse_prefix(0, 10).unwrap(), expected);
        // Asking for more than the genome has yields all of it, sorted.
        assert_eq!(db.dense_sparse_prefix(1, 1000).unwrap(), tiny);

        let mut loosened = open(std::io::Cursor::new(build_with_mode(
            &sketches,
            SmallGenomeMode::Loosen,
        )))
        .unwrap();
        let before = sorted_hits(&loosened.screen(&sample_from(&[(tiny[0], 7)])));
        let kmers_before = loosened.screen_kmers(1);
        loosened
            .enable_small_genome_screen(SPARSE_TARGET_MIN_DEFAULT)
            .unwrap();
        // `big` already cleared the floor, so only the 40-k-mer `tiny` is picked
        // up -- and its reconstructed set is the one already stored.
        let small = loosened.small_screen.as_ref().unwrap();
        assert_eq!(small.num_genomes(), 1);
        assert_eq!(small.num_kmers(), tiny.len());
        assert_eq!(loosened.screen_kmers(1), kmers_before);
        assert_eq!(
            sorted_hits(&loosened.screen(&sample_from(&[(tiny[0], 7)]))),
            before
        );
    }
}
