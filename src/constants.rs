pub const EM_ABUND_CUTOFF: f64 = 0.01;
pub const PAIR_REGEX: &str = r"(.+)(_?1|_?2)(\..+)";
pub const CUTOFF_PVALUE:f64 = 0.9999999999;
pub const SAMPLE_SIZE_CUTOFF: usize = 25;
pub const MEDIAN_ANI_THRESHOLD: f64 = 2.;
pub const QUERY_FILE_SUFFIX: &str = ".syldb";
pub const SAMPLE_FILE_SUFFIX: &str = ".sylsp";
pub const QUERY_FILE_SUFFIX_VALID : [&str;2] = [QUERY_FILE_SUFFIX, ".sylqueries"];
pub const SAMPLE_FILE_SUFFIX_VALID : [&str;2] = [SAMPLE_FILE_SUFFIX, ".sylsample"];
pub const MIN_ANI_DEF: f64 = 0.9;
pub const MIN_ANI_P_DEF: f64 = 0.95;
pub const MAX_MEDIAN_FOR_MEAN_FINAL_EST: f64 = 15.;
pub const DEREP_PROFILE_ANI: f64 = 0.975;
pub const MAX_DEDUP_COUNT: u32 = 4;
pub const MAX_DEDUP_LEN: usize = 10000000;
pub const DEFAULT_FPR: f64 = 0.0001;
/// Base seed for algorithms whose output must be repeatable across runs.
pub const DEFAULT_RNG_SEED: u64 = 7;
pub const MED_KMER_FOR_ID_EST: f64 = 3.;
pub const SCREEN_C_DEFAULT: usize = 3000;
pub const SCREEN_MIN_ANI_DEFAULT: f64 = 85.;
/// Two-stage seekable database: a small bincoded sparse (screen) index plus
/// per-genome Golomb-Rice compressed dense blocks loaded on demand.
pub const TWO_STAGE_DB_SUFFIX: &str = ".syl2db";
/// Below this (compressed, for .gz inputs) file size, the multi-threaded
/// producer/consumer read-sketching pipeline isn't worth its setup cost;
/// fall back to the legacy single-threaded-per-file path regardless of the
/// available thread budget.
pub const MIN_BYTES_FOR_SKETCH_PIPELINE: u64 = 4_000_000;
/// Default reads per batch handed from the sketching pipeline's I/O thread to
/// worker threads.
pub const DEFAULT_SKETCH_BATCH_SIZE: usize = 2000;
/// Default in-flight batches buffered per worker thread between the I/O
/// thread and worker threads in the sketching pipeline.
pub const DEFAULT_SKETCH_CHANNEL_DEPTH: usize = 2;
/// Default maximum accumulated sequence bytes per sketching-pipeline batch --
/// a safety cap alongside DEFAULT_SKETCH_BATCH_SIZE, whichever is hit first.
/// A no-op for typical short-read data (2000 x ~150bp is far below this);
/// bounds long-read (Nanopore/PacBio) batches, which a fixed record count
/// alone could otherwise let grow to many tens of megabytes.
pub const DEFAULT_SKETCH_BATCH_MAX_BYTES: usize = 4_000_000;
/// Default minimum stage-1 sparse/screen k-mers per genome in a `.syl2db`
/// (see `write_two_stage_db`'s adaptive floor).
pub const SPARSE_TARGET_MIN_DEFAULT: usize = 50;
/// Bits per key in the band-1 bloom pre-filter of a value-banded screen index
/// (`--small-genome-screen band`). Two hash functions, so ~2% false-positive
/// rate: enough that sample k-mers in band 1's value range almost never reach
/// the (cache-missing) binary search.
pub const BAND_BLOOM_BITS_PER_KEY: usize = 8;
