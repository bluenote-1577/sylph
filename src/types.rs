use std::collections::HashMap;
// bytecheck can be used to validate your data if you want
use std::hash::{BuildHasherDefault, Hasher};
use std::collections::HashSet;
use serde::{Deserialize, Serialize};

pub type Kmer = u64;
pub const BYTE_TO_SEQ: [u8; 256] = [
    0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

#[inline]
pub fn mm_hash(bytes: &[u8]) -> usize {
    let mut key = usize::from_ne_bytes(bytes.try_into().unwrap()) as usize;
    key = !key.wrapping_add(key << 21); // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8); // key * 265
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4); // key * 21
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    return key;
}

#[inline]
pub fn mm_hash64(kmer: u64) -> u64 {
    let mut key = kmer;
    key = !key.wrapping_add(key << 21); // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8); // key * 265
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4); // key * 21
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    return key;
}

pub struct MMHasher {
    hash: usize,
}

impl Hasher for MMHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        self.hash = mm_hash(bytes);
    }
    #[inline]
    fn finish(&self) -> u64 {
        self.hash as u64
    }
}

impl Default for MMHasher {
    #[inline]
    fn default() -> MMHasher {
        MMHasher { hash: 0 }
    }
}

//Implement minimap2 hashing, will test later.
pub type MMBuildHasher = BuildHasherDefault<MMHasher>;
pub type MMHashMap<K, V> = HashMap<K, V, MMBuildHasher>;
pub type MMHashSet<K> = HashSet<K, MMBuildHasher>;

#[derive(Default, Deserialize, Serialize, Debug, PartialEq)]
pub struct SequencesSketch{
    pub kmer_counts: MMHashMap<Kmer, u32>,
    pub c: usize,
    pub k: usize,
    pub file_name: String
}

//Encoding kmer_counts as vec speeds up serialize/deserialize by
//a magnitude. 
#[derive(Default, Deserialize, Serialize, Debug, PartialEq)]
pub struct SequencesSketchEncode{
    pub kmer_counts: Vec<(Kmer, u32)>,
    pub c: usize,
    pub k: usize,
    pub file_name: String
}

impl SequencesSketchEncode{
    pub fn new(sketch: SequencesSketch) -> SequencesSketchEncode{
        let mut vec_map = Vec::with_capacity(sketch.kmer_counts.len());
        for (key,val) in sketch.kmer_counts.into_iter(){
            vec_map.push((key,val));
        }
        return SequencesSketchEncode{kmer_counts: vec_map, file_name: sketch.file_name, c: sketch.c, k: sketch.k};
    }
}

impl SequencesSketch{
    pub fn new(file_name: String, c: usize, k: usize) -> SequencesSketch{
        return SequencesSketch{kmer_counts : HashMap::default(), file_name, c, k}
    }
    pub fn from_enc(sketch: SequencesSketchEncode) -> SequencesSketch{
        let mut new_map = MMHashMap::default();
        new_map.reserve(sketch.kmer_counts.len());
        for item in sketch.kmer_counts.into_iter(){
            new_map.insert(item.0, item.1);
        }
        return SequencesSketch{kmer_counts: new_map, file_name: sketch.file_name, c: sketch.c, k: sketch.k};
    }
}

#[derive(Deserialize, Serialize, Debug, PartialEq)]
#[derive(Default, Clone)]
pub struct GenomeSketch{
    pub genome_kmers: Vec<Kmer>,
    pub file_name: String,
    pub first_contig_name: String,
    pub c: usize,
    pub k: usize,
}

#[derive(Deserialize, Serialize, Debug, PartialEq)]
#[derive(Default, Clone)]
pub struct MultGenomeSketch{
    pub genome_kmer_counts: Vec<Vec<Kmer>>,
    pub file_names: Vec<String>,
    pub contig_names: Vec<String>,
    pub c: usize,
    pub k: usize,
}

#[derive(Deserialize, Serialize, Debug, PartialEq)]
pub struct AniResult<'a>{
    pub ani: f64,
    pub cov_cor_ani: Option<f64>,
    pub seq_name: &'a str,
    pub gn_name: &'a str,
    pub contig_name: &'a str,
    pub mean_cov: f64,
    pub median_cov: f64,
    pub containment_index: (usize,usize),
    pub lambda: Option<f64>,
}
