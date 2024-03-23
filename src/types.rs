//Various byte-tables and hashing methods are taken from miniprot by Heng Li. Attached below is their license:
//The MIT License

// **** miniprot LICENSE ***
//Copyright (c) 2022-     Dana-Farber Cancer Institute
//
//Permission is hereby granted, free of charge, to any person obtaining
//a copy of this software and associated documentation files (the
//"Software"), to deal in the Software without restriction, including
//without limitation the rights to use, copy, modify, merge, publish,
//distribute, sublicense, and/or sell copies of the Software, and to
//permit persons to whom the Software is furnished to do so, subject to
//the following conditions:
//
//The above copyright notice and this permission notice shall be
//included in all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
//EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
//NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
//BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
//ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.
//******************************


use std::collections::HashMap;

// bytecheck can be used to validate your data if you want
use std::hash::{BuildHasherDefault, Hasher};
use std::collections::HashSet;
use smallvec::SmallVec;
use serde::{Deserialize, Serialize};
use fxhash::FxHashMap;

#[derive(Deserialize, Serialize, Debug, PartialEq)]
pub enum AdjustStatus {
    Lambda(f64),
    Low,
    High,
}

impl Default for AdjustStatus {
    fn default() -> Self {AdjustStatus::Low }
}

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
    key = (!key).wrapping_add(key << 21); // key = (key << 21) - key - 1;
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
    pub kmer_counts: FxHashMap<Kmer, u32>,
    pub c: usize,
    pub k: usize,
    pub file_name: String,
    pub sample_name: Option<String>,
    pub paired: bool,
    pub mean_read_length: f64,
}

//Encoding kmer_counts as vec speeds up serialize/deserialize by
//a magnitude. 
#[derive(Default, Deserialize, Serialize, Debug, PartialEq)]
pub struct SequencesSketchEncode{
    pub kmer_counts: Vec<(Kmer, u32)>,
    pub c: usize,
    pub k: usize,
    pub file_name: String,
    pub sample_name: Option<String>,
    pub paired: bool,
    pub mean_read_length: f64,
}

impl SequencesSketchEncode{
    pub fn new(sketch: SequencesSketch) -> SequencesSketchEncode{
        let mut vec_map = Vec::with_capacity(sketch.kmer_counts.len());
        for (key,val) in sketch.kmer_counts.into_iter(){
            vec_map.push((key,val));
        }
        return SequencesSketchEncode{kmer_counts: vec_map, file_name: sketch.file_name, c: sketch.c, k: sketch.k, paired: sketch.paired, mean_read_length: sketch.mean_read_length, sample_name: sketch.sample_name };
    }
}

impl SequencesSketch{
    pub fn new(file_name: String, c: usize, k: usize, paired: bool, sample_name: Option<String>, mean_read_length: f64) -> SequencesSketch{
        return SequencesSketch{kmer_counts : HashMap::default(), file_name, c, k, paired, sample_name, mean_read_length}
    }
    pub fn from_enc(sketch: SequencesSketchEncode) -> SequencesSketch{
        let mut new_map = FxHashMap::default();
        new_map.reserve(sketch.kmer_counts.len());
        for item in sketch.kmer_counts.into_iter(){
            new_map.insert(item.0, item.1);
        }
        return SequencesSketch{kmer_counts: new_map, file_name: sketch.file_name, c: sketch.c, k: sketch.k, paired: sketch.paired, mean_read_length: sketch.mean_read_length, sample_name: sketch.sample_name};
    }
}

#[derive(Deserialize, Serialize, Debug, PartialEq, Hash, PartialOrd, Eq, Ord, Default, Clone)]
pub struct GenomeSketch{
    pub genome_kmers: Vec<Kmer>,
    pub pseudotax_tracked_nonused_kmers: Option<Vec<Kmer>>,
    pub file_name: String,
    pub first_contig_name: String,
    pub c: usize,
    pub k: usize,
    pub gn_size: usize,
    pub min_spacing: usize,
}

#[derive(Deserialize, Serialize, Debug, PartialEq)]
#[derive(Default, Clone)]
pub struct MultGenomeSketch{
    pub genome_kmer_index: Vec<(Kmer,SmallVec<[u32;1]>)>,
    pub file_names: Vec<String>,
    pub contig_names: Vec<String>,
    pub c: usize,
    pub k: usize,
}

#[derive(Debug, PartialEq)]
pub struct AniResult<'a>{
    pub naive_ani: f64,
    pub final_est_ani: f64,
    pub final_est_cov: f64,
    pub seq_name: String,
    pub gn_name: &'a str,
    pub contig_name: &'a str,
    pub mean_cov: f64,
    pub median_cov: f64,
    pub containment_index: (usize,usize),
    pub lambda: AdjustStatus,
    pub ani_ci: (Option<f64>,Option<f64>),
    pub lambda_ci: (Option<f64>,Option<f64>),
    pub genome_sketch: &'a GenomeSketch,
    pub rel_abund: Option<f64>,
    pub seq_abund: Option<f64>,
    pub kmers_lost: Option<usize>,

}
