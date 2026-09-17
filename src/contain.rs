use crate::cmdline::*;
use std::path::Path;
use std::io::prelude::*;
use std::io;
use std::io::BufWriter;
use fxhash::FxHashMap;
use crate::constants::*;
use crate::inference::*;
use crate::sketch::*;
use crate::twostage_db::TwoStageDb;
use crate::types::*;
use log::*;
use rayon::prelude::*;
use statrs::distribution::{DiscreteCDF, Poisson};
use std::fs::File;
use std::io::BufReader;
use std::sync::Mutex;

fn print_ani_result(ani_result: &AniResult, pseudotax: bool, writer: &mut Box<dyn Write + Send>) {
    let print_final_ani = format!("{:.2}", f64::min(ani_result.final_est_ani * 100., 100.));
    let lambda_print;
    if let AdjustStatus::Lambda(lambda) = ani_result.lambda {
        lambda_print = format!("{:.3}", lambda);
    } else if ani_result.lambda == AdjustStatus::High {
        lambda_print = format!("HIGH");
    } else {
        lambda_print = format!("LOW");
    }
    let low_ani = ani_result.ani_ci.0;
    let high_ani = ani_result.ani_ci.1;
    let low_lambda = ani_result.lambda_ci.0;
    let high_lambda = ani_result.lambda_ci.1;

    let ci_ani;
    if low_ani.is_none() || high_ani.is_none() {
        ci_ani = "NA-NA".to_string();
    } else {
        ci_ani = format!(
            "{:.2}-{:.2}",
            low_ani.unwrap() * 100.,
            high_ani.unwrap() * 100.
        );
    }

    let ci_lambda;
    if low_lambda.is_none() || high_lambda.is_none() {
        ci_lambda = "NA-NA".to_string();
    } else {
        ci_lambda = format!("{:.2}-{:.2}", low_lambda.unwrap(), high_lambda.unwrap());
    }


    //"Sample_file\tQuery_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\tEff_cov\tANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tNaive_ANI\tContig_name",

    if !pseudotax{
        writeln!(writer, 
            "{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{:.0}\t{:.3}\t{}/{}\t{:.2}\t{}",
            ani_result.seq_name,
            ani_result.gn_name,
            print_final_ani,
            ani_result.final_est_cov,
            ci_ani,
            lambda_print,
            ci_lambda,
            ani_result.median_cov,
            ani_result.mean_cov,
            ani_result.containment_index.0,
            ani_result.containment_index.1,
            ani_result.naive_ani * 100.,
            ani_result.contig_name,
        ).expect("Error writing to file");
    }
    else{
        writeln!(writer,
            "{}\t{}\t{:.4}\t{:.4}\t{}\t{:.3}\t{}\t{}\t{}\t{:.0}\t{:.3}\t{}/{}\t{:.2}\t{}\t{}",
            ani_result.seq_name,
            ani_result.gn_name,
            ani_result.rel_abund.unwrap(),
            ani_result.seq_abund.unwrap(),
            print_final_ani,
            ani_result.final_est_cov,
            ci_ani,
            lambda_print,
            ci_lambda,
            ani_result.median_cov,
            ani_result.mean_cov,
            ani_result.containment_index.0,
            ani_result.containment_index.1,
            ani_result.naive_ani * 100.,
            ani_result.kmers_lost.unwrap(),
            ani_result.contig_name,
        ).expect("Error writing to file");

    }
}

fn get_chunks(indices: &Vec<usize>, steps: usize) -> Vec<Vec<usize>>{
    let mut start = 0;
    let mut end = steps;
    let len = indices.len();
    let mut return_chunks = vec![];

    while start < len {
        if end > len {
            end = len;
        }

        let chunk: Vec<usize> = (start..end).collect();
        start = end;
        end += steps;
        return_chunks.push(chunk);
    }
    return_chunks
}

fn minimum_ani_fraction(args: &ContainArgs) -> f64 {
    if let Some(minimum_ani) = args.minimum_ani {
        minimum_ani / 100.0
    } else if args.pseudotax {
        MIN_ANI_P_DEF
    } else {
        MIN_ANI_DEF
    }
}

/// One `--screen-dump` row: genome id, source file, contig name, matched and
/// total screen k-mers, naive and adjusted ANI, median coverage.
type ScreenDumpRow = (usize, String, String, usize, usize, f64, f64, f64);

/// Two-stage stage 1 + 2 against one `.syl2db`: screen `sequence_sketch`
/// against the database's pooled `screen_index` (a single inverted pass over
/// the sample at the sparse, per-database `screen_c`), then decode and return
/// *dense* sketches for only the genomes that pass. The returned set is one
/// database's contribution to the combined active-genome list the (expensive)
/// dense profiling pass runs against.
fn compute_dense_survivors(
    args: &ContainArgs,
    db: &TwoStageDb,
    sequence_sketch: &SequencesSketch,
) -> Vec<GenomeSketch> {
    let start = std::time::Instant::now();
    let screen_index = &db.screen_index;
    // Stage 1: cheap, permissive screen (query-like settings, no CIs).
    let mut screen_args = args.clone();
    screen_args.pseudotax = false;
    screen_args.minimum_ani = Some(args.screen_ani);
    screen_args.no_ci = true;
    // The screen scores against sketches sub-sampled to `db.screen_c`, so the
    // dense `-M/--min-number-kmers` floor would over-reject here: a genome has
    // only ~length/screen_c sparse k-mers, so the dense floor M corresponds to a
    // length of M*screen_c, silently dropping smaller genomes (viruses,
    // plasmids, short contigs) that single-stage profiling would report. Scale
    // the floor to the screen resolution so a genome that could clear the dense
    // floor also clears the screen; genomes truly below the floor are still
    // rejected at the dense stage. `effective_screen_c` is the *densest* rate
    // anything in this database is screened at (see `write_two_stage_db`'s floor
    // and `--screen-small-genomes`), so this scaling automatically becomes less
    // aggressive when a small genome forced a denser rate. `--min-contain` is
    // deliberately deferred to the decoded dense sketch below: applying an
    // absolute dense evidence floor to this much sparser screen could discard a
    // valid dense candidate.
    let effective_screen_c = db.effective_screen_c();
    screen_args.min_number_kmers = args.min_number_kmers * db.c as f64 / effective_screen_c as f64;

    // Stage 1 = "Path B": one inverted pass over the sample produces, per genome,
    // the same matched-coverage multiset the per-genome `get_stats` loop would
    // collect; feeding it to the same `finalize_stats` + `min_number_kmers`
    // checks reproduces the survivor set exactly, in O(sample) rather than
    // O(reference) work. (See experiments/7_mphf_screen_again.)
    let hits: Vec<(u32, Vec<u32>)> = db.screen(sequence_sketch).into_iter().collect();
    let gather_s = start.elapsed().as_secs_f64();
    // Optional per-survivor dump. The genome id and contig name are what make a
    // row identifiable: a database built from multi-record FASTA (a viral vOTU
    // set, say) has one source file for millions of genomes.
    let dump: Mutex<Vec<ScreenDumpRow>> = Mutex::new(vec![]);
    let mut survivors: Vec<usize> = hits
        .into_par_iter()
        .filter_map(|(g, covs)| {
            let g = g as usize;
            let n_kmers = db.screen_kmers(g as u32) as usize;
            // Mirror get_stats: reject genomes below the (scaled) k-mer floor.
            if (n_kmers as f64) < screen_args.min_number_kmers {
                return None;
            }
            let contain_count = covs.len();
            // finalize_stats applies the screen min-ANI gate (returns None below it).
            let fin = finalize_stats(&screen_args, db.k, n_kmers, contain_count, covs, None)?;
            if args.screen_dump.is_some() {
                dump.lock().unwrap().push((
                    g,
                    db.genome_file_name(g as u32).to_string(),
                    db.genome_first_contig_name(g as u32).to_string(),
                    contain_count,
                    n_kmers,
                    fin.naive_ani * 100.,
                    fin.final_est_ani * 100.,
                    fin.median_cov,
                ));
            }
            Some(g)
        })
        .collect();
    survivors.sort_unstable();
    let stage1_s = start.elapsed().as_secs_f64();
    log::info!(
        "{}: stage-1 screen (c={}, min-ANI {}) kept {} / {} candidate genomes",
        sequence_sketch.file_name, effective_screen_c, args.screen_ani, survivors.len(), screen_index.num_genomes()
    );
    if let Some(path) = &args.screen_dump {
        let mut f = BufWriter::new(File::create(path).expect("could not create --screen-dump file"));
        writeln!(f, "Genome_id\tGenome_file\tContig_name\tscreen_matched_kmers\tscreen_total_kmers\tnaive_ani\tscreen_adjusted_ani\tscreen_median_cov").unwrap();
        // Sorted, so two runs' dumps are directly comparable: the survivor set is
        // the stage-1 output, and the parallel filter above collects it in
        // whatever order threads finish.
        let mut rows = dump.into_inner().unwrap();
        rows.sort_unstable_by_key(|r| r.0);
        for (g, file, contig, m, t, na, ea, mc) in rows {
            writeln!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{}",
                g, file, contig, m, t, na, ea, mc
            )
            .unwrap();
        }
        log::info!("Wrote stage-1 screen dump to {}", path);
    }

    // Stage 2: decode each survivor's dense block into a short-lived sketch,
    // run the (prefetch-friendly) pass-1 profiling on it, and keep it only if
    // it passes -- so the discarded majority is freed immediately and never
    // cached, keeping peak RAM proportional to the genomes that survive rather
    // than to the (much larger) screen-survivor set.
    let dense: Mutex<Vec<GenomeSketch>> = Mutex::new(vec![]);
    survivors.par_iter().for_each(|i| {
        match db.decode_dense(*i as u32) {
            Ok(g) => {
                // Keep only genomes that pass pass-1 profiling; drop (free) the rest.
                if get_stats(args, &g, sequence_sketch, None, false).is_some() {
                    dense.lock().unwrap().push(g);
                }
            }
            Err(e) => warn!("Could not decode dense block for genome index {}: {}", i, e),
        }
    });
    let dense = dense.into_inner().unwrap();
    log::info!(
        "{}: stage-2 dense profiling (c={}) against {} genomes",
        sequence_sketch.file_name, db.c, dense.len()
    );
    // One machine-readable line per (sample, database) for benchmark harnesses
    // (debug-level: the format is for parsers, not for people). `gather` is the
    // inverted pass over the sample plus (if enabled) the reader-side
    // small-genome intersection; `stage1` adds the per-survivor ANI
    // finalisation; `stage2` is the dense decode + pass-1 profiling.
    log::debug!(
        "bench_screen sample={} db={} sample_kmers={} screen_c={} effective_screen_c={} \
         gather_s={:.4} stage1_s={:.4} survivors={} stage2_s={:.4} dense_kept={} genomes={}",
        sequence_sketch.file_name,
        db.path,
        sequence_sketch.kmer_counts.len(),
        db.screen_c,
        effective_screen_c,
        gather_s,
        stage1_s,
        survivors.len(),
        start.elapsed().as_secs_f64() - stage1_s,
        dense.len(),
        db.len(),
    );
    dense
}

pub fn contain(mut args: ContainArgs, pseudotax_in: bool) {

    if pseudotax_in{
        args.pseudotax = true;
    }

    let level;
    if args.trace {
        level = log::LevelFilter::Trace;
    } else if args.debug {
        level = log::LevelFilter::Debug;
    }
    else{
        level = log::LevelFilter::Info;
    }
    
    simple_logger::SimpleLogger::new()
        .with_level(level)
        .init()
        .unwrap();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    let out_writer = match args.out_file_name {
        Some(ref x) => {
            let path = Path::new(&x);
            Box::new(BufWriter::new(File::create(&path).unwrap())) as Box<dyn Write + Send>
        }
        None => Box::new(BufWriter::new(io::stdout())) as Box<dyn Write + Send>,
    };

    if args.estimate_read_counts{
        args.estimate_unknown = true;
        log::info!("--estimate-read-counts detected, also enabling -u. Sequence_abundance column will be set to estimated read counts, not abundance. This is still experimental.");
    }

    log::info!("Obtaining sketches...");
    let mut genome_sketch_files = vec![];
    let mut genome_files = vec![];
    let mut two_stage_db_files = vec![];
    let mut read_sketch_files = vec![];
    let mut read_files = vec![];

    let mut all_files = args.files.clone();
    all_files.extend(args.databases.clone());

    if let Some(ref newline_file) = args.file_list{
        let file = File::open(newline_file).unwrap();
        let reader = BufReader::new(file);
        for line in reader.lines() {
            all_files.push(line.unwrap());
        }

    }


    for file in all_files.iter(){

        let mut genome_sketch_good_suffix = false;
        for suff in QUERY_FILE_SUFFIX_VALID{
            if file.ends_with(suff){
                genome_sketch_good_suffix = true;
                break
            }
        }

        let mut sample_sketch_good_suffix = false;
        for suff in SAMPLE_FILE_SUFFIX_VALID{
            if file.ends_with(suff){
                sample_sketch_good_suffix = true;
                break
            }
        }

        if file.ends_with(TWO_STAGE_DB_SUFFIX){
            two_stage_db_files.push(file);
        } else if genome_sketch_good_suffix{
            genome_sketch_files.push(file);
        } else if sample_sketch_good_suffix{
            read_sketch_files.push(file);
        } else if is_fasta(&file) {
            genome_files.push(file);
        } else if is_fastq(&file) {
            read_files.push(vec![file]);
        } else {
            warn!(
                "{} file extension is not a sketch or a fasta/fastq file.",
                &file
            );
        }
    }

    if args.first_pair.len() != args.second_pair.len() {
        error!("Different number of paired sequences (-1, -2) for sketching. Exiting.");
        std::process::exit(1);
    }

    // zip together the first and second pair files, push them to read_files
    for (first, second) in args.first_pair.iter().zip(args.second_pair.iter()) {
        read_files.push(vec![first,second]);
    }

    for read in args.reads.iter() {
        read_files.push(vec![read]);
    }

    if genome_sketch_files.is_empty() && genome_files.is_empty() && two_stage_db_files.is_empty(){
        log::error!("No genome files found; see sylph query/profile -h for help. Exiting");
        std::process::exit(1);
    }

    if read_sketch_files.is_empty() && read_files.is_empty(){
        log::error!("No read files found; see sylph query/profile -h for help. Exiting");
        std::process::exit(1);
    }

    // Any number of .syl2db files can be combined with each other and with
    // plain .syldb/fasta genomes in one call: each .syl2db is screened and
    // decoded independently per sample (see the per-sample loop below), and
    // its survivors are unioned with the (always fully-loaded, unscreened)
    // plain genomes before profiling/reassignment.
    let two_stage_dbs: Vec<TwoStageDb> = two_stage_db_files.iter().map(|f| {
        log::info!("Opening two-stage database {} (loading stage-1 sparse index)...", f);
        let open_start = std::time::Instant::now();
        let mut db = crate::twostage_db::open_file(f)
            .unwrap_or_else(|e| panic!("{} is not a valid two-stage database: {}", f, e));
        let index_s = open_start.elapsed().as_secs_f64();
        // Reader-side small-genome screen: only touches genomes the *stored*
        // index left short, and reconstructs exactly the set a densified build
        // would have stored, so it changes nothing on a database already built
        // with --small-genome-screen loosen/band and composes safely with either.
        if args.screen_small_genomes > 0 {
            db.enable_small_genome_screen(args.screen_small_genomes)
                .unwrap_or_else(|e| panic!("Could not read back dense blocks of small genomes in {}: {}", f, e));
            let s = db.small_screen.as_ref().unwrap();
            log::info!(
                "{}: screening {} genome(s) with < {} stored stage-1 k-mers directly, from {} \
                 k-mers read back out of their dense blocks",
                f, s.num_genomes(), args.screen_small_genomes, s.num_kmers()
            );
        }
        let small = db.small_screen.as_ref();
        log::debug!(
            "bench_open db={} genomes={} screen_c={} effective_screen_c={} band_keys={} \
             band_owners={} small_genomes={} small_kmers={} index_s={:.4} open_s={:.4}",
            f,
            db.len(),
            db.screen_c,
            db.effective_screen_c(),
            db.screen_index.band.as_ref().map_or(0, |b| b.num_keys()),
            db.screen_index.band.as_ref().map_or(0, |b| b.num_owners()),
            small.map_or(0, |s| s.num_genomes()),
            small.map_or(0, |s| s.num_kmers()),
            index_s,
            open_start.elapsed().as_secs_f64(),
        );
        db
    }).collect();
    for db in &two_stage_dbs {
        if db.is_empty() {
            log::error!("Two-stage database contains no genomes. Exiting");
            std::process::exit(1);
        }
    }

    let genome_sketches = get_genome_sketches(&args, &genome_sketch_files, &genome_files);
    log::info!("Finished obtaining genome sketches.");

    if two_stage_dbs.is_empty() && genome_sketches.is_empty() {
        log::error!("No genome sketches found; see sylph query/profile -h for help. Exiting");
        std::process::exit(1);
    }
    if !genome_sketches.is_empty()
        && genome_sketches.first().unwrap().pseudotax_tracked_nonused_kmers.is_none() && args.pseudotax {
        log::error!("Attempting profiling, but *.syldb was sketched with the --disable-profiling option. Exiting");
        std::process::exit(1);
    }
    // (No equivalent pseudotax-completeness check is needed for two_stage_dbs:
    // db-convert already requires 100% pseudotax coverage at conversion time.)

    // Duplicate genomes loaded across multiple sources (e.g. the same genome
    // present in both a .syl2db and a plain .syldb) collide as identical
    // GenomeSketch entries downstream in reassignment/dereplication, causing
    // nondeterministic output (whichever copy a racing thread inserts first
    // "wins"). Warn loudly rather than silently producing inconsistent runs.
    {
        let mut seen: FxHashMap<(&str, &str), usize> = FxHashMap::default();
        for db in &two_stage_dbs {
            for g in 0..db.len() {
                *seen.entry((db.genome_file_name(g as u32), db.genome_first_contig_name(g as u32))).or_insert(0) += 1;
            }
        }
        for gs in &genome_sketches {
            *seen.entry((gs.file_name.as_str(), gs.first_contig_name.as_str())).or_insert(0) += 1;
        }
        let dups: Vec<&(&str, &str)> = seen.iter().filter(|(_, &n)| n > 1).map(|(k, _)| k).collect();
        if !dups.is_empty() {
            log::warn!(
                "{} genome(s) appear more than once across the combined input sources \
                 (same file name + contig name loaded from two different databases): {:?}. \
                 This can cause nondeterministic reassignment/dereplication output -- consider \
                 removing the duplicate(s) from one of the input databases.",
                dups.len(), dups
            );
        }
    }

    // k must match exactly across every loaded source.
    let mut ks: Vec<usize> = two_stage_dbs.iter().map(|db| db.k).collect();
    if !genome_sketches.is_empty() { ks.push(genome_sketches[0].k); }
    ks.dedup();
    if ks.len() > 1 {
        log::error!("Inconsistent -k across combined genome sources ({:?}). Exiting", ks);
        std::process::exit(1);
    }
    let db_k = ks[0];

    // The sample sketch must not be sparser than any active genome's own rate
    // (get_stats hard-exits otherwise): use the minimum -c across every loaded
    // source (two-stage dense rate and plain -c alike).
    let effective_genome_c = two_stage_dbs.iter().map(|db| db.c)
        .chain(genome_sketches.iter().map(|g| g.c))
        .min()
        .expect("at least one genome source must be present (checked above)");

    log::info!(
        "Loaded {} two-stage database(s) ({} genomes total) and {} plain genome sketch(es).",
        two_stage_dbs.len(),
        two_stage_dbs.iter().map(|d| d.len()).sum::<usize>(),
        genome_sketches.len(),
    );

    let num_raw_read_files = read_files.len();
    let step;
    if let Some(sample_threads) = args.sample_threads{
        if sample_threads > 0{
            step = sample_threads;
        }
        else{
            step = 1;
        }
    }
    else{
        // Historically `profile` (pseudotax) reserved fewer samples-in-flight
        // than `query` via the `threads/3 + 1` floor, leaving headroom for the
        // per-sample reassignment pass. That floor was inert whenever
        // `num_raw_read_files >= threads/3 + 1` (the common case), since with
        // single-threaded-per-file sketching `step` only mattered via
        // `min(step, num_raw_read_files)` in `get_chunks` below -- but it
        // becomes an active bottleneck now that `step` also sizes
        // `threads_per_file` for the multi-threaded sketch pipeline: with few
        // files and many threads, inflating `step` starves each file's
        // pipeline of threads for no benefit (nested rayon parallelism in the
        // reassignment pass doesn't need reserved threads; it work-steals from
        // the same global pool). Match `query`'s formula.
        step = usize::max(1, usize::min(num_raw_read_files, args.threads))
    }

    // Threads available per file's own sketching pipeline, given `step` files
    // are processed concurrently out of the shared `args.threads` budget.
    let threads_per_file = (args.threads / step.max(1)).max(1);

    let read_sketch_files_as_vec = read_sketch_files.clone().into_iter().map(|x| vec![x]).collect::<Vec<Vec<&String>>>();
    read_files.extend(read_sketch_files_as_vec);
    let sequence_index_vec = (0..read_files.len()).collect::<Vec<usize>>();
    let out_writer:Mutex<Box<dyn Write + Send>> = Mutex::new(out_writer);

    let chunks = get_chunks(&sequence_index_vec, step);

    print_header(args.pseudotax,&mut *out_writer.lock().unwrap(), args.estimate_unknown);
    chunks.into_iter().for_each(|chunk| {
        chunk.into_par_iter().for_each(|j|{
            let is_sketch = j >= read_files.len() - read_sketch_files.len();
            let sequence_sketch = get_seq_sketch(&args, &read_files[j], is_sketch, effective_genome_c, db_k, threads_per_file);
            if sequence_sketch.is_some(){
                let first_read_file = read_files[j][0];
                let sequence_sketch = sequence_sketch.unwrap();

                let kmer_id_opt;
                if args.seq_id.is_some(){
                    kmer_id_opt = Some((args.seq_id.unwrap()/100.).powf(sequence_sketch.k as f64));
                }
                else{
                    kmer_id_opt = get_kmer_identity(&sequence_sketch, args.estimate_unknown);
                    // None unless -u/--estimate-unknown was given, in which case
                    // there is no identity to report. Unwrapping it here aborted
                    // every `--debug` run without `-u`.
                    if let Some(kmer_id) = kmer_id_opt{
                        log::debug!("{} has estimated identity {:.3}.", &first_read_file, kmer_id.powf(1./sequence_sketch.k as f64) * 100.);
                    }
                }

                // Screen each two-stage database, then union its survivors with
                // the (always fully-loaded, unscreened) plain genomes into one
                // combined active-genome list for this sample.
                let mut dense_local: Vec<GenomeSketch> = Vec::new();
                for db in &two_stage_dbs {
                    dense_local.extend(compute_dense_survivors(&args, db, &sequence_sketch));
                }
                let active_sketches: Vec<&GenomeSketch> =
                    dense_local.iter().chain(genome_sketches.iter()).collect();

                let stats_vec_seq: Mutex<Vec<AniResult>> = Mutex::new(vec![]);
                active_sketches.par_iter().for_each(|genome_sketch| {
                    let res = get_stats(&args, genome_sketch, &sequence_sketch, None, args.log_reassignments);
                    if let Some(res) = res {
                        stats_vec_seq.lock().unwrap().push(res);
                    }
                });

                let mut stats_vec_seq = stats_vec_seq.into_inner().unwrap();
                estimate_true_cov(&mut stats_vec_seq, kmer_id_opt, args.estimate_unknown, sequence_sketch.mean_read_length, sequence_sketch.k);

                if args.pseudotax{
                    log::info!("{} taxonomic profiling; reassigning k-mers for {} genomes...", &first_read_file, stats_vec_seq.len());
                    let winner_map = winner_table(&stats_vec_seq, args.log_reassignments);
                    let remaining_genomes = stats_vec_seq.iter().map(|x| x.genome_sketch).collect::<Vec<&GenomeSketch>>();
                    let stats_vec_seq_2 = Mutex::new(vec![]);
                    remaining_genomes.into_par_iter().for_each(|genome_sketch|{
                        let res = get_stats(&args, &genome_sketch, &sequence_sketch, Some(&winner_map), args.log_reassignments);
                        if res.is_some() {
                            stats_vec_seq_2.lock().unwrap().push(res.unwrap());
                        }
                    });
                    stats_vec_seq = derep_if_reassign_threshold(&stats_vec_seq, stats_vec_seq_2.into_inner().unwrap(), args.redundant_ani, sequence_sketch.k);
                    //stats_vec_seq = stats_vec_seq_2.into_inner().unwrap();
                    estimate_true_cov(&mut stats_vec_seq, kmer_id_opt, args.estimate_unknown, sequence_sketch.mean_read_length, sequence_sketch.k);
                    log::info!("{} has {} genomes passing profiling threshold. ", &first_read_file, stats_vec_seq.len());

                    let mut bases_explained = 1.;
                    if args.estimate_unknown{
                        bases_explained = estimate_covered_bases(&stats_vec_seq, &sequence_sketch, sequence_sketch.mean_read_length, sequence_sketch.k);
                        log::info!("{} has {:.2}% of reads detected in database by profile", &first_read_file, bases_explained * 100.);
                        
                    }

                    let total_cov = stats_vec_seq.iter().map(|x| x.final_est_cov).sum::<f64>();
                    let total_seq_cov = stats_vec_seq.iter().map(|x| x.final_est_cov * x.genome_sketch.gn_size as f64).sum::<f64>();
                    for thing in stats_vec_seq.iter_mut(){
                        thing.rel_abund = Some(thing.final_est_cov/total_cov * 100.);
                    }
                    for thing in stats_vec_seq.iter_mut(){
                        if args.estimate_read_counts{
                            thing.seq_abund = Some((thing.final_est_cov * thing.genome_sketch.gn_size as f64 / sequence_sketch.mean_read_length * bases_explained).round());
                        }
                        else{
                            let seq_abund = thing.final_est_cov * thing.genome_sketch.gn_size as f64 / total_seq_cov * 100. * bases_explained;
                            thing.seq_abund = Some(seq_abund);
                        }
                    }
                }

                if args.pseudotax{
                    stats_vec_seq.sort_by(|x,y| y.rel_abund.unwrap().partial_cmp(&x.rel_abund.unwrap()).unwrap());
                }
                else{
                    stats_vec_seq.sort_by(|x,y| y.final_est_ani.partial_cmp(&x.final_est_ani).unwrap());
                }
                
                let mut out_writer = out_writer.lock().unwrap();
                for res in stats_vec_seq{
                    print_ani_result(&res, args.pseudotax, &mut *out_writer);
                }
            }
            if read_files[j].len() > 1{
                log::info!("Finished paired sample {}.", &read_files[j][0]);
            }
            else{
                log::info!("Finished sample {}.", &read_files[j][0]);
            }
        });
    });

    log::info!("sylph finished.");
}

fn derep_if_reassign_threshold<'a>(results_old: &Vec<AniResult>, results_new: Vec<AniResult<'a>>, ani_thresh: f64, k: usize) -> Vec<AniResult<'a>>{
    let ani_thresh = ani_thresh/100.;

    let mut gn_sketch_to_contain = FxHashMap::default();
    for result in results_old.iter(){
        gn_sketch_to_contain.insert(result.genome_sketch, result);
    }

    let threshold = f64::powf(ani_thresh, k as f64);
    let mut return_vec = vec![];
    for result in results_new.into_iter(){
        let old_res = &gn_sketch_to_contain[result.genome_sketch];
        let num_kmer_reassign = (old_res.containment_index.0 - result.containment_index.0) as f64;
        let reass_thresh = threshold * result.containment_index.1 as f64;
        if num_kmer_reassign < reass_thresh{
            return_vec.push(result);
        }
        else{
            log::debug!("genome {} had num k-mers reassigned = {}, threshold was {}, removing.", result.gn_name, num_kmer_reassign, reass_thresh);
        }
    }
    return return_vec;
}

fn estimate_true_cov(results: &mut Vec<AniResult>, kmer_id_opt: Option<f64>, 
                     estimate_unknown: bool, read_length: f64, k: usize){
    let mut multiplier = 1.;
    if estimate_unknown{
        multiplier = read_length / (read_length - k as f64 + 1.);
    }
    if estimate_unknown && kmer_id_opt.is_some(){
        let id = kmer_id_opt.unwrap();
        for res in results.iter_mut(){
            res.final_est_cov = res.final_est_cov / id * multiplier ;
        }
    }
}

fn estimate_covered_bases(results: &Vec<AniResult>, sequence_sketch: &SequencesSketch, read_length: f64, k: usize) -> f64{
    let multiplier = read_length / (read_length - (k as f64) + 1.);

    let mut num_covered_bases = 0.;
    for res in results.iter(){
        num_covered_bases += (res.genome_sketch.gn_size as f64) * res.final_est_cov
    }
    let mut num_total_counts = 0;
    for count in sequence_sketch.kmer_counts.values(){
        num_total_counts += *count as usize;
    }
    let num_tentative_bases = sequence_sketch.c * num_total_counts;
    let num_tentative_bases = num_tentative_bases as f64 * multiplier;
    if num_tentative_bases == 0.{
        return 0.;
    }
    return f64::min(num_covered_bases as f64 / num_tentative_bases, 1.);
}

fn winner_table<'a>(results : &'a Vec<AniResult>, log_reassign: bool) -> FxHashMap<Kmer, (f64,&'a GenomeSketch, bool)> {
    let mut kmer_to_genome_map : FxHashMap<_,_> = FxHashMap::default();
    for res in results.iter(){
        //let gn_sketch = &genome_sketches[res.genome_sketch_index];
        let gn_sketch = res.genome_sketch;
        for kmer in gn_sketch.genome_kmers.iter(){
            let v = kmer_to_genome_map.entry(*kmer).or_insert((res.final_est_ani, res.genome_sketch, false));
            if res.final_est_ani > v.0{
                *v = (res.final_est_ani, gn_sketch, true);
            }
        }
        
        if gn_sketch.pseudotax_tracked_nonused_kmers.is_some(){
            for kmer in gn_sketch.pseudotax_tracked_nonused_kmers.as_ref().unwrap().iter(){
                let v = kmer_to_genome_map.entry(*kmer).or_insert((res.final_est_ani, res.genome_sketch, false));
                if res.final_est_ani > v.0{
                    *v = (res.final_est_ani, gn_sketch, true);
                }
            }
        }
    }

    //log reassigned kmers
    if log_reassign{
        log::info!("------------- Logging k-mer reassignments -----------------");
        let mut sketch_to_index = FxHashMap::default();
        for (i,res) in results.iter().enumerate(){
            log::info!("Index\t{}\t{}\t{}", i, res.genome_sketch.file_name, res.genome_sketch.first_contig_name);
            sketch_to_index.insert(res.genome_sketch, i);
        }
        (0..results.len()).into_par_iter().for_each(|i|{
            let res = &results[i];
            let mut reassign_edge_map = FxHashMap::default();
            for kmer in res.genome_sketch.genome_kmers.iter(){
                let value = kmer_to_genome_map[kmer].1;
                if value != res.genome_sketch{
                    let edge_count = reassign_edge_map.entry((sketch_to_index[value],i)).or_insert(0);
                    *edge_count += 1;
                }
            }
            for (key,val) in reassign_edge_map{
                if val > 10{
                    log::info!("{}->{}\t{}\tkmers reassigned", key.0, key.1, val);
                }
            }
        });
    }

    return kmer_to_genome_map;
}

fn print_header(pseudotax: bool, writer: &mut Box<dyn Write + Send>, estimate_unknown: bool) {
    if !pseudotax{
        writeln!(writer,
            //"Sample_file\tQuery_file\tAdjusted_ANI\tNaive_ANI\tANI_5-95_percentile\tEff_cov\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tContig_name",
            "Sample_file\tGenome_file\tAdjusted_ANI\tEff_cov\tANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tNaive_ANI\tContig_name",
            ).expect("Error writing to file.");
    }
    else{
        let cov_head;
        if estimate_unknown{
            cov_head = "True_cov";
        }
        else{
            cov_head = "Eff_cov";
        }
        writeln!(writer,
            "Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\t{}\tANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tNaive_ANI\tkmers_reassigned\tContig_name", cov_head
            ).expect("Error writing to file.");
    }
}

fn get_genome_sketches(
    args: &ContainArgs,
    genome_sketch_files: &Vec<&String>,
    genome_files: &Vec<&String>,
) -> Vec<GenomeSketch> {
    let mut lowest_genome_c = None;
    let mut current_k = None;

    let genome_sketches = Mutex::new(vec![]);

    for genome_sketch_file in genome_sketch_files {
        let file = File::open(genome_sketch_file).expect(&format!("The sketch `{}` could not be opened. Exiting", genome_sketch_file));
        let genome_reader = BufReader::with_capacity(10_000_000, file);
        let genome_sketches_vec: Vec<GenomeSketch> = bincode::deserialize_from(genome_reader)
            .expect(&format!(
                "The sketch `{}` is not a valid sketch. Perhaps it is an older, incompatible version ",
                &genome_sketch_file
            ));
        if genome_sketches_vec.is_empty() {
            continue;
        }
        let c = genome_sketches_vec.first().unwrap().c;
        let k = genome_sketches_vec.first().unwrap().k;
        if lowest_genome_c.is_none() {
            lowest_genome_c = Some(c);
        } else if lowest_genome_c.unwrap() < c {
            lowest_genome_c = Some(c);
        }
        if current_k.is_none() {
            current_k = Some(genome_sketches_vec.first().unwrap().k);
        } else if current_k.unwrap() != k {
            error!("Query sketches have inconsistent -k. Exiting.");
            std::process::exit(1);
        }
        genome_sketches.lock().unwrap().extend(genome_sketches_vec);
    }

    genome_files.into_par_iter().for_each(|genome_file|{
        if lowest_genome_c.is_some() && lowest_genome_c.unwrap() < args.c{
            error!("Value of -c for contain is {} -- greater than the smallest value of -c for a genome sketch {}. Continuing without sketching.", args.c, lowest_genome_c.unwrap());
        }
        else if current_k.is_some() && current_k.unwrap() != args.k{
            error!("-k {} is not equal to -k {} found in sketches. Continuing without sketching.", args.k, current_k.unwrap());
        }
        else {
            if args.individual{
            let indiv_gn_sketches = sketch_genome_individual(args.c, args.k, genome_file, args.min_spacing_kmer, args.pseudotax);
                genome_sketches.lock().unwrap().extend(indiv_gn_sketches);

            }
            else{
                let genome_sketch_opt = sketch_genome(args.c, args.k, &genome_file, args.min_spacing_kmer, args.pseudotax);
                if genome_sketch_opt.is_some() {
                    genome_sketches.lock().unwrap().push(genome_sketch_opt.unwrap());
                }
            }
        }
    });

    return genome_sketches.into_inner().unwrap();
}

fn get_seq_sketch(
    args: &ContainArgs,
    read_file: &Vec<&String>,
    is_sketch_file: bool,
    genome_c: usize,
    genome_k: usize,
    threads: usize,
) -> Option<SequencesSketch> {
    if is_sketch_file {
        let read_file = read_file[0];
        let read_sketch_file = read_file;
        let file = File::open(read_sketch_file).expect(&format!(
            "The sketch `{}` could not be opened",
            &read_sketch_file
        ));
        let read_reader = BufReader::with_capacity(10_000_000, file);
        let read_sketch: SequencesSketch = bincode::deserialize_from(read_reader).expect(
            &format!("The sketch `{}` is not a valid sketch. Perhaps it is an older incompatible version ", read_sketch_file),
        );
        if read_sketch.c > genome_c {
            error!("{} is a pre-sketched sample at -c {}; this is sparser than required by the loaded genome sketch(es) (-c {}). A pre-sketched sample can't be re-sketched denser after the fact -- re-sketch it at -c <= {} first, or use a raw fastq input instead. Continuing without sketching.", read_file, read_sketch.c, genome_c, genome_c);
            return None;
        }
        else if read_sketch.c < genome_c{
            info!("{} is a pre-sketched sample at -c {}, denser than required by the loaded genome sketch(es) (-c {}) -- using it as-is.", read_file, read_sketch.c, genome_c);
        }

        return Some(read_sketch);
    } else {
        // Raw reads are sketched fresh here, so unlike the pre-sketched-file
        // branch above there's no unrecoverable case: if the requested -c is
        // sparser than what the loaded genome(s) need, just sketch at the
        // denser rate genome_c instead -- never sparser than requested, and
        // always safe/correct (a sample can't be sparser than any genome it's
        // compared against).
        let effective_c = args.c.min(genome_c);
        if args.c > genome_c{
            info!("{} value of -c for reads is {}; this is sparser than required by the loaded genome sketch(es) (-c {}). Sketching at the denser -c {} instead.", read_file[0], args.c, genome_c, effective_c);
        }
        if genome_k != args.k {
            error!(
                "{} -k {} is not equal to -k {} found in sketches. Continuing without sketching.",
                read_file[0], args.k, genome_k
            );
            return None;
        } else {
            let pipeline_params = crate::parallel_sketch::PipelineParams {
                batch_records: args.sketch_batch_size,
                batch_max_bytes: args.sketch_batch_max_bytes,
                channel_depth: args.sketch_channel_depth,
                num_shards: args.sketch_shards,
            };
            if read_file.len() == 1{
                let read_sketch_opt = if crate::parallel_sketch::should_use_pipeline(args.no_sketch_pipeline, threads, &read_file[0]) {
                    crate::parallel_sketch::sketch_sequences_needle_parallel(&read_file[0], effective_c, args.k, None, false, threads, &pipeline_params)
                } else {
                    sketch_sequences_needle(&read_file[0], effective_c, args.k, None, false)
                };
                return read_sketch_opt;
            }
            else if read_file.len() == 2{
                let read_sketch_opt = if crate::parallel_sketch::should_use_pipeline(args.no_sketch_pipeline, threads, &read_file[0]) {
                    crate::parallel_sketch::sketch_pair_sequences_parallel(&read_file[0], &read_file[1], effective_c, args.k, None, false, DEFAULT_FPR, threads, &pipeline_params)
                } else {
                    sketch_pair_sequences(&read_file[0], &read_file[1], effective_c, args.k, None, false, DEFAULT_FPR)
                };
                return read_sketch_opt;
            }
            else{
                panic!("Internal Error: read_file has length {}. Something went wrong...", read_file.len());
            }
        }
    }
}

fn get_stats<'a>(
    args: &ContainArgs,
    genome_sketch: &'a GenomeSketch,
    sequence_sketch: &SequencesSketch,
    winner_map: Option<&FxHashMap<Kmer, (f64,& GenomeSketch, bool)>>,
    log_reassign: bool
) -> Option<AniResult<'a>> {
    if genome_sketch.k != sequence_sketch.k {
        log::error!(
            "k parameter for reads {} != k parameter for genome {}",
            sequence_sketch.k,
            genome_sketch.k
        );
        std::process::exit(1);
    }
    if genome_sketch.c < sequence_sketch.c {
        log::error!(
            "c parameter for reads {} > c parameter for genome {}",
            sequence_sketch.c,
            genome_sketch.c
        );
        std::process::exit(1);
    }
    let mut contain_count = 0;
    let mut covs = vec![];
    let gn_kmers = &genome_sketch.genome_kmers;
    if (gn_kmers.len() as f64) < args.min_number_kmers {
        return None
    }

    let mut kmers_lost_count = 0;
    for kmer in gn_kmers.iter() {
        if sequence_sketch.kmer_counts.contains_key(kmer) {
            if sequence_sketch.kmer_counts[kmer] == 0{
                continue
            }
            if winner_map.is_some(){
                let map = &winner_map.unwrap();
                if map[kmer].1 != genome_sketch{
                    kmers_lost_count += 1;
                    continue
                }
                contain_count += 1;
                covs.push(sequence_sketch.kmer_counts[kmer]);

            }
            else{
                contain_count += 1;
                covs.push(sequence_sketch.kmer_counts[kmer]);
            }
        }
    }

    let n_kmers = gn_kmers.len();
    if contain_count < args.min_contain {
        log::debug!(
            "Discarding {}/{}: {} contained k-mers is below the minimum evidence floor {}",
            genome_sketch.file_name,
            genome_sketch.first_contig_name,
            contain_count,
            args.min_contain,
        );
        return None;
    }
    let reassign_log = if winner_map.is_some() && log_reassign {
        Some((
            genome_sketch.file_name.as_str(),
            genome_sketch.first_contig_name.as_str(),
            kmers_lost_count,
        ))
    } else {
        None
    };
    let fin = finalize_stats(args, genome_sketch.k, n_kmers, contain_count, covs, reassign_log)?;

    let seq_name = if let Some(sample) = &sequence_sketch.sample_name {
        sample.clone()
    } else {
        sequence_sketch.file_name.clone()
    };
    let kmers_lost = if winner_map.is_some() {
        Some(kmers_lost_count)
    } else {
        None
    };

    Some(AniResult {
        naive_ani: fin.naive_ani,
        final_est_ani: fin.final_est_ani,
        final_est_cov: fin.final_est_cov,
        seq_name,
        gn_name: genome_sketch.file_name.as_str(),
        contig_name: genome_sketch.first_contig_name.as_str(),
        mean_cov: fin.mean_cov,
        median_cov: fin.median_cov,
        containment_index: (contain_count, n_kmers),
        lambda: fin.lambda,
        ani_ci: fin.ani_ci,
        lambda_ci: fin.lambda_ci,
        genome_sketch,
        rel_abund: None,
        seq_abund: None,
        kmers_lost,
    })
}

/// Scalar outputs of the coverage-correction + ANI estimation.
struct Finalized {
    naive_ani: f64,
    final_est_ani: f64,
    final_est_cov: f64,
    median_cov: f64,
    mean_cov: f64,
    lambda: AdjustStatus,
    ani_ci: (Option<f64>, Option<f64>),
    lambda_ci: (Option<f64>, Option<f64>),
}

/// Coverage-correction + ANI math shared by `get_stats` (materialized genome)
/// and the streaming two-stage pass-1 (which never builds a genome `Vec`).
/// Consumes the matched coverage counts `covs` and the genome k-mer count
/// `n_kmers`; applies the minimum-ANI gate (returns None below it). When
/// `reassign_log` is set, logs a drop during pseudotax reassignment.
fn finalize_stats(
    args: &ContainArgs,
    k: usize,
    n_kmers: usize,
    contain_count: usize,
    mut covs: Vec<u32>,
    reassign_log: Option<(&str, &str, usize)>,
) -> Option<Finalized> {
    if covs.is_empty() {
        return None;
    }
    let naive_ani = f64::powf(contain_count as f64 / n_kmers as f64, 1. / k as f64);
    covs.sort();
    let median_cov = covs[covs.len() / 2] as f64;
    let pois = Poisson::new(median_cov).unwrap();
    let mut max_cov = f64::MAX;
    if median_cov < 30. {
        for i in covs.len() / 2..covs.len() {
            let cov = covs[i];
            if pois.cdf(cov.into()) < CUTOFF_PVALUE {
                max_cov = cov as f64;
            } else {
                break;
            }
        }
    }

    let mut full_covs = vec![0; n_kmers - contain_count];
    for cov in covs.iter() {
        if (*cov as f64) <= max_cov {
            full_covs.push(*cov);
        }
    }
    let mean_cov = full_covs.iter().sum::<u32>() as f64 / full_covs.len() as f64;
    let geq1_mean_cov = full_covs.iter().sum::<u32>() as f64 / covs.len() as f64;

    let use_lambda;
    if median_cov > MEDIAN_ANI_THRESHOLD {
        use_lambda = AdjustStatus::High
    } else {
        let test_lambda;
        if args.ratio {
            test_lambda = ratio_lambda(&full_covs, args.min_count_correct)
        } else if args.mme {
            test_lambda = mme_lambda(&full_covs)
        } else if args.nb {
            test_lambda = binary_search_lambda(&full_covs)
        } else if args.mle {
            test_lambda = mle_zip(&full_covs, k as f64)
        } else {
            test_lambda = ratio_lambda(&full_covs, args.min_count_correct)
        };
        if test_lambda.is_none() {
            use_lambda = AdjustStatus::Low
        } else {
            use_lambda = AdjustStatus::Lambda(test_lambda.unwrap());
        }
    }

    let final_est_cov;
    if let AdjustStatus::Lambda(lam) = use_lambda {
        final_est_cov = lam
    } else if median_cov < MAX_MEDIAN_FOR_MEAN_FINAL_EST {
        final_est_cov = geq1_mean_cov;
    } else if args.mean_coverage {
        final_est_cov = geq1_mean_cov;
    } else {
        final_est_cov = median_cov;
    }

    let opt_lambda;
    if use_lambda == AdjustStatus::Low || use_lambda == AdjustStatus::High {
        opt_lambda = None
    } else {
        opt_lambda = Some(final_est_cov)
    };

    let opt_est_ani = ani_from_lambda(opt_lambda, mean_cov, k as f64, &full_covs);

    let final_est_ani;
    if opt_lambda.is_none() || opt_est_ani.is_none() || args.no_adj {
        final_est_ani = naive_ani;
    } else {
        final_est_ani = opt_est_ani.unwrap();
    }

    let min_ani = minimum_ani_fraction(args);
    if final_est_ani < min_ani {
        if let Some((gn, ctg, lost)) = reassign_log {
            log::info!(
                "Genome/contig {}/{} has ANI = {} < {} after reassigning {} k-mers ({} contained k-mers after reassign)",
                gn, ctg, final_est_ani * 100., min_ani * 100., lost, contain_count
            );
        }
        return None;
    }

    let (mut low_ani, mut high_ani, mut low_lambda, mut high_lambda) = (None, None, None, None);
    if !args.no_ci && opt_lambda.is_some() {
        let bootstrap = bootstrap_interval(&full_covs, k as f64, args);
        low_ani = bootstrap.0;
        high_ani = bootstrap.1;
        low_lambda = bootstrap.2;
        high_lambda = bootstrap.3;
    }

    Some(Finalized {
        naive_ani,
        final_est_ani,
        final_est_cov,
        median_cov,
        mean_cov: geq1_mean_cov,
        lambda: use_lambda,
        ani_ci: (low_ani, high_ani),
        lambda_ci: (low_lambda, high_lambda),
    })
}


fn ani_from_lambda(lambda: Option<f64>, _mean: f64, k: f64, full_cov: &[u32]) -> Option<f64> {
    if lambda.is_none() {
        return None;
    }
    let mut contain_count = 0;
    let mut _zero_count = 0;
    for x in full_cov {
        if *x != 0 {
            contain_count += 1;
        } else {
            _zero_count += 1;
        }
    }

    let lambda = lambda.unwrap();
    let adj_index =
        contain_count as f64 / (1. - f64::exp(-lambda)) / full_cov.len() as f64;
    let ret_ani;
    //let ani = f64::powf(1. - pi, 1./k);
    let ani = f64::powf(adj_index, 1. / k);
    if ani < 0. || ani.is_nan() {
        ret_ani = None;
    } else {
        if ani > 1. {
            ret_ani = Some(ani)
        } else {
            ret_ani = Some(ani);
        }
    }
    return ret_ani;
}

fn bootstrap_interval(
    covs_full: &Vec<u32>,
    k: f64,
    args: &ContainArgs,
) -> (Option<f64>, Option<f64>, Option<f64>, Option<f64>) {
    // A local generator avoids mutating thread-local global RNG state when
    // genomes are finalized concurrently through rayon.
    let mut rng = fastrand::Rng::with_seed(DEFAULT_RNG_SEED);
    let num_samp = covs_full.len();
    let iters = 100;
    let mut res_ani = vec![];
    let mut res_lambda = vec![];

    for _ in 0..iters {
        let mut rand_vec = vec![];
        rand_vec.reserve(num_samp);
        for _ in 0..num_samp {
            rand_vec.push(covs_full[rng.usize(..covs_full.len())]);
        }
        let lambda;
        if args.ratio {
            lambda = ratio_lambda(&rand_vec, args.min_count_correct);
        } else if args.mme {
            lambda = mme_lambda(&rand_vec);
        } else if args.nb {
            lambda = binary_search_lambda(&rand_vec);
        } else if args.mle {
            lambda = mle_zip(&rand_vec, k);
        } else {
            lambda = ratio_lambda(&rand_vec,args.min_count_correct);
        }
        let ani = ani_from_lambda(lambda, mean(&rand_vec).unwrap().into(), k, &rand_vec);
        if ani.is_some() && lambda.is_some() {
            if !ani.unwrap().is_nan() && !lambda.unwrap().is_nan() {
                res_ani.push(ani);
                res_lambda.push(lambda);
            }
        }
    }
    res_ani.sort_by(|x, y| x.partial_cmp(y).unwrap());
    res_lambda.sort_by(|x, y| x.partial_cmp(y).unwrap());
    if res_ani.len() < 50 {
        return (None, None, None, None);
    }
    let suc = res_ani.len();
    let low_ani = res_ani[suc * 5 / 100 - 1];
    let high_ani = res_ani[suc * 95 / 100 - 1];
    let low_lambda = res_lambda[suc * 5 / 100 - 1];
    let high_lambda = res_lambda[suc * 95 / 100 - 1];

    return (low_ani, high_ani, low_lambda, high_lambda);
}


fn get_kmer_identity(seq_sketch: &SequencesSketch, estimate_unknown: bool) -> Option<f64>{

    if !estimate_unknown{
        return None
    }

    let mut median = 0;
    let mut mov_avg_median = 0.;
    let mut n = 1.;
    for count in seq_sketch.kmer_counts.values(){
        if *count > 1{
            if *count > median{
                median += 1;
            }
            else{
                median -= 1;
            }
            mov_avg_median += median as f64;
            n += 1.;
        }
    }

    mov_avg_median /= n;
    log::debug!("Estimated continuous median k-mer count for {} is {:.3}", &seq_sketch.file_name, mov_avg_median);
    
    let mut num_1s = 0;
    let mut num_not1s = 0;
    for count in seq_sketch.kmer_counts.values(){
        if *count == 1{
            num_1s += 1;
        }
        else{
            num_not1s += *count;
        }
    }
    //0.1 so no div by 0 error
    let eps = num_not1s as f64 / (num_not1s as f64 + num_1s as f64 + 0.1);
    //dbg!("Automatic id est, 1-to-2 ratio, 2-to-3", eps.powf(1./31.), num_1s as f64 / num_2s as f64, two_to_three);

    if mov_avg_median < MED_KMER_FOR_ID_EST && seq_sketch.mean_read_length < 400.{
        log::info!("{} short-read sample has high diversity compared to sequencing depth (approx. avg depth < 3). Using 99.5% as read accuracy estimate instead of automatic detection for --estimate-unknown.", &seq_sketch.file_name);
        return Some(0.995f64.powf(seq_sketch.k as f64));
    }

    if eps < 1.{
        return Some(eps)
    }
    else{
        return Some(1.)
    }
}

#[cfg(test)]
mod min_number_kmers_tests {
    use super::*;
    use clap::Parser;

    fn parse_contain_args(subcommand: &str, extra_args: &[&str]) -> ContainArgs {
        let mut argv = vec!["sylph", subcommand];
        argv.extend_from_slice(extra_args);
        argv.extend(["database.syldb", "sample.sylsp"]);

        let cli = Cli::try_parse_from(argv).unwrap();
        match cli.mode {
            Mode::Profile(mut args) => {
                args.pseudotax = true;
                args
            }
            Mode::Query(args) => args,
            _ => panic!("expected profile or query arguments"),
        }
    }

    fn synthetic_genome_passes(args: &ContainArgs, genome_kmers: usize, contained: usize) -> bool {
        let genome = GenomeSketch {
            genome_kmers: (0..genome_kmers as u64).collect(),
            pseudotax_tracked_nonused_kmers: Some(vec![]),
            file_name: "genome".to_string(),
            first_contig_name: "contig".to_string(),
            c: 200,
            k: 31,
            gn_size: 1_000,
            min_spacing: 0,
        };
        let mut sequence = SequencesSketch::new("sample".to_string(), 200, 31, false, None, 150.0);
        sequence
            .kmer_counts
            .extend((0..contained as u64).map(|kmer| (kmer, 1)));

        get_stats(args, &genome, &sequence, None, false).is_some()
    }

    #[test]
    fn min_number_kmers_defaults_to_10() {
        let args = parse_contain_args("profile", &[]);

        assert_eq!(args.min_number_kmers, 10.0);
    }

    #[test]
    fn explicit_min_number_kmers_overrides_default() {
        let args = parse_contain_args("profile", &["-M", "32.5"]);

        assert_eq!(args.min_number_kmers, 32.5);
    }

    #[test]
    fn min_number_kmers_discards_small_genome_sketches() {
        let mut args = parse_contain_args("profile", &[]);
        args.minimum_ani = Some(0.0);
        args.no_ci = true;

        assert!(!synthetic_genome_passes(&args, 9, 9));
    }

    #[test]
    fn min_contain_defaults_to_7_for_all_genome_sizes() {
        let mut args = parse_contain_args("profile", &[]);
        assert_eq!(args.min_contain, 7);

        args.minimum_ani = Some(0.0);
        args.no_ci = true;

        assert!(!synthetic_genome_passes(&args, 10, 6));
        assert!(synthetic_genome_passes(&args, 10, 7));
        assert!(!synthetic_genome_passes(&args, 500, 6));
        assert!(synthetic_genome_passes(&args, 500, 7));
    }

    #[test]
    fn explicit_min_contain_overrides_default() {
        let mut args = parse_contain_args("profile", &["--min-contain", "3"]);
        assert_eq!(args.min_contain, 3);

        args.minimum_ani = Some(0.0);
        args.no_ci = true;
        assert!(!synthetic_genome_passes(&args, 500, 2));
        assert!(synthetic_genome_passes(&args, 500, 3));
    }
}
