use assert_cmd::prelude::*; // Add methods on commands
use std::str;
use std::fs;
use std::path::Path;
use serial_test::serial;
use std::process::Command; // Run programs

fn fresh(){
    let _ = fs::remove_dir_all("./tests/results/test_sketch_dir");
}

#[serial]
#[test]
fn test_sketch_commands() {
   let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("test_files/e.coli-EC590.fasta.gz")
        .arg("test_files/e.coli-K12.fasta.gz")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("-o")
        .arg("tests/results/test_sketch_dir/db")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("profile")
        .arg("./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp")
        .arg("./tests/results/test_sketch_dir/db.syldb")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("profile")
        .arg("-l")
        .arg("./test_files/list.txt")
        .assert();
    assert.success().code(0);


    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("profile")
        .arg("./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("profile")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("-i")
        .arg("-m")
        .arg("90")
        .assert();
    assert.success().code(0);

    let mut cmd= Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-1")
        .arg("./test_files/t1.fq")
        .arg("-2")
        .arg("./test_files/t2.fq")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);
    assert!(Path::new("./tests/results/test_sketch_dir/t1.fq.paired.sylsp").exists(), "Output file was not created");
    fresh();

    let mut cmd= Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("--l1")
        .arg("./test_files/pair_list1.txt")
        .arg("--l2")
        .arg("./test_files/pair_list2.txt")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);
    assert!(Path::new("./tests/results/test_sketch_dir/t1.fq.paired.sylsp").exists(), "Output file was not created");

    fresh();
    let mut cmd= Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-g")
        .arg("./test_files/t1.fq")
        .arg("-r")
        .arg("./test_files/t2.fq")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir/testdb")
        .assert();
    assert.success().code(0);
    assert!(Path::new("./tests/results/test_sketch_dir/t2.fq.sylsp").exists(), "Output file was not created");
    assert!(Path::new("./tests/results/test_sketch_dir/testdb.syldb").exists(), "Output file was not created");
}

#[serial]
#[test]
fn test_profile_vs_query(){
    fresh();

    let mut output = Command::cargo_bin("sylph").unwrap();
    let output = output
        .arg("profile")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .output()
        .expect("Output failed");
    let stdout = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    dbg!(stdout.matches('\n').count());
    assert!(stdout.matches('\n').count() == 2);

    let mut output = Command::cargo_bin("sylph").unwrap();
    let output = output
        .arg("query")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("./test_files/e.coli-K12.fasta.gz")
        .output()
        .expect("Output failed");
    let stdout = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    dbg!(stdout.matches('\n').count());
    println!("{}",stdout);
    assert!(stdout.matches('\n').count() == 4);
}

#[serial]
#[test]
fn test_sketch_list(){
    fresh();
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-r")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir/db")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);
    assert!(Path::new("./tests/results/test_sketch_dir/e.coli-EC590.fasta.gz.sylsp").exists(), "Output file was not created");
    assert!(Path::new("./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp").exists(), "Output file was not created");
    assert!(!Path::new("./tests/results/test_sketch_dir/db.syldb").exists(), "Output file was created");
    fresh();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-g")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir/db")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);
    assert!(!Path::new("./tests/results/test_sketch_dir/e.coli-EC590.fasta.gz.sylsp").exists(), "Output file was created");
    assert!(!Path::new("./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp").exists(), "Output file was created");
    assert!(Path::new("./tests/results/test_sketch_dir/db.syldb").exists(), "Output file was not created");
    fresh();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("--gl")
        .arg("test_files/list.txt")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir/db")
        .assert();
    assert.success().code(0);
    assert!(Path::new("./tests/results/test_sketch_dir/db.syldb").exists(), "Output file was not created");
    fresh();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("--rl")
        .arg("test_files/list.txt")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir/db")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);
    assert!(!Path::new("./tests/results/test_sketch_dir/db.syldb").exists(), "Output file was not created");
    assert!(Path::new("./tests/results/test_sketch_dir/e.coli-EC590.fasta.gz.sylsp").exists(), "Output file was not created");
    assert!(Path::new("./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp").exists(), "Output file was not created");
    fresh();

}
#[serial]
#[test]
fn test_profile_disabling(){
    fresh();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-g")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir/db")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .arg("--disable-profiling")
        .assert();
    assert.success().code(0);

    let mut output = Command::cargo_bin("sylph").unwrap();
    let assert = output
        .arg("profile")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("./tests/results/test_sketch_dir/db.syldb")
        .assert();
    assert.failure().code(1);

    let mut output = Command::cargo_bin("sylph").unwrap();
    let assert = output
        .arg("query")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("./tests/results/test_sketch_dir/db.syldb")
        .assert();
    assert.success().code(0);

    fresh();
}
#[serial]
#[test]
fn test_sketch_fasta_fastq_concord(){
    fresh();
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("-o")
        .arg("./tests/results/test_sketch_dir/db")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);

    let mut output = Command::cargo_bin("sylph").unwrap();
    let out1 = output
        .arg("profile")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("./tests/results/test_sketch_dir/db.syldb")
        .output()
        .expect("Fail");

    let mut output = Command::cargo_bin("sylph").unwrap();
    let out2 = output
        .arg("profile")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .output()
        .expect("Fail");

    let mut output = Command::cargo_bin("sylph").unwrap();
    let out3 = output
        .arg("profile")
        .arg("./tests/results/test_sketch_dir/o157_reads.fastq.gz.sylsp")
        .arg("./tests/results/test_sketch_dir/db.syldb")
        .output()
        .expect("Fail");

    let stdout1 = str::from_utf8(&out1.stdout).expect("Output was not valid UTF-8");
    let stdout2 = str::from_utf8(&out2.stdout).expect("Output was not valid UTF-8");
    let stdout3 = str::from_utf8(&out3.stdout).expect("Output was not valid UTF-8");

    assert!(stdout1 == stdout2);
    assert!(stdout1 == stdout3);
    assert!(stdout2 == stdout3);

    fresh();
}
#[serial]
#[test]
fn test_sample_names(){
    fresh();
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-1")
        .arg("test_files/t1.fq")
        .arg("-2")
        .arg("test_files/t2.fq")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .arg("--lS")
        .arg("./test_files/single_sample.txt")
        .assert();
    assert.success().code(0);
    assert!(Path::new("./tests/results/test_sketch_dir/SAMPLE_TEST.paired.sylsp").exists(), "Output file was not created");
    fresh();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("test_files/t1.fq")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .arg("--lS")
        .arg("./test_files/sample_list.txt")
        .assert();
    assert.success().code(0);
    assert!(Path::new("./tests/results/test_sketch_dir/S1.sylsp").exists(), "Output file was not created");
    assert!(Path::new("./tests/results/test_sketch_dir/S2.sylsp").exists(), "Output file was not created");

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd
        .arg("profile")
        .arg("./tests/results/test_sketch_dir/S2.sylsp")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .output().unwrap();
    let stdout = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    dbg!(&stdout);
    assert!(stdout.contains("S2"));
    assert!(!stdout.contains("o157_reads"));

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-1")
        .arg("test_files/t1.fq")
        .arg("-2")
        .arg("test_files/t2.fq")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .arg("-S")
        .arg("SAMPLE_TEST_S")
        .assert();
    assert.success().code(0);
    assert!(Path::new("./tests/results/test_sketch_dir/SAMPLE_TEST_S.paired.sylsp").exists(), "Output file was not created, -S");

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-1")
        .arg("test_files/t1.fq")
        .arg("test_files/t1.fq")
        .arg("-2")
        .arg("test_files/t2.fq")
        .arg("test_files/t2.fq")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .arg("-S")
        .arg("SAMPLE_TEST_S")
        .arg("SAMPLE_TEST_S1")
        .assert();
    assert.success().code(0);
    assert!(Path::new("./tests/results/test_sketch_dir/SAMPLE_TEST_S1.paired.sylsp").exists(), "Output file was not created, -S");

    fresh();
}
#[serial]
#[test]
fn test_fpr(){
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-1")
        .arg("test_files/t1.fq")
        .arg("-2")
        .arg("test_files/t2.fq")
        .arg("-d ")
        .arg("./tests/results/test_sketch_dir")
        .arg("0")
        .assert();
    assert.success().code(0);
    fresh();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-1")
        .arg("test_files/t1.fq")
        .arg("-2")
        .arg("test_files/t2.fq")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .arg("--fpr")
        .arg("0.001")
        .assert();
    assert.success().code(0);
    fresh();
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-1")
        .arg("test_files/t1.fq")
        .arg("-2")
        .arg("test_files/t2.fq")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .arg("--fpr")
        .arg("2")
        .assert();
    assert.failure().code(1);
    fresh();

}
#[serial]
#[test]
fn test_raw_inputs_profile_simple(){
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("profile")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("-1")
        .arg("test_files/t1.fq")
        .arg("-2")
        .arg("test_files/t2.fq")
        .assert();
    assert.success().code(0);
    fresh();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("profile")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("-1")
        .arg("test_files/t1.fq")
        .assert();
    assert.failure().code(1);
    fresh();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("profile")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("-1")
        .arg("test_files/k12_R1.fq")
        .arg("test_files/t1.fq")
        .arg("-2")
        .arg("test_files/k12_R2.fq")
        .arg("test_files/t1.fq")
        .assert();
    assert.success().code(0);
    fresh();
    
}

#[serial]
#[test]
fn test_estimate_read_counts(){
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd
        .arg("profile")
        .arg("--estimate-read-counts")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("./test_files/o157_reads.fastq.gz")
        .output()
        .expect("output failed");
    let stdout_1 = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    let mut lines = stdout_1.lines();
    dbg!(stdout_1);
    lines.next();
    let output = lines.next().unwrap();
    let split : Vec<&str> = output.split('\t').collect();
    assert!(split[3].parse::<f64>().unwrap() > 1000.0);

    fresh();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd
        .arg("profile")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("./test_files/o157_reads.fastq.gz")
        .output()
        .expect("output failed");
    let stdout_1 = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    let mut lines = stdout_1.lines();
    dbg!(stdout_1);
    lines.next();
    let output = lines.next().unwrap();
    let split : Vec<&str> = output.split('\t').collect();
    assert!(split[3].parse::<f64>().unwrap() < 101.00);

    fresh();

}

#[serial]
#[test]
fn test_raw_inputs_profile_with_sketch(){
    
    let mut output = Command::cargo_bin("sylph").unwrap();
    let output = output
        .arg("profile")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("-1")
        .arg("test_files/k12_R1.fq")
        .arg("-2")
        .arg("test_files/k12_R2.fq")
        .output()
        .expect("Output failed");
    let stdout_1 = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-1")
        .arg("test_files/k12_R1.fq")
        .arg("-2")
        .arg("test_files/k12_R2.fq")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);

    let mut output = Command::cargo_bin("sylph").unwrap();
    let output = output
        .arg("profile")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./tests/results/test_sketch_dir/k12_R1.fq.paired.sylsp")
        .output()
        .expect("Output failed");
    let stdout_2 = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");

    assert!(stdout_1 == stdout_2);
}

#[serial]
#[test]
fn test_inspect(){
   let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("test_files/e.coli-EC590.fasta.gz")
        .arg("test_files/e.coli-K12.fasta.gz")
        .arg("test_files/o157_reads.fastq.gz")
        .arg("-o")
        .arg("tests/results/test_sketch_dir/db")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let assert = cmd
        .arg("sketch")
        .arg("-1")
        .arg("test_files/k12_R1.fq")
        .arg("-2")
        .arg("test_files/k12_R2.fq")
        .arg("-d")
        .arg("./tests/results/test_sketch_dir")
        .assert();
    assert.success().code(0);

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd
        .arg("inspect")
        .arg("./tests/results/test_sketch_dir/k12_R1.fq.paired.sylsp")
        .output()
        .expect("Output failed");

    let stdout = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    assert!(stdout.contains("k12_R1.fq"));

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd
        .arg("inspect")
        .arg("./tests/results/test_sketch_dir/db.syldb")
        .output()
        .expect("Output failed");
    let stdout = str::from_utf8(&output.stdout).expect("Output was not valid UTF-8");
    assert!(stdout.contains("e.coli-EC590.fasta.gz"));
    assert!(stdout.contains("e.coli-K12.fasta.gz"));

}

#[serial]
#[test]
fn test_two_stage_db_convert_and_profile(){
    fresh();
    let dir = "./tests/results/two_stage_db";
    let _ = fs::remove_dir_all(dir);
    fs::create_dir_all(dir).unwrap();

    // Dense (-c 50) database carrying profiling k-mers.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("./test_files/e.coli-K12.fasta.gz")
        .arg("-o").arg(format!("{}/db_c50", dir))
        .assert().success().code(0);

    // Dense (-c 50) read sample.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("-d").arg(dir)
        .assert().success().code(0);

    let dense_db = format!("{}/db_c50.syldb", dir);
    let sample = format!("{}/o157_reads.fastq.gz.sylsp", dir);

    // Convert the dense db into a two-stage seekable database: dense blocks at
    // c=50, sparse stage-1 screen index at c=200.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("convert-db-two-screen")
        .arg(&dense_db)
        .arg("--screen-c").arg("200")
        .arg("-o").arg(format!("{}/db2", dir))
        .assert().success().code(0);
    let two_stage_db = format!("{}/db2.syl2db", dir);
    assert!(Path::new(&two_stage_db).exists(), "convert-db-two-screen did not produce a .syl2db");
    // The two-stage db should be no larger than the dense .syldb it came from
    // (dense blocks are Golomb-Rice compressed; only the small sparse index adds).
    let dense_sz = fs::metadata(&dense_db).unwrap().len();
    let two_sz = fs::metadata(&two_stage_db).unwrap().len();
    assert!(two_sz < dense_sz, "compressed two-stage db ({} B) not smaller than dense db ({} B)", two_sz, dense_sz);

    // Profile against the .syl2db directly -- two-stage is auto-detected from
    // the file extension, no flag needed: stage 1 screens via the sparse
    // index, stage 2 decodes only the screened genomes' dense blocks.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd.arg("profile")
        .arg(&two_stage_db).arg(&sample)
        .output().expect("Output failed");
    assert!(output.status.success());
    let from_db2 = str::from_utf8(&output.stdout).expect("not UTF-8").to_string();
    assert!(from_db2.contains("e.coli-o157.fasta.gz"));

    // The detected genome set must equal a plain single-stage dense profile.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd.arg("profile")
        .arg(&dense_db).arg(&sample)
        .output().expect("Output failed");
    let single = str::from_utf8(&output.stdout).expect("not UTF-8").to_string();

    let detected = |tsv: &str| -> Vec<String> {
        let mut v: Vec<String> = tsv.lines().skip(1)
            .filter_map(|l| l.split('\t').nth(1).map(|s| s.to_string()))
            .collect();
        v.sort();
        v
    };
    assert_eq!(detected(&from_db2), detected(&single),
        "two-stage .syl2db and single-stage detected different genome sets");

    // `query` (containment-only, no reassignment) also auto-detects a .syl2db.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd.arg("query").arg(&two_stage_db).arg(&sample)
        .output().expect("Output failed");
    assert!(output.status.success(), "query against a .syl2db should succeed");
    let query_out = str::from_utf8(&output.stdout).expect("not UTF-8").to_string();
    assert!(query_out.contains("e.coli-o157.fasta.gz"));
}

#[serial]
#[test]
fn test_two_stage_individual_records(){
    fresh();
    let dir = "./tests/results/two_stage_indiv";
    let _ = fs::remove_dir_all(dir);
    fs::create_dir_all(dir).unwrap();

    // Dense (-c 50) database built with --individual-records: e.coli-o157 has two
    // records, so multiple database entries share one file name -- the case that
    // must be preserved per record by convert-db-two-screen (and rejected by the densify
    // fallback).
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50").arg("-i")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("./test_files/e.coli-K12.fasta.gz")
        .arg("-o").arg(format!("{}/db_c50", dir))
        .assert().success().code(0);

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("-d").arg(dir)
        .assert().success().code(0);

    let dense_db = format!("{}/db_c50.syldb", dir);
    let sample = format!("{}/o157_reads.fastq.gz.sylsp", dir);

    // Convert to a two-stage db (per-record blocks are written individually).
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("convert-db-two-screen").arg(&dense_db)
        .arg("--screen-c").arg("200")
        .arg("-o").arg(format!("{}/db2", dir))
        .assert().success().code(0);
    let two_stage_db = format!("{}/db2.syl2db", dir);

    // Per-record key = genome_file (col 2) + contig name (last col).
    let detected = |tsv: &str| -> Vec<String> {
        let mut v: Vec<String> = tsv.lines().skip(1)
            .filter_map(|l| {
                let cols: Vec<&str> = l.split('\t').collect();
                if cols.len() < 2 { return None; }
                Some(format!("{}\t{}", cols[1], cols[cols.len() - 1]))
            })
            .collect();
        v.sort();
        v
    };

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let two = cmd.arg("profile").arg(&two_stage_db).arg(&sample)
        .output().expect("Output failed");
    assert!(two.status.success());
    let two = str::from_utf8(&two.stdout).expect("not UTF-8").to_string();
    assert!(two.contains("e.coli-o157.fasta.gz"));

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let single = cmd.arg("profile").arg(&dense_db).arg(&sample)
        .output().expect("Output failed");
    let single = str::from_utf8(&single.stdout).expect("not UTF-8").to_string();

    // convert-db-two-screen + auto-detected two-stage must reproduce single-stage
    // per-record detections (no collapsing/merging of records sharing a file
    // name).
    assert_eq!(detected(&two), detected(&single),
        "two-stage .syl2db lost or merged individual records vs single-stage");
}

/// Combining a `.syl2db` and a plain `.syldb` in one `profile` call must union
/// their genomes (screen survivors from the two-stage db + all genomes from the
/// plain db) and detect the same genomes as a single dense database containing
/// everything.
#[serial]
#[test]
fn test_two_stage_mixed_sources(){
    fresh();
    let dir = "./tests/results/two_stage_mixed";
    let _ = fs::remove_dir_all(dir);
    fs::create_dir_all(dir).unwrap();

    // o157 goes into a two-stage db; K12 and EC590 stay as a plain dense db.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("-o").arg(format!("{}/db_o157", dir))
        .assert().success().code(0);
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50")
        .arg("./test_files/e.coli-K12.fasta.gz")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("-o").arg(format!("{}/db_plain", dir))
        .assert().success().code(0);
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("-d").arg(dir)
        .assert().success().code(0);

    let sample = format!("{}/o157_reads.fastq.gz.sylsp", dir);

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("convert-db-two-screen").arg(format!("{}/db_o157.syldb", dir))
        .arg("--screen-c").arg("200")
        .arg("-o").arg(format!("{}/db_o157_2", dir))
        .assert().success().code(0);
    let two_stage_db = format!("{}/db_o157_2.syl2db", dir);
    let plain_db = format!("{}/db_plain.syldb", dir);

    // Mixed: one .syl2db + one plain .syldb in the same profile call.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let mixed = cmd.arg("profile").arg(&two_stage_db).arg(&plain_db).arg(&sample)
        .output().expect("Output failed");
    assert!(mixed.status.success());
    let mixed = str::from_utf8(&mixed.stdout).expect("not UTF-8").to_string();

    // Reference: a single dense .syldb containing all three genomes together.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("./test_files/e.coli-K12.fasta.gz")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("-o").arg(format!("{}/db_all", dir))
        .assert().success().code(0);
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let single = cmd.arg("profile").arg(format!("{}/db_all.syldb", dir)).arg(&sample)
        .output().expect("Output failed");
    let single = str::from_utf8(&single.stdout).expect("not UTF-8").to_string();

    let detected = |tsv: &str| -> Vec<String> {
        let mut v: Vec<String> = tsv.lines().skip(1)
            .filter_map(|l| l.split('\t').nth(1).map(|s| s.to_string()))
            .collect();
        v.sort();
        v
    };
    assert_eq!(detected(&mixed), detected(&single),
        "mixed .syl2db + .syldb profile detected a different genome set than a combined single-stage database");
}

/// The three `--small-genome-screen` modes are different *storage* for the same
/// screen, so they must profile identically: `loosen` (widen the pooled index),
/// `band` (extra keys in a separate value band) and `none` + reader-side
/// `--screen-small-genomes`. A deliberately coarse `--screen-c` with a high
/// `--min-sparse-kmers` puts even these multi-Mbp genomes below the floor, which
/// is what makes them exercise the mixed-rate paths.
#[serial]
#[test]
fn test_two_stage_small_genome_screen_modes_agree(){
    fresh();
    let dir = "./tests/results/two_stage_screen_modes";
    let _ = fs::remove_dir_all(dir);
    fs::create_dir_all(dir).unwrap();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("./test_files/e.coli-K12.fasta.gz")
        .arg("./test_files/e.coli-EC590.fasta.gz")
        .arg("-o").arg(format!("{}/db", dir))
        .assert().success().code(0);
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("-d").arg(dir)
        .assert().success().code(0);
    let sample = format!("{}/o157_reads.fastq.gz.sylsp", dir);

    let floor = "4000";
    let convert = |mode: &str| -> String {
        let out = format!("{}/db_{}", dir, mode);
        let mut cmd = Command::cargo_bin("sylph").unwrap();
        cmd.arg("convert-db-two-screen").arg(format!("{}/db.syldb", dir))
            .arg("--screen-c").arg("20000")
            .arg("--min-sparse-kmers").arg(floor)
            .arg("--small-genome-screen").arg(mode)
            .arg("-o").arg(&out)
            .assert().success().code(0);
        format!("{}.syl2db", out)
    };
    let (loosen_db, band_db, none_db) = (convert("loosen"), convert("band"), convert("none"));

    // Not applying the floor must leave the file smaller (no extra keys at all),
    // and banding must not cost more than loosening the whole index.
    let size = |p: &str| fs::metadata(p).unwrap().len();
    assert!(size(&none_db) < size(&band_db), "the floor must add keys somewhere");
    assert!(size(&band_db) <= size(&loosen_db),
        "banding the extra keys should not cost more than widening the pooled index");

    // Each profile also dumps its stage-1 survivors: the three modes store the
    // same screen k-mer set, so that dump must be identical, not merely lead to
    // the same profile.
    let profile = |db: &str, tag: &str, extra: &[&str]| -> (String, String) {
        let dump = format!("{}/screen_{}.tsv", dir, tag);
        let mut cmd = Command::cargo_bin("sylph").unwrap();
        let out = cmd.arg("profile").arg(db).arg(&sample).args(extra)
            .arg("--screen-dump").arg(&dump)
            .output().expect("Output failed");
        assert!(out.status.success());
        (
            str::from_utf8(&out.stdout).expect("not UTF-8").to_string(),
            fs::read_to_string(&dump).unwrap(),
        )
    };
    let (loosen, loosen_screen) = profile(&loosen_db, "loosen", &[]);
    assert!(loosen.contains("e.coli-o157.fasta.gz"));
    assert!(
        loosen_screen.starts_with("Genome_id\tGenome_file\tContig_name\t")
            && loosen_screen.lines().count() > 1,
        "expected a screen dump with identifiable rows, got: {}",
        loosen_screen
    );
    let (band, band_screen) = profile(&band_db, "band", &[]);
    assert_eq!(band, loosen, "band mode profiled differently from loosen mode");
    assert_eq!(band_screen, loosen_screen, "band mode passed different stage-1 survivors");
    let (reader, reader_screen) = profile(&none_db, "reader", &["--screen-small-genomes", floor]);
    assert_eq!(reader, loosen,
        "reader-side small-genome screening profiled differently from loosen mode");
    assert_eq!(reader_screen, loosen_screen,
        "reader-side small-genome screening passed different stage-1 survivors");
    // Enabling the reader-side screen on an already-densified database must also
    // be a no-op rather than double-counting.
    let (both, both_screen) = profile(&loosen_db, "loosen_reader", &["--screen-small-genomes", floor]);
    assert_eq!(both, loosen,
        "reader-side screening changed a database that was already densified at build time");
    assert_eq!(both_screen, loosen_screen,
        "reader-side screening changed the survivors of an already-densified database");
}

/// A genome with very few total dense k-mers (a short contig/virus) must
/// trigger the convert-db-two-screen warning, and the resulting .syl2db must still be
/// valid and profile the other (normal-sized) genomes correctly.
#[serial]
#[test]
fn test_two_stage_small_genome_warning(){
    fresh();
    let dir = "./tests/results/two_stage_small_genome";
    let _ = fs::remove_dir_all(dir);
    fs::create_dir_all(dir).unwrap();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("./test_files/tiny_virus.fasta")
        .arg("-o").arg(format!("{}/db", dir))
        .assert().success().code(0);
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("-d").arg(dir)
        .assert().success().code(0);

    let sample = format!("{}/o157_reads.fastq.gz.sylsp", dir);

    // --min-contain 0 disables the (unrelated) default dense-k-mer-count
    // exclusion, so tiny_virus is kept in the db and hits the adaptive-floor
    // sparse-screen path below instead of being dropped outright.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let convert = cmd.arg("convert-db-two-screen").arg(format!("{}/db.syldb", dir))
        .arg("--screen-c").arg("200")
        .arg("--min-contain").arg("0")
        .arg("-o").arg(format!("{}/db2", dir))
        .output().expect("Output failed");
    assert!(convert.status.success());
    let stderr = str::from_utf8(&convert.stderr).expect("not UTF-8").to_string();
    assert!(
        stderr.contains("tiny_virus") && stderr.contains("stage-1 screen k-mers"),
        "expected a small-genome warning mentioning tiny_virus and its screen k-mer count, got: {}",
        stderr
    );

    // The db is still valid and profiles the normal-sized genome correctly.
    let two_stage_db = format!("{}/db2.syl2db", dir);
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd.arg("profile").arg(&two_stage_db).arg(&sample)
        .output().expect("Output failed");
    assert!(output.status.success());
    let out = str::from_utf8(&output.stdout).expect("not UTF-8").to_string();
    assert!(out.contains("e.coli-o157.fasta.gz"));
}

/// By default (`--min-contain 7`), genomes with fewer dense k-mers than the
/// `profile`/`query` hit threshold are dropped from the two-stage db entirely,
/// with a warning -- rather than kept only to silently never be reported.
#[serial]
#[test]
fn test_two_stage_min_contain_excludes_tiny_genome(){
    fresh();
    let dir = "./tests/results/two_stage_min_contain";
    let _ = fs::remove_dir_all(dir);
    fs::create_dir_all(dir).unwrap();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("./test_files/tiny_virus.fasta")
        .arg("-o").arg(format!("{}/db", dir))
        .assert().success().code(0);

    // Default --min-contain (7): tiny_virus has only 4 dense k-mers at -c 50,
    // so it's excluded with a warning rather than silently kept-but-unreachable.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let convert = cmd.arg("convert-db-two-screen").arg(format!("{}/db.syldb", dir))
        .arg("--screen-c").arg("200")
        .arg("-o").arg(format!("{}/db2", dir))
        .output().expect("Output failed");
    assert!(convert.status.success());
    let stderr = str::from_utf8(&convert.stderr).expect("not UTF-8").to_string();
    assert!(
        stderr.contains("tiny_virus") && stderr.contains("--min-contain=7") && stderr.contains("excluded"),
        "expected a --min-contain exclusion warning mentioning tiny_virus, got: {}",
        stderr
    );

    let two_stage_db = format!("{}/db2.syl2db", dir);
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd.arg("inspect").arg(&two_stage_db)
        .output().expect("Output failed");
    assert!(output.status.success());
    let out = str::from_utf8(&output.stdout).expect("not UTF-8").to_string();
    assert!(out.contains("e.coli-o157"), "surviving genome should still be present: {}", out);
    assert!(!out.contains("tiny_virus"), "excluded genome should not be present: {}", out);
}

/// The same genome loaded from two different sources in one profile call
/// (e.g. present in both a .syl2db and a plain .syldb) must trigger a
/// duplicate-genome warning rather than silently producing nondeterministic
/// reassignment output.
#[serial]
#[test]
fn test_two_stage_duplicate_genome_warning(){
    fresh();
    let dir = "./tests/results/two_stage_dup";
    let _ = fs::remove_dir_all(dir);
    fs::create_dir_all(dir).unwrap();

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50")
        .arg("./test_files/e.coli-o157.fasta.gz")
        .arg("-o").arg(format!("{}/db", dir))
        .assert().success().code(0);
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("sketch").arg("-c").arg("50")
        .arg("./test_files/o157_reads.fastq.gz")
        .arg("-d").arg(dir)
        .assert().success().code(0);

    let sample = format!("{}/o157_reads.fastq.gz.sylsp", dir);
    let plain_db = format!("{}/db.syldb", dir);

    let mut cmd = Command::cargo_bin("sylph").unwrap();
    cmd.arg("convert-db-two-screen").arg(&plain_db)
        .arg("--screen-c").arg("200")
        .arg("-o").arg(format!("{}/db2", dir))
        .assert().success().code(0);
    let two_stage_db = format!("{}/db2.syl2db", dir);

    // Same genome loaded via both the .syl2db and the plain .syldb it came from.
    let mut cmd = Command::cargo_bin("sylph").unwrap();
    let output = cmd.arg("profile").arg(&two_stage_db).arg(&plain_db).arg(&sample)
        .output().expect("Output failed");
    assert!(output.status.success(), "duplicate genomes should warn, not fail");
    let stderr = str::from_utf8(&output.stderr).expect("not UTF-8").to_string();
    assert!(
        stderr.contains("appear more than once"),
        "expected a duplicate-genome warning, got: {}",
        stderr
    );
}
