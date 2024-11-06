use assert_cmd::prelude::*; // Add methods on commands
use std::str;
use std::fs;
use std::path::Path;
use serial_test::serial;
use std::process::Command; // Run programs

fn fresh(){
    Command::new("rm")
        .arg("-r")
        .args(["./tests/results/test_sketch_dir"])
        .spawn();
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
