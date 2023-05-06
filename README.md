# sylph -  abundanced-corrected minhash for ANI querying against shotgun metagenomes

## Introduction

TODO

##  Install

#### Option 1: Build from source

Requirements:
1. [rust](https://www.rust-lang.org/tools/install) programming language and associated tools such as cargo are required and assumed to be in PATH.
2. A c compiler (e.g. GCC)
3. make

Building takes a few minutes (depending on # of cores).

```sh
git clone https://github.com/bluenote-1577/sylph
cd skani

# If default rust install directory is ~/.cargo
cargo install --path . --root ~/.cargo
sylph TODO

# If ~/.cargo doesn't exist use below commands instead
#cargo build --release
#./target/release/sylph TODO
```

#### Option 2: Pre-built x86-64 linux statically compiled executable

We offer a pre-built statically compiled executable for x86-64 linux systems. That is, if you're on a x86-64 linux system, you can just download the binary and run it without installing anything. 

```sh
wget https://github.com/bluenote-1577/sylph/releases/download/latest/sylph
chmod +x sylph
./sylph -h
```

Note: the binary is compiled with a different set of libraries (musl instead of glibc), possibly impacting performance (slightly). Probably not a huge deal. 

## Quick start

```sh
TODO

```

## Tutorials and manuals

TODO

## Output

TODO

## Citation

TODO
