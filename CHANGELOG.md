## sylph v0.2.0 release: pseudotax improved - 2023-09-XX

### BREAKING
- Sylph's *.sylqueries are no longer compatible with older versions of sylph (< v0.2). Files will need to be resketched. 

### Major
- Fixed a major bug for the `--pseudotax` option that required redesigning file formats. Please use `--enable-pseudotax` when using using `contain --pseudotax` from now on.
- `--pseudotax` option gives relative abundances now. We are gaining some confidence that this approach gives a rough, but surprisingly decent taxonomic classification.  
- Changed how `Eff_cov` is calculated. We just use the median coverage now, except when we apply coverage-adjustment 

### Minor
- Fixed command line ambiguity for sketching outputs. `-s` has been replaced with `-d` for `sylph sketch`.


## sylph v0.1.0 release - 2023-09-03

### Major

- Added `--pseudotax` option, similar to the `-w` option in mash screen, where k-mers are assigned to the highest ANI genome so redundancy is removed. The output is a very rough taxonomic classification of the sample. 

### Minor

- Some fixes and parameter changes from the v0.0.x releases. 
