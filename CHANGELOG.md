## sylph v0.3.0 release: first class support for pseudotax, now called "profile" - 2023-10-01

Continuing development of sylph taxonomic profiling. 

### BREAKING

- `--pseudotax` option in previous version is now a new command called `profile`.
- Databases are enabled for profiling by default. 
- Changed file suffices to `syldb` and `sylsp`.

### Major
- Default parameter changes. --min-spacing is set to 30 now. 
- Made profiling faster with some algorithmic tweaks. 
- Coverage calculated slightly differently
- Many small software changes with respect to threading and outputs

## sylph v0.2.0 release: pseudotax improved - 2023-09-19

### BREAKING
- Sylph's *.sylqueries are no longer compatible with older versions of sylph (< v0.2). Files will need to be resketched. 

### Major
- Fixed a major bug for the `--pseudotax` option that required redesigning file formats. Please use `--enable-pseudotax` when using using `contain --pseudotax` from now on.
- `--pseudotax` option gives relative abundances now. We are gaining some confidence that this approach gives a rough, but surprisingly decent taxonomic classification.  
- Changed how `Eff_cov` is calculated. We just use the median coverage now, except when we apply coverage-adjustment 

### Minor
- Fixed command line ambiguity for sketching outputs. `-s` has been replaced with `-d` for `sylph sketch`.
- Sylph outputs the results after processing every sample, instead of batching results, now


## sylph v0.1.0 release - 2023-09-03

### Major

- Added `--pseudotax` option, similar to the `-w` option in mash screen, where k-mers are assigned to the highest ANI genome so redundancy is removed. The output is a very rough taxonomic classification of the sample. 

### Minor

- Some fixes and parameter changes from the v0.0.x releases. 
