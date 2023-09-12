### Manual

A whole-genome aligner, based on the algorithm of SibeliaZ-LCB, but operating on a variation graph instead of a compacted de Bruijn graph.
We have provided two alternative ways of performing the MSA in each collinear block - one based on SPOA and another one using a modification of the library poapy.

The project can be used with the following command.
```
<path_to_wga> -i <path_to_input>
```
where the path_to_wga is the path to script wga_scripts/calculate_wga.py 
and path_to_input is the relative path to the input .gfa file from the source directory of the project.
The following not required flags can be used.
 - -a - Abundance pruning parameter, default to 150. Vertices occurring more than a times are not considered as collinear walk seeds.
 - -b - Maximal bubble size (for each walk), default to 200 residues.
 - -m - Minimal collinear walk length, default to 50.
 - -s - Seed sorting mode, default to 'no'. Possible values:
                        'no' (random vertex order), 'nr_occurrences' (most occurrences first),
                        'length' (longest labels first).
 - --align - Use to align sequences within blocks. Default.
 - --no-align - Use in order not to align sequences within blocks. Tool aligns by default.
 - --match - Match score in alignment (used if --align is True). Default to 5.
 - --mismatch - Mismatch penalty in alignment (used if --align is True). Default to -4.
 - --gap - Gap penalty in alignment (used if --align is True). Default to -8.
 - --align_mode - Alignment mode --- 'poapy' or 'spoa'. Default to 'poapy'.

After the code is ran, the following folders will be created in the source directory:
 - blocks --- contains .gff files with block coordinates,
 - maf --- contains .maf files with WGA,
 - genome_idx_to_name --- contains dictionaries translating indices of genomes in the variation graph to their names, 
 - vertex_name_to_idx --- contains dictionaries translating names of vertices to their indoces in the graph.

References

Minkin, I., Medvedev, P. Scalable multiple whole-genome alignment and locally collinear block construction with SibeliaZ.
Nat Commun 11, 6327 (2020). https://doi.org/10.1038/s41467-020-19777-8

https://github.com/ljdursi/poapy

Vaser R, Sović I, Nagarajan N, Šikić M. Fast and accurate de novo genome assembly from long uncorrected reads. Genome Res. 2017 May;27(5):737-746. doi: 10.1101/gr.214270.116. Epub 2017 Jan 18. PMID: 28100585; PMCID: PMC5411768.
