The project can be used with the following command:
<path_to_wga> -i <path_to_input>
where the path_to_wga is the path to script wga_scripts/calculate_wga.py 
and path_to_input is the relative path to the input .gfa file from the source directory of the project.
The following not required flags can be used:
    - -a - Abundance pruning parameter, default to 150. 
        Vertices occurring more than a times are not considered as collinear walk seeds.
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