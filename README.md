# DLattice
A system for spatiotemporal data analysis.

## Background
The system was developed as part of a PhD Research. The related thesis can be found at http://mtc-m16d.sid.inpe.br/col/sid.inpe.br/mtc-m19/2011/11.15.22.02/doc/publicacao.pdf

### How to use
The spatiotemporal data to be analysed must be presented in a gnl file, where the first line contains the metadata and the next line contains the values as a sequence. For example, a time series composed by 1000 measures should be stored in a gnl file where the first line would be: 0 0 0 1000; and the second line would contain the 1000 measures. A 2-dimensional 64x64 image sequence composed by 1000 images should be stored in a gnl file where the first line would be: 0 64 64 1000, and the second line would contain the values of pixels linearly arranged, rowsXcolumns, for each image in sequence.

Compile: >> g++ -o <filename> -I<path_to_NR_folder> gnlsys.cpp

Run: >> ./gnlsys <path_to_folder_containing_gnl_files_related_to_one_single_data_system> 

### Output
The output is an xml file containing the analysis results.
