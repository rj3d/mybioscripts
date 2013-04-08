mybioscripts
============
A collection of various bioinformatics scripts.

cachedwriter
------------
A simple utility to reduce the frequency disk writes.
Useful when doing multiple file writes with a lot of logic in between writes.

fastarator
------------
A simple utility for reading fasta files.
Can either return an list of all of the title, sequence pairs or iterate over them. 

getORFs
------------
A rough script to find all ORFs of a given size window in an input set of sequences.
Currently all paramaters for 
Will be updated and made more user friendly.

getSpacers
------------
Finds all spacers (Cas9 targeting sequences) in an input set of sequences that will
not have off-targets within the set.

parallelBlast
------------
A rough script to decompose a single blast job into multiple blast jobs that run in parallel.
Currently only runs blastn and parameters must be hardcoded.
Will be updated and made more user friendly.

parallelPMM
------------
Given an input set of protein sequences runs multiple patmatmotifs jobs to find motifs.
Generates a summary and distribution file of the runs at the end.
