# introduction of the api with python to AliTV

This script is mainly implement functions to generate a suitable json file which could fit to AliTV.

To achieve the goal, there have three major parts within this library. Not all parts are compulsory.

1. alignment parts: you must have some genomes which may have a large of numbers. It should be nucleotide? the protein parts would implement later.
2. Order part: you must pass a tree to the file which decide the order of these genomes. This script could not build the tree for you. But you could also disable the tree visualization which leave the tree part blank. And use some user-defined order or something else.
3. Annotation parts: you must pass some protein id or information to this script or use some kegg annotation tables.

## The usage of the AliTV

This script has some special functions which not in perl-interface.
It could accept the order of the genomes, the order will be used to stepwise align the sequence. And the order would used to taken as the order which finally print out on the AliTV.

##
`python`


### alignment ways
For now, it only support `blast`.

Maybe later it will implement modules about some other alignment software.
It will generate a temporary directory which stodge all alignment results on it. It would be not deleted unless you pass the `force` parameter.

### The annotation table for annotation/color of the genes or some self-defined regions.

```
The provide table should follow below formats
1. no header
2. separator is tab(\t)
3. each line should be look like
   genome name {TAB} contig id {TAB} start {TAB} end {TAB} annotation name (like some genes)
   e.g.

```

## extra_bin for generate some required file to this library.
1. format_anno_table.py
It could accept a file like `example_data/format_anno_table/name2genes.txt`. You could get the name from anywhere including hmmer search results/self curation/experimental evidence.
And the other one needed is the directory which stodge all the genbank files which should contain identical sequences to the fasta file you passed to blast/alignment. Genbank files should end with suffix `.gbk`.
It could help you to convert the requested genes into an annotation table. Or you should retrieve the corresponding start,end,contig information with corresponding gene ID from the truncated/complete genbank file.

2. truncate_genome_from_target.py


## Potential disasters

1. To perform blast or other alignment, it need to deal with the file name. And if there are any **'_to_'** in the orignal file name. It might raise errors. Details could be seen in the **nearly line 79** of `toolkit/alignment.py`
2.


