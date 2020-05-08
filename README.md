# Introduction to the api with python to AliTV

This project implements functions to generate a validated json file for visualization using AliTV.

To achieve the goal, there have three major parts within this module. Not all parts are compulsory. In general, you could following the pipelines to achieve following purposes: 
> 1. truncate the genomes into suitable regions based on target genes or others.
> 2. stepwise align above truncated genomic sequences instead of pairwise alignment which could save the computational costs. (the order to stepwise alignment could follow a tree or a list of genome IDs.)
> 3. Input a file with information of genes need to be annotated. It could retrieve correct positional information from truncated genomic files and annotated them correctly into the visualization using AliTV.


## Installation

To install the requirement of this project, you just simply run 
`pip install -r requirement.txt`

## Usage

This script has some special functions which not in perl-interface.
It could accept the order of the genomes, the order will be used to stepwise align the sequence. And the order would used to taken as the order which finally print out on the AliTV.

##
`python`


### * Alignment
For now, it only support `blast`. The default parameters for blast operation are `-evalue 1e-3 -outfmt 6 `

Maybe it will implement modules using other alignment software in the future.
It will generate a temporary directory which deposits all alignment results on it. It would be not deleted for now.

### * The annotation table for annotation/color of the genes or some self-defined regions.

The provide table should follow below formats

> 1. no header
> 2. The separator is tab  `\t`
> 3. each line should be look like genome name {TAB} contig id {TAB} start {TAB} end {TAB} annotation name (like some genes)
   e.g.


## Extra binary executions (extra_bin) for generating suitable files to this library.
1. **format_anno_table.py**

It could help you to convert the requested genes into an annotation table. Or you should retrieve the corresponding start,end,contig information with corresponding gene ID from the truncated/complete genbank file by yourself.
> Required files/parameters:
> 
> `name2gene`: It could accept a file similar to `example_data/format_anno_table/name2genes.txt`. You could get the required genes in anyways including hmmer /self designation/Experimental evidence.
>
> `input directory`: The directory comprises of all the genbank files which should be identical to the fasta file you passed to blast/alignment. Genbank files should end with suffix `gbk`. If the suffix of your files are `gbff` or others, you could pass corresponding suffix to the parameters `-s`. 
> \
> 

2. **truncate_genome_from_target.py**

> Required files/parameters:
> 
> `target_files`: It should contain information about the target genes you want to focus for each genomes. The example file is deposited at `example_data/truncate_genome_from_target/gene.infile`. For each genomes, you could pass multiple genes, but finally it only consider the left most and the right most genes within each contig/chromosome.
> `input directory`: The directory should contains `genbank` files of each genomes. The default suffix is `gbk`. If the suffix of your files are `gbff` or others, you could pass corresponding suffix to the parameters `-s`
> 
> `output directory`: Auto generate if it doesn't exist.
> 
> `fuzzy_match`: This parameters controls **the ways mapping** the requested genes name to the actual name within the genbank files. If you pass this parameter, you just need to pass the unique parts of the name. For example, you should pass `000123_123` to indicate the 123th gene of genome 000123. If you turn on the fuzzy_match, you just need to pass `_123` to indicate the 123th gene. It would be useful when you know the exact positions of the requested genes.
> 
> `num_p`: This parameter controls the length of the truncated regions. Default is `50p`. It will truncate the region covering 100 CDS(left+right) in total started from requested genes. You could also use `bp` to extract the region like `15e3bp` (15000bp). But using `bp` would not generate a exact length of region due to the consideration to the completeness of CDS. 



## Potential disasters

1. To perform blast or other alignment, it need to deal with the file name. And if there are any **'_to_'** in the orignal file name. It might raise errors. Details could be seen in the **nearly line 79** of `toolkit/alignment.py`
2.


