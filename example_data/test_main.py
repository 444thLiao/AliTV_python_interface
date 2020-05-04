import sys
from os.path import dirname, join

example_dir = dirname(__file__)
lib_dir = dirname(dirname(__file__))
sys.path.insert(0, lib_dir)
from main import main

if __name__ == '__main__':
    main(genome_list=None,
         tree_file=join(example_dir, "species.newick"),
         indir=join(example_dir, "split_gbk", "fna_dir"),
         odir=join(example_dir, "ali_odir"),
         annotation_table=join(example_dir, "annotation.tab"),
         enable_stepwise=True,
         alignment_ways="blastn",
         parallel=0,
         only_align=False,
         force=False,
         suffix='fna',
         exact=False)

# python3 ~/software/AliTV_python_interface/extra_bin/format_anno_table.py -i ../gene_trees/name2genes.txt -indir ./split_gbk/ -o ./annotation.tab
