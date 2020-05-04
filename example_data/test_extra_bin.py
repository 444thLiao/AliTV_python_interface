import sys
from os.path import dirname, join

example_dir = dirname(__file__)
lib_dir = dirname(dirname(__file__))
sys.path.insert(0, lib_dir)


if __name__ == '__main__':
    print("testing main from extra_bin/format_anno_table.py    ")
    from extra_bin.format_anno_table import main as format_anno_table

    format_anno_table(infile=join(example_dir,"format_anno_table",'name2genes.txt'),
                      indir=join(example_dir,"split_gbk"),
                      ofile=join(example_dir,"format_anno_table","example_out.tab"))
    # python3 ./extra_bin/format_anno_table.py -i ./example_data/format_anno_table/name2genes.txt -indir ./example_data/split_gbk/ -o ./example_data/format_anno_table/example_out.tab

    print("testing main from extra_bin/truncate_genome_from_target.py    ")
    from extra_bin.truncate_genome_from_target import main as truncate_genome_from_target

    truncate_genome_from_target(infile=,
                                indir=,
                                odir=,
                                suffix=,
                                fuzzy_match=,
                                num_p=,
                                force=,)
    # python3 ./extra_bin/format_anno_table.py -i ./example_data/format_anno_table/name2genes.txt -indir ./example_data/split_gbk/ -o ./example_data/format_anno_table/example_out.tab

