"""
This script is mainly for quickly mapping locus id into detailed information

Accept a file with gene locus and its annotated name AND a directory which stodge genbank files which could be used to map.
"""

import sys
from glob import glob
from os.path import *

import click
from tqdm import tqdm

sys.path.insert(0, dirname(dirname(__file__)))
from extra_bin.truncate_genome_from_target import get_all_CDS_from_gbk


def parsed_infile(infile,indir):
    locus2gene = {}
    locus2genome = {}
    gbk2gbk_obj = {}
    f = open(infile)
    for row in tqdm(f):
        row = row.strip('\n')
        objs = row.split('\t')
        if len(objs) == 2:
            gene, locus  = objs
            locus2gene[locus] = gene

        elif len(objs) == 3:
            gene, locus , gbk_file = objs

            locus2gene[locus] = gene
            if not exists(gbk_file):
                gbk_file = join(indir,gbk_file)
            if not exists(gbk_file):
                tqdm.write(f"wrong input gbk name {gbk_file}")
            gbk_name = basename(gbk_file).rpartition('.')[0]
            locus2genome[locus] = gbk_name
            gbk2gbk_obj[gbk_name] = get_all_CDS_from_gbk(gbk_file)[0]

    return locus2gene, locus2genome, gbk2gbk_obj


def core(locus2gene, locus2genome, gbk2gbk_obj,
         ofile):
    rows = []
    for locus, gene in locus2gene.items():
        genome = locus2genome[locus]
        gbk_obj = gbk2gbk_obj.get(genome,{}).get(locus)
        if gbk_obj is None:
            continue
        contig_name = gbk_obj['contig_name']
        start = str(int(gbk_obj['start']))
        end = str(int(gbk_obj['end']))
        rows.append('\t'.join([genome,
                               contig_name, start, end,
                               gene]))
    with open(ofile, 'w') as f1:
        f1.write("\n".join(rows))


def find_f(genomes, indir, suffix='gbk'):
    genome2path = {}
    for genome in genomes:
        files = glob(join(indir, f"{genome}.{suffix}"))
        if len(files) >= 2:
            print("detect multiple files... may get wrong")
        genome2path[genome] = get_all_CDS_from_gbk(files[0])
    return genome2path


def main(infile, indir, ofile, suffix='gbk'):
    locus2gene, locus2genome, gbk2gbk_obj = parsed_infile(infile,indir)
    if not gbk2gbk_obj and indir is not None:
        tqdm.write(f"search in {indir}")
        genomes = set([g for l, g in locus2genome.items()])
        gbk2gbk_obj = find_f(genomes, indir, suffix=suffix)
    core(locus2gene, locus2genome, gbk2gbk_obj, ofile)


@click.command()
@click.option("-i", "infile")
@click.option("-indir", "indir", default=None)
@click.option("-o", "ofile")
@click.option("-s", "suffix", default='gbk')
def cli(infile, indir, ofile, suffix):
    suffix = suffix.strip('.')
    main(infile, indir, ofile, suffix)


if __name__ == '__main__':
    cli()

# python3 ~/software/AliTV_python_interface/extra_bin/format_anno_table.py -i ../gene_trees/name2genes.txt -indir ./split_gbk/ -o ./annotation.tab
