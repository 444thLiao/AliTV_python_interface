"""
This script is for alignment (stepwise or pairwise)

"""
import itertools
import multiprocessing as mp
import os
from os.path import *
from subprocess import check_call

from tqdm import tqdm

from toolkit import process_path, gbk2fna


def run(cmd):
    check_call(cmd,
               shell=True,
               stdout=open('/dev/null','w'),
               stderr=open('/dev/null','w'))


def gbk2fna_path(f):
    return join(dirname(f),
                'fna',
                basename(f).replace('gbk', 'fna'))


def align_unit(f1, f2, ofile, method='blastn',
               force=True, parallel=0, extra_params=''):
    if method == 'blastn':
        if force or not exists(ofile):
            blastn_path = os.popen('which blastn').read().strip('\n')
            cmd = f"{blastn_path} -query {f1} -subject {f2} -evalue 1e-3 -outfmt 6 -out {ofile} " + extra_params
            if parallel == 0:
                os.system(cmd)
            else:
                return cmd
        return ''
    # elif method is None:
    #     pass
    # else:
    #     #todo: implement other software
    #     pass


def alignment_batch(genomes,
                    names,
                    odir,
                    alignment_ways='blastn',
                    parallel=0,
                    force=True,
                    how='stepwise'):
    new_genomes = []
    for g in genomes:
        # if genbank file was passed, then it will convert it into fna first.
        if g.endswith('gbk') and not exists(gbk2fna_path(g)):
            g = gbk2fna(g,
                        gbk2fna_path(g))
        elif exists(gbk2fna_path(g)):
            g = gbk2fna_path(g)
        new_genomes.append(g)
    genomes = new_genomes[::] # use [::] to avoid potential overwritten errors

    g2name = dict(zip(genomes, names))
    odir = process_path(odir)
    if not exists(odir):
        os.makedirs(odir, exist_ok=1)
    ali_record_file = join(odir, 'align_record.txt')
    f1 = open(ali_record_file, 'w')

    if how == 'stepwise':
        iter_objs = list(zip(genomes[:-1], genomes[1:]))
    elif how == 'pairwise':
        # if the way to perform alignment is pairwise, it will generate a permutation set of given genomes
        # to ABCDE genomes, it will generate both (A,B) and (B,A). But it will not generate a set which align to itself.
        iter_objs = list(itertools.permutations(genomes, 2))
    else:
        # not implement yet
        return
    iter_objs = tqdm(iter_objs) if parallel == 0 else iter_objs
    tqdm.write(f"start running the '{how}' alignment " + "based on designated order." if how == 'stepwise' else ' .')
    cmds = []
    for file1, file2 in iter_objs:
        name1 = g2name[file1]
        name2 = g2name[file2]
        ofile = join(odir, name1 + '_to_' + name2 + '.aliout')
        f1.write('\t'.join([file1, file2, name1, name2, ofile, '\n']))
        # in some special case there is '_to_' in file name... it might happen... when it happen, you should need to know
        # it records some file
        cmds.append(align_unit(file1, file2, ofile,
                               method=alignment_ways,
                               force=force,
                               parallel=parallel))
        # if parallel is setting 0, it mean run one by one.
    if parallel != 0:
        parallel = len(cmds) if parallel == -1 else int(parallel)
        with mp.Pool(processes=parallel) as tp:
            r = list(tqdm(tp.imap(run, cmds), total=len(cmds)))
    f1.close()

    return genomes
