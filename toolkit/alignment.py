import itertools
import multiprocessing as mp
import os
from os.path import *
from subprocess import check_call

from tqdm import tqdm

from toolkit import process_path


def run(cmd):
    check_call(cmd,
               shell=True)


def align_unit(f1, f2, ofile, method='blastn',
               force=True, parallel=0):
    if method == 'blastn':
        if force or not exists(ofile):
            cmd = f"blastn -query {f1} -subject {f2} -evalue 1e-3 -outfmt 6 -out {ofile}"
            if parallel == 0:
                os.system(cmd)
            else:
                return cmd
    # elif method is None:
    #     pass
    # else:
    #     #todo
    #     pass


def alignment_batch(genomes,
                    names,
                    odir,
                    alignment_ways='blastn',
                    parallel=0,
                    force=True,
                    how='stepwise'):
    g2name = dict(zip(genomes,names))
    odir = process_path(odir)
    if not exists(odir):
        os.makedirs(odir, exist_ok=1)
    ali_record_file = join(odir, 'align_record.txt')
    f1 = open(ali_record_file, 'w')

    if how == 'stepwise':
        iter_objs = list(zip(genomes[:-1], genomes[1:]))
    elif how == 'pairwise':
        iter_objs = list(itertools.combinations(genomes, 2))
    else:
        # not implement yet
        return
    iter_objs = tqdm(iter_objs) if parallel == 0 else iter_objs
    tqdm.write(f"start running the {how} alignment based on designated order.")

    cmds = []
    for file1, file2 in iter_objs:
        name1 = g2name[file1]
        name2 = g2name[file2]
        ofile = join(odir, name1 + '_to_' + name2 + '.aliout')
        f1.write('\t'.join([file1, file2, name1, name2, ofile]))
        # in case there is '_to_' in file name... it might happen...
        # it records some file
        cmds.append(align_unit(file1, file2, ofile,
                               method=alignment_ways,
                               force=force,
                               parallel=parallel))
        # if parallel is default 0, it mean run one by one.
    if parallel != 0:
        parallel = len(cmds) if parallel == -1 else int(parallel)
        with mp.Pool(processes=parallel) as tp:
            r = list(tqdm(tp.imap(run, cmds), total=len(cmds)))
    f1.close()
