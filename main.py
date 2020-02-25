#################################################################################
####  This script is mainly written as a interface to the AliTV v1.06
####
####  written by tianhua liao on 2020/02/17
#################################################################################
from os.path import join

import click

from toolkit import get_files_from_dir, alignment_batch


def main(genome_list=None,
         tree_file=None,
         indir=None,
         disable_tree=False,
         odir=None,
         enable_stepwise=True,
         force=True,
         alignment_ways='blastn',
         parallel=0,
         ):
    # prepare IO
    ali_odir = join(odir, 'tmp_ali_out')

    # get a list of (name,file_path)
    if genome_list is None and tree_file is None and enable_stepwise:
        print('Warning, no designated order is given but enable stepwise alignment. ')
        order_files = get_files_from_dir(genome_list,
                                         indir=indir,
                                         how='from file')
    elif genome_list is not None:
        order_files = get_files_from_dir(genome_list,
                                         indir=indir,
                                         how='from file')
    elif tree_file is not None:
        order_files = get_files_from_dir(tree_file,
                                         indir=indir,
                                         how='from tree')
    else:
        raise Exception('wrong')
    # run alignment and get info from blast infor
    names = [_[0] for _ in order_files]
    files = [_[1] for _ in order_files]
    alignment_batch(files,
                    names=names,
                    odir=ali_odir,
                    alignment_ways=alignment_ways,
                    how='stepwise' if enable_stepwise else None,
                    force=force,
                    parallel=parallel)
    blast_info = 'a'
    # run alignment


@click.command()
def cli():
    pass


if __name__ == '__main__':
    cli()
