#################################################################################
####  This script is mainly written as a interface to the AliTV v1.06
####
####  written by tianhua liao on 2020/02/17
#################################################################################
from os.path import join

import click
from tqdm import tqdm

from default_setting import json_obj_default
from toolkit import get_files_from_dir, alignment_batch, gname2seq_num, get_link_info, get_chrome_info, read_annotation_table, nwk2json, modify_json_from_config, IO_json




def main(genome_list=None,
         tree_file=None,
         config_file=None,
         indir='./',
         odir='./',
         annotation_table=None,
         enable_stepwise=True,
         force=True,
         alignment_ways='blastn',
         parallel=0,
         only_align=False,
         suffix='fna',
         pairwise=False,
         exact=False,
         ):
    # prepare the IO
    ali_odir = join(odir, 'tmp_ali_out')
    json_txt = None  # TREE
    # get a list of (name,file_path)
    if genome_list is None and tree_file is None and enable_stepwise:
        # no genome list file which indicate order
        # no tree file which also indicate order
        # enable stepwise is meaningless
        print('Warning, no designated order is given but enable stepwise alignment. ')
        order_files = get_files_from_dir(file=genome_list,
                                         indir=indir,
                                         suffix=suffix, )
        # it will return list of files with suffix in indir
    elif genome_list is not None:
        order_files = get_files_from_dir(file=genome_list,
                                         indir=indir,
                                         how='from file')
    elif tree_file is not None:
        order_files = get_files_from_dir(file=tree_file,
                                         indir=indir,
                                         how='from tree')
        order_names = [_[0] for _ in order_files]
        json_txt = nwk2json(tree_file, odir, subset_names=order_names)
    else:
        raise Exception('wrong')
    # run alignment and get info from blast output
    order_names = [_[0] for _ in order_files]
    files = [_[1] for _ in order_files]
    tqdm.write('running alignment')
    ali_how = 'pairwise' if pairwise else 'stepwise'
    used_seq = alignment_batch(genomes=files,
                               names=order_names,
                               odir=ali_odir,
                               alignment_ways=alignment_ways,
                               how=ali_how,
                               force=force,
                               parallel=parallel)
    # if only run the alignment program, it could stop here.
    # it will generate output directory and stodge the alignment results in it.
    if only_align:
        return
    # get information from blast results and raw gbk files
    name2seq = gname2seq_num(used_seq)
    link_dict, linkfea_dict = get_link_info(ali_odir,
                                            name2seq)

    chr_dict = get_chrome_info(files,
                               name2seq, )
    # if not annotations need to add
    json_obj = modify_json_from_config(json_obj_default, config_file).copy()
    json_obj['filters']['karyo']['chromosomes'] = {}
    for seq_name in name2seq.values():
        _dict = {seq_name: {"reverse": False,
                            "visible": True}}
        json_obj['filters']['karyo']['chromosomes'].update(_dict)

    json_obj['filters']['karyo']['genome_order'] = order_names
    json_obj['filters']['karyo']['order'] = list(name2seq.values())
    # data part
    json_obj['data']['features']['link'] = {}
    json_obj['data']['features']['link'].update(linkfea_dict)
    json_obj['data']['links'] = {}
    json_obj['data']['links'].update(link_dict)
    if json_txt is not None:
        json_obj['data']['tree'] = json_txt
    json_obj['data']['karyo']['chromosomes'] = chr_dict

    if annotation_table is not None:
        conf_fea, data_fea = read_annotation_table(annotation_table, name2seq)
        json_obj['data']['features'].update(data_fea)
        json_obj["conf"]['features']['supportedFeatures'].update(conf_fea)

    # auto conf part
    json_obj['conf']['graphicalParameters']['canvasHeight'] = len(order_names) * 120
    json_obj['conf']['graphicalParameters']['tickDistance'] = 10000

    IO_json(file=join(odir, 'aliTV_input.json'),
            way='w',
            json_obj=json_obj)


@click.command()
@click.option("-gl", "genome_list", default=None)
@click.option("-tf", "tree_file", default=None)
@click.option("-indir", "indir", default="./")
@click.option("-odir", "odir", default="./")
@click.option("-at", "annotation_table", default=None)
@click.option("-no-stepwise", "enable_stepwise", is_flag=True, default=True)
@click.option("-pairwise", "pairwise", is_flag=True, default=False)
@click.option("-align", "alignment_ways", default='blastn')
@click.option("-p", "parallel", default=0)
@click.option("-only_align", "only_align", is_flag=True, default=False)
@click.option("-s", "suffix", default='fna')
@click.option("-f", "force", is_flag=True, default=False)
@click.option("-exact", "exact", is_flag=True, default=False)
def cli(genome_list, tree_file, indir, odir, annotation_table, enable_stepwise, pairwise, force, alignment_ways, parallel, only_align, suffix, exact):
    if pairwise:
        enable_stepwise = False

    main(genome_list=genome_list,
         tree_file=tree_file,
         indir=indir,
         odir=odir,
         annotation_table=annotation_table,
         enable_stepwise=enable_stepwise,
         pairwise=pairwise,
         alignment_ways=alignment_ways,
         parallel=parallel,
         only_align=only_align,
         force=force,
         suffix=suffix,
         exact=exact)


if __name__ == '__main__':
    cli()

# python3 ~/software/AliTV_python_interface/main.py -tf ./species.newick -indir ./split_gbk/fna_dir -odir ./ali_odir/ -at ./annotation.tab
