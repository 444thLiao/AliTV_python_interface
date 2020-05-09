#################################################################################
####  This script is mainly written as a interface to the AliTV v1.06
####
####  written by tianhua liao on 2020/02/17
#################################################################################
from os.path import join

import click
from tqdm import tqdm

from setting import json_obj_default
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
    ali_how = 'pairwise' if not enable_stepwise else 'stepwise'
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
    # sto the name2seq
    with open(join(ali_odir, '..', 'genome_name2seq_num.txt'), 'w') as f1:
        f1.write('\n'.join([f"{k}\t{v}"
                            for k, v in name2seq.items()]))

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
@click.option("-gl", "genome_list", default=None, help='A file comprises of order genome names, without suffix. Default None')
@click.option("-tf", "tree_file", default=None, help="Newick formatted tree file. The order of the leaves will be taken as the order of finall visualization")
@click.option("-c", "config_file", default=None,
              help="configuration file. Default None. You could pass some preset parameters with it. Or you could change parameters after the generation first json using `extra_bin/change_parameters.py `")
@click.option("-indir", "indir", default="./", help="input directory which contain fasta files for alignment")
@click.option("-odir", "odir", default="./", help="output directory which contains temporary alignment files and final json.")
@click.option("-at", "annotation_table", default=None,
              help="annotation file which comprises of information of annotated genes. The required format could refer to the output of `./extra_bin/format_anno_table.py` ")
@click.option("-no-stepwise", "enable_stepwise", is_flag=True, default=True, help="pairwise alignment")
@click.option("-align", "alignment_ways", default='blastn', help="alignment ways. Default is blastn for nucleotide sequences.")
@click.option("-p", "parallel", default=0, help="the number of tasks you want to parallel run. Default is 0. It means no parallel running. ")
@click.option("-only_align", "only_align", is_flag=True, default=False,
              help="only perform the alignment instead of generation of json file. It would stop after the alignment completed. ")
@click.option("-s", "suffix", default='fna', help="The suffix of nucleotide sequence files under `indir` ")
@click.option("-f", "force", is_flag=True, default=False, help="Re-run or not if the alignment results have been existed.")
def cli(genome_list, tree_file, indir, odir, config_file,
        annotation_table, enable_stepwise, force, alignment_ways,
        parallel, only_align, suffix):
    main(genome_list=genome_list,
         config_file=config_file,
         tree_file=tree_file,
         indir=indir,
         odir=odir,
         annotation_table=annotation_table,
         enable_stepwise=enable_stepwise,
         alignment_ways=alignment_ways,
         parallel=parallel,
         only_align=only_align,
         force=force,
         suffix=suffix)


if __name__ == '__main__':
    cli()

# python3 ~/software/AliTV_python_interface/main.py -tf ./species.newick -indir ./split_gbk/fna_dir -odir ./ali_odir/ -at ./annotation.tab
