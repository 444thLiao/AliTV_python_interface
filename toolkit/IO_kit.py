import os
from collections import defaultdict
from glob import glob
from os.path import expanduser, abspath, exists, join, basename, dirname

import pandas as pd
from Bio import SeqIO
from ete3 import Tree,TreeNode
from tqdm import tqdm

default_colors = ['#2E91E5', '#E15F99', '#1CA71C', '#FB0D0D', '#DA16FF', '#222A2A', '#B68100', '#750D86', '#EB663B', '#511CFB', '#00A08B', '#FB00D1', '#FC0080', '#B2828D', '#6C7C32', '#778AAE', '#862A16', '#A777F1', '#620042', '#1616A7', '#DA60CA', '#6C4516', '#0D2A63', '#AF0038']


def process_path(path):
    if path is None:
        return None
    if not '/' in path:
        path = './' + path
    if path.startswith('~'):
        path = expanduser(path)
    if path.startswith('.'):
        path = abspath(path)
    return path


def read_tree(nwk_file):
    for f in [0, 1, 3, 5, 8]:
        try:
            t = Tree(nwk_file, format=f)
            return t
        except:
            pass


def gbk2fna(gbk, fna):
    fna = process_path(fna)
    if not exists(dirname(fna)):
        os.makedirs(dirname(fna))
    records = SeqIO.parse(gbk, format='genbank')
    with open(fna, 'w') as f1:
        SeqIO.write(records, f1, format='fasta-2line')
    return fna


def gbk2faa(gbk, faa):
    faa_handle = open(faa, 'w')
    for seq_record in SeqIO.parse(gbk, format="genbank"):
        for seq_feature in seq_record.features:
            if seq_feature.type == "CDS":
                assert len(seq_feature.qualifiers['translation']) == 1
                faa_handle.write(">%s from %s\n%s\n" % (
                    seq_feature.qualifiers['locus_tag'][0],
                    seq_record.name,
                    seq_feature.qualifiers['translation'][0]))
    faa_handle.close()


def get_protein_info_from_gbk(gbk, protein):
    for seq_record in SeqIO.parse(gbk, format="genbank"):
        for seq_feature in seq_record.features:
            if seq_feature.type == "CDS" and protein in seq_feature.qualifiers['locus_tag']:
                s, e = seq_feature.location.start.real, seq_feature.location.end.real
                s = str(s)
                e = str(e)
                return seq_record.id, s, e


def get_files_from_dir(file=None,
                       indir=None,
                       suffix='fna',
                       how='from file'):
    """

    :param file:
    :param indir:
    :param suffix:
    :param how:
    :return:
    """
    indir = process_path(indir)
    if not exists(indir):
        raise IOError(f"indir {indir} doesn't exists")
    if file is None:
        file_list = glob(join(indir,f'*.{suffix}'))
        file_list = sorted(file_list)
        file_list = list(zip([basename(f).rpartition('.')[0] for f in file_list],
                             file_list))
        return file_list,[_[0] for _ in file_list]
    file_list = []
    if how == 'from file':
        names = list(open(file).readlines())
    elif how == 'from tree':
        t = read_tree(file)
        names = t.get_leaf_names()
    else:
        return
    for name in names:
        f = glob(join(indir, f'{name}*'))
        if len(f) == 1:
            file_path = process_path(f[0])
            file_list.append((name, file_path))
        elif len(f) > 1:
            print(f"Warning! Duplicated files for id {name}. It will be passed")
            pass
        else:
            pass
    return file_list, [_[0] for _ in file_list]


def to_name2seq_num(genomes):
    # in fna format
    name2seq = {}
    c = 0
    for g in genomes:
        rn = basename(g).replace('.fna', '')
        records = SeqIO.parse(g, format='fasta')
        for record in records:
            c_id = record.id
            name2seq[rn + '_' + c_id] = f"seq{c}"
            c += 1
    return name2seq


def get_link_info(indir, name2seq, suffix='aliout'):
    """
    The file which stodge link info is for blastn (format 6)
    If the
    :param indir:
    :param name2seq:
    :return: return two objects. first one is using the raw genome name. the latter is using translated features name.
    """
    link_dict = defaultdict(lambda: defaultdict(dict))
    linkfea_dict = defaultdict(dict)
    count = 1

    for blastout in glob(join(indir, f'*.{suffix}')):
        g1, g2 = basename(blastout).split('_to_')
        g2 = g2.replace(f'.{suffix}', '')

        for row in open(blastout):
            rows = row.split('\t')
            qid, sid, iden, qstart, qend, sstart, send = [rows[_]
                                                          for _ in [0, 1, 2, 6, 7, 8, 9]]
            fea1 = 'linkfeature' + '{:0>8}'.format(len(linkfea_dict) + 1)
            linkfea_dict[fea1]["end"] = int(qend)
            linkfea_dict[fea1]["start"] = int(qstart)
            linkfea_dict[fea1]["karyo"] = name2seq[g1 + '_' + qid]
            fea2 = 'linkfeature' + '{:0>8}'.format(len(linkfea_dict) + 1)
            linkfea_dict[fea2]["end"] = int(send)
            linkfea_dict[fea2]["start"] = int(sstart)
            linkfea_dict[fea2]["karyo"] = name2seq[g2 + '_' + sid]

            _dict = {'link' + '{:0>8}'.format(count):
                         {'identity': float(iden),
                          'source': fea1,
                          'target': fea2}
                     }
            link_dict[g1][g2].update(_dict)
            link_dict[g2][g1].update(_dict)
            count += 1
    return link_dict, linkfea_dict


def get_chrome_info(genome_files, name2seq,stodge_seq=False):
    _dict = {}
    for gf in tqdm(genome_files):
        if gf.endswith('gbk'):
            records = SeqIO.parse(gf, format='genbank')
            for record in records:
                gn = basename(gf).replace('.gbk', '')
                id_name = name2seq[gn + '_' + record.id]
                _dict[id_name] = {}
                _dict[id_name]['genome_id'] = gn
                _dict[id_name]['length'] = len(record)
                _dict[id_name]['name'] = id_name
                _dict[id_name]['seq'] = '' if not stodge_seq else str(record.seq)
        else:
            try:
                records = SeqIO.parse(gf, format='fasta')
            except:
                raise Exception('unknown input... not genbank with gbk as suffix. not fasta file')
            for record in records:
                gn = basename(gf).rpartition('.')[0]
                id_name = name2seq[gn + '_' + record.id]
                _dict[id_name] = {}
                _dict[id_name]['genome_id'] = gn
                _dict[id_name]['length'] = len(record)
                _dict[id_name]['name'] = id_name
                _dict[id_name]['seq'] = '' if not stodge_seq else str(record.seq)
    return _dict


def read_annotation_table(f, name2seq):
    """
    read from a table with annotation information
    the table should like
    1. no header
    2. separator is tab(\t)
    3. should be
       genome name {\t} contig id {\t} start {\t} end {\t} annotation name (like some genes)
    :param f:
    :param name2seq:
    :return:
    """
    conf_fea = defaultdict(dict)
    data_fea = {}
    df = pd.read_csv(f, sep='\t', header=None)
    gb = df.groupby(4)
    for idx, gene in enumerate(list(gb.groups.keys())):
        conf_fea[gene]['color'] = default_colors[idx]
        conf_fea[gene]['form'] = "rect"
        conf_fea[gene]['height'] = "30"
        conf_fea[gene]['visible'] = True
    for gene, sub_df_index in gb.groups.items():
        sub_df = df.loc[sub_df_index,:]
        data_fea[gene] = []
        for idx, row in sub_df.iterrows():
            _dict = dict(start=row[2],
                         end=row[3],
                         name=gene,
                         karyo=name2seq[f"{row[0]}_{row[1]}"])
            data_fea[gene].append(_dict)
    return conf_fea, data_fea

def deep_scan(tree):
    if tree.is_leaf():
        return [{'children':[{"name":tree.name}]}]
    else:
        return_v = []
        child = tree.children
        for c in child:
            return_v.append({"children":deep_scan(c)})
        return return_v

def nwk2json(treefile,odir,subset_names=[]):
    # following the alitv way...
    t = read_tree(treefile)
    if subset_names:
        t.prune(subset_names)
        text = t.write()
        if not exists(odir):
            os.makedirs(odir)
        with open(join(odir,'used.newick'),'w') as f1:
            f1.write(text)
    json_v = deep_scan(t)
    return json_v
