import os
from collections import defaultdict
from glob import glob
from os.path import expanduser, abspath, exists, join, basename, dirname

import pandas as pd
from Bio import SeqIO
from ete3 import Tree
from tqdm import tqdm

default_colors = ['#2E91E5', '#E15F99', '#1CA71C', '#FB0D0D', '#DA16FF', '#222A2A', '#B68100', '#750D86',
                  '#EB663B', '#511CFB', '#00A08B', '#FB00D1', '#FC0080', '#B2828D',
                  '#6C7C32', '#778AAE', '#862A16', '#A777F1', '#620042', '#1616A7', '#DA60CA', '#6C4516', '#0D2A63', '#AF0038']


def process_path(path):
    if path is None:
        return None
    if not '/' in path and ':' not in path:
        # : for windows, / for linux
        path = './' + path
    if path.startswith('~'):
        path = expanduser(path)
    if path.startswith('.'):
        path = abspath(path)
    return path


def read_tree(nwk_file):
    """
    0	flexible with support values	((D:0.723274,F:0.567784)1.000000:0.067192,(B:0.279326,H:0.756049)1.000000:0.807788);
    1	flexible with internal node names	((D:0.723274,F:0.567784)E:0.067192,(B:0.279326,H:0.756049)B:0.807788);
    2	all branches + leaf names + internal supports	((D:0.723274,F:0.567784)1.000000:0.067192,(B:0.279326,H:0.756049)1.000000:0.807788);
    3	all branches + all names	((D:0.723274,F:0.567784)E:0.067192,(B:0.279326,H:0.756049)B:0.807788);
    4	leaf branches + leaf names	((D:0.723274,F:0.567784),(B:0.279326,H:0.756049));
    5	internal and leaf branches + leaf names	((D:0.723274,F:0.567784):0.067192,(B:0.279326,H:0.756049):0.807788);
    6	internal branches + leaf names	((D,F):0.067192,(B,H):0.807788);
    7	leaf branches + all names	((D:0.723274,F:0.567784)E,(B:0.279326,H:0.756049)B);
    8	all names	((D,F)E,(B,H)B);
    9	leaf names	((D,F),(B,H));
    100	topology only	((,),(,));
    :param nwk_file:
    :return:
    """
    for f in [0, 1, 3, 5, 8]:
        try:
            t = Tree(nwk_file, format=f)
            # ignore the distance, only get the topology
            nt = Tree(t.write(format=8), format=8)
            return nt
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
                if "translation" not in seq_feature.qualifiers:
                    continue
                if len(seq_feature.qualifiers['translation']) != 1:
                    continue
                faa_handle.write(">%s from %s\n%s\n" % (
                    seq_feature.qualifiers['locus_tag'][0],
                    seq_record.name,
                    seq_feature.qualifiers['translation'][0]))
    faa_handle.close()


def get_protein_info_from_gbk(gbk, protein):
    for seq_record in SeqIO.parse(gbk, format="genbank"):
        for seq_feature in seq_record.features:
            if seq_feature.type == "CDS" and protein in seq_feature.qualifiers['locus_tag']:
                start, end = seq_feature.location.start.real, seq_feature.location.end.real
                start = str(start)
                end = str(end)
                return seq_record.id, start, end


def get_files_from_dir(file=None,
                       indir=None,
                       suffix='fna',
                       how='from file'):
    """

    :param file: a file list all genomes need to be processed. the order is important
    :param indir:
    :param suffix:
    :param how:
    :return:
    return object is a list of tuple which stodge from the genome name to the path of file.
    Not using dict is because the name could be duplicated.

    """
    indir = process_path(indir)
    if not exists(indir):
        raise IOError(f"indir '{indir}' doesn't exist")
    if file is None:
        file_list = glob(join(indir,f'*.{suffix}'))
        file_list = sorted(file_list)
        file_list = list(zip([basename(f).rpartition('.')[0] for f in file_list],
                             file_list))
        return file_list,[_[0] for _ in file_list]
    file_list = []
    if how == 'from file':
        names = list(open(file).readlines())
        names = [_ for _ in names if _]
    elif how == 'from tree':
        t = read_tree(file)
        names = t.get_leaf_names()
        tqdm.write(f"In given tree, {len(names)} were presented ")
    else:
        return
    for name in names:
        f = glob(join(indir, f'{name}.*'))
        if len(f) == 1:
            file_path = process_path(f[0])
            file_list.append((name, file_path))
        elif len(f) > 1:
            print(f"Warning! Duplicated files for id {name}. It will be passed. please check the genome or tree file carefully.")
            pass
        else:
            continue
    tqdm.write(f"Finally, {len(file_list)} were retrieved ")
    return file_list


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
    If the format is not belong to the format 6, it might have some silent errors.
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


def get_chrome_info(genome_files, name2seq, stodge_seq=False):
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
    the table should follow below format
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
        sub_df = df.loc[sub_df_index, :]
        data_fea[gene] = []
        for idx, row in sub_df.iterrows():
            format_names = [f"{row[0]}_{row[1]}",
                            f"{row[0]}",
                            f"{row[0]}_{row[1].split('.')[0]}"]

            _t = [format_name
                  for format_name in format_names
                  if format_name in name2seq]
            if _t:
                format_name = _t[0]
                _dict = dict(start=row[2],
                             end=row[3],
                             name=gene,
                             karyo=name2seq[format_name])
            data_fea[gene].append(_dict)
    return conf_fea, data_fea


def deep_scan(current_node, root=None):
    """
    It implement a function which could generate a suitable tree hierarchical format to AliTV.

    :param current_node: it should be the node/root which from etetools.
    :param root: normally, it should not be used.
    :return:
    """
    # aligned leafs
    if root is None:
        root = current_node
    _, max_depth = root.get_farthest_leaf()
    if current_node.is_leaf():
        ad = max_depth - root.get_distance(current_node)
        # abbrev. across depth
        _dict = {"name": current_node.name}
        for _ in range(int(ad)):
            _dict = {"children": [_dict]}
        return _dict
    else:
        return_v = []
        child = current_node.children
        for c in child:
            _re = deep_scan(c, root=root)
            if type(_re) == dict:
                _re = [_re]
            return_v.append({"children": _re})
        return return_v


def nwk2json(treefile, odir, subset_names=[]):
    """
    With above deep_scan function, it could convert a newick plaintext file into a json format which is suitable for AliTV.
    :param treefile: path/ete3.Tree object
    :param odir: output directory instead of output file... it just want to fix the output file name.
    :param subset_names: If you pass parts of leafs to it, it will help you to truncated.
    :return:
    """
    # following the alitv way...
    t = read_tree(treefile)
    if subset_names:
        t.prune(subset_names)
        text = t.write()
        if not exists(odir):
            os.makedirs(odir)
        with open(join(odir, 'used.newick'), 'w') as f1:
            f1.write(text)
    json_v = {"children": deep_scan(t)}
    return json_v
