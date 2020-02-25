from collections import defaultdict
from glob import glob
from os.path import expanduser, abspath, exists, join, basename

from ete3 import Tree


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
    t = Tree(nwk_file, )
    return t


def get_files_from_dir(file,
                       indir=None,
                       how='from file'):
    indir = process_path(indir)
    if not exists(indir):
        raise IOError(f"indir {indir} doesn't exists")
    files_list = []
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
            files_list.append((name, file_path))
        elif len(f) > 1:
            print(f"Warning! Duplicated files for id {name}. ")
            pass
        else:
            pass
    return files_list


def get_link_info(indir, name2seq):
    link_dict = defaultdict(lambda: defaultdict(dict))
    linkfea_dict = defaultdict(dict)
    count = 1
    for blastout in glob(join(indir, '*.aliout')):
        g1, g2 = basename(blastout).split('_to_')
        g2 = g2.replace('.aliout', '')
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
