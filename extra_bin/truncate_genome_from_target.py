import io
import os
from collections import defaultdict
from glob import glob
from os.path import *

import click
import gffutils
from Bio import SeqIO
from tqdm import tqdm


def get_information(infile):
    # format should be like
    # genome name\t target locus ID \n
    # the target locus ID could be multiple but separate each one with comma(,)
    # if you turn on the fuzzy match ( you could simple pass part of the name of locus )
    # But be careful for that ...
    # it only useful for the data with exact number but different prefix of locus.
    rows = [row.split('\t') for row in open(infile).read().split('\n') if row]
    genome2target_locus = {genome: [locus for locus in mlocus.split(',')]
                           for genome, mlocus in rows
                           }
    return genome2target_locus


# read gff part (redundant part)
def read_gff(gff_file, id_spec):
    tmp_dir = join(os.environ.get('HOME'), '.tmp')
    if not exists(tmp_dir):
        os.makedirs(tmp_dir)
    odb_file = join(tmp_dir, basename(gff_file) + '.db')
    if not exists(odb_file):
        db = gffutils.create_db(gff_file,
                                odb_file,
                                id_spec=id_spec,
                                force=True,
                                sort_attribute_values=True,
                                merge_strategy='merge')
    else:
        db = gffutils.FeatureDB(odb_file, keep_order=True, sort_attribute_values=True)
    return db


def get_all_gene_pos(genome_file, CDS_names=None, id_spec='ID'):
    """
    major part for reading gff file.
    1. giving a gff file.
    2. iterating all features and retrieve its id. (normally is locus ID)
    3. return gene2pos info like {locus1: {start:1,end:100,strand:'-',previous end:0}}
    3.5. and also record the contig order
    """
    if not exists(genome_file):
        tqdm.write('wrong file path[%s]... pass it' % genome_file)
        return None, None
    db = read_gff(genome_file, id_spec)
    if CDS_names is None:
        all_CDS = [_.id for _ in db.all_features()]
    else:
        all_CDS = CDS_names
    gene2pos = defaultdict(dict)
    order_contig_list = []
    last_end = None
    contig_list = []
    contig_names = []
    contig_name = ''
    for fea_id in all_CDS:
        fea = db[fea_id]
        if 'gene' not in fea.id:
            gene2pos[fea_id]['contig_name'] = fea.seqid

            if fea.seqid != contig_name:
                if contig_list:
                    order_contig_list.append(tuple(contig_list))
                    contig_names.append(fea.seqid)
                contig_list = []
                contig_name = fea.seqid
            gene2pos[fea_id]['strand'] = fea.strand
            gene2pos[fea_id]['start'] = fea.start
            gene2pos[fea_id]['end'] = fea.end
            gene2pos[fea_id]['previous end'] = last_end
            last_end = fea.end
            contig_list.append(fea_id)
    # get final contig
    order_contig_list.append(tuple(contig_list))
    contig_names.append(fea.seqid)
    contig2order_fea = dict(zip(contig_names,
                                tuple(order_contig_list)))
    return gene2pos, order_contig_list


# read gbk part
def read_gbk(gbk_file):
    # provide a better/robust way to read genbank file
    rows = open(gbk_file).readlines()
    rows = [process_locus_name_of_gbk(_)
            if '_length_' in _ and _.startswith('LOCUS') else _
            for _ in rows]
    records = SeqIO.parse(io.StringIO(''.join(rows)), format='genbank')
    records = list(records)
    return records


def process_locus_name_of_gbk(row):
    infos = row.split(' ')
    infos = [_ for _ in infos if _]
    if len(infos) <= 2:
        # print(infos)
        pass
    try:
        length = infos[1].split('_length_')[1].split('_')[0]
    except:
        raise Exception('')
        # import pdb;pdb.set_trace()
    name = 'LOCUS' + ' ' * 7 + infos[1].split('_length')[0] + f' {length} bp  ' + '  DNA  linear  20-Jan-2020'
    return name + '\n'


def get_all_CDS_from_gbk(gbk_file, tag='locus_tag'):
    """
    major part for reading gbk file.
    1. giving a gbk file.
    2. iterating all features and retrieve its id. (normally is locus ID)
    3. return gene2pos info like {locus1: {start:1,end:100,strand:'-',previous end:0}}
    3.5. and also record the contig order
    """
    if not exists(gbk_file):
        tqdm.write('wrong file path[%s]... pass it' % gbk_file)
        return None, None
    contigs = read_gbk(gbk_file)
    order_contig_list = []

    gene2pos = defaultdict(dict)
    last_end = None
    for contig in contigs:
        contig_list = []
        contig_name = contig.id
        all_cds = [fea for fea in contig.features if fea.type == 'CDS']
        for fea in all_cds:
            if tag in fea.qualifiers:
                fea_id = fea.qualifiers[tag][0]
            else:
                fea_id = fea.qualifiers['locus_tag'][0]
            gene2pos[fea_id]['contig_name'] = contig_name
            gene2pos[fea_id]['strand'] = '+' if fea.strand == 1 else '-'
            gene2pos[fea_id]['start'] = int(fea.location.start)
            gene2pos[fea_id]['end'] = int(fea.location.end)
            gene2pos[fea_id]['previous end'] = last_end
            last_end = int(fea.location.end)
            contig_list.append(fea_id)
        order_contig_list.append(tuple(contig_list))
    contig2order_fea = dict(zip([contig.id for contig in contigs],
                                tuple(order_contig_list)))
    gene2pos = dict(gene2pos)
    return gene2pos, contig2order_fea


# start to split
def split_gbk(genome_files, odir,
              target_gene_dict, num_p,
              suffix, force=False, fuzzy_match=False):
    if not exists(odir):
        os.makedirs(odir)
    for genome_file in tqdm(genome_files):
        genome_name = basename(genome_file).rpartition('.')[0]
        ofile = join(odir, f'{genome_name}.{suffix}')


        # get target genes
        target_genes = target_gene_dict.get(genome_name, None)
        if target_genes is None:
            tqdm.write(f"{genome_name} doesn't have target_gene, please check it")
            continue

        # read genomic information of proteins/CDS
        if exists(genome_file) and (not exists(ofile) or force):
            records = read_gbk(genome_file)
            if genome_file.endswith('.gbk'):
                genome2gene_info, contig2order_fea = get_all_CDS_from_gbk(genome_file)
                g2info = genome2gene_info
            else:
                # only use genbank. not implement others
                raise Exception('')
        elif exists(ofile) and not force:
            tqdm.write("output file existed, it will bypass it unless you use force")
            continue
        else:
            tqdm.write("no validated genome", genome_file)
            continue

        if fuzzy_match:
            # if using fuzzy match
            # renew the target_genes
            target_genes = [right_CDS
                            for _target_gene in target_genes
                            for right_CDS in g2info
                            if _target_gene in right_CDS]
        else:
            # confirm all retrieved target_genes are placed in g2info.
            assert all([True if _target_gene in g2info else False
                        for _target_gene in target_genes])

        # get centre position of target genes, if there are multiple genes,

        contig2pos_remained = defaultdict(list)
        for gene in target_genes:
            pos1 = g2info[gene]['start']
            pos2 = g2info[gene]['end']
            contig = g2info[gene]['contig_name']
            contig2pos_remained[contig].append((int(pos1), int(pos2), gene))
        # below is the part to extract designated regions
        if num_p[1] == 'bp':
            expand_len = num_p[0]
            # filter out some contig which extend too long, like expand over 2 times of expand_len
            contigs_list = []
            for _contig, pos_list in list(contig2pos_remained.items()):
                min_pos = min(pos_list, key=lambda x: x[0])[0]
                max_pos = max(pos_list, key=lambda x: x[1])[1]

                if (max_pos - min_pos) >= 2 * expand_len:
                    contig2pos_remained.pop(_contig)
                    continue

                start = int(min_pos - expand_len)
                start = 0 if start < 0 else start
                end = int(max_pos + expand_len)

                contig_obj = [record
                              for record in records
                              if record.id == _contig]
                if contig_obj:
                    contig_obj = contig_obj[0]
                else:
                    continue
                subset_contig = contig_obj[start:end]
                contigs_list.append(subset_contig)

            if contigs_list:
                with open(ofile, 'w') as f1:
                    SeqIO.write(contigs_list, f1, format='genbank')

        elif num_p[1] == 'CDS':
            expand_CDS = num_p[0]
            contigs_list = []
            for contig, info in list(contig2pos_remained.items()):
                left_most_cds = min(info, key=lambda x: x[0])[-1]
                right_most_cds = max(info, key=lambda x: x[0])[-1]

                order_fea = contig2order_fea[contig]
                left_most_idx = order_fea.index(left_most_cds)
                right_most_idx = order_fea.index(right_most_cds)
                left_bound_idx = 0 if left_most_idx - expand_CDS <= 0 else left_most_idx - expand_CDS
                right_bound_idx = len(order_fea) - 1 if right_most_idx + expand_CDS >= len(order_fea) else right_most_idx + expand_CDS

                pos_end = g2info[order_fea[right_bound_idx]]['end']
                pos_start = g2info[order_fea[left_bound_idx]]['start']

                contig_obj = [record
                              for record in records
                              if record.id == contig]
                if contig_obj:
                    contig_obj = contig_obj[0]
                else:
                    continue
                subset_contig = contig_obj[pos_start:pos_end]
                contigs_list.append(subset_contig)
            with open(ofile, 'w') as f1:
                SeqIO.write(contigs_list, f1, format='genbank')
        else:
            raise Exception('unknown num_p parameters')
    # convert them into fasta in order to blast further
    # it will create sub directory called fna under the output directory.
    for gbk in tqdm(glob(join(odir, '*.gbk'))):
        a = SeqIO.parse(gbk, format='genbank')
        base_gbk = basename(gbk)
        if not exists(join(odir, 'fna')):
            os.makedirs(join(odir, 'fna'))
        with open(join(odir,
                       'fna',
                       base_gbk.replace('.gbk', '.fna')), 'w') as f1:
            SeqIO.write(a, f1, format='fasta')


def parse_nump(num_p):
    num_p = num_p.lower() if type(num_p) == str else None
    if num_p.endswith('bp'):
        return (int(float(num_p.replace('bp', ''))),
                'bp')
    elif num_p.endswith('p'):
        return (int(num_p.strip('p')),
                'CDS')
    else:
        raise Exception('unknown num_p parameters,must end with bp or p.')


def main(infile, indir, odir, suffix, fuzzy_match, num_p, force):
    all_gbk = glob(join(indir, f'*.{suffix}'))
    genome2target_locus = get_information(infile)
    num_p = parse_nump(num_p)
    split_gbk(genome_files=all_gbk,
              odir=odir,
              suffix=suffix,
              target_gene_dict=genome2target_locus,
              num_p=num_p,
              force=force,
              fuzzy_match=fuzzy_match)


@click.command()
@click.option("-i", "infile")
@click.option("-indir", "indir")
@click.option("-odir", "odir")
@click.option("-s", "suffix", default='gbk')
@click.option("-r", "num_p",
              help='maximum number of locus you want to truncate. Normally it just a approximate number. '
                   'You could use 50p to mention 50 protein/CDS or 15e3bp to mention 15000bp ',
              default='50p')
@click.option("--fuzzy-match", "fuzzy_match", is_flag=True, default=False)
@click.option("-f", "force", is_flag=True, default=False)
def cli(infile, indir, odir, suffix, fuzzy_match, num_p, force):
    main(infile, indir, odir, suffix, fuzzy_match, num_p, force)


if __name__ == '__main__':
    cli()

    # python3 ~/software/AliTV_python_interface/extra_bin/truncate_genome_from_target.py -i ./test.infile -indir ./raw_gbk -odir ./split_gbk -f -r 50p
