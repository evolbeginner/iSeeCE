#! /bin/env python

import sys
import re
import os
import getopt


########################################################################
pair_files = []
blast_files = []
gff_files = []
gene_objsh = {}
distance = 100000
number_of_spacers = 10
pairsh = {}
blast_pairsh = {}
pair_sep = "-"
e_value_cutoff = 1e-10
aln_length_cutoff = 0
identity_cutoff = 0
de_duplicate_clustersh = {}
gff_type = None
include_list_file = None
attr = "ID|Parent"
gene_regexp = None
 

genes_included = {}
genes_by_chr = {}


########################################################################
class Gene_from_gff():
    def __init__(self):
        pass

def show_help():
    script_basename=os.path.basename(sys.argv[0])
    print "Usage of " + script_basename
    sys.exit()


def read_pair_files(pair_files, sep='-'):
    pairsh={}
    for pair_file in pair_files:
        in_fh = open(pair_file, 'r')
        for line in in_fh.readlines():
            line=line.strip('\r\n')
            pairsh[line] = ''
    return(pairsh)


def read_blast_files(blast_files, e_value_cutoff=1e-10, gene_regexp=None, aln_length_cutoff=0, identity=0):
    tmp_pairsh = {}
    for blast_file in blast_files:
        in_fh = open(blast_file, 'r')
        for line in in_fh.readlines():
            line=line.strip('\r\n')
            line_array = line.split("\t")
            gene1, gene2, identity, aln_length, e_value = line_array[0], line_array[1], line_array[2], line_array[3], line_array[10]
            if gene_regexp:
                m = re.search(gene_regexp, gene1)
                if m:
                    gene1 = m.group(1)
                m = re.search(gene_regexp, gene2)
                if m:
                    gene2 = m.group(1)

            aln_length = int(aln_length)
            identity = float(identity)
            if gene1 == gene2:
                continue
            if float(e_value) > 1e-10:
                if float(e_value) <= e_value_cutoff:
                    if aln_length >= aln_length_cutoff and identity >= identity_cutoff:
                        pass
                    else:
                        continue
                else:
                    continue
            tmp_pairsh["-".join(sorted([gene1,gene2]))] = ''
        in_fh.close()
    return(tmp_pairsh)


def read_gff(gff_files, genes_included, attr, type="gff3"):
    gene_objsh = {}
    for gff_file in gff_files:
        in_fh = open(gff_file,'r')
        for line in in_fh:
            #chr1	At1NC000020	2497	2816	-1
            if re.search("^#", line):
                continue
            line = line.strip('\r\n')
            line_array = line.split('\t')
            if type == "saf":
                chr, gene, start, end = line_array[0:4]
            elif type == "gff3":
                chr, start, end = line_array[0], line_array[3], line_array[4]
                m = re.search(attr + "=([^;]+)", line_array[-1])
                if not m:
                    continue
                gene = m.groups(1)[0]
            if genes_included and not gene in genes_included:
                continue
            m = re.search("[<>]([0-9]+)", str(start))
            if m:
                start = m.group(1)
            m = re.search("[<>]([0-9]+)", str(end))
            if m:
                end = m.group(1)
            start=int(start)
            end=int(end)
            gene_obj=Gene_from_gff()
            gene_obj.chr = chr
            gene_obj.end = end
            gene_obj.start = start
            gene_objsh[gene] = gene_obj # put every gene_obj into gene_objsh
        in_fh.close
    return(gene_objsh)


def parse_pairsh(pairsh, sep='-'):
    gene_relations = {}
    for k,v in pairsh.iteritems():
        list = k.split(sep)
        for index,gene in enumerate(list):
            if not gene in gene_relations:
                gene_relations[gene] = {}
            gene_relations[gene][list[abs(1-index)]] = ''
    return(gene_relations)


def get_tandem_relations(seq_objsh, gene_relations, distance, genes_by_chr, number_of_spacers):
    tandem_relations = {}
    for gene, paired_genesh in seq_objsh.iteritems():
        gene_chr = seq_objsh[gene].chr
        gene_start = seq_objsh[gene].start
        if not gene in gene_relations:
            continue
        #for paired_gene in gene_relations[gene].keys():
        #    print "-".join(sorted([gene,paired_gene]))
        for paired_gene in gene_relations[gene].keys():
            if not paired_gene in seq_objsh:
                continue
            if seq_objsh[paired_gene].chr != gene_chr:
                continue
            if abs(seq_objsh[paired_gene].start - gene_start) > distance:
                continue
            if abs(genes_by_chr[gene_chr].index(gene) - genes_by_chr[gene_chr].index(paired_gene)) > number_of_spacers:
                continue
            for i in [gene, paired_gene]:
                if not i in tandem_relations:
                    tandem_relations[i] = {}
            tandem_relations[gene][paired_gene] = ''
            tandem_relations[paired_gene][gene] = ''
    return(tandem_relations)


def iteratively_find_genes_in_cluster(gene_hash, genes_in_cluster):
    for paired_gene in gene_hash.keys():
        for paired_gene2 in tandem_relations[paired_gene].keys():
            if not paired_gene2 in genes_in_cluster:
                genes_in_cluster.append(paired_gene2)
                iteratively_find_genes_in_cluster(tandem_relations[paired_gene2], genes_in_cluster)


def get_tandem_clusters(tandem_relations):
    tandem_clusters = []
    for gene, hash in tandem_relations.iteritems():
        genes_in_cluster = [gene] + tandem_relations[gene].keys()
        iteratively_find_genes_in_cluster(tandem_relations[gene], genes_in_cluster)
        tandem_clusters.append(genes_in_cluster)
    return(tandem_clusters)


def read_include_list(include_list_file):
    genes = {}
    in_fh = open(include_list_file, "r")
    for line in in_fh.readlines():
        line = line.strip('\n\r')
        genes[line.split("\t")[0]] = ""
    in_fh.close
    return(genes)


########################################################################
try:
    opts, args = getopt.getopt(
        sys.argv[1:],
        "i:n:d:g:p:e:",
        ["input=","in=","blast8=","gff=","pair=","pair_file=","number_of_spacers=","distance=","sep=","e_value=","aln_length=","length=","identity=","type=","include_list=", "list_included=", "attr=", "gene_regexp="],
    )

except getopt.GetoptError:
    print "Illegal params!"
    show_help()
for op, value in opts:
    if op == "-p" or op == "--pair" or op == "--pair_file":
        pair_files.append(value)
    if op == "-i" or op == "--in" or op == '--input' or op == "--blast8":
        blast_files.append(os.path.expanduser(value))
    if op == "-g" or op == "--gff":
        gff_files.append(value)
    if op == "-n" or op == "--number_of_spacers":
        number_of_spacers = int(value)
    if op == "-d" or op == "--distance":
        distance = int(value)
    if op == "--sep":
        pair_sep = value
    if op == "--e_value" or op == "-e":
        e_value_cutoff = float(value)
    if op == "--aln_length" or op == "--length":
        aln_length_cutoff = int(value)
    if op == "--identity":
        identity_cutoff = float(value)
    if op == "--type":
        gff_type = value
    if op == "--include_list" or op == "--list_included":
        include_list_file = value
    if op == "--attr":
        attr = value
    if op == "--gene_regexp":
        gene_regexp = value

if not gff_files:
    print "gff must be given! Exiting ......"
    sys.exit()


########################################################################
if include_list_file:
    genes_included = read_include_list(include_list_file)

pairsh = read_pair_files(pair_files, pair_sep)

blast_pairsh = read_blast_files(blast_files, e_value_cutoff, gene_regexp, aln_length_cutoff, identity_cutoff)

pairsh.update(blast_pairsh)

gene_relations = parse_pairsh(pairsh, pair_sep)

gene_objsh = read_gff(gff_files, genes_included, attr, gff_type)

for gene, gene_obj in gene_objsh.iteritems():
    if not gene_obj.chr in genes_by_chr:
        genes_by_chr[gene_obj.chr] = []
    genes_by_chr[gene_obj.chr].append(gene)

for chr, genes in genes_by_chr.iteritems():
    genes.sort(key=lambda i: gene_objsh[i].start)

tandem_relations = get_tandem_relations(gene_objsh, gene_relations, distance, genes_by_chr, number_of_spacers)

tandem_clusters = get_tandem_clusters(tandem_relations)

for cluster in tandem_clusters:
    de_duplicate_clustersh["\t".join(sorted(cluster))] = ''

for k,v in de_duplicate_clustersh.iteritems():
    print k 


