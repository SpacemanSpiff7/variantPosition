""" This script takes a file with regions of interest and maps gene names to those regions if they exist """
import sys
import os
import re
import pandas as pd


def convert(chrom):
    split = re.findall(r"\d+", chrom)
    # only interested in chrom# if it doesn't contain any other numbers
    # other numbers correspond to different parts and should be treted as strings
    if len(split) == 1:
        return split[0]
    return chrom


class Gene:
    def __init__(self, chrom, start, stop, names):
        self.chr = convert(chrom)
        self.start = start
        self.stop = stop
        self.names = names
        self.isSingleChrom = isinstance(self.chr, int)

    def __add__(self, other):
        if isinstance(other, Gene):
            if self.chr != other.chr:
                return self.chr + other.chr
            else:
                return self.start + other.start
        else:
            return -100

    def __sub__(self, other):
        if isinstance(other, Gene):
            if self.chr != other.chr:
                return self.chr - other.chr
            else:
                return self.start - other.start
        else:
            return -100


# this file has a list of regions in the genome that are  
# separated by a large number of base pairs (bp specified by user)
in_file = pd.read_csv(sys.argv[1], header=None, sep='\t')
# the file should have these 3 tab-separated columns with no headers
in_file.columns = ["chromosome", "start", "stop"]


# this is the file that conatins all of the genes and their positions
# genefile = pd.read_csv('/Users/simonelongo/Desktop/QuinlanGroup/variantPosition/sept5chr22/genes.chr22.bed', sep='\t')
def import_genes(path):
    gene_dict = dict()
    temp = pd.read_csv(path, sep='\t')
    # import each gene by chromosome
    for i, r in temp.iterrows():
        new_gene = Gene(r["#chrom"], r["txStart"], r["txEnd"], r["name2"])
        gene_dict.update()
        if len(gene_dict[new_gene.chr]) == 0:
            gene_dict[new_gene.chr] = [new_gene]
        else :
            gene_dict[new_gene.chr].append(new_gene)
    del temp  # clear from memory, no longer needed

    # sort by starting position
    for key in gene_dict:
        gene_dict[key].sort()

    return gene_dict


# this is a dictionary with a list of positions mapped to each unique chromosome key
genefile = import_genes('/Users/simonelongo/Desktop/QuinlanGroup/variantPosition/sept5chr22/genes.chr22.bed')


def get_gene_name(current_row):
    chr_num = convert(current_row["chromosome"])
    add_index = genefile[chr_num].searchsorted(current_row["start"]) - 1
    gene_in_zone = genefile[chr_num][add_index]
    if gene_in_zone.stop >= current_row["stop"]:
        return gene_in_zone.names
    return ""


for index, row in in_file.iterrows():
    gene = get_gene_name(row)
    if len(gene) > 0:
        row["gene"] = gene

filename = os.path.splitext(sys.argv[1])[0] + "_named_genes.bed"
in_file.to_csv(filename, index=False, sep='\t')

# class GeneList:
#     def __init__(self):
#         self.all_genes = []
#
#     def add_gene(self, gene):
#         if gene.isSingleChrom:
#             self.all_genes.append(gene)
