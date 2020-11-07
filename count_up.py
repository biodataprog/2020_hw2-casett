#!/usr/bin/env python3

# this is a python script template
# this next line will download the file using curl

gff="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz"
fasta="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz"

import os,gzip,itertools,csv,re

# this is code which will parse FASTA files
# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'

def aspairs(f):
    seq_id = ''
    sequence = ''
    for header,group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence



if not os.path.exists(gff):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz")

if not os.path.exists(fasta):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz")


#add counters
num_genes = 0
total_length = 0
    
with gzip.open(gff,"rt") as fh:
    # now add code to process this
    gff = csv.reader(fh,delimiter="\t")
    for row in gff:
        if row[0].startswith("#"):
            continue
        
        #if the feature is a gene, add 1 
        #also, get the length of the gene, end - start
        if row[2] == "gene":
            num_genes += 1
            total_length += (int(row[4]) - int(row[3]))
    
    print("The total number of genes is:", num_genes, "\nThe total length of genes is:", total_length)


 
with gzip.open(fasta,"rt") as f:
    seqs = dict(aspairs(f))

#define counters            
total_length_genome = 0

for k,v in seqs.items():
    #for each sequence v in seq.items - add length to total genome
    total_length_genome += len(v)

print("The total length of the genome is:", total_length_genome)
print("The % of the genome that is coding is:", (total_length / total_length_genome) * 100)
