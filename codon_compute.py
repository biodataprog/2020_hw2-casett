#!/usr/bin/env python3

import os, gzip, itertools

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

url1="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2/cds/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
url2="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/cds/Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"
file1="Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
file2="Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"
sp1 = "S. enterica"
sp2 = "M. tuberculosis"

if not os.path.exists(file1):
    os.system("curl -O %s"%(url1))

if not os.path.exists(file2):
    os.system("curl -O %s"%(url2))

with gzip.open(file1,"rt") as fh:
    seqs = aspairs(fh)
    
    num_genes = 0
    total_len = 0
    gc = 0
    codon_count = {}
    num_codons = 0
    
    for seq in seqs:
        seqname  = seq[0]
        seqstring= seq[1]
        
        num_genes += 1
        total_len += len(seqstring)
        
        #codon counting probably needs to be its own helper function
        #but a problem for another day
        codon = ""
        for base in seqstring:
            codon += base
            
            if base in ["G", "C"]:
                gc += 1
            
            if len(codon) == 3:
                if codon in codon_count:
                    codon_count[codon] += 1
                else:
                    codon_count[codon] = 1
                        
                codon = ""
    
        
    for codon in codon_count:
        num_codons += codon_count[codon]
        #print(codon, codon_count[codon]) 

    print("The total number of genes for", sp1, ":", num_genes)
    print("The total genome length for", sp1, ":",total_len)
    print("The % GC content for", sp1, "is:", gc/total_len*100)
    print("The total number of codons for", sp1, ":", num_codons)


with gzip.open(file2,"rt") as fh:
    seqs = aspairs(fh)
    
    num_genes = 0
    total_len = 0
    gc = 0
    codon_count_sp2 = {}
    num_codons = 0
    
    for seq in seqs:
        seqname  = seq[0]
        seqstring= seq[1]
        
        num_genes += 1
        total_len += len(seqstring)
        
        
        codon = ""
        for base in seqstring:
            codon += base
            
            if base in ["G", "C"]:
                gc += 1
            
            if len(codon) == 3:
                if codon in codon_count_sp2:
                    codon_count_sp2[codon] += 1
                else:
                    codon_count_sp2[codon] = 1
                        
                codon = ""
    
        
    for codon in codon_count_sp2:
        num_codons += codon_count_sp2[codon]
        #print(codon, codon_count[codon]) 

    print("The total number of genes for", sp2, ":", num_genes)
    print("The total genome length for", sp2, ":",total_len)
    print("The % GC content for", sp2, "is:", gc/total_len*100)
    print("The total number of codons for", sp2, ":", num_codons)


print("Codon \t", sp1, "\t", sp2)
for codon in codon_count:
    if codon in codon_count_sp2:
        print(codon, "\t", codon_count[codon], "\t", codon_count_sp2[codon])
    else:
        print(codon, "\t", codon_count[codon], "\t", 0)
for codon in codon_count_sp2:
    if codon not in codon_count:
        print(codon, "\t", 0, "\t", codon_count_sp2[codon])
        
        
