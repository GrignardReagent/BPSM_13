#!/bin/python3

import os
import pandas as pd

### 1. Processing tabular data

# Using os and doing it the old school way
flygenedata = open('data.csv') # This line has to be ran every time or the for loop won't work.
for geneline in flygenedata:
    # Define our variables from the split() list elements
    gene = geneline.split(",")
    species = gene[0]
    geneseqs = gene[1].upper()
    seqlengths = len(gene[1])
    genenames = gene[2]
    expression_level = int(gene[3])
    # Print out the gene names for all genes from the species Drosophila melanogaster or Drosophila simulans.
    if "melanogaster" in species or "simulans" in species:
        print("melanogaster or simulans: ", species, genenames)
    # Print out the gene names for all genes that are between 90 and 110 bases long.
    if seqlengths > 90 and seqlengths < 110:
        print("gene length between 90 and 110: ", species, genenames)
    # Print out the gene names for all genes whose AT content is less than 0.5 and whose expression level is greater than 200.
    atcontent = (geneseqs.count("A")+geneseqs.count("T"))/seqlengths
    if atcontent < 0.5 and expression_level > 200:
        print("AT content is less than 0.5 and whose expression level is greater than 200: ", species, genenames)
    # Print out the gene names for all genes whose name begins with "k" or "h" except those belonging to Drosophila melanogaster.
    if genenames.startswith("k") or genenames.startswith("h") and ("Drosophila melanogaster" not in species):
        print("Start with k or h, except those belong to melanogaster: ", species, genenames)
    # For each gene, print out a message giving the gene name and saying whether its AT content is high (greater than 0.65), low (less than 0.45) or medium (between 0.45 and 0.65).
    atstatus = ''
    if atcontent > 0.65:
        atstatus = "high"
    if atcontent < 0.45:
        atstatus = 'low'
    if atcontent >= 0.45 and atcontent <= 0.65:
        atstatus = 'medium'

    print(genenames, ' has a', atstatus, 'AT content')


### 2. Kmer counting
sequencein="atatatatatcgcgtatatatacgactatatgcattaattatagcatatcgatatatatatcgatattatatcgcattatacgcgcgtaattatatc"
possible_kmer_sizes = list(range(2,7))
kmersfound_minimum = 3
for window in possible_kmer_sizes:
    kmersfound = []
    kmerrange = list(range(0,len(sequencein)))
    for startingbase in kmerrange:
        if (startingbase + window) < len(sequencein) + 1:
            seqout = sequencein[startingbase:startingbase+window]
            kmersfound = kmersfound + [seqout]
    nonredundantset = list(set(kmersfound))
    for kmerfreq in nonredundantset:
        if kmersfound.count(kmerfreq) > kmersfound_minimum:
            print('Lots! ', kmerfreq, ' ', str(kmersfound.count(kmerfreq)))
        else:
            print(str(kmerfreq), ' ', str(kmersfound.count(kmerfreq)))


sequencein = input("Please enter the raw DNA sequence you want to analyse\n").upper()
print("Thanks for entering the sequence", sequencein, "(which is",len(sequencein), "characters long)")

### 3.
# INPUT
seqs = ['ATTGTACGG', 'AATGAACCG', 'AATGAACCC', 'AATGGGAAT']
# PROCESS
for_range = list(range(0,3))
for i in for_range:
    query_1 = list(seqs[i])
    for j in for_range:
        # Go along each base in each of the two sequences to compare if identical or not; score +1 if it is
        count = 0
        query_2 = list(seqs[j])
        for base in list(range(0,len(query_1))):
            # print("Index", str(base), ":", str(query_1[base]), ",", str(query_2[base]))
            if query_1[base] == query_2[base]:
                count += 1
        # OUTPUT
        # DONT include the self-comparisons in the summary output
        if seqs[i] != seqs[j]:
            print(count," identities between ", seqs[i], " and ", seqs[j])
            print("\t", 100*(count/len(seqs[i]), " percent similarity between ", seqs[i], " and ", seqs[j]))

