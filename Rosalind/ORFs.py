#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 19:37:20 2021

@author: Sid
"""

# open rosalind_orf.txt file 
with open('/Users/Sid/Documents/2020-21/BCB 5250/Rosalind DataSets/rosalind_orf.txt',"r") as f:
    f.readline()
    dna=''
    for line in f:
        dna+=line.strip()
    dna=list(dna)

# create empty list to hold reverse strand 
reverse_complement=[]
# iterate thru input dna seq
for i in range(len(dna)):
    ORF=len(dna)-1-i # set current ORF 
    if dna[ORF]=="T": # scan DNA for T
        dna[ORF]="U" # replace with U 
        reverse_complement.append("A") # add A to reverse complement
    else:
        if dna[ORF] == "A": # scan DNA for A 
            reverse_complement.append("U") # add U to reverse complement 
        else:
            if dna[ORF]=="C": # scan DNA for C 
                reverse_complement.append("G") # add G to reverse complement 
            else:
                reverse_complement.append("C")

rna=''.join(dna)
reverse_complement=''.join(reverse_complement)

# define codon table 
codon={"UUU": "F", "CUU": "L", "AUU": "I", "GUU": "V", "UUC": "F", "CUC": "L", "AUC": "I", "GUC": "V", "UUA": "L",
            "CUA": "L", "AUA": "I", "GUA": "V", "UUG": "L", "CUG": "L", "AUG": "M", "GUG": "V", "UCU": "S", "CCU": "P",
            "ACU": "T", "GCU": "A", "UCC": "S", "CCC": "P", "ACC": "T", "GCC": "A", "UCA": "S", "CCA": "P", "ACA": "T",
            "GCA": "A", "UCG": "S", "CCG": "P", "ACG": "T", "GCG": "A", "UAU": "Y", "CAU": "H", "AAU": "N", "GAU": "D",
            "UAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D", "UAA": "Stop", "CAA": "Q", "AAA": "K", "GAA": "E",
            "UAG": "Stop", "CAG": "Q", "AAG": "K", "GAG": "E", "UGU": "C", "CGU": "R", "AGU": "S", "GGU": "G",
            "UGC": "C", "CGC": "R", "AGC": "S", "GGC": "G", "UGA": "Stop", "CGA": "R", "AGA": "R", "GGA": "G",
            "UGG": "W", "CGG": "R", "AGG": "R", "GGG": "G"}
# create list to hold final answer
polypeptides=[]
start_codon = 'AUG'
# iterate across RNA string 
for i in range(len(rna)-2):
    if rna[i:i+3]==start_codon: # search for start codon 
        j=i # set j == i 
        prot='' # create empty string to hold translated protein 
        letter='AUG' 
        # iterate across substring until stop codon reached
        while codon[letter]!="Stop":
            prot+=codon[letter] # add translated a.a. to protein 
            j+=3 # move to next codon (+3 spots )
            if j>len(rna)-3: # check for end of string (or <3 spaces till end)
                break
            letter=rna[j:j+3] # set next a.a.
        # check for stop codon and for uniqueness of translated protein 
        if (codon[letter]=="Stop") and (prot not in polypeptides):
            # add protein to list iff unique and translation terminated 
            polypeptides.append(prot)
            
# iterate across reverse complement string 
for i in range(len(reverse_complement)-2):
    # repeat same process from lines 50-65
    if reverse_complement[i:i+3]=='AUG':
        j=i
        prot=''
        letter='AUG'
        while codon[letter]!="Stop":
            prot+=codon[letter]
            j+=3
            if j>len(reverse_complement)-3:
                break
            letter=reverse_complement[j:j+3]
        if codon[letter]=="Stop" and prot not in polypeptides:
            polypeptides.append(prot)
# print final list of polypeptides 
for i in polypeptides:
    print(i)
