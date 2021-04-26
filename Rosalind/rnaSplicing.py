#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 23:29:55 2021

@author: Sid
"""

# import file 
file = open('/Users/Sid/Documents/2020-21/BCB 5250/Rosalind DataSets/rosalind_splc.txt', "r")

codon={"UUU": "F", "CUU": "L", "AUU": "I", "GUU": "V", "UUC": "F", "CUC": "L", "AUC": "I", "GUC": "V", "UUA": "L",
            "CUA": "L", "AUA": "I", "GUA": "V", "UUG": "L", "CUG": "L", "AUG": "M", "GUG": "V", "UCU": "S", "CCU": "P",
            "ACU": "T", "GCU": "A", "UCC": "S", "CCC": "P", "ACC": "T", "GCC": "A", "UCA": "S", "CCA": "P", "ACA": "T",
            "GCA": "A", "UCG": "S", "CCG": "P", "ACG": "T", "GCG": "A", "UAU": "Y", "CAU": "H", "AAU": "N", "GAU": "D",
            "UAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D", "UAA": "Stop", "CAA": "Q", "AAA": "K", "GAA": "E",
            "UAG": "Stop", "CAG": "Q", "AAG": "K", "GAG": "E", "UGU": "C", "CGU": "R", "AGU": "S", "GGU": "G",
            "UGC": "C", "CGC": "R", "AGC": "S", "GGC": "G", "UGA": "Stop", "CGA": "R", "AGA": "R", "GGA": "G",
            "UGG": "W", "CGG": "R", "AGG": "R", "GGG": "G"}

def extract_sequence(file):
    sequences = [] # create list to hold sequences 
    result = ""       
    
    # iterate across each line 
    for line in file:
        # find new seq by '>' char
        if ">" in line:
            sequences.append(result)
            result = ""
        else:
            # check for '/n' char
            if "\n" in line:
                result += line[:len(line) - 1]
            else:
                result += line
    # add final result to list of sequences
    sequences.append(result)
    sequences.remove('')
    return sequences

# run above function 
sequences = extract_sequence(file)

# store dna string 
dna = sequences[0]
# store introns in reverse length order
introns = sequences[1:]
introns.sort(reverse = True)

# iterate across list of introns 
for intron in introns:
    # remove them from dna once aligned
    dna = dna.replace(intron,"")
# replace T with U (dna has no U)
dna = dna.replace("T","U")

# print final answer
for i in range(0, len(dna) - 3, 3):
    print(codon[dna[i: i +3]],end='')
