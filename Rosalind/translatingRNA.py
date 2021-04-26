#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 11:31:45 2021

@author: Sid
"""

def protein (rna):
    codon = ""
    protein = "" # create empty string for protein 
    for i in range(0, len(rna) - 3, 3):
        codon = rna[i:i+3]
        if (codon == "UUU" or codon == "UUC"):
            protein = protein + "F"
        elif (codon == "UUA" or codon == "UUG"):
            protein = protein + "L"
        elif (codon == "UCU" or codon == "UCC" or codon == "UCA" or codon == "UCG"):
            protein = protein + "S"
        elif (codon == "UAU" or codon == "UAC"):
            protein = protein + "Y"
        elif (codon == "UAA" or codon == "UAG" or codon == "UGA"):
            break # stop codon 
        elif (codon == "UGU" or codon == "UGC"):
            protein = protein + "C"
        elif (codon == "UGG"):
            protein = protein + "W"
        elif (codon == "CUU" or codon == "CUC" or codon == "CUA" or codon == "CUG"):
            protein = protein + "L"
        elif (codon == "CCU" or codon == "CCC" or codon == "CCA" or codon == "CCG"):
            protein = protein + "P"
        elif (codon == "CAU" or codon == "CAC"):
            protein = protein + "H"
        elif (codon == "CAA" or codon == "CAG"):
            protein = protein + "Q"
        elif (codon == "CGU" or codon == "CGC" or codon == "CGA" or codon == "CGG" or codon == "AGA" or codon == "AGG"):
            protein = protein + "R"
        elif (codon == "AUU" or codon == "AUC" or codon == "AUA"):
            protein = protein + "I"
        elif (codon == "AUG"):
            protein = protein + "M"
        elif (codon == "ACU" or codon == "ACC" or codon == "ACA" or codon == "ACG"):
            protein = protein + "T"
        elif (codon == "AAU" or codon == "AAC"):
            protein = protein + "N"
        elif (codon == "AAA" or codon == "AAG"):
            protein = protein + "K"
        elif (codon == "AGU" or codon == "AGC"):
            protein = protein + "S"
        elif (codon == "GUU" or codon == "GUC" or codon == "GUA" or codon == "GUG"):
            protein = protein + "V"
        elif (codon == "GCU" or codon == "GCC" or codon == "GCA" or codon == "GCG"):
            protein = protein + "A"
        elif (codon == "GAU" or codon == "GAC"):
            protein = protein + "D"
        elif (codon == "GAA" or codon == "GAG"):
            protein = protein + "E"
        elif (codon == "GGU" or codon == "GGC" or codon == "GGA" or codon == "GGG"):
            protein = protein + "G"
    return protein

rnaSeq = open(r'/Users/Sid/Downloads/rosalind_prot.txt').read()
aaSeq = protein(rnaSeq)
print("Translated protein sequence is: ", aaSeq)
