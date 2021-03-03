#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 22:33:26 2021

@author: Sid
"""

# open txt file containing DNA sequences 
dna = open(r'/Users/Sid/Documents/2020-21/BCB 5250/Rosalind DataSets/rosalind_rna.txt').read() 
# replace each instance of 'T' with 'U' 
rna = dna.replace("T", "U") 
# print original and new strings
print("Original DNA strand is: ", dna)
print("RNA sequence for associated DNA strand is: ", rna)

