#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 11 19:10:22 2021

@author: Sid
"""

from Bio import SeqIO           
           
# copied from prior Rosalind assignments to open files
sequences = []
handle = open('/Users/Sid/Documents/2020-21/BCB 5250/Rosalind DataSets/rosalind_sseq.txt', 'r')
for record in SeqIO.parse(handle, 'fasta'):
    sequences.append(str(record.seq))
handle.close()
target_seq = sequences[0]
subseq = sequences[1]

index = 0 # initialize var to hold index        
matching_indices = [] # initialize list to hold matching indices of substring
                          
for i in range(len(subseq)): # iterate across target seq
    for j in range(index, len(target_seq)): # iterate across target seq by incr index 
        index += 1 # increment index                     
        if len(matching_indices) < len(subseq): # check there is space on string
            if subseq[i] == target_seq[j]: # when target seq ID's
                matching_indices.append(index) # add subseq match index to list 
                break
print(*matching_indices, sep=' ')    
