#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 20:44:57 2021

@author: Sid
"""

from Bio import SeqIO
""" Note: Instead of comparing all seqs against each other, we simply 
    need to compare the longest (a) and shortest (b) strings. The longest 
    common substring between a, b will be true for all strings in between."""

sequences = [] # empty list to hold sequence                 
handle = open('/Users/Sid/Documents/2020-21/BCB 5250/Rosalind DataSets/rosalind_lcsm.txt', 'r') # open file   
for record in SeqIO.parse(handle, 'fasta'): # iterate thru file 
    sequence = [] # empty list                   
    seq = '' # create empty str for current sequence             
    for nt in record.seq: #           
        seq += nt # add nt to current sequence        
    sequences.append(seq) # add current sequence to list   
handle.close() # close file 

sorted_seqs = sorted(sequences, key=len) # sort list of seqs by length 
shortest_seq = sorted_seqs[0] # store shortest seq
comp_seq = sorted_seqs[1:] # store longest seq 
common_substring = '' # create empty str for current sequence 
for i in range(len(shortest_seq)): # iterate len of shortest seq 
    for j in range(i, len(shortest_seq)): # iterate across shortest seq 
        m = shortest_seq[i:j + 1] # add one more nt to curr position 
        found = False # boolean to indicate substr identification 
        for sequ in comp_seq:
            if m in sequ: # check for substr being present in longest seq     
                found = True             
            else: # if substr not found      
                found = False # keep false      
                break                    
        if found and len(m) > len(common_substring): # check if new substr is longer 
            common_substring = m # if so 
print(common_substring)    
