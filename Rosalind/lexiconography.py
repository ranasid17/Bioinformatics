#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 15:34:16 2021

@author: Sid
"""

import itertools # python lib for permutations, combinations 

# copied from previous Rosalind submission 
## code to open file and store data
with open('/Users/Sid/Documents/2020-21/BCB 5250/Rosalind DataSets/rosalind_lexf.txt') as f:
    inputs = f.read().split()
    basePairs = inputs[:-1]
    n = int(inputs[-1]) # store string length (last value given in dataset)

# calc all permutations of given chars for n number of repeats 
permuateBases = itertools.product(basePairs, repeat = n) 
# create empty list to store final soln 
output = []
# iterate across list of permutations 
for i, j in enumerate(list(permuateBases)):
    # create empty string for current pair 
    curr_permutation = ''
    # iterate across jth pair 
    for elt in j:
        # add current enumeration to list of permutations
        curr_permutation += str(elt)
    # output current permutation
    output.append(curr_permutation)
# Note: Regular print output does not work for submissio 
# print(output)
# must print each item in own line for correct submission 
for i in output:
    print(i, end="\n")