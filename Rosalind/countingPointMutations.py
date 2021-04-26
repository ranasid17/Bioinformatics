#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 18:33:04 2021

@author: Sid
"""

# Read input sequences from file and store in list 
inputSequences = []
with open('/Users/Sid/Downloads/rosalind_hamm.txt') as f:
    inputSequences = f.readlines()

# Function to calculate Hamming Distance of two sequences
# Input: 
    # 1) Input array: has two elts, each a DNA seq 
# Output: 
    # 1) Integer value of hamming distance 
def hammingDistance(inputArray): 
    strandOne = inputArray[0]
    strandTwo = inputArray[1]
    
    mismatch = 0 # initialize count for num mismatches
    for i in range(0, len(strandOne)): # iter thru strandOne
        if strandOne[i] != strandTwo[i]: # when a mismatch found
            mismatch = mismatch + 1 # incr total mistamtches 
    return mismatch


test = hammingDistance(inputSequences)
