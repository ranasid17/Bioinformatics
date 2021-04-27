#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 09:48:49 2021

@author: Sid
"""

data = [] # empty array to hold input sequence 
with open('/Users/Sid/Downloads/rosalind_lgis.txt', 'r') as f:  
    for line in f:                      
        for j in line.split():         
            data.append(int(j)) 
perm = data[1:]


def longestIncreasingSubseq(seq):
    P = [None] * len(seq)
    M = [None] * len(seq)
    # There must exist at least one elt in the list 
    L = 1
    # Therefore we know there exists an increasing subseq w length == 1
    M[0] = 0
    # iter across seq from second elt to end (second elt bc already stored first)
    for i in range(1, len(seq)):
        # apply binary search s.t. largest j less than/eq to L providing 
        lo = 0
        hi = L        
        # binary search algo will ignore upper bound so apply manual check
        if seq[M[hi - 1]] < seq[i]:
            j = hi
        # meanwhile lower bound should be set to end of search.
        else:
            # binary search algo 
            while hi - lo > 1:
                mid = (hi + lo) // 2
                if seq[M[mid - 1]] < seq[i]:
                    lo = mid
                else:
                    hi = mid

            j = lo
        P[i] = M[j - 1]
        # reconstruct search results 
        if j == L or seq[i] < seq[M[j]]:
            M[j] = i
            L = max(L, j + 1)

    result = []
    pos = M[L - 1]
    for k in range(L):
        result.append(seq[pos])
        pos = P[pos]
    # return in correct order
    return (result[::-1])

# this is very similar to LIS function
def longestDecreasingSubseq(seq):
    P = [None] * len(seq)
    M = [None] * len(seq)
    
    L = 1 
    M[0] = 0
    for i in range(1, len(seq)):
        lo = 0
        hi = L
        # use M[hi-1] > seq[i] rather than < 
        if seq[M[hi - 1]] > seq[i]:
            j = hi
        else:
            while hi - lo > 1:
                mid = (hi + lo) // 2
                # use M[mid-1] > seq[i] rather than < 
                if seq[M[mid - 1]] > seq[i]:
                    lo = mid
                else:
                    hi = mid

            j = lo
        P[i] = M[j - 1]
        # use > instead of < again
        if j == L or seq[i] > seq[M[j]]:
            M[j] = i
            L = max(L, j + 1)

    result = []
    pos = M[L - 1]
    for k in range(L):
        result.append(seq[pos])
        pos = P[pos]

    return (result[::-1])

incr = longestIncreasingSubseq(perm)
decr = longestDecreasingSubseq(perm)

# print(*incr)
# print(*decr)
