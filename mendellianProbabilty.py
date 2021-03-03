#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 22:33:26 2021

@author: Sid
"""

# open file containing three ints and store in list 
def extractIntegers(): 
     ints = []
     with open(r'/Users/Sid/Downloads/rosalind_iprb.txt') as f:
         for line in f:
             data = line.split()
             ints.append(int(data[0]))
             ints.append(int(data[1]))
             ints.append(int(data[2]))
     return ints
# method to calculate probability carrying dom allele 
def homozygousRecProbability(ints=[]): 
    k = 16 # num homozygous dominant orgs
    m = 22 # num heterozygous orgs
    n = 17 # num homozygous recesssive orgs 

    # calculate probabilities of how to obtain homozygous rec phenotpye 
    probHet_Het = (m/(k+m+n))*((m-1)/(k+m+n-1))*(1/4) 
    probHet_HomRec = (m/(k+m+n))*(n/(k+m+n-1))*(1/2) 
    probHomRec_Het = (n/(k+m+n))*(m/(k+m+n-1))*(1/2) 
    probHomRec_HomRec = (n/(k+m+n))*((n-1)/(k+m+n-1))
    # probability of homozygous dom allele is (1 - sum(above))
    probHomDom = 1-(probHet_Het+probHet_HomRec+probHomRec_Het+probHomRec_HomRec)
    return probHomDom