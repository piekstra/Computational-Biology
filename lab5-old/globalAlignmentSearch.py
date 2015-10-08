#!/usr/bin/python

# global alignment brute force search implementation
# for each k from 0 to L (where L is the min length of the two strings)
#      Create all subsequences of size k for each string
#           Align letters from each subsequence to calculate alignment score
#           All non-chosen letters match with "-" for alignment score
#
# author: Tammy VanDeGrift
# CS 423, Lab 5
# Fall 2015

import itertools # to use combinations method
import time      # to calculate running times

# Determine the score of the optimal global alignment of S and T
def globalAlignmentSearch(S, T):
    
    # Scoring system for gaps, matches, and mismatches
    MATCH = 5
    MISMATCH = -4
    GAP = -6

    # variables to keep track of best alignment found so far
    bestScore = GAP * (len(S) + len(T)) # aligning both strings with dashes
    bestAlignmentS = []
    bestAlignmentT = []

    # find length of shorter string to use for subset selection size
    L = min(len(S), len(T))

    # create all possible subsets of size K for K = 0 to L
    for k in range(L+1):
        # to see progress of program execution
        print("size of chosen subset: " + str(k))
        
        # create list of all subsets
        setsOfS = [x for x in itertools.combinations(S, k)]
        setsOfT = [y for y in itertools.combinations(T, k)]
        
        # go through each subset of S and match letters to each
        # subset of T
        for i in range(len(setsOfS)):
            for j in range(len(setsOfT)):
                # calculate score aligning each character from each set
                subseqS = setsOfS[i]
                subseqT = setsOfT[j]
                alignScore = 0
                for m in range(k):
                    if(subseqS[m] == subseqT[m]):
                        alignScore += MATCH
                    else:
                        alignScore += MISMATCH
                # now add the gaps for the letters of
                # S that are not part of the subset
                alignScore += GAP * (len(S) - k)
                # now add the gaps for the letters of
                # T that are not part of the subset
                alignScore += GAP * (len(T) - k)

                # determine if better alignment is found
                if(alignScore > bestScore):
                    bestScore = alignScore
                    bestAlignmentS = subseqS
                    bestAlignmentT = subseqT
                    
    # for debugging purposes, print the two subsequences that are aligned
    # uncomment if you want to see these printed
    # print(bestAlignmentS)
    # print(bestAlignmentT)
    return bestScore


###########################################
s = "AAAGGTAACGTACAT"
t = "ATCAATCGCGACGAG"
t0 = time.clock()
optimalScore = globalAlignmentSearch(s, t)
print("Global alignment score via search: " + str(optimalScore))
t1 = time.clock()

print("Total time: " + str(t1-t0))
    
