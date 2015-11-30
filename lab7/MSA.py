#!/usr/bin/python

################################################
# Lab7 CS423 
# Shaun Stice, Caleb Piekstra
# Fall 2015
################################################

import random, math


####################################################################
# Create a 2D table with the given number of rows and columns
# and fills all entries with value given as a parameter
# (function completed for you)
####################################################################
def createTable(numRows, numCols, value):
     table = []
     row = 0
     # create 2D table initialized with value
     while (row < numRows):
          table.append([])    
          col = 0
          while (col < numCols):
               table[row].append(value)
               col = col + 1
          row = row + 1
     return table


######################################################################
# Determines the score of the optimal global alignment of two sequences
# Taken from lab 5
# Just calculates the score (not the alignment; do not need to
# print out tables; just returns the optimal score)
######################################################################
def globalAlignmentScore(s1, s2):
     # add a blank string as padding
     s1 = ' ' + s1
     s2 = ' ' + s2
     
     # Scoring system
     MATCH = 5
     MISMATCH = -4
     GAP = -6

     # set table size
     NUM_ROWS = len(s2)
     NUM_COLS = len(s1)

     # Create table and fill it with zeros
     c = createTable(NUM_ROWS, NUM_COLS, 0)

     # Create table for getting back the optimal alignment, fill table with "F"
     # Suggest you use "D", "L", and "T" for diagonal, left, and top
     d = createTable(NUM_ROWS, NUM_COLS, "F")
	
     # complete this part of the function
     # implement dynamic programming algorithm for global alignment here
     # fill in entries in cost table and direction table
     for i in range(1, NUM_ROWS):
          c[i][0] = GAP*i
     for j in range(1, NUM_COLS):
          c[0][j] = GAP*j
     for i in range(1, NUM_ROWS):
          for j in range(1, NUM_COLS):
               if s1[j] == s2[i]:
                    m = MATCH
               else:
                    m = MISMATCH
               left = c[i][j-1] + GAP
               top = c[i-1][j] + GAP
               diag = c[i-1][j-1] + m
               c[i][j] = max(left,top,diag)
               if c[i][j] == diag:
                    d[i][j] = "D"
               elif c[i][j] == top:
                    d[i][j] = "T"
               else:                    
                    d[i][j] = "L"
                    
     # return optimal score (lower right-hand cell in table]
     return c[NUM_ROWS-1][NUM_COLS-1]


#################################################################
# Generate a random DNA sequence with given length and nucleotide
# probabilities
# complete this function
# A previous lab exercise may be helpful
#################################################################
def generateRandomSequence(length, A_prob, C_prob, G_prob, T_prob):
    seq = ""
    for i in range(0, length):
        ran = random.random()
        if ran < A_prob:
            seq += "A"
        elif ran < A_prob + T_prob:
            seq += "T"
        elif ran < A_prob + T_prob + C_prob:
            seq += "C"
        else:
            seq += "G"
    return seq


##################################################################
# Genenerate a random DNA sequence comparable in nucleotide
# distribution as the given input sequence
# complete this function
##################################################################
def generateComparableRandomSequence(s):
    # calculate nucleotide frequencies and call generateRandomSequence function
    totalLen = len(s)
    sList = list(s)
    percentA = sList.count('A')/totalLen
    percentC = sList.count('C')/totalLen
    percentT = sList.count('T')/totalLen
    percentG = sList.count('G')/totalLen

    # to create the random DNA sequence similar in nucleotide composition and same length
    return generateRandomSequence(totalLen, percentA, percentC, percentT, percentG) 


####################################################################
# For the two input sequences s1 and s2, calculate the distance score D where
# D is definted as 100.0*(-ln(S_norm)). 
#
# S_norm = (S_global - S_rand) / (S_iden - S_rand)
#
# S_global is the optimal global alignment score between s1 and s2.
#
# S_iden is the average of the global alignment scores of s1 aligned
# with s1 and s2 aligned with s2.
#
# S_rand is the average of 1000 global alignment scores between
# sequences similar in composition as s1 and s2.
# complete this function
###################################################################
def calculateDistanceScore(s1, s2):
    sGlobal = globalAlignmentScore(s1, s2)
    sIden = (globalAlignmentScore(s1, s1) + globalAlignmentScore(s2, s2))/2.0

    sRandSum = 0
    for _ in range(1000):
        s1Rand = generateComparableRandomSequence(s1)
        s2Rand = generateComparableRandomSequence(s2)
        sRandSum += globalAlignmentScore(s1Rand, s2Rand)
    sRand = sRandSum/1000
    sNorm = (sGlobal - sRand) / (sIden - sRand)    
    D = 100.0 * (-(math.log(sNorm)))
    return D


########################################################
### End of functions ###################################
########################################################

# DNA sequences
seq1 = "CGATAGTGCTATATCTAGCGCCGTCTAGATGCATTATACGATATCG"
seq2 = "AACGACATGGCTCGTGCTATTACGCGCGAATATCC"
seq3 = "ATAGTGCTATACTCGTGCTATTCTAGATGCCGCGATATAT"
seq4 = "GGATAGGCTATATCTAGCGCGTCTAGATGCATTTACGATATC"
seq5 = "TACGACATGCGCTCGTGCATATTAGCGCGCGATATATCG"

# calculate global alignment scores for all ten pairs
print("Global alignment scores")
print("seq1 and seq2: " + str(globalAlignmentScore(seq1, seq2)))
print("seq1 and seq3: " + str(globalAlignmentScore(seq1, seq3)))
print("seq1 and seq4: " + str(globalAlignmentScore(seq1, seq4)))
print("seq1 and seq5: " + str(globalAlignmentScore(seq1, seq5)))
print("seq2 and seq3: " + str(globalAlignmentScore(seq2, seq3)))
print("seq2 and seq4: " + str(globalAlignmentScore(seq2, seq4)))
print("seq2 and seq5: " + str(globalAlignmentScore(seq2, seq5)))
print("seq3 and seq4: " + str(globalAlignmentScore(seq3, seq4)))
print("seq3 and seq5: " + str(globalAlignmentScore(seq3, seq5)))
print("seq4 and seq5: " + str(globalAlignmentScore(seq4, seq5)))

# calculate distance score for all ten pairs
print("")
print("Distance scores")
print("seq1 and seq2: " + str(calculateDistanceScore(seq1, seq2)))
print("seq1 and seq3: " + str(calculateDistanceScore(seq1, seq3)))
print("seq1 and seq4: " + str(calculateDistanceScore(seq1, seq4)))
print("seq1 and seq5: " + str(calculateDistanceScore(seq1, seq5)))
print("seq2 and seq3: " + str(calculateDistanceScore(seq2, seq3)))
print("seq2 and seq4: " + str(calculateDistanceScore(seq2, seq4)))
print("seq2 and seq5: " + str(calculateDistanceScore(seq2, seq5)))
print("seq3 and seq4: " + str(calculateDistanceScore(seq3, seq4)))
print("seq3 and seq5: " + str(calculateDistanceScore(seq3, seq5)))
print("seq4 and seq5: " + str(calculateDistanceScore(seq4, seq5)))


