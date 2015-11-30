#!/usr/bin/python

##################################################
# CS 423 Lab 6 random alignment experiment
# Taylor Spooner, Caleb Piekstra
##################################################
import random, math


# prints histogram (list of ints) to outputFile
def printHistogram(histogram, outputFile):
    with open(outputFile, 'w') as out:
        for idx, el in enumerate(histogram):
            out.write("%d\t%d\n" % (idx, el))
    return

    
# create random DNA sequence of length n with each nucleotide randomly
# selected at rate .25
def createRandomSequence(n):
    seq = ""
    for i in range(0, n):
        ran = random.random()
        if ran < 0.25:
            seq += "A"
        elif ran < 0.50:
            seq += "T"
        elif ran < 0.75:
            seq += "C"
        else:
            seq += "G"
    return seq
        
    
# runs an experiment with numberOfAlignments pairs of random
# DNA sequences
def alignManyRandomSequences(numberOfAlignments):
    hist = [0] * 150 # python code for creating a list of 150 0’s
    for i in range(0, numberOfAlignments):
        s1 = createRandomSequence(30)
        s2 = createRandomSequence(30)
        score = localAlignmentScore(s1, s2)
        hist[score] += 1
    printHistogram(hist, "manyRandom.txt")
    return # delete stub

    
# you may copy/paste function from lab 5
# just calculate the optimal local alignment score (do not print tables
# and perform the actual alignment) of s1 and s2
def localAlignmentScore(s1, s2):
    # add a blank string as padding
    s1 = ' ' + s1
    s2 = ' ' + s2

    # Scoring system
    MATCH = 5
    MISMATCH = -4
    GAP = -6

    # Max value of alignment found
    maxValue = 0

    # To keep track of the location of the maximum alignment
    maxRowPos = 0     
    maxColPos = 0

    # set table size
    NUM_ROWS = len(s2)
    NUM_COLS = len(s1)

    # Create table and fill it with zeros
    c = createTable(NUM_ROWS, NUM_COLS, 0)

    # Creates table for getting back the optimal alignment, fill table with "F"
    # uses "D", "L", and "T" for diagonal, left, and top
    d = createTable(NUM_ROWS, NUM_COLS, "F")

    # The above automatically sets the gap penalties in left column and top row

    # implements dynamic programming algorithm for local alignment
    # fills in entries in cost table and direction table

    # update table to show alignment penalties
    for i in range(1, NUM_ROWS):
      for j in range(1, NUM_COLS):

           # set value of match/mismatch
           m = MATCH if s1[j] == s2[i] else MISMATCH

           # calculate costs based on direction
           left = c[i][j-1] + GAP
           top = c[i-1][j] + GAP
           diag = c[i-1][j-1] + m

           # set cost
           cost = max(0, left,top,diag)
           c[i][j] = cost

           # set direction based on cost
           # Precedence: F > D > T > L               
           if cost == 0:
                d[i][j] = "F"
           elif cost == diag:
                d[i][j] = "D"
           elif cost == top:
                d[i][j] = "T"
           else:
                d[i][j] = "L"

           # update max value and its position
           if cost > maxValue:
                maxValue = cost
                maxRowPos = i
                maxColPos = j
    return maxValue


# Create a 2D table with the given number of rows and columns
# and fills all entries with value
# (function completed for you)
def createTable(NUM_ROWS, NUM_COLS, value):
    table = []
    row = 0
    while (row < NUM_ROWS):
        table.append([])
        col = 0
        while (col < NUM_COLS):
            table[row].append(value)
            col = col + 1
        row = row + 1
    return table

########################################################
### End of functions ###################################
########################################################

alignManyRandomSequences(1000)

