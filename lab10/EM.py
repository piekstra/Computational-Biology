#!/usr/bin/python

# Lab10 CS423 fall 2015
# Expectation-Maximization Algorithm for Finding Motifs
# Sara Perkins, Caleb Piekstra

import random, math

# Suppose *sequence* is a DNA sequence, and suppose
# A1,C1,G1,T1,A2,C2,G2,T2,A3,C3,G3,T3,A4,C4,G4,T4
# represent a model for a motif of 4 nucleotides:
#
#               1       2       3       4
#       A       A1      A2      A3      A4
#       C       C1      C2      C3      C4
#       G       G1      G2      G3      G4
#       T       T1      T2      T3      T4
#
# This function calculates and returns the 4-mer
# in *sequence* that best matches the motif model.
# 
# Tammy apologizes for the long list of parameters, but
# this is the clearest way for you to know which frequency
# is which, instead of putting the data in a 2D list
# complete this function
def getBestMatchInSequence(sequence, A1, C1, G1, T1, A2, C2, G2, T2, A3, C3, G3, T3, A4, C4, G4, T4):
    A = [A1, A2, A3, A4]
    C = [C1, C2, C3, C4]
    G = [G1, G2, G3, G4]
    T = [T1, T2, T3, T4]

	# holds the best sequence
    bestseq = ""
	# holds the score of the best sequence
    bestscore = 0.0

	# loop through the sequence and determine which 4-mer has the best score
	# by multiplying the probabilties of each nuc ocurring in it's position
    for i in range(0, len(sequence)-3):
        subseq = sequence[i:i+4]
        score = 1.0
        for idx, nuc in enumerate(subseq):
            if nuc == 'A': score *= A[idx]
            elif nuc == 'C': score *= C[idx]
            elif nuc == 'G': score *= G[idx]
            else: score *= T[idx]
		# keep track of the highest scoring sequence
        if score > bestscore:
            bestseq = subseq
            bestscore = score
    
	# return the best found sequence or emptystring if none (shouldn't happen)
    return bestseq 

# Calculates and returns the frequency of the given
# nucleotide, nt, at given position in four sequences.
# If frequency is 0, return .01 (so info content calculation is
# well-defined)
# complete this function
def getFrequencyOfNucleotideAtPosition(nt, position, seq1, seq2, seq3, seq4):
    sequences = [seq1, seq2, seq3, seq4]
    count = 0
	# loop through and count the occurrences of the nuc
    for seq in sequences:
        if seq[position] == nt:
            count += 1
    frequency = count/4.0
	# Default to 0.001 so that math.log works if the 
	# sequence has a frequency of 0
    if frequency == 0.0:
        return 0.001
    return frequency

# Calculates and returns the information content of
# a 4-mer motif profile, given the 16 values in the profile
# Assumes background frequency of 25% for each nucleotide
# complete this function
def calcInfoContent(A1, C1, G1, T1, A2, C2, G2, T2, A3, C3, G3, T3, A4, C4, G4, T4):
    prof = [[A1, C1, G1, T1], [A2, C2, G2, T2], [A3, C3, G3, T3], [A4, C4, G4, T4]]
	# calculates I for each position (1-4) and then sums all of those I's 
    return sum([sum([weight * math.log(weight/.25, 2) for weight in col]) for col in prof])

#########################################################
# Prints the specified motif model to the screen.
#               1       2       3       4
#       A       A1      A2      A3      A4
#       C       C1      C2      C3      C4
#       G       G1      G2      G3      G4
#       T       T1      T2      T3      T4
#
# (function completed for you)
#########################################################
def printMotif(A1, C1, G1, T1, A2, C2, G2, T2, A3, C3, G3, T3, A4, C4, G4, T4):
    print ("\t 1 \t 2 \t 3 \t 4")
    print ("A\t" + str(A1) + "\t" + str(A2) + "\t" + str(A3) + "\t" + str(A4))
    print ("C\t" + str(C1) + "\t" + str(C2) + "\t" + str(C3) + "\t" + str(C4))
    print ("G\t" + str(G1) + "\t" + str(G2) + "\t" + str(G3) + "\t" + str(G4))
    print ("T\t" + str(T1) + "\t" + str(T2) + "\t" + str(T3) + "\t" + str(T4) + "\n")

##########################################################
# Uses expectation-maximization algorithm for finding
# the best motif in the four sequences
#
# (function completed for you)
##########################################################
def findMotif(seq1, seq2, seq3, seq4):

    # remember the motif instances from the previous iteration so we know 
    # when algorithm converges
    old_instance1 = ""
    old_instance2 = ""
    old_instance3 = ""
    old_instance4 = ""

    # randomly choose starts of 4-mers from each sequence  
    randomStart1 = random.randint(0, len(seq1)-4)
    randomStart2 = random.randint(0, len(seq2)-4)
    randomStart3 = random.randint(0, len(seq3)-4)
    randomStart4 = random.randint(0, len(seq4)-4)
    
    # create the random 4-mers from each sequence
    instance1 = seq1[randomStart1:randomStart1+4]  
    instance2 = seq2[randomStart2:randomStart2+4]
    instance3 = seq3[randomStart3:randomStart3+4]
    instance4 = seq4[randomStart4:randomStart4+4]
    
    # repeat two steps of EM until convergence.
    while(old_instance1 != instance1 or old_instance2 != instance2 or old_instance3 != instance3 or old_instance4 != instance4):
        
        # calculate a motif model for the 4 motif instances (step 1 of EM)
        A1 = getFrequencyOfNucleotideAtPosition("A", 0, instance1, instance2, instance3, instance4)
        C1 = getFrequencyOfNucleotideAtPosition("C", 0, instance1, instance2, instance3, instance4)
        G1 = getFrequencyOfNucleotideAtPosition("G", 0, instance1, instance2, instance3, instance4)
        T1 = getFrequencyOfNucleotideAtPosition("T", 0, instance1, instance2, instance3, instance4)
        
        A2 = getFrequencyOfNucleotideAtPosition("A", 1, instance1, instance2, instance3, instance4)
        C2 = getFrequencyOfNucleotideAtPosition("C", 1, instance1, instance2, instance3, instance4)
        G2 = getFrequencyOfNucleotideAtPosition("G", 1, instance1, instance2, instance3, instance4)
        T2 = getFrequencyOfNucleotideAtPosition("T", 1, instance1, instance2, instance3, instance4)
        
        A3 = getFrequencyOfNucleotideAtPosition("A", 2, instance1, instance2, instance3, instance4)
        C3 = getFrequencyOfNucleotideAtPosition("C", 2, instance1, instance2, instance3, instance4)
        G3 = getFrequencyOfNucleotideAtPosition("G", 2, instance1, instance2, instance3, instance4)
        T3 = getFrequencyOfNucleotideAtPosition("T", 2, instance1, instance2, instance3, instance4)
        
        A4 = getFrequencyOfNucleotideAtPosition("A", 3, instance1, instance2, instance3, instance4)
        C4 = getFrequencyOfNucleotideAtPosition("C", 3, instance1, instance2, instance3, instance4)
        G4 = getFrequencyOfNucleotideAtPosition("G", 3, instance1, instance2, instance3, instance4)
        T4 = getFrequencyOfNucleotideAtPosition("T", 3, instance1, instance2, instance3, instance4)

        # print the motif model to the screen (comment out when running findBestMotif)
        #printMotif(A1, C1, G1, T1, A2, C2, G2, T2, A3, C3, G3, T3, A4, C4, G4, T4)

        # re-assign old instances as the current instances to determine convergence
        old_instance1 = instance1
        old_instance2 = instance2
        old_instance3 = instance3
        old_instance4 = instance4
        
        # find best match in each sequence (step 2 of EM)     
        instance1 = getBestMatchInSequence(seq1, A1, C1, G1, T1, A2, C2, G2, T2, A3, C3, G3, T3, A4, C4, G4, T4)
        instance2 = getBestMatchInSequence(seq2, A1, C1, G1, T1, A2, C2, G2, T2, A3, C3, G3, T3, A4, C4, G4, T4)
        instance3 = getBestMatchInSequence(seq3, A1, C1, G1, T1, A2, C2, G2, T2, A3, C3, G3, T3, A4, C4, G4, T4)
        instance4 = getBestMatchInSequence(seq4, A1, C1, G1, T1, A2, C2, G2, T2, A3, C3, G3, T3, A4, C4, G4, T4)
        
        # print out each of the 4 sequences and instances (comment out when running findBestMotif)         
##        print (seq1 + "\t" + instance1)
##        print (seq2 + "\t" + instance2)
##        print (seq3 + "\t" + instance3)
##        print (seq4 + "\t" + instance4 + "\n")

    return [A1, C1, G1, T1, A2, C2, G2, T2, A3, C3, G3, T3, A4, C4, G4, T4]

##################################################################
# Runs findMotif 10000 times and returns the motif
# with the highest information content
# When calling this function, be sure to comment out
# printing in the findMotif function
# (function completed for you)
##################################################################
def findBestMotif(seq1, seq2, seq3, seq4):
    bestInfoContent = 0.0
    bestModel = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]    #stored 2D table as a 1-dimensional list
    for i in range(10000):
        m = findMotif(seq1, seq2, seq3, seq4)
        infoContent = calcInfoContent(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11], m[12], m[13], m[14], m[15])

        # keep best motif found so far
        if infoContent > bestInfoContent:
                bestModel = m
                bestInfoContent = infoContent
	
    # display information about best motif
    m = bestModel
    # print best motif
    printMotif(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11], m[12], m[13], m[14], m[15])

    # find best match in each sequence     
    instance1 = getBestMatchInSequence(seq1, m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11], m[12], m[13], m[14], m[15])
    instance2 = getBestMatchInSequence(seq2, m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11], m[12], m[13], m[14], m[15])
    instance3 = getBestMatchInSequence(seq3, m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11], m[12], m[13], m[14], m[15])
    instance4 = getBestMatchInSequence(seq4, m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11], m[12], m[13], m[14], m[15])
            
    # print out each of the 4 sequences and instances         
    print (seq1 + "\t" + instance1)
    print (seq2 + "\t" + instance2)
    print (seq3 + "\t" + instance3)
    print (seq4 + "\t" + instance4 + "\n")

    # print info content
    print ("Information content: " + str(bestInfoContent))

######################################
# Run functions here #
######################################

# search for 4-mer (length 4) motif in the following 4 sequences
seq1 = "GTATACGATGTCTAGTATCAGCGGCATTAG"
seq2 = "TAGCTGTACGTAGCGGCTTTAGCTGCAT"
seq3 = "GACAGTCAGCGTTAGCTATATGCT"
seq4 = "GCAGCAGTTGAGCAGCGATGATTTATCG"

findBestMotif(seq1, seq2, seq3, seq4)
# findBestMotif(seq1, seq2, seq3, seq4)	# be sure to comment out printing when running findBestMotif

