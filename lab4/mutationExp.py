#!/usr/bin/python

##########################################################
# CS423 Lab 4, fall 2015
# 
# Caleb Piekstra
# Connor Haas
##########################################################

import random
import math

## generateRandomChar
#
# Returns G C T or A randomly
#
def generateRandomChar():
    randomNum = random.random()
    if (randomNum < 0.25):
        return "G"
    elif (randomNum < 0.50):
        return "C"
    elif (randomNum < 0.75):
        return "A"
    return "T"

## generateRandomSequence
#
# Given a length for the sequence to generate, 'N',
# returns a sequence of that length containing
# random bases
#
def generateRandomSequence(N):
    count = 0
    sequence = ""
    while (count < N):
        sequence = sequence + generateRandomChar()
        count = count + 1
  
    return sequence   

## mutate
#
# Given a string sequence, picks a random
# nucleotide in the sequence and changes it
# to one of the three other bases, randomly
#
def mutate(seq):
    seq = list(seq)
    nucForC = "ATG"
    nucForG = "ATC"
    nucForT = "ACG"
    nucForA = "TCG"
    
    randomChosen = int(random.random()*len(seq))
    mutatedChar = seq[randomChosen]
    mutationDecide = random.randint(0,2)
    if(mutatedChar == "C"):
        seq[randomChosen] = nucForC[mutationDecide]
    if(mutatedChar == "T"):
        seq[randomChosen] = nucForT[mutationDecide]
    if(mutatedChar == "G"):
        seq[randomChosen] = nucForG[mutationDecide]
    if(mutatedChar == "A"):
        seq[randomChosen] = nucForA[mutationDecide]
        
    return "".join(seq)	# stub (remove when you complete the code)

## hammingDistance
#
# Given two strings, calculates and returns the
# hamming distance between them. (Count of differences)
#
def hammingDistance(s1, s2):
    countBad = 0
    if(len(s1) != len(s2)):
        return -1
    for idx, nuc in enumerate(s1):
        if nuc != s2[idx]:
            countBad += 1
       
    return countBad	


## jukesCantor
#
# Given a hamming distance and the length of a sequence,
# calculates and returns the jukes-cantor correction
#
def jukesCantor(H, L):        
    return (L) * (-3.0/4.0) * math.log(1 - (4.0/3.0)*(H/L))

## experiment
#
# Given a sequence, number of mutations to perform on
# that sequence, and an output file:
#   randomly mutates a single nucleotide in the sequence
#   as many times as is specified by expLength
# After each mutation, the hamming distance and jukes-
# cantor correction are printed out to the outputFile
#
def experiment(seq, expLength, outputFile):
    with open(outputFile, 'w') as out:     #######this line is awesome for file IO
        mutatedSeq = seq
        for i in range(0, expLength):
            mutatedSeq = mutate(mutatedSeq)
            Ham = hammingDistance(mutatedSeq, seq)
            Cheese = jukesCantor(Ham, len(mutatedSeq))
            out.write("%d\t%0.4f\r\n" % (Ham,Cheese))
    return	

########################################################
### End of functions ###################################
########################################################

# Run mutation experiment
s = generateRandomSequence(1000)
experiment(s, 1000, "mutationResultsTest.txt")
