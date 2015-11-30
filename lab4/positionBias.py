#!/usr/bin/python

##########################################################
# CS423 Lab 4, fall 2015
# 
# Authors: Caleb Piekstra, Connor Haas
##########################################################

########################################################
# convertFileToSequence -- converts FASTA file to string
# from lab 1 
########################################################
def convertFileToSequence(filename):

    # read in file 
    file = open(filename)

    # read in first line
    header = file.readline()
    if (header[0] == '>'):
        print("in FASTA format")
    else:
        print("invalid format")
        return
  
    # read in rest of file
    sequence = file.read()
    
    # close file
    file.close()

    # remove all newline characters
    sequence = sequence.replace("\n", "")
    sequence = sequence.replace("\r", "") # needed to get rid of some newline characters
    return sequence

            
########################################################
# convertFileToSequence -- converts FASTA file to string
# from lab 1 
########################################################        
def codonPosition(inputFile1, inputFile2):
    # convert the contents of each file to a string
    sequence1 = convertFileToSequence(inputFile1)
    sequence2 = convertFileToSequence(inputFile2)
    
    # determine the total number of codons
    numCodons1 = int(len(sequence1) / 3)
    numCodons2 = int(len(sequence2) / 3)
    
    # total number of codons
    totalNumCodons = numCodons1 if numCodons1 < numCodons2 else numCodons2
    # get the shortest sequence
    shortestSeq = sequence1 if len(sequence1) < len(sequence2) else sequence2
        
    # Make these floats initially so that when we 
    # calculate a percentage they are in float form
    totalFirstDiffs = 0.0
    totalSecondDiffs = 0.0
    totalThirdDiffs = 0.0

    # Loop through and compare the first, second, and third
    # positions in the codons of the two sequcnes, counting
    # the number of differences for each position
    for i in range (0, len(shortestSeq), 3):
        # abort early if the loop reaches a codon that isn't len 3
        if (len(shortestSeq[i:i+3]) < 3):
            break

                # if the first position in the two codons is different,
                # increment the total number of differences in position 1
        if sequence1[i] != sequence2[i]:
            totalFirstDiffs += 1
            
                # if the second position in the two codons is different,
                # increment the total number of differences in position 2
        if sequence1[i+1] != sequence2[i+1]:
            totalSecondDiffs += 1
            
                # if the third position in the two codons is different,
                # increment the total number of differences in position 3
        if sequence1[i+2] != sequence2[i+2]:
            totalThirdDiffs += 1

    # calculate the percent frequency of differences for each position
    percentFirstDiffFrequency = (totalFirstDiffs / totalNumCodons) * 100
    percentSecondDiffFrequency = (totalSecondDiffs / totalNumCodons) * 100
    percentThirdDiffFrequency = (totalThirdDiffs / totalNumCodons) * 100

    # print out the results
    print ("Position 1 Different: %0.2f%%" % percentFirstDiffFrequency)
    print ("Position 2 Different: %0.2f%%" % percentSecondDiffFrequency)
    print ("Position 3 Different: %0.2f%%" % percentThirdDiffFrequency)

# Run the function using the seq1.txt and seq2.txt files (for the write-up)
codonPosition("seq1.txt", "seq2.txt")    
