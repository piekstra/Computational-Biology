#!/usr/bin/python

# Starter Code
# BIO/CS 423 Lab 1
# Fall 2015
# Caleb Piekstra, 

####################################################################
# Exercise 0
# convertFileToSequence - takes a FASTA file and returns the
# sequence as a string
#####################################################################
def convertFileToSequence(filename):
    # read in file
    file = open(filename)

    # read in first line
    header = file.readline()
    if (header[0] == '>'):
        print("in FASTA format")
    else:
        print("invalid format")
        file.close()
        return
  
    # read in rest of file
    sequence = file.read()
    
    # close file
    file.close()

    # remove all return and newline characters
    sequence = sequence.replace("\r", "")
    sequence = sequence.replace("\n", "")
    return sequence

##### Code to test and run
seq = convertFileToSequence("sequence.txt")
print(seq)



