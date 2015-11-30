#!/usr/bin/python

# Starter Code
# BIO/CS 423 Lab 1
# Fall 2015
# Caleb Piekstra, Shaun Stice

####################################################################
# Exercise 0
# convertFileToSequence - takes a FASTA file and returns the
# sequence as a string
#####################################################################
import re
newline = '\r\n'

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

def calcGC(filename):
    sequence = convertFileToSequence(filename)
    countGC = sequence.count("G") + sequence.count("C")
    percentage = (float(countGC) / float(len(sequence))) * 100
    print (str(round(percentage, 2)) + "% GC")

def findAmbChars(filename):
    sequence = convertFileToSequence(filename)
    Goodchar = "ACGT"
    for i in range(0, len(sequence)):
         if sequence[i] not in Goodchar or sequence[i] not in Goodchar.lower():
             print(str(i)+":"+sequence[i])

def printCodons(filename):
    sequence = convertFileToSequence(filename)
    for i in range(0, len(sequence), 3):
        print (sequence[i:i+3])

def createFASTA(string,filename):
    with open(filename,'w') as file:
        file.write('> '+newline)
        for i in range(0, len(string),60):
            file.write(string[i:i+60]+newline)

def createReverseComplement(strSequence):
    return strSequence.replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c").upper()[::-1]

def findSubsequence(sequence,subseq,filename):
    with open(filename,'a') as file:        
        file.write(filename+newline)
        file.write(subseq+newline)
        startingLocations = [match.start() for match in re.finditer(subseq, sequence)]
        for startingLocation in startingLocations:
            file.write(str(startingLocation)+newline)
    

def findSubsequenceReverseComp(sequence, subseq, filename):
    with open(filename,'w') as file:        
        file.write("Reverse Complement" + newline)
        findSubsequence(createReverseComplement(sequence), subseq, filename)
                          
##### Code to test and run
yeastFile="P:\CS423\lab1\yeastGenome.txt"
humanFile="P:\CS423\lab1\human22.txt"
ecoliFile="P:\CS423\lab1\eColiGenome.txt"
sequence="P:\CS423\lab1\sequence.txt"
ambigFile="P:\CS423\lab1\\findAmbig.txt"
rad55="P:\CS423\lab1\\rad55.txt"

#1-a
#calcGC(yeastFile)
#1-b
#calcGC(ecoliFile)
#1-c
#calcGC(humanFile)

#2
#findAmbChars(ambigFile)

#3
#printCodons(rad55)

#4
#import string, random
#chars = string.ascii_uppercase
#createFASTA((''.join([chars[int(random.random()*len(chars))] for i in range(0, 80)])), 'P:\CS423\lab1\\random.txt')

#5-a
#file="P:\\CS423\\lab1\\chr05.txt"
#seq = convertFileToSequence(file)
#findSubsequence(seq, 'TATAAA', "P:\\CS423\\lab1\\chr05_tata_subs.txt")
#5-b
#file="P:\\CS423\\lab1\\eColiGenome.txt"
#seq = convertFileToSequence(file)
#findSubsequence(seq, 'TATAAA', "P:\\CS423\\lab1\\ecoli_tata_subs.txt")

#6-a
#file="P:\\CS423\\lab1\\chr05.txt"
#seq = convertFileToSequence(file)
#findSubsequenceReverseComp(seq, 'TATAAA', "P:\\CS423\\lab1\\chr05_tata_reverse_subs.txt")
#6-b
#file="P:\\CS423\\lab1\\eColiGenome.txt"
#seq = convertFileToSequence(file)
#findSubsequenceReverseComp(seq, 'TATAAA', "P:\\CS423\\lab1\\ecoli_tata_reverse_subs.txt")

#7-a
#seq = convertFileToSequence(yeastFile)
#findSubsequence(seq, 'TATACA', "P:\\CS423\\lab1\\chr05_tata_subs_tataca.txt")
#findSubsequence(seq, 'TATAGA', "P:\\CS423\\lab1\\chr05_tata_subs_tataga.txt")
#findSubsequence(seq, 'TATATA', "P:\\CS423\\lab1\\chr05_tata_subs_tatata.txt")
#7-b
#findSubsequenceReverseComp(seq, 'TATACA', "P:\\CS423\\lab1\\chr05_tata_reverse_subs_tataca.txt")
#findSubsequenceReverseComp(seq, 'TATAGA', "P:\\CS423\\lab1\\chr05_tata_reverse_subs_tataga.txt")
#findSubsequenceReverseComp(seq, 'TATATA', "P:\\CS423\\lab1\\chr05_tata_reverse_subs_tatata.txt")

#8
seq = convertFileToSequence(yeastFile)
findSubsequence(seq, 'CCAAT', "P:\\CS423\\lab1\\chr05_tata_subs_ccaat.txt")


#seq = convertFileToSequence("P:\CS423\lab1\yeastGenome.txt")
#print(seq)
