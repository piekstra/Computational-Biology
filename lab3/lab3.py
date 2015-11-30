#!/usr/bin/python

#########################
# Lab 3, BIO/CS423
# Names:
#########################

# library for functions for random number generation
import random
import re

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
            
####################################################################
# createFASTA - takes a string and the name of a text file and a
# file header and creates a FASTA formatted file containing a
# header line and thestring, printed out 60 characters per line 
# last line may have fewer than 60)
####################################################################
def createFASTA(sequence, filename, fileheader):
    # Open the file for writing
    with open(filename,'w') as file:
        # Write a FASTA header to the file
        file.write('> ' + fileheader + newline)

        # Write the sequence to the file, 60 characters per line
        for i in range(0, len(string),60):
            file.write(sequence[i:i+60]+newline)

def generateRandomDNA(seqLength, filename):
    sequence = ""
    for i in range (0, seqLength):
        randomNum = random.random()
        if randomNum < 0.31:
            sequence += 'A'
        elif randomNum < 0.62:
            sequence += 'T'
        elif randomNum < 0.81:
            sequence += 'C'
        else:
            sequence += 'G'
    createFASTA(sequence, 'Random DNA Sequence', 'randomDNA.txt')

#for exercise 2
#code for help in translating codons to amino acids
#to save you some typing
#use Format->Uncomment out region to use this code
#note that out is the name of the file object used for writing to
#an output file
def translate(inputfile, outputfile):
    inputSequence = convertFileToSequence(inputfile)
    with open(outputfile, 'w') as out:
        out.write("> " + inputfile + " translated to amino acids\r\n")
        for i in range (0, len(inputSequence), 3):
            codon = inputSequence[i:i+3]
            if i % 60 == 0:
                out.write("\r\n")
            if(codon == "TTT" or codon =="TTC" ):
                out.write("F")
            elif(codon == "TTA" or codon =="TTG" or codon =="CTT" or codon =="CTC" or codon == "CTA" or codon =="CTG"):
                out.write( "L")
            elif(codon == "ATT" or codon =="ATC" or codon =="ATA"):
                out.write( "I")
            elif(codon == "GTT" or codon =="GTC" or codon =="GTA"or codon =="GTG"):
                out.write( "V")
            elif(codon == "TTA" or codon =="TTG" or codon =="CTT"or codon =="CTC"or codon =="CTA" or codon =="CTG"):
                out.write( "L")
            elif(codon == "ATG" ):
                out.write("M")
            elif(codon == "TCT" or codon =="TCC" or codon =="TCA" or codon =="TCG"):
                out.write( "S")
            elif(codon == "CCT" or codon =="CCC" or codon =="CCA" or codon =="CCG"):
                out.write( "P")
            elif(codon == "ACT" or codon =="ACC" or codon =="ACA" or codon =="ACG"):
                out.write( "T" )
            elif(codon == "GCT" or codon =="GCC" or codon =="GCA" or codon =="GCG"):
                out.write( "A")
            elif(codon == "TAT" or codon =="TAC"):
                out.write( "Y")
            elif(codon == "TAA" or codon =="TAG" or codon =="TGA"):
                out.write( "X")
            elif(codon == "CAT" or codon =="CAC"):
                out.write( "H")
            elif(codon == "CAA" or codon =="CAG"):
                out.write( "Q")
            elif(codon == "AAT" or codon =="AAC" ):
                out.write( "N")
            elif(codon == "AAA" or codon =="AAG" ):
                out.write( "K")
            elif(codon == "GAT" or codon =="GAC"):
                out.write("D")
            elif(codon == "GAA" or codon =="GAG"):
                out.write("E")
            elif(codon == "TGT" or codon =="TGC"):
                out.write("C")
            elif(codon == "TGG"):
                out.write( "W" )
            elif(codon == "CGT" or codon =="CGC" or codon =="CGA" or codon =="CGG"):
                out.write( "R")
            elif(codon == "AGT" or codon =="AGC"):
                out.write( "S")
            elif(codon == "AGA" or codon =="AGG" ):
                out.write( "R")
            elif(codon == "GGT" or codon =="GGC" or codon =="GGA" or codon =="GGG"):
                out.write("G")
                
def findORFs(inputfile, outputfile):
    sequence = convertFileToSequence(inputfile)
    # m is a MatchObject
    startCodonIndexes = [m.start() for m in re.finditer('ATG', sequence)]
    stopCodonIndexes =  [m.start() for m in re.finditer('TAA', sequence)] +
                        [m.start() for m in re.finditer('TAG', sequence)] +
                        [m.start() for m in re.finditer('TGA', sequence)]
    stopCodonIndexes.sort()
    with open(outputfile, 'w') as out:
        for stopIdx in stopCodonIndexes:
            for startIdx in startCodonIndexes:
                if startIdx < stopIdx:
                    if (stopIdx - startIdx) % 3 == 0:
                        out.write('> ORF ' + str(startIdx) + '\r\n')
                        out.write(sequence[startIdx:stopIdx+3] + '\r\n')
                    
                
