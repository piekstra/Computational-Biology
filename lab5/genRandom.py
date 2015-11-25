# Caleb Piekstra and Sara Perkins
# CS423
# Lab 5: genRandom

#!/usr/bin/python

import random

#############
# genRandom -- generates a random sequence of input length
# with a choosing rate provided in the form of percent GC and AT
# content desired, writes the sequnece to a file
#############
def genRandom(lenStr, percentAT, percentCG, outputFileName):
    # open file, set fasta header
    with open(outputFileName, 'w') as out:
        out.write("> random sequence\r\n")

        # write the sequence to the file
        for i in range(0, lenStr):
            ran = random.random()

            # choose the nucleotide based on the given desired content
            # and random number produced
            if ran < percentAT/2:
                out.write("A")
            elif ran < percentAT:
                out.write("T")
            elif ran < percentAT + (percentCG / 2):
                out.write("C")
            else:
                out.write("G")

######## code for yeast and fruitfly needed for writeup                
#genRandom(1404, 1 - 0.38, 0.38, "yeast_random.txt")
#genRandom(1368, 1 - 0.42, 0.42, "fruit_fly_random.txt")

        
