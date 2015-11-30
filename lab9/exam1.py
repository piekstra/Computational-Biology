########################################
# CS 423 Midterm exam 1, problem 1
# Name: Caleb Piekstra
#
########################################

import random

validNucs = "ATCG"
    
def calcDinucleotide(seq, nuc, outFile):
    if nuc not in validNucs and nuc not in validNucs.lower():
        print ("Error: not a nucleotide")
        return
    
    dinucleotidePairs = {nuc+let:0 for let in validNucs}
    dinucCount = 0
    
    for i in range(0, len(seq)):
        dinuc = seq[i:i+2]
        if len(dinuc) < 2:
            break
        dinucCount += 1
        if dinuc in dinucleotidePairs.keys():
            dinucleotidePairs[dinuc] += 1
            
    with open(outFile, 'w') as out:
        for pair in dinucleotidePairs.keys():
            out.write("%s:\t%d\n" % (pair, dinucleotidePairs[pair]))
        out.write("Total dinuc: %d\n" % (dinucCount))
    return


def genRandom(N):
    seq = ""
    for i in range(0, N):
        ran = random.random()
        if ran < 0.4:
            seq += "A"
        elif ran < 0.7:
            seq += "T"
        elif ran < 0.9:
            seq += "C"
        else:
            seq += "G"
    return seq





################### end of functions #####
# testing for calcDinucleotide
calcDinucleotide("ATAACTTGACATTAGAAACTT", "A", "A_di.txt")
calcDinucleotide("ATAACTTGACATTAGAAACTT", "C", "C_di.txt")
calcDinucleotide("ATAACTTGACATTAGAAACTT", "G", "G_di.txt")
calcDinucleotide("ATAACTTGACATTAGAAACTT", "T", "T_di.txt")


# testing for genRandom
s = genRandom(100)
print(s)
print(s.count("A"))
print(s.count("T"))
print(s.count("C"))
print(s.count("G"))
