#!/usr/bin/python

import random

def genRandom(lenStr, percentAT, percentCG, outputFileName):
    with open(outputFileName, 'w') as out:
        out.write("> random sequence\r\n")
        for i in range(0, lenStr):
            ran = random.random()
            if ran < percentAT/2:
                out.write("A")
            elif ran < percentAT:
                out.write("T")
            elif ran < percentAT + (percentCG / 2):
                out.write("C")
            else:
                out.write("G")
                
genRandom(1404, 1 - 0.38, 0.38, "yeast_random.txt")
genRandom(1368, 1 - 0.42, 0.42, "fruit_fly_random.txt")

        