# Caleb Piekstra and Sara Perkins
# CS423
# Lab 5: Timing Experiment

import globalAlignment, random
#import localAlignment   ## comment out once you have this module working
import time      # to calculate running times

#######
# randomSeq
# generates randomly composed DNA sequences (25%A, 25%C, 25%T, 25%G)
# of length N and returns it as a string
# returns the sequence
#######
def randomSeq(N):
    seq = ""
    for i in range(0, N):
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

########
# Timing experiment
# creates two random sequences of same length
# calculates time to find the optimal score
# writes the results to a file
# repeats for sequences of length 50 to 5000 stepping by 50
########
with open("longTimingExperiment.txt", 'w') as out:
    for N in range(50, 5001, 50):
        # generate two random sequences
        ranSeq1 = randomSeq(N)
        ranSeq2 = randomSeq(N)

        # start clock and find global score
        t0 = time.clock()
        optimalScore = globalAlignment.globalAlignmentScore(ranSeq1, ranSeq2)
        t1 = time.clock()
        
        # print the time results to the file
        totalTime = str(t1-t0)
        print("Seq Len: %d\tTotal time: %s" % (N, totalTime))
        out.write("%d\t%s\n" % (N, totalTime))
