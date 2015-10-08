import globalAlignment, random
#import localAlignment   ## comment out once you have this module working
import time      # to calculate running times

# generates randomly composed DNA sequences (25%A, 25%C, 25%T, 25%G)
# of length N and returns it as a string
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

    
with open("longTimingExperiment.txt", 'w') as out:
    for N in range(50, 5001, 50):
        ranSeq1 = randomSeq(N)
        ranSeq2 = randomSeq(N)
        t0 = time.clock()
        optimalScore = globalAlignment.globalAlignmentScore(ranSeq1, ranSeq2)
        t1 = time.clock()

        totalTime = str(t1-t0)
        print("Seq Len: %d\tTotal time: %s" % (N, totalTime))
        out.write("%d\t%s\n" % (N, totalTime))