#!/usr/bin/python
############################################################
# CS 423 lab 6
# (put your names here)
# Finding high-scoring pairs of words; first step in BLASTp
#############################################################

# readMatrix
# read data from file and store into score matrix
# returns score matrix as 2D array and the amino acid table order
# (completed for you)
def readMatrix(inputFile):
    # create 2D scoring matrix 24 x 24 (20 AAs, plus non-determinates)
    score = [[0 for x in range(24)] for x in range(24)]
    
    f = open(inputFile)

    # remove top six lines of comments plus row of AAs from BLOSUM 62 file
    for i in range(6):
        f.readline()

    # create array for the order of amino acids
    legend = f.readline().split()
    
    # assign values in score
    row = 0
    for line in f:
        # parse the line
        items = line.split()  #splits on whitespace
        # first value in item is the amino acid letter, so we want to
        # ignore it by starting index at 1
        for i in range(1,len(items)):
            score[row][i-1] = int(items[i])
        row+=1
    # return both the 2D table of scores and the 1D legend of the AA order
    return [score, legend]

# determine match score of word1 and word2 (length 3 words)
# given [score, legend] data structure
# for example, in BLOSUM 62, the match score between "AHK" and "BHK"
# is 11
# (complete this function)
def matchScore(word1, word2, mat):
    # note: if the lengths of word1 and word2 are not 3, then
    # return 0
	if len(word1) != 3 or len(word2) != 3:
		return 0
    
	scoreTable = mat[0]
	legend = mat[1]
	score = 0
	
    # note: mat is a list of 2 items: the 2D score array and
    # the 1D legend of amino acid characters
	for idx, letter in enumerate(word1):
		# figure out what column/row the letters are in using the legend
		w1Idx = legend.index(letter)
		w2Idx = legend.index(word2[idx])
		score += scoreTable[w1Idx][w2Idx]
	
	return score
    
# compute all high-scoring words given the input word
# scoring threshold T and the [score, legend] data structure
# print all high scoring words to the output file
# in format given in lab handout
def highScoringWords(word, T, mat, outputFile):
	#aminoAcids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X']
	aminoAcids = mat[1][:-1]
	with open(outputFile, 'w') as out:
		out.write("High scoring words for %s, Threshold: %d\n" % (word, T))
		for letter1 in aminoAcids:
			for letter2 in aminoAcids:
				for letter3 in aminoAcids:
					word2 = letter1+letter2+letter3
					score = matchScore(word, word2, mat)
					if score > T:
						out.write("%s\t%d\n" % (word2, score))
	return


######## end of functions ##############

sc = readMatrix("BLOSUM62.txt")

v = matchScore("AHK", "BHK", sc)
print(v)

highScoringWords("AHK", 11, sc, "AHK_results.txt")
