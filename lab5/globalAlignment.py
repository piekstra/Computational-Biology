#!/usr/bin/python
# PUT YOUR NAME HERE
# CS 423, Lab 5
# fall 2015

######################################################################
# Determine the score of the optimal global alignment of two strings
# students: complete this function
######################################################################
def globalAlignmentScore(s1, s2):
     # add a blank string as padding
     s1 = ' ' + s1
     s2 = ' ' + s2
     # Scoring system
     MATCH = 5
     MISMATCH = -4
     GAP = -6

     # set table size
     NUM_ROWS = len(s2)
     NUM_COLS = len(s1)

     # Create table and fill it with zeros
     c = createTable(NUM_ROWS, NUM_COLS, 0)

     # Create table for getting back the optimal alignment, fill table with "F"
     # Suggest you use "D", "L", and "T" for diagonal, left, and top
     d = createTable(NUM_ROWS, NUM_COLS, "F")
	
     # complete this part of the function
     # implement dynamic programming algorithm for global alignment here
     # fill in entries in cost table and direction table
     for i in range(1, NUM_ROWS):
          c[i][0] = GAP*i
##          d[i][0] = "T"
     for j in range(1, NUM_COLS):
          c[0][j] = GAP*j
##          d[0][j] = "L"
     for i in range(1, NUM_ROWS):
          for j in range(1, NUM_COLS):
               if s1[j] == s2[i]:
                    m = MATCH
               else:
                    m = MISMATCH
               left = c[i][j-1] + GAP
               top = c[i-1][j] + GAP
               diag = c[i-1][j-1] + m
               c[i][j] = max(left,top,diag)
               if c[i][j] == diag:
                    d[i][j] = "D"
               elif c[i][j] == top:
                    d[i][j] = "T"
               else:                    
                    d[i][j] = "L"
                    
     #
     #
     #
	
     # Print out table (only useful for small tables - used for debugging)
     # Comment out when you are satisfied that the algorithm is working
     printTable(c, "costs.txt")
     printTable(d, "directions.txt")

     # find optimal alignment
     align(d, s1, s2, "alignment.txt")

     # return optimal score (lower right-hand cell in table]
     return c[NUM_ROWS-1][NUM_COLS-1]

####################################################################
# Create a 2D table with the given number of rows and columns
# and fills all entries with value given as a parameter
# (function completed for you)
####################################################################
def createTable(numRows, numCols, value):
     table = []
     row = 0
     # create 2D table initialized with value
     while (row < numRows):
          table.append([])    
          col = 0
          while (col < numCols):
               table[row].append(value)
               col = col + 1
          row = row + 1
     return table

##################################################################
# Print 2D table to file (only useful for small tables for short
# strings)
# Should have tabs between the values on each row
# Useful function for debugging purposes
# Complete this function
##################################################################
def printTable(table, filename):
     NUM_ROWS = len(table)
     NUM_COLS = len(table[0])
     with open(filename, 'w') as out:
          for i in range(1, NUM_ROWS):
               row = ""
               for j in range(1, NUM_COLS):
                    row += str(table[i][j]) + '\t'
##               print (row+'\n')
               out.write(row+'\n')
     return
	

################################################################
# Reconstruct the optimal alignment and print the alignment
# to a file. Because the sequences can be long, print the 
# alignment 50 characters on one line, the other string of 50 characters
# on the next line, and then skip one line, as follows:
# AATT--GGCTATGCT--C-G-TTACGCA-TTACT-AA-TCCGGTC-AGGC
# AAATATGG---TGCTGGCTGCTT---CAGTTA-TGAACTCC---CCAGGC
#
# TATGGGTGCTATGCTCG--T--TACG-CA
# TCAT--TGG---TGCTGGCTGCTT--ACA
#
# Complete this function
# direction is a 2D table, seq1 and seq2 are the original DNA
# sequences to align, and filename is the name of the output file
###############################################################
def align(direction, s1, s2, filename):
     # complete this function - find optimal alignment
     
     
     row = len(direction)
     col = len(direction[0])
     alignS1 = ""
     alignS2 = ""

     while(row >= 0 and col >= 0):
          d = direction[row][col]
          if d == "D":
               row -= 1
               col -= 1
          elif d == "T":
               row -= 1
          elif d == "L":
               col -= 1
          else:
               break
     	
     # print the strings to the output file, 50 characters at a time
	out = open(filename, 'w')
     return


### End of functions ###################################


###################################################
### Testing #######################################
###################################################

# Calculate global alignment score of two sequences
s = "AGCGTCTA"
t = "TGCATCTCG"
optimalScore = globalAlignmentScore(s, t)

print(s)
print(t)
print("Global alignment score: " + str(optimalScore))



