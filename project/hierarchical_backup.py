import random
import math
import sys

##
# Class for performing hierarchical clustering to generate trees.
# Uses microarray data to output a tree where each leaf is a gene.
#
class Hierarchical:

# Constructor
#
    def __init__ (self, microarray, verbose=False):
        self.verbose = verbose
        self.mdata = microarray
        self.numGenes = len(microarray)
        self.numExps = len(microarray[0])
        self.distTable = self.initTable()

    ##
    # calcDist
    # 
    # Method for calculating distances between two genes based on
    # microarray experiment data. Takes two genes' experiment results
    # as input and returns the distance between them as a float
    #
    def calcDist (self, initial, dest):
        sumSquares = 0
        for i in range(len(initial)):
            sumSquares += abs(dest[i] - initial[i])**2
        return math.sqrt(sumSquares)

    ##
    # initTable
    # 
    # Creates the initial distance table. Calculates the distances
    # between each pair of genes and places it in the appropriate entry
    # in the table.
    #
    def initTable(self):
        # Size of square table containing distances and gene group labels
        tableDim = self.numGenes+1

        # Initialize table as a table of placeholder values
        distTable = [['X' for x in range(tableDim)] for x in range(tableDim)]

        # Fill in each of the gene group labels
        for i in range(self.numGenes):
            label = "G" + str(i+1)
            distTable[0][i] = label
            distTable[i+1][self.numGenes] = label

        # Place the distances between pairs of genes in the appropriate
        # table cells
        for i in range(self.numGenes-1):
            for j in range(i+1,self.numGenes):
                distTable[i+1][j] = self.calcDist(self.mdata[i],self.mdata[j])
                
        return distTable

    ##
    # getMinDistance
    #
    # returns the coordinates in the table of the minimum distance in
    # the table
    #
    def getMinDistance(self):
        
        tableDim = len(self.distTable)
        
        # Arbitrarily large so anything will be smaller
        minVal = sys.maxsize

        coords = (-1,-1)

        # Loop over all distance values in the array
        for i in range(len(self.distTable)-2):
            for j in range(i+1,len(self.distTable[0])-1):
                if self.distTable[i+1][j] < minVal:
                    minVal = self.distTable[i+1][j]
                    coords = (i+1,j)
        return coords

    ##
    # mergeClusters
    #
    # Combines two clusters in the table to create a new cluster.
    # Creates a new distance table with the merged clusters
    #
    def mergeClusters(self):
        # Size of a table 1 smaller than the current distance table
        tableDim = len(self.distTable)-1

        # Size of old distance table
        oldDim = tableDim + 1

        # Coordinates of minimum distance
        minimum = self.getMinDistance()
        minCol = minimum[1]
        minRow = minimum[0]

        # Initialize table as a table of placeholder values
        newTable = [['X' for x in range(tableDim)] for x in range(tableDim)]

        # Label rows and columns of new table
        col = 0
        row = 1
        for i in range(oldDim):
            if i != minCol and i != minRow - 1:
                label = self.distTable[0][i]
                newTable[0][col] = label
                newTable[row][tableDim-1] = label
                col += 1
                row += 1
        label = '{' + self.distTable[minRow][oldDim-1] + ',' + self.distTable[0][minCol] + '}'
        newTable[0][tableDim-2] = newTable[tableDim-1][tableDim-1] = label

        row = col = 1
        # Iterate over rows
        for i in range(1,oldDim-1):
            if i != minRow and i != minCol + 1:
                col = row
                # Iterate over columns
                for j in range(i,oldDim-1):
                    if j != minCol and j != minRow - 1:
                        newTable[row][col] = self.distTable[i][j]
                        col += 1
                row += 1

        clustDist1 = []
        clustDist2 = []

        # Get distances from one of our clusters to the others
        for i in range(1, minCol+1):
            if i != minRow:
                clustDist1.append(self.distTable[i][minCol])
        for i in range(minCol+1, oldDim-1):
            if i != minRow-1:
                clustDist1.append(self.distTable[minCol+1][i])

        # Get distances from the other cluster to the others
        for i in range(1, minRow):
            if i != minCol + 1:
                clustDist2.append(self.distTable[i][minRow-1])
        for i in range(minRow, oldDim-1):
            if i != minCol:
                clustDist2.append(self.distTable[minRow][i])

        for i in range(1,tableDim-1):
            newTable[i][tableDim-2] = (clustDist1[i-1] + clustDist2[i-1])/2
            
        self.distTable = newTable
        return

    
    ##
    # hierarchicalCluster
    #
    # Performs hierarchical clustering
    #
    def hierarchicalCluster(self):
        while len(self.distTable) > 2:
            print (self.distTable[0][0])
            self.mergeClusters()
            
                
            #for i in range(len(self.distTable)):
                #print(self.distTable[i])
            #print("==========================================================================")
        return self.distTable[0][0]

# TESTS
if __name__ == "__main__":
    exp1 = []
    exp2 = []
    for i in range(12):
        exp1.append(random.randint(0,10))
        exp2.append(random.randint(0,10))
        
    exps = (exp1, exp2)

    microarrayData = list(zip(*exps))
    print(microarrayData)
    h = Hierarchical(microarrayData, False)
    #print(h.mdata)
    #print("==========================================================================")
    print(h.hierarchicalCluster())
