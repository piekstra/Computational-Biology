import math

#########################################################################
# QT (class)
#
# Author: Triton Pitassi
#
# Last Modified: 12/2/2015
#
# Description: This class contains the necessary funtions to run the
# QT Clustering algorithm. 
#########################################################################
class QT:

    def __init__(self, verbose=False):
        self.verbose = verbose
        self.divider = "================================================================"

    #########################################################################
    # QTClustering
    #
    # Description:
    #       This is the main function for implementing the QT Clustering
    # algorithm. It uses the indicies of the genes in order to perform the
    # clustering. Creates and uses a distance table and a diameter truth
    # table in order to determine the clusters. It loops through the set S
    # containing all indicies until S is empty then returns the final
    # clusters that were found
    #       
    # Parameters:
    #       self - the object pointer
    #       data - a list of lists where each inner list contains the
    #              experiment data for one gene
    #       diameter - the maximum diameter of each cluster
    #
    # Return:
    #       a list of lists where each inner list is a cluster
    #########################################################################
    def QTClustering(self, data, diameter):

        # number of genes
        numGenes = len(data)

        # make set S
        setS = self.createSetS(numGenes)

        # initialize list of final clusters
        finalClusters = []

        # create distance table
        distTable = self.createDistanceTable(self.createTable(numGenes, numGenes, 0), data, numGenes)

        # create table to tell if genes are within diameter d
        diameterTable = self.createDiameterTable(self.createTable(numGenes, numGenes, True), distTable, diameter, numGenes)

        if self.verbose:
            print("Total number of genes:", numGenes, "\n")
            print("Distance Table")
            self.printDistTable(distTable)
            print("\t")
            print("Diameter Table")
            self.printDiameterTable(diameterTable)
            print("\n" + self.divider + "\n")

        # while there are still items left to cluster
        while setS:
            
            # find each gene's personal cluster
            clusters = self.findClusters(diameterTable, distTable, setS)

            # determine the biggest cluster
            biggestCluster = self.findBiggestCluster(clusters)

            # add list to list of final clusters
            finalClusters.append(biggestCluster)

            if self.verbose:
                print("The set of genes left to cluster:", setS, "\n")
                
                print("The clusters for this iteration are:")
                for cluster in clusters:
                    print(str(clusters[0]) + ":" + str(cluster))

                print("\nThe biggest cluster for this iteration is:", biggestCluster, "\n")
                print(self.divider, "\n") 
            # remove genes that have been clustered
            setS = self.removeClusteredGenes(setS, biggestCluster)

        # return the final clusters
        return [sorted(cluster) for cluster in finalClusters]


    #########################################################################
    # createSetS
    #
    # Description: creates a list of the indicies of each gene to use in
    # the QT algorithm
    #
    # Parameters:
    #       self - the object pointer
    #       numGenes - the number of genes to cluster
    #
    # Return:
    #       a list of ints representing gene indicies
    #########################################################################
    def createSetS(self, numGenes):

        # initialize the list
        setS = []

        # put the indicies in setS
        for i in range(0, numGenes):

            setS.append(i)

        return setS


    #########################################################################
    # createDistanceTable
    #
    # Description: make a table that holds the distances between each pair
    # of genes. The indicies represent genes, so the distance at [i][j] is
    # the same as the distance at [j][i].
    #
    # Parameters:
    #       self - the object pointer
    #       distTable - the table to fill
    #       data - the data set represented by a list of lists
    #       numGenes - the number of genes to cluster
    #
    # Return:
    #       a list of lists representing a table of distances
    #########################################################################
    def createDistanceTable(self, distTable, data, numGenes):

        #find all of the distances
        for i in range(0, numGenes):
            for j in range(0, numGenes):

                # calculate the distance
                distTable[i][j] = self.dist(data[i], data[j])

        # return the distance table
        return distTable


    #########################################################################
    # createDiameterTable
    #
    # Description: make a table that holds a boolean telling whether the
    # distance between a pair of genes is less than or equal to the
    # diameter i.e. if they could possibly be in the same cluster. The
    # indicies represent genes, so the boolean at [i][j] is the same as the
    # distance at [j][i].
    #
    # Parameters:
    #       self - the object pointer
    #       diameterTable - the table to fill
    #       distTable - the table that contains distances between pairs of
    #                   genes
    #       diameter - the diameter of each cluster
    #       numGenes - the number of genes to cluster
    #
    # Return:
    #       a list of lists representing a table of booleans
    #########################################################################
    def createDiameterTable(self, diameterTable, distTable, diameter, numGenes):

        for i in range(0, numGenes):
            for j in range(0, numGenes):

                # if the distance between the genes is greater than the diameter, set to False
                if distTable[i][j] > diameter:
                    
                    diameterTable[i][j] = False

        # return the diameter table
        return diameterTable


    #########################################################################
    # findClusters
    #
    # Description: finds each gene's personal clusters. Starts with the
    # closest gene and adds the next closest until the distance exceeds the
    # diameter.
    #
    # Parameters:
    #       self - the object pointer
    #       diameterTable - the table telling whether the distance between
    #                       two genes is less than the diameter
    #       distTable - the table that contains distances between pairs of
    #                   genes
    #       setS - a list of indicies representing the genes to cluster
    #
    # Return:
    #       a list of lists where each inner list is a cluster
    #########################################################################
    def findClusters(self, diameterTable, distTable, setS):

        # initialize list of clusters
        clusters = []

        # for every gene index left in S
        for i in range(0, len(setS)):

            # copy S to a list so that we can remove items from it without changing S
            tempSetS = list(setS)

            # remove the gene
            tempSetS.remove(setS[i])
            clusters.append([setS[i]])

            # initialize boolean to keep track of whether we are done finding this gene's personal cluster
            fitInCluster = True

            # while there is a gene that fits in the diameter of the cluster
            while fitInCluster:

                # if there aren't any genes left to possibly put in a cluster,
                # then we're done with this cluster
                if not tempSetS:

                    break

                # initialize variables for finding the closest gene
                closestGene = -1
                closestGeneDist = float("inf")

                # find closest gene
                for gene in tempSetS:

                    # if we haven't found a gene yet or this gene is closer than the closest we've found, make it the closest gene
                    if closestGene == -1 or distTable[setS[i]][gene] < closestGeneDist:

                        closestGene = gene
                        closestGeneDist = distTable[setS[i]][gene]

                # determine if the closest gene fits within the diameter of the cluster
                for gene in clusters[i]:

                    # if the distance is greater than the diameter, this gene doesn't fit
                    if not diameterTable[gene][closestGene]:

                        fitInCluster = False

                        # since we know the gene doesn't fit, we don't have to check the rest of the genes in the cluster
                        break

                # if the next closest gene fits, put it in the cluster and remove it from our temporary S
                if fitInCluster:

                    clusters[i].append(closestGene)
                    tempSetS.remove(closestGene)

        # return a list with each gene's personal cluster
        return clusters


    #########################################################################
    # findBiggestCluster
    #
    # Description: finds the biggest cluster in a group of clusters
    #
    # Parameters:
    #       self - the object pointer
    #       clusters - a list of each gene's personal cluster (only includes
    #                  genes left to cluster)
    #
    # Return:
    #       a list representing a cluster
    #########################################################################
    def findBiggestCluster(self, clusters):

        # initialize variables for finding biggest cluster
        biggestCluster = []
        lenBiggestCluster = 0

        # for each cluster, determine if it is larger than the biggest cluster we've found so far
        for cluster in clusters:

            if len(cluster) > lenBiggestCluster:

                lenBiggestCluster = len(cluster)
                biggestCluster = cluster

        # return the biggest cluster found
        return biggestCluster


    #########################################################################
    # removeClusteredGenes
    #
    # Description: removes genes from set S if they have been put in the
    # group of final clusters
    #
    # Parameters:
    #       self - the object pointer
    #       setS - a list of indicies representing the genes to cluster
    #       biggestCluster - the biggest cluster found in this iteration
    #
    # Return:
    #       set S with the genes in the biggest cluster removed
    #########################################################################
    # we've found the biggest cluster for this iteration so remove those gene indicies from S
    def removeClusteredGenes(self, setS, biggestCluster):

        for geneIdx in biggestCluster:

            setS.remove(geneIdx)

        return setS
            

    #########################################################################
    # createTable
    #
    # Description: Create a 2D table with the given number of rows and columns
    # and fills all entries with value given as a parameter
    #
    # Parameters:
    #       self - the object pointer
    #       numRows - the number of rows to put in the table
    #       numCols - the number of columns to put in the table
    #       value - the temporary value to store in each cell
    #
    # Return:
    #       a table with with the specified number of rows and columns
    #########################################################################
    def createTable(self, numRows, numCols, value):

        return [list([value]*numCols) for _ in range(numRows)]


    #########################################################################
    # printDistTable
    #
    # Description: print the distance table to the console
    #
    # Parameters:
    #       self - the object pointer
    #       table - the distance table to print
    #
    # Return:
    #       void
    #########################################################################
    def printDistTable(self, table):
        # get table dimensions
        NUM_ROWS = len(table)
        NUM_COLS = len(table[0])

        # print top row
        print(" \t", end="")
        for i in range(0, NUM_ROWS):
            print(repr(i).ljust(1), "\t", end="")

        print("\t")

        # print table contents
        for i in range(0, NUM_ROWS):
            print(i, "\t", end="")
            for j in range(0, NUM_COLS):
                print(repr(float("{0:.2f}".format(table[i][j]))).ljust(1), "\t", end="")

            print("\t")

        return


    #########################################################################
    # printDiameterTable
    #
    # Description: print the diameter table to the console
    #
    # Parameters:
    #       self - the object pointer
    #       table - the diameter table to print
    #
    # Return:
    #       void
    #########################################################################
    def printDiameterTable(self, table):
        # get table dimensions
        NUM_ROWS = len(table)
        NUM_COLS = len(table[0])

        # print top row
        print(" \t", end="")
        for i in range(0, NUM_ROWS):
            print(repr(i).ljust(2), "\t", end="")

        print("\t")

        # print table contents
        for i in range(0, NUM_ROWS):
            print(i, "\t", end="")
            for j in range(0, NUM_COLS):
                print(repr(table[i][j]).ljust(2), "\t", end="")

            print("\t")

        return

    #########################################################################
    # dist
    #
    # Description: Calculates the distance between two gene vectors.
    #
    # Parameters:
    #       self - the object pointer
    #       gene1 - the first gene vector
    #       gene2 - the second gene vector
    #
    # Return:
    #       the distance between the two gene vectors
    #########################################################################
    def dist(self, gene1, gene2):

        # initialize distance
        distance = 0

        # calculate (x2 - x1)^2 for each term
        for i in range(0, len(gene1)):
            
            distance += (gene1[i] - gene2[i])**2

        # calculate square root and return
        return math.sqrt(distance)


if __name__ == "__main__":
    
    # UNIT TESTS
    gene1 = [2, 5, 1]
    gene2 = [6, 3, 2]
    gene3 = [1, 1, 3]
    gene4 = [5, 3, 6]
    gene5 = [2, 3, 3]
    gene6 = [4, 4, 2]

    genes = []
    genes.append(gene1)
    genes.append(gene2)
    genes.append(gene3)
    genes.append(gene4)
    genes.append(gene5)
    genes.append(gene6)

    qt = QT(True)
    print("The final clusters are:", qt.QTClustering(genes, 5))


