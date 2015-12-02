import random
# Graphing imports
import matplotlib.pyplot as plt
import matplotlib.colors as colors

## 
#
# Description: Runs the KMeans clustering algorithm
#
# K-means Clustering: Given N items with known distances between 
#   items and K clusters, assigns the N items into one of K clusters
#   such that the total distance from each item to its cluster
#   is minimized.
#
class KMeans:
    
    ##
    # constructor
    #
    # Description: Initializes instance variables used for the algorithm
    #
    # Parameters:
    #   self - The object pointer
    #   microarray - A 2D list where each sublist is a gene containing
    #       its expression data
    #   k - The upper bound on the number of initial clusters
    #   verbose - Optional parameter to print out information
    #
    def __init__ (self, verbose=False):
        # whether or not to print extra information to console
        self.verbose = verbose
    
    
    ##    
    # center
    #
    # Description: Given a cluster of genes, calculates and
    #   returns the center of the cluster
    #
    # Parameters:
    #   cluster - A cluster of genes
    #
    # Returns:
    #   The calculated center of the cluster
    #
    def center (self, cluster):
        # get the length of the cluster
        listLen = float(len(cluster))
        # return a list containing the center of the cluster
        return [sum(exp)/listLen for exp in zip(*cluster)]
    
    
    ##    
    # sqDist
    #
    # Description:
    #
    # Parameters:
    #
    #   self - the object pointer
    #   initial - 
    #
    # Returns:
    #   The square distance between the initial gene and dest
    #   gene
    #
    def sqDist (self, initial, dest):
        # holds the square distance between two genes
        sumSquares = 0
        
        # loop for the length of the initial gene
        for i in range(len(initial)):
            # square the distance between each expression in
            # the destination and initial genes and add it 
            # to the total sum
            sumSquares += abs(dest[i] - initial[i])**2
            
        # return the square distance
        return sumSquares
    
    
    ##    
    # newClusters
    #
    # Description:
    #
    # Parameters:
    #
    # Returns:
    #
    def newClusters(self, len=None):
        return [[] for _ in range(len)]
        

    ##    
    # plot
    #
    # Description: Thickens
    #
    # Parameters:
    #
    # Returns:
    #
    def plot(self, clusters, figure):
        # only allow plotting of clusters with
        # genes with 2 expressions (2D graph)
        if self.numExps != 2:
            return
        # remove the gene label from each gene
        refinedClusters = [[gene[0] for gene in cluster] for cluster in clusters]
        # modify clusters to have all x's in one tuple and all
        # y's in the other tuple (seperate exp1 and exp2)
        xyClusters = [list(zip(*cluster)) for cluster in refinedClusters]
        clusterColors = colors.cnames.values()[:len(xyClusters)]
        plt.figure(figure)
        for idx, cluster in enumerate(xyClusters):
            plt.scatter(cluster[0], cluster[1], color=clusterColors[idx])
        
        
    ##    
    # kmeans
    #
    # Description: 
    #
    # Parameters:
    #
    # Returns:
    #        
    def kmeans (self, mdata, k):
        # the number of genes in the microarray
        numGenes = len(mdata)
        # the number of expressions per gene
        numExps = len(mdata[0])
        
        # 
        if self.verbose:
            print "K-means Clustering:\n\tk = %d\n\tnum genes: %d\n\tnum expressions: %d" % (k, numGenes, numExps)
            
        clusters = self.newClusters(k)
        # for each item, randomly assign it to one of the k clusters
        for geneIdx, gene in enumerate(mdata):
            clusters[random.randint(0,k-1)].append((gene, geneIdx))
        # remove any empty clusters
        clusters = filter(None, clusters)
        counter = 0
        while True:
            counter += 1
            if numExps == 2:
                self.plot(clusters, counter)
            newClusters = self.newClusters(len(clusters))
            centers = [self.center(zip(*cluster)[0]) for cluster in clusters]
            if self.verbose:
                print "\ncluster arrangement %d:" % (counter)
                print "Gene%s\t%s\t cluster assignment" % (' '.join(['  ' for _ in range(2,numExps)]), '\t'.join(["sq dist to center C%d" % i for i in range(1, len(centers)+1)]))
            for geneIdx, gene in enumerate(mdata):
                distancesToCenters = [self.sqDist(gene, center) for center in centers]
                assignedClusterIdx = min(xrange(len(distancesToCenters)), key=distancesToCenters.__getitem__)
                newClusters[assignedClusterIdx].append((gene, geneIdx))
                if self.verbose:
                    print "%s\t%s\t\t\t %s" % (gene, '\t\t\t'.join(["%0.2f" % dist for dist in distancesToCenters]), ("C%d" % (assignedClusterIdx+1)))
            # remove any empty clusters
            newClusters = filter(None, newClusters)
            if newClusters == clusters:
                return clusters
            else:
                clusters = newClusters
        
        # display the data in figures if there are 2 expressions 
        # (2D plot!)
        if numExps == 2:
            plt.show()


# Test code for directly running the file
if __name__ == "__main__":                
    # preprocessing (input)    
    exp1 = [1, 6, 5, 2, 6, 0, 1, 1, 2, 4, 5, 6, 7]
    exp2 = [4, 2, 3, 5, 4, 5, 3, 1, 2, 4, 5, 6, 7]
    exps = (exp1, exp2)

    # preliminary data
    microarrayData = list(zip(*exps))
    k = 4

    # initialize the kmeans object
    kmeans = KMeans(k, verbose=True)
    
    # run the algorithm!
    finalClusters = kmeans.kmeans(microarrayData)
    
    # print out the results
    print "\nFinal set of gene clusters:"
    for clusterIdx, cluster in enumerate(finalClusters):
        print "\tCluster %d: %s" % (clusterIdx+1, ["gene" + str(idx+1) for gene, idx in cluster])
    print ""
    if kmeans.numExps == 2:
        plt.show()
