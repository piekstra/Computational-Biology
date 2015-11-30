import random
import matplotlib.pyplot as plt
import matplotlib.colors as colors

class KMeans:
    
    def __init__ (self, microarray, k, verbose=False):
        self.verbose = verbose
        self.mdata = microarray
        self.numGenes = len(microarray)
        self.numExps = len(microarray[0])
        self.k = k
        if self.verbose:
            print "K-means Clustering:\n\tk = %d\n\tnum genes: %d\n\tnum expressions: %d" % (self.k, self.numGenes, self.numExps)
        
    def center (self, cluster):
        listLen = float(len(cluster))
        return [sum(exp)/listLen for exp in zip(*cluster)]

    def sqDist (self, initial, dest):
        sumSquares = 0
        for i in range(len(initial)):
            sumSquares += abs(dest[i] - initial[i])**2
        return sumSquares
    
    def newClusters(self, len=None):
        if not len:
            len = self.k
        return [[] for _ in range(len)]
    
    def cluster (self):
        ""
    
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
            
    def kmeans (self):
        clusters = self.newClusters()
        # for each item, randomly assign it to one of the k clusters
        for geneIdx, gene in enumerate(self.mdata):
            clusters[random.randint(0,self.k-1)].append((gene, geneIdx))
        # remove any empty clusters
        clusters = filter(None, clusters)
        counter = 0
        while True:
            counter += 1
            if self.numExps == 2:
                self.plot(clusters, counter)
            newClusters = self.newClusters(len(clusters))
            centers = [self.center(zip(*cluster)[0]) for cluster in clusters]
            if self.verbose:
                print "\ncluster arrangement %d:" % (counter)
                print "Gene%s\t%s\t cluster assignment" % (' '.join(['  ' for _ in range(2,self.numExps)]), '\t'.join(["sq dist to center C%d" % i for i in range(1, len(centers)+1)]))
            for geneIdx, gene in enumerate(self.mdata):
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

if __name__ == "__main__":                
    # preprocessing (input)    
    exp1 = [1, 6, 5, 2, 6, 0, 1, 1, 2, 4, 5, 6, 7]
    exp2 = [4, 2, 3, 5, 4, 5, 3, 1, 2, 4, 5, 6, 7]
    exps = (exp1, exp2)

    # preliminary data
    microarrayData = list(zip(*exps))
    print microarrayData
    sys.exit(0)
    numExps = len(exps)
    k = 8

    kmeans = KMeans(microarrayData, k, verbose=True)
    finalClusters = kmeans.kmeans()

    print "\nFinal set of gene clusters:"
    for clusterIdx, cluster in enumerate(finalClusters):
        print "\tCluster %d: %s" % (clusterIdx+1, ["gene" + str(idx+1) for gene, idx in cluster])
    print ""
    if kmeans.numExps == 2:
        plt.show()
