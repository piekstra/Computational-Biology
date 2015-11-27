import random

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
            newClusters = self.newClusters(len(clusters))
            centers = [self.center(zip(*cluster)[0]) for cluster in clusters]
            if self.verbose:
                print "\ncluster arrangement %d:" % (counter)
                print "Gene%s\t%s\t cluster assignment" % (' '.join(['  ' for _ in range(2,self.numExps)]), '\t'.join(["sq dist to center C%d" % i for i in range(1, len(centers)+1)]))
            for geneIdx, gene in enumerate(microarrayData):
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

                
# preprocessing (input)    
exp1 = [1, 6, 5, 2, 6, 0, 1, 1, 2, 4, 5, 6, 7]
exp2 = [4, 2, 3, 5, 4, 5, 3, 1, 2, 4, 5, 6, 7]
exps = (exp1, exp2, exp1, exp1, exp2, exp1, exp1)

# preliminary data
microarrayData = [gene for gene in zip(*exps)]
numExps = len(exps)
k = 4

kmeans = KMeans(microarrayData, k, verbose=True)
finalClusters = kmeans.kmeans()

print "\nFinal set of gene clusters:"
for clusterIdx, cluster in enumerate(finalClusters):
    print "\tCluster %d: %s" % (clusterIdx+1, ["gene" + str(idx+1) for gene, idx in cluster])
print ""
