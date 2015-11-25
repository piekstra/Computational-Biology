import random

class KMeans:
    def center (self, cluster):
        listLen = float(len(cluster))
        return [sum(exp)/listLen for exp in zip(*cluster)]

    def sqDist (self, initial, dest):
        sumSquares = 0
        for i in range(len(initial)):
            sumSquares += abs(dest[i] - initial[i])**2
        return sumSquares
    
kmeans = KMeans()

# preprocessing (input)    
exp1 = [1, 6, 5, 2, 6, 0, 1]
exp2 = [4, 2, 3, 5, 4, 5, 3]
exps = (exp1, exp2)

# preliminary data
microarrayData = [merged for merged in zip(*exps)]
numExps = len(exps)

k = 4
clusters = [[] for _ in range(k)]
# for each item, randomly assign it to one of the k clusters
for geneIdx, gene in enumerate(microarrayData):
    clusters[random.randint(0,k-1)].append((gene, geneIdx))
clusters = [cluster for cluster in clusters if not(not cluster)]


## main    
print "K-means Clustering:\n\tk = %d\n\tnum genes: %d\n\tnum expressions: %d" % (k, len(microarrayData), numExps)
# Unit test!
# c1Idxs = [0, 1, 4, 5]
# c2Idxs = [2,3,6]
# c1 = [(microarrayData[i], i) for i in c1Idxs]
# c2 = [(microarrayData[i], i) for i in c2Idxs]
# clusters = [c1, c2]

counter = 0
newClusterAssignments = True
while newClusterAssignments:
    counter += 1
    newClusters = [[] for _ in range(k)]
    centers = [kmeans.center(zip(*cluster)[0]) for cluster in clusters]
    print "\ncluster arrangement %d:" % (counter)
    print "Gene\t%s\t cluster assignment" % ('\t'.join(["sq dist to center C%d" % i for i in range(1, len(centers)+1)]))
    for geneIdx, gene in enumerate(microarrayData):
        distancesToCenters = [kmeans.sqDist(gene, center) for center in centers]
        assignedClusterIdx = min(xrange(len(distancesToCenters)), key=distancesToCenters.__getitem__)
        newClusters[assignedClusterIdx].append((gene, geneIdx))
        print "%s\t%s\t\t\t %s" % (gene, '\t\t\t'.join(["%0.2f" % dist for dist in distancesToCenters]), ("C%d" % (assignedClusterIdx+1)))
    # remove empty clusters
    newClusters = [cluster for cluster in newClusters if cluster]
    if newClusters == clusters:
        newClusterAssignments = False
    else:
        clusters = newClusters

print "\nFinal set of gene clusters:"
for clusterIdx, cluster in enumerate(clusters):
    print "\tCluster %d: %s" % (clusterIdx+1, ["gene" + str(idx+1) for gene, idx in cluster])
print ""
