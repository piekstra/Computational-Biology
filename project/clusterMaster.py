# graphing imports!
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from csvReader import CSVReader
# k-means clustering
from k_means import KMeans
# QT clustering
from QT import QT

# Print lots of stuff
VERBOSE = True

inputFile = "microarraydata.csv"
csvReader = CSVReader()

microarrayData = csvReader.read(inputFile)
if VERBOSE:
    print ("\nParsed Microarray Data:\n")
    for idx, gene in enumerate(microarrayData):
        print ("Gene %d:" % (idx+1), gene)

# k-means algorithm!
k = 6
kmeans = KMeans(verbose=False)
kMeansFinalClusters = kmeans.kmeans(microarrayData, k)
print ("\n\nFinal sets of gene clusters:\n")
print ("Using K-Means Algorithm:")
for clusterIdx, cluster in enumerate(kMeansFinalClusters):
    print ("\tCluster %d: {%s}" % (clusterIdx+1, ', '.join([str(idx+1) for gene, idx in cluster])))
print ("")

# QT algorithm!
diameter = 6
qt = QT(verbose=False)
qtFinalClusters = qt.QTClustering(microarrayData, diameter)
print ("Using QT Algorithm:")
for clusterIdx, cluster in enumerate(qtFinalClusters):
    print ("\tCluster %d: {%s}" % (clusterIdx+1, ', '.join([str(geneNum+1) for geneNum in cluster])))
print ("")
