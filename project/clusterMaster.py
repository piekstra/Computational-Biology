# graphing imports!
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from csvReader import CSVReader
# k-means clustering
from k_means import KMeans
# QT clustering
from QT import QT
# Hierarchical clustering
from hierarchical import Hierarchical

# Print lots of stuff
VERBOSE = False

# the input file
inputFile = "microarraydata.csv"
# A CSV file reader
csvReader = CSVReader()

# get the microarray data from the csv file
microarrayData = csvReader.read(inputFile)

# print out the data in the microarray if in VERBOSE mode
if VERBOSE:
    print ("\nParsed Microarray Data:\n")
    for idx, gene in enumerate(microarrayData):
        print ("Gene %d:" % (idx+1), gene)

## k-means algorithm!
# set the k-value (max potential clusters)
k = 6

# holds a reference to a KMeans object 
kmeans = KMeans(verbose=VERBOSE)
# get the clusters determined by the algorithm
kMeansFinalClusters = kmeans.kmeans(microarrayData, k)

# print out the clusters
print ("\n\nFinal sets of gene clusters:\n")
print ("Using K-Means Algorithm:")
for clusterIdx, cluster in enumerate(kMeansFinalClusters):
    print ("\tCluster %d: {%s}" % (clusterIdx+1, ', '.join([str(idx+1) for gene, idx in cluster])))
print ("")

## QT algorithm!
# set the diameter for the QT algorithm
diameter = 6

# holds a reference to a QT object 
qt = QT(verbose=VERBOSE)
# get the clusters determined by the algorithm
qtFinalClusters = qt.QTClustering(microarrayData, diameter)

# print out the clusters
print ("Using QT Algorithm:")
for clusterIdx, cluster in enumerate(qtFinalClusters):
    print ("\tCluster %d: {%s}" % (clusterIdx+1, ', '.join([str(geneNum+1) for geneNum in cluster])))
print ("")

## Hierarchical Clustering algorithm!

# holds a reference to a QT object 
hc = Hierarchical(microarrayData, VERBOSE)
# get the clusters determined by the algorithm
hcFinalClusters = hc.hierarchicalCluster()

# print out the clusters
print ("Using Hierarchical Clustering Algorithm:")
print ("Tree: %s" % hcFinalClusters)
