# graphing imports!
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import re

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
inputFile = "ALL-AML-TRANSPOSED.csv"
# the output file
outputFile = "results.txt"
# A CSV file reader
csvReader = CSVReader()

# get the microarray data from the csv file
microarrayData = csvReader.read(inputFile)
microarrayLabels = csvReader.getLabels(inputFile)
print ("File %s parsed succesfully!\n\tRows:\t\t%d\n\tColumns:\t%d" % (inputFile, len(microarrayData), len(microarrayData[0])))
print ("\nLabels: {%s}" % (', '.join(microarrayLabels)))

## k-means algorithm!
# set the k-value (max potential clusters)
k = 3
# holds a reference to a KMeans object 
kmeans = KMeans(verbose=VERBOSE)
# get the clusters determined by the algorithm
kMeansFinalClusters = kmeans.kmeans(microarrayData, k)

## QT algorithm!
# set the diameter for the QT algorithm
diameter = 100000# holds a reference to a QT object 
qt = QT(verbose=VERBOSE)
# get the clusters determined by the algorithm
qtFinalClusters = qt.QTClustering(microarrayData, diameter)

## Hierarchical Clustering algorithm!
# holds a reference to a QT object 
hc = Hierarchical(microarrayData, VERBOSE)
# get the clusters determined by the algorithm
hcFinalClusters = hc.hierarchicalCluster()

# print out the data in the microarray if in VERBOSE mode
if VERBOSE:
    print ("\nParsed Microarray Data:\n")
    for idx, gene in enumerate(microarrayData):
        print ("Gene %d:" % (idx+1), gene)

# print out the k-means clusters
print ("\n\nFinal sets of gene clusters:\n")
print ("Using K-Means Algorithm:")
for clusterIdx, cluster in enumerate(kMeansFinalClusters):
    print ("\tCluster %d: {%s}" % (clusterIdx+1, ', '.join([microarrayLabels[idx] for gene, idx in cluster])))
print ("")

# print out the QT clusters
print ("Using QT Algorithm:")
for clusterIdx, cluster in enumerate(qtFinalClusters):
    print ("\tCluster %d: {%s}" % (clusterIdx+1, ', '.join([microarrayLabels[geneNum] for geneNum in cluster])))
print ("")

# print out the hierarchical tree
print ("Using Hierarchical Clustering Algorithm:")
# establish a regular expression to match the genes
geneRegex = re.compile('G\d+')
# get the gene indexes
geneNumbers = [int(gene[1:]) - 1 for gene in geneRegex.findall(hcFinalClusters)]
# use the gene indexes to get the labels for the genes
# and replace the genes with their labels
labeledGenes = re.sub(geneRegex, lambda match: microarrayLabels[geneNumbers.pop(0)], hcFinalClusters)
# print out the labeled 'tree'
print ("Tree: %s" % (labeledGenes))

with open(outputFile, 'w') as out:
    # print out the k-means clusters
    out.write("\n\nProgram results for gene clusters:\n\n")
    out.write("\nUsing K-Means Algorithm with k = %d:\n" % (k))
    for clusterIdx, cluster in enumerate(kMeansFinalClusters):
        out.write("\tCluster %d: {%s}\n" % (clusterIdx+1, ', '.join([microarrayLabels[idx] for gene, idx in cluster])))
    out.write("")

    # print out the QT clusters
    out.write("\nUsing QT Algorithm with diameter %d:\n" % (diameter))
    for clusterIdx, cluster in enumerate(qtFinalClusters):
        out.write("\tCluster %d: {%s}\n" % (clusterIdx+1, ', '.join([microarrayLabels[geneNum] for geneNum in cluster])))
    out.write("")

    # print out the hierarchical tree
    out.write("\nUsing Hierarchical Clustering Algorithm:\n")
    out.write("Tree: %s\n" % hcFinalClusters)
    
