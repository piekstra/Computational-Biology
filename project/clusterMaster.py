# graphing imports!
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# clustering
from csvReader import CSVReader
from k_means import KMeans

inputFile = "microarraydata.csv"
k = 8


csvReader = CSVReader()
microarrayData = csvReader.read(inputFile)
print microarrayData

kmeans = KMeans(microarrayData, k, verbose=True)
finalClusters = kmeans.kmeans()

print "\nFinal set of gene clusters:"
for clusterIdx, cluster in enumerate(finalClusters):
    print "\tCluster %d: %s" % (clusterIdx+1, ["gene" + str(idx+1) for gene, idx in cluster])
print ""
if kmeans.numExps == 2:
    plt.show()