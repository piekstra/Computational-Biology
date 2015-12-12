System Requirements:
	Python 3.X
	matplotlib python package
	possibly numpy package

To install missing packages run the following:
	pip install matplotlib
	pip install numpy
Note: pip is python's package manager 

To run the program, simply run clusterMaster.py using Python 3.X

Example:
	$ python3 clusterMaster.py

To edit the k value for k-means clustering and the diameter for QT clustering
you have to edit the 'k' and 'diameter' variable in clusterMaster.py

The clustering results as well as a tree of the genes is printed out to the
console and written to the file 'results.txt'

To modify the microarray data, edit the 'microarraydata.csv' file.
Note that the program expects that each column will be labeled with
an expression header and each row's first column will contain the 
gene name. The values contained in the rest of the file are comma-
seperated integers.


