## 
#
# Author: Caleb Piekstra
#
# Last Modified: 12/11/2015
#
# Description: Parases a CSV file that contains microarray data
#   in the form of genes on rows and expression values on columns
#
class CSVReader:
    
    
    ##    
    # read
    #
    # Description: Given a csv file name, opens the file
    #   and parses out the microarray data
    #
    # Parameters:
    #   fileName - The name of the file to read
    #   verbose - Optional parameter to print out 
    #       extra information
    #
    # Returns:
    #   The microarray data that was parsed from the file
    #
    def read(self, fileName, verbose=False):
        # create an array to hold the microarray data
        microArrayData = []
        
        # open the file for reading
        with open(fileName, 'r') as csv:
            # go through each line
            for lineNo, line in enumerate(csv):
                # skip the experiment column headers
                if lineNo == 0:
                    if verbose:
                        print (filter(None, line.rstrip().split(',')))
                else:
                    # gather the gene data skipping the first column (gene label)
                    gene = line.rstrip().split(',')[1:]
                    if verbose: 
                        print (gene)
                    # add the gene to the microarray data
                    microArrayData.append(list(map(int, gene)))
                    
        # return the microarray data in the form of a 2D list
        # where each inner list is a gene's expression data
        return microArrayData
    
    
    ##    
    # getLabels
    #
    # Description: Given a csv file name, opens the file
    #   and parses out the labels
    #
    # Parameters:
    #   fileName - The name of the file to read
    #   verbose - Optional parameter to print out 
    #       extra information
    #
    # Returns:
    #   The labels
    #
    def getLabels(self, fileName, verbose=False):
        # create an array to hold the labels
        labels = []
        
        # open the file for reading
        with open(fileName, 'r') as csv:
            # return the labels (column 0)
            labels = [line.rstrip().split(',')[0] for line in csv][1:]
                    
        # return the labels
        return labels


# Test code for directly running the file        
if __name__ == "__main__":         
    # initialize the reader
    csv = CSVReader()
    # read the data and printo out the results
    microArrayData = csv.read("microarraydata.csv", verbose=True)
