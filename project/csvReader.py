class CSVReader:
    def read(self, fileName, verbose=False):
        microArrayData = []
        with open(fileName, 'r') as csv:
            for lineNo, line in enumerate(csv):
                if lineNo == 0:
                    if verbose:
                        print (filter(None, line.rstrip().split(',')))
                else:
                    gene = line.rstrip().split(',')[1:]
                    if verbose: 
                        print (gene)
                    microArrayData.append(list(map(int, gene)))
        return microArrayData

if __name__ == "__main__":         
    csv = CSVReader()
    microArrayData = csv.read("microarraydata.csv", verbose=True)