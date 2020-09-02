
[a, b] = size(GHResultsProcessing)

GHResultsProcessing(b+1).deviceName = fileName;
GHResultsProcessing(b+1).Jsc = [simulationResults.Jsc]
GHResultsProcessing(b+1).PCE = [simulationResults.PCE]
GHResultsProcessing(b+1).G = [simulationResults.growthfactor]