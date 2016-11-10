#==================================================================================================
# Sarah Middleton
# Junhyong Kim Lab
# University of Pennsylvania
# Last update: Nov 2016
#==================================================================================================
# step2_classification.py
# Given a .scoremat file (generated by step1_threading.py), outputs fold classifications for each
# sequence using a nearest neighbor (NN) classifier. Sequences with a NN distance of <= 17.5 are
# marked as "high" confidence. All predictions are output to the .fold_preds.txt file. 
# 
# Notes:
# - Requires Numpy and scikit-learn
#
# Usage:
#    python step2_classification.py SCOREMAT
#  
# Usage examples: 
#    python ~/src/step2_classification.py ~/demo/demo.scoremat
#==================================================================================================
import sys, os
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KNeighborsClassifier
from optparse import OptionParser


#------------------------------------------------
# Read/set up variables
#------------------------------------------------

# parse args
usageMsg = "Usage: %prog SCOREMAT [options]"
parser = OptionParser(usage=usageMsg)
parser.add_option("--cpus", action="store", type='int', default=1, dest="MAX_CPU", help="Maximum number of CPUs to use. Default is [%default].")

# read/process args
(opts, args) = parser.parse_args()
if len(args) == 1:
	testDataFile = args[0]
else:
	print(">> Missing input file. Use -h for help. Exiting.")
	sys.exit()
scriptDir = os.path.dirname(sys.argv[0])

# parameters
THRESH = 17.5
K = 1
CPU = opts.MAX_CPU


#------------------------------------------------
# Create classifier
#------------------------------------------------
print("")
print("Reading data...")

# read training data
trainDataFile = "%s/train/scope20_all.scoremat" % scriptDir
train_idList = []
train_data = []
ins = open(trainDataFile, 'r')
train_header = ins.readline().rstrip('\r\n').split()

for line in ins:
    lineParts = line.strip('\r\n').split()
    train_idList.append(lineParts[0])
    train_data.append([float(x) for x in lineParts[1:]])
	
ins.close()

# get fold assignment of each training example
train_foldIDs = []
train_foldCodes = []
foldToCode = {}
codeToFold = {}
count = 0

for seqID in train_idList:
    classif, protID = seqID.split("_", 1)
    classifParts = classif.split(".")
    
    foldID = "%s.%s" % (classifParts[0], classifParts[1])
    if foldID not in foldToCode:
        count += 1
        foldToCode[foldID] = count
        codeToFold[count] = foldID
    
    train_foldIDs.append(foldID)
    train_foldCodes.append(foldToCode[foldID])

# z-scale data
zScaler = StandardScaler()
train_data_scaled = zScaler.fit_transform(train_data)

# create KNN classifier
clf_knn = KNeighborsClassifier(K, n_jobs=CPU)
clf_knn.fit(train_data_scaled, train_foldCodes)


#------------------------------------------------
# Read in query data
#------------------------------------------------
test_idList = []
test_data = []
ins = open(testDataFile, 'r')
test_header = ins.readline().rstrip('\r\n').split()
for line in ins:
    lineParts = line.strip('\r\n').split()
    test_idList.append(lineParts[0])
    test_data.append([float(x) for x in lineParts[1:]])
ins.close()

# z-scale data using same parameters as training set
test_data_scaled = zScaler.transform(test_data)

# classify using 1NN
print("Classifying...")
predicted = clf_knn.predict(test_data_scaled)
distToClosest, _ = clf_knn.kneighbors(test_data_scaled, n_neighbors=1)
print("Finished.")
print("")

# separate into "classified" and "not classified"
distToClosest2 = np.array([x[0] for x in distToClosest])
belowThresh = distToClosest2 < THRESH
print("%s of %s (%.2f%%) classified with high confidence (NN Dist <= 17.5)" % (sum(belowThresh), len(belowThresh), (float(sum(belowThresh)) / len(belowThresh) * 100)))
print("")

# print results
outfile_preds = "%s.fold_preds.txt" % testDataFile

outs = open(outfile_preds, 'w')
outs.write("SeqID\tPredFold\tNNDistance\tConfidence\n")

for i in range(len(test_idList)):
    seqId = test_idList[i]
    predFold = codeToFold[predicted[i]]
    nndist = distToClosest2[i]
    conf = "Low"
    if nndist <= THRESH:
        conf = "High"
        
    outStr = "%s\t%s\t%s\t%s\n" % (seqId, predFold, nndist, conf)
    outs.write(outStr)

outs.close()

print("Fold predictions printed to %s" % outfile_preds)
print("")
