#!/usr/bin/env python

"""
A script to work out the stability of clusters from a set of clusterings of data.

NOTE: It expects to be pointed at a folder, containing folders of text files, each with the genes
contained in that cluster (such as output by PythonScripts\GeneEnrichmentAnalysis\auxScripts\gephiClustersToList.py)

It then does the following:
For a chosen clustering (though depending on run time I guess could do for all) :
	Per module/cluster:
		find overlap numbers with all other modules/clusters per other clustering
		Store the maximal overlap number, along with p-value and increment sums for overlapping genes

		Calculate the median overlap and MAD (median absolute deviation) for cluster.
		Calculate the median p-Values and MAD for cluster

##----------------------------------------IMPORTANT----------------------------------------##
Note that precentages are with regards to the clusters from selectClusterFolder.  So you could get 100% overlap 
but this only means that the genes in that cluster are always found together in a single cluster, not that they
are the ONLY genes in the cluster, indeed they could be in a much larger cluster in other clusterings.

##----------------------------------------SIDE NOTE----------------------------------------##
I wanted to work out overlaps of all pairs of genes across clusters.  This is a stupid idea, for istance 
for a cluster with 4000 genes there are 7,998,000 possible combinations... starts to become a mad number

"""
import sys
import os
import re
import glob
import math
import numpy as np
from scipy import special
import scipy.stats as stats
from collections import defaultdict


workFolder = r"/PATH/PGCNA_normB_0.3/results/PGCNA/EPG3/LEIDENALG/BEST" # change path
selectClusterFolder = r"Clust000" # change cluster number
geneListsFolder = r"ClustersLists"

##----------------------------------------Change less often----------------------------------------##
LOG_BASE = 10
MAX_LOG = 100

outFolder = "ClstStability"


#  Relating to outputing overlapping genes
# CLUST_RE = re.compile("(clust\d+)", re.IGNORECASE)
CLUST_RE = re.compile("^([^\.]+)", re.IGNORECASE)

########################################################################################
outFolder = os.path.join(workFolder, outFolder, selectClusterFolder)

if not os.path.exists(outFolder):
    os.makedirs(outFolder)
###########################################################################################
##----------------------------------------Methods----------------------------------------##
###########################################################################################


def runMain(finishT="Finished!"):
    # global # variable.  Declare each global variable and pass new values/defaults
    pass
    main(finishT=finishT)


def percentile(N, percent, key=lambda x: x):
    """
    Find the percentile of a list of values.

    @parameter N - is a list of values. Note N MUST BE already sorted.
    @parameter percent - a float value from 0.0 to 1.0.
    @parameter key - optional key function to compute value from each element of N.

    @return - the percentile of the values
    """
    N = sorted(N)

    if not N:
        return None
    k = (len(N) - 1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (k - f)
    d1 = key(N[int(c)]) * (c - k)
    return d0 + d1


def mad(arr):
    """Median Absolute Deviation: a "Robust" version of standard deviation.
    Indices variabililty of the sample.
    https://en.wikipedia.org/wiki/Median_absolute_deviation
    """
    # arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    arr = np.array(arr)  # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))


def make_key_naturalSort():
    """
    A factory function: creates a key function to use in sort.
    Sort data naturally
    """

    def nSort(s):
        def convert(text): return int(text) if text.isdigit() else text
        def alphanum_key(key): return [convert(c)
                                       for c in re.split("([0-9]+)", key)]

        return alphanum_key(s)

    return nSort


def hyper(totalPopSize, popOneSize, drawSize, popOneInDraw, maxZ=37.67):
    """
    Calculates probabilities and exact Z-score from a hypergeometric distribution
    see: http://en.wikipedia.org/wiki/Hypergeometric_distribution

    User has to provide 4 numbers:

    totalPopSize : size of popOne + popTwo
    popOneSize: size of just popOne
    drawSize: number of entities drawn without replacement
    popOneInDraw: number of popOne entities in draw
    """

    def zScoreFromPval(pValue):
        """
        For an p-value returns a corresponding Zscore
        """
        z = -stats.norm.ppf(pValue / 2)
        if math.isinf(z):
            return maxZ
        else:
            return z

    phyperSeen = popOneInDraw - 1

    pOne = stats.hypergeom.sf(phyperSeen, totalPopSize, popOneSize, drawSize)
    pTwo = stats.hypergeom.cdf(
        popOneInDraw, totalPopSize, popOneSize, drawSize)

    finalPVal = pOne

    #  Calculate Mean
    mean = (drawSize * popOneSize) / float(totalPopSize)
    variance = (
        drawSize
        * ((totalPopSize - drawSize) / float(totalPopSize - 1))
        * (popOneSize / float(totalPopSize))
        * (1 - (popOneSize / float(totalPopSize)))
    )

    standardDev = math.sqrt(variance)

    if pTwo < pOne:
        finalPVal = pTwo
        zScore = -zScoreFromPval(finalPVal)
    else:
        zScore = zScoreFromPval(finalPVal)

    return [finalPVal, mean, standardDev, zScore]


def getSignificance(numGenes, numOne, numTwo, overlapNum):

    exactVals = hyper(numGenes, numOne, numTwo, overlapNum)
    exactMean, exactSD, exactZscore = exactVals[1:]
    phyperSeen = overlapNum - 1

    pOne = stats.hypergeom.sf(phyperSeen, numGenes, numOne, numTwo)
    pTwo = stats.hypergeom.cdf(phyperSeen, numGenes, numOne, numTwo)
    finalPVal = pOne
    if pTwo < pOne:
        finalPVal = pTwo

    return exactMean, exactSD, exactZscore, finalPVal


def overlapSig(f1, clust1, f2, gPC, natSort, information, comparePos, geneSums, clustRE="Clust(\d+)"):

    totalGeneDict = {}

    #  First pass find all unique genes (background)

    for clust in gPC[f1]:
        for gene in gPC[f1][clust]:
            totalGeneDict[gene] = 1
    for clust in gPC[f2]:
        for gene in gPC[f2][clust]:
            totalGeneDict[gene] = 1

    #  Create numpy arrays to store information
    yNum = len(gPC[f2].keys())
    overlaps = np.zeros((yNum), dtype=np.int32)
    clusterSizes2 = np.zeros((yNum), dtype=np.int32)
    overlapGenes = defaultdict(list)
    clusterOrder = []

    clusterSizes1 = len(gPC[f1][clust1])

    for j, clust2 in enumerate(sorted(gPC[f2], key=natSort)):
        clusterSizes2[j] = len(gPC[f2][clust2])
        clusterOrder.append(clust2)
        for gene in gPC[f1][clust1]:
            if gene in gPC[f2][clust2]:
                overlaps[j] += 1
                overlapGenes[clust1 + ":" + clust2].append(gene)

    #  Get info on maximal overlapping cluster
    maxOverlapIndex = np.argmax(overlaps)
    maxOverlap = max(overlaps)
    maxCluster = clusterOrder[maxOverlapIndex]

    # From perspective of selectClusterFolder cluter
    overlapPercent = 100 * (maxOverlap / float(clusterSizes1))

    #  Calculate significances, for maximal overlap
    mean, sD, zscore, pVal = getSignificance(
        len(
            totalGeneDict), clusterSizes1, clusterSizes2[maxOverlapIndex], maxOverlap
    )

    #  Store information
    information[comparePos] = [maxOverlap,
                               overlapPercent, pVal, zscore, clusterSizes1]

    for g in overlapGenes[clust1 + ":" + maxCluster]:
        geneSums[g] += 1

    ########################################################################################


def processGeneSets(wF, gLF, outF, scF, clustRE=None):

    natSort = make_key_naturalSort()

    def keyOne():
        return defaultdict(dict)

    genesPerCluster = defaultdict(keyOne)

    print("\nStoring gene sets:")
    for foldP in sorted(glob.glob(os.path.join(wF, gLF, "*")), key=natSort):
        foldN = os.path.basename(foldP)
        if not os.path.isdir(foldP):
            continue
        print("\t...", foldN)
        for path in sorted(glob.glob(os.path.join(foldP, "*.txt")), key=natSort):
            clustN = os.path.basename(path)
            # print("\t\t...",clustN)

            for line in open(path):
                line = line.rstrip().upper()
                genesPerCluster[foldN][clustN][line] = 1

    print("\nWorking out stability for cluster")
    #  Create numpy array to store information
    compareNumber = len(genesPerCluster) - 1  # Don't include comparison folder
    perClusterInfo = {}

    for clust in sorted(genesPerCluster[scF], key=natSort):
        information = np.zeros((compareNumber, 5))
        comparePos = 0
        geneSums = defaultdict(int)

        print("\t...", clust)

        for fold1 in sorted(genesPerCluster, key=natSort):
            #  Don't compare with self..
            if fold1 == scF:
                continue

            print("\t...\t", scF, "-vs-", fold1)
            overlapSig(scF, clust, fold1, genesPerCluster, natSort,
                       information, comparePos, geneSums, clustRE=clustRE)
            comparePos += 1

        #  Store info
        perClusterInfo[clust] = [information, geneSums]

    print("\n\nWorking out final numbers:")
    outFile = open(os.path.join(outFolder, "#clusterInfo.txt"), "w")  #
    # print("Cluster\tClusterSize\tMedian_Overlap\tMAD_Overlap\tMedian_Percent\tMAD_Percent\tMedian_Pval\tMAD_Pval\tMedian_Zscore\tMAD_ZScore",file=outFile)
    print(
        "Module\tModuleSize\tMedian_Overlap\tMAD_Overlap\tModuleStability\tMAD_Percent\tMedian_Pval\tMAD_Pval\tMedian_Zscore\tMAD_ZScore",
        file=outFile,
    )
    #  Calculate averages across clusterings
    for clust in sorted(perClusterInfo, key=natSort):
        print("\t...", clust)
        information, geneSums = perClusterInfo[clust]

        medianOverlap = percentile(information.T[0], 0.5)
        madOverlap = mad(information.T[0])
        medianPercent = percentile(information.T[1], 0.5)
        madPercent = mad(information.T[1])
        medianPval = percentile(information.T[2], 0.5)
        madPval = mad(information.T[2])
        medianZscore = percentile(information.T[3], 0.5)
        madZscore = mad(information.T[3])

        #  Finally print out
        print(
            clust,
            information.T[4][0],
            medianOverlap,
            madOverlap,
            medianPercent,
            madPercent,
            medianPval,
            madPval,
            medianZscore,
            madZscore,
            sep="\t",
            file=outFile,
        )

        #  Create file for gene overlaps
        geneFile = open(os.path.join(outFolder, clust + "_geneSums.txt"), "w")
        print("Gene\tOverlapSum\tPercentage", file=geneFile)

        for gene in sorted(genesPerCluster[scF][clust], key=natSort):
            if gene in geneSums:
                percentage = 100 * (geneSums[gene] / float(compareNumber))
                print(gene, geneSums[gene], percentage,
                      sep="\t", file=geneFile)
            else:
                print(gene, "0", "0", sep="\t", file=geneFile)

    outFile.close()


########################################################################################
##----------------------------------------MAIN----------------------------------------##
########################################################################################


def main(finishT="Finished!"):

    processGeneSets(workFolder, geneListsFolder, outFolder,
                    selectClusterFolder, clustRE=CLUST_RE)

    print("\n", finishT, sep="")


##----------------------------------------------------------------------------------------##
##----------------------------------------__main__----------------------------------------##
if __name__ == "__main__":
    main()