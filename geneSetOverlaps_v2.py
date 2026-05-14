#!/usr/bin/env python
"""
For a folder containing subfolders of gene lists, compare
each list in one folder against each list in the other folder and
generate a heatmap of significance.

This is much like the WGCNA Fisher's exact method.  But hopefully nicely
automated...

Written by Matthew Care who was bored of using WGCNA and that horrible R thingy.

NOTE: The version for Gene Signatures is geneSignatureEnrichmentOverlaps.py
NOTE2: When you have used both geneSignatureEnrichmentOverlaps.py and geneSetOverlaps.py
        then use geneSetOverlapsVsgseOverlaps.py to compare/contrast

##----------------------------------------IMPORTANT----------------------------------------##
NOTE: this calculates correct p-values for the overlaps seen **EVEN** when the two data-sets contain different genes

The p-values are calculated using the stats.hypergeom.cdf :

stats.hypergeom.cdf(k,M,n,N), where k = overlap, M = total genes, n = total type 1, N = draw size (total type 2)

Thus for example if two data-sets have 11,224 non-redundant genes and we want to calculate the overlap of two modules:

Module 1: 74 genes
Module 2: 602 genes
Overlap : 9 genes (for the calculation we want to find prob of >= observed, but cdf/phyper(lower.tail=FALSE) calculate > observed, thus need to use **observed - 1**)

=== stats.hypergeom.cdf(8,11224,74,602) = 0.01701453

Equivalent in R is the phyper function phyper(q, m, n, k), where q = overlap, m = total type 1, n = total black (totalGenes - total type 1), k = draw size (total type 2)

phyper(8,74,11150,602,lower.tail = FALSE) == phyper(8,602,10622,74,lower.tail=FALSE) i.e. it doesn't matter which is used as the draw size.

So what we're asking is, given that we draw 602 genes (module 2), and 9 of those overlap with module 1 (n=74) what's the probability given that there are 11,224 genes in total (or 11150 that aren't module 1)

NOTE2: Can use General\\drawHistogramFrom2DMatrix.py to draw histograms of all -log10 p-values to visualise overlap significance range.

# MERGING PDFS
Merge pdfs in linux with the two commands
pdfunite *.pdf all.pdf
pdfxup -x 6 -y 5 -m 0 -o Wrapped.pdf all.pdf
"""

import glob
import math
import os
import re
import sys
from collections import defaultdict

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from scipy import special

workFolder = r"PATH/Data/" # change path

geneListsFolder = r"geneListsNamed"
outFolder = "#OL_000_001_v2_labelled" # change numbers to those of selected clusterings

OVERWRITE_EXISTING = True  # Useful to just repopulate the changed/deleted files --> set to False

SPLIT_GENE = True  # If True will use first gene only from set.
SPLIT_GENE_BY = "#"

########################################################################################
LOG_BASE = 10
MAX_LOG = 50
LINE_WIDTH = 0.1  # 0.3 Width of line to draw grid
PLOT_WIDTH = 25  # 25
PLOT_HEIGHT = 25  # 25

ALTER_AXIS_LABEL_SIZES = True  # Useful for putting in presentations; along with PLOT_WIDTH=20 #baseline is false
AXIS_LABEL_SIZE = 50 # 100 matt baseline. This just does x and y, not labelling the module names. 

PLOT_COLOUR_MAP = "jet"  # e.g jet,rainbow, bwr, Reds etc. From http://matplotlib.org/examples/color/colormaps_reference.html
# If True will add in each box the number of overlapping genes.
PLOT_OVERLAP_NUMBERS = True
OUTPUT_GENE_LISTS = False

#  fontdict = {'fontsize': rcParams['axes.titlesize'],
#  'fontweight': rcParams['axes.titleweight'],
#  'verticalalignment': 'baseline',
#  'horizontalalignment': loc}
#  CHANGE LESS OFTEN
########################################################################################
#  Make values below VAL_MIN appear as white
VAL_MIN = 2  # 1.3; Log value 0.05 would be 1.301029996 for LOG_BASE 10 and 0.01 would be 2
myCmap = plt.get_cmap(PLOT_COLOUR_MAP).copy()
myCmap.set_under("white")
#  Make text appear correctly
#  A number between 0 --4 , to set threshold of when text appears white given a colour background.  > more likely to be white
TEXT_COLOUR_THRESHOLD = 2.2  # Default 2.2 for jet cmap"

#  Relating to outputing overlapping genes
# CLUST_RE = re.compile(r"(clust\d+)", re.IGNORECASE)
CLUST_RE = re.compile(r"(M\d+).*", re.IGNORECASE)
CLUST_RE = re.compile(r"(.*)", re.IGNORECASE)
# CLUST_RE = re.compile(r"(n\d+).*", re.IGNORECASE)
# CLUST_RE = re.compile("^([^\.]+)", re.IGNORECASE)
OUT_GENE_PREFIX = "OLgenes"
OUT_GENE_PREFIX = "OLgenes"


########################################################################################
if SPLIT_GENE:
    outFolder = outFolder + "_splitG"
outFolder = outFolder + "_ValMin" + str(VAL_MIN)

outFolder = os.path.join(workFolder, outFolder)
if not os.path.exists(outFolder):
    os.makedirs(outFolder)

###########################################################################################
##----------------------------------------Methods----------------------------------------##
###########################################################################################


def runMain(finishT="Finished!"):
    # global # variable.  Declare each global variable and pass new values/defaults
    pass
    main(finishT=finishT)


def make_key_naturalSort():
    """
    A factory function: creates a key function to use in sort.
    Sort data naturally
    """

    def nSort(s):
        def convert(text):
            return int(text) if text.isdigit() else text

        def alphanum_key(key):
            return [convert(c) for c in re.split("([0-9]+)", key)]

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
    pTwo = stats.hypergeom.cdf(popOneInDraw, totalPopSize, popOneSize, drawSize)

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
    pTwo = stats.hypergeom.cdf(overlapNum, numGenes, numOne, numTwo)
    finalPVal = pOne
    if pTwo < pOne:
        finalPVal = pTwo

    return exactMean, exactSD, exactZscore, finalPVal


def overlapSig(f1, f2, gPC, natSort, outF, outFGenes="OLgenes", clustRE=r"Clust(\d+)"):
    totalGeneDict = {}

    global VAL_MIN

    #  First pass find all unique genes (background)

    for clust in gPC[f1]:
        for gene in gPC[f1][clust]:
            totalGeneDict[gene] = 1
    for clust in gPC[f2]:
        for gene in gPC[f2][clust]:
            totalGeneDict[gene] = 1

    print("\t\t...Total unique genes:", len(totalGeneDict))

    #  Print total gene list
    totalGenesF = open(os.path.join(outF, f1 + "_vs_" + f2 + "_Genes.txt"), "w")
    for g in sorted(totalGeneDict, key=natSort):
        print(g, file=totalGenesF)
    totalGenesF.close()

    #  Create numpy arrays to store information
    xHeader = sorted(gPC[f1].keys(), key=natSort)
    yHeader = sorted(gPC[f2].keys(), key=natSort)
    xNum = len(gPC[f1].keys())
    yNum = len(gPC[f2].keys())

    print("\t\t...Clusters:", xNum, "by", yNum)

    logPVals = np.ones((xNum, yNum))
    enrichedLogPVals = np.ones((xNum, yNum))
    zScores = np.zeros((xNum, yNum))
    overlaps = np.zeros((xNum, yNum), dtype=np.int32)
    clusterSizes1 = np.zeros((xNum), dtype=np.int32)
    clusterSizes2 = np.zeros((yNum), dtype=np.int32)

    overlapGenes = defaultdict(list)

    print("\t...\tCalculate overlaps")
    for i, clust1 in enumerate(sorted(gPC[f1], key=natSort)):
        clusterSizes1[i] = len(gPC[f1][clust1])
        for j, clust2 in enumerate(sorted(gPC[f2], key=natSort)):
            clusterSizes2[j] = len(gPC[f2][clust2])

            for gene in gPC[f1][clust1]:
                if gene in gPC[f2][clust2]:
                    overlaps[i][j] += 1
                    overlapGenes[clust1 + ":" + clust2].append(gene)

            #  Calculate significance
            mean, sD, zscore, pVal = getSignificance(len(totalGeneDict), clusterSizes1[i], clusterSizes2[j], overlaps[i][j])
            try:
                tempLog = -math.log(pVal, LOG_BASE)
                if tempLog > MAX_LOG:
                    tempLog = MAX_LOG
                logPVals[i][j] = tempLog

                if zscore > 0:
                    enrichedLogPVals[i][j] = tempLog

            except ValueError:
                logPVals[i][j] = MAX_LOG
                if zscore > 0:
                    enrichedLogPVals[i][j] = MAX_LOG
            zScores[i][j] = zscore

    #  Make final headers
    for i, fName in enumerate(xHeader):
        xHeader[i] = xHeader[i] + ": " + str(clusterSizes1[i])
    for i, fName in enumerate(yHeader):
        yHeader[i] = yHeader[i] + ": " + str(clusterSizes2[i])

    ##----------------------------------------OUTPUT RESULTS----------------------------------------##
    # Zscores
    zScoreFile = open(os.path.join(outF, f1 + "_vs_" + f2 + "_Zscores.txt"), "w")
    print("File1/File2", "\t".join(yHeader), sep="\t", file=zScoreFile)
    for i in range(zScores.shape[0]):
        print(xHeader[i], "\t".join(map(str, zScores[i])), sep="\t", file=zScoreFile)
    zScoreFile.close()

    # Overlaps
    overlapFile = open(os.path.join(outF, f1 + "_vs_" + f2 + "_Overlaps.txt"), "w")
    print("File1/File2", "\t".join(yHeader), sep="\t", file=overlapFile)
    for i in range(overlaps.shape[0]):
        print(xHeader[i], "\t".join(map(str, overlaps[i])), sep="\t", file=overlapFile)
    overlapFile.close()

    # Enriched log vals
    logPvalFile = open(
        os.path.join(outF, f1 + "_vs_" + f2 + "_Log" + str(LOG_BASE) + "values.txt"),
        "w",
    )
    print("File1/File2", "\t".join(yHeader), sep="\t", file=logPvalFile)
    for i in range(enrichedLogPVals.shape[0]):
        print(
            xHeader[i],
            "\t".join(map(str, enrichedLogPVals[i])),
            sep="\t",
            file=logPvalFile,
        )
    logPvalFile.close()

    if OUTPUT_GENE_LISTS:
        outFold2 = os.path.join(outF, outFGenes + f1 + "_vs_" + f2)
        if not os.path.exists(outFold2):
            os.makedirs(outFold2)
        # Overlapping genes
        for inf in overlapGenes:
            clust1, clust2 = inf.replace(" ", "").split(":")

            try:
                finalClust1 = re.search(clustRE, clust1).groups()[0]
                finalClust2 = re.search(clustRE, clust2).groups()[0]
            except AttributeError:
                print("Somethings wrong with clustRE for this set of clusters:", inf)
                sys.exit()

            outGeneFile = open(os.path.join(outFold2, finalClust1 + "-vs-" + finalClust2 + ".txt"), "w")
            for gene in sorted(overlapGenes[inf], key=natSort):
                print(gene, file=outGeneFile)
            outGeneFile.close()

    ##----------------------------------------MATPLOTLIB PLOTS----------------------------------------##
    plt.figure(figsize=(PLOT_WIDTH, PLOT_HEIGHT))
    ax = plt.gca()

    def show_values(pc2, main, fmt="%.0f", **kw):
        pc2.update_scalarmappable()
        main.update_scalarmappable()
        ax = pc2.axes

        # Get the array data and flatten it
        pc2_array = pc2.get_array()
        if pc2_array.ndim > 1:
            pc2_array = pc2_array.flatten()

        for i, (p, color) in enumerate(zip(pc2.get_paths(), main.get_facecolors())):
            x, y = p.vertices[:-2, :].mean(0)

            #  Sort text colour so we can see it
            if np.all(sum(color) > TEXT_COLOUR_THRESHOLD):
                color = (0.0, 0.0, 0.0)
            else:
                color = (1.0, 1.0, 1.0)

            if PLOT_OVERLAP_NUMBERS:
                # Get the value for this specific cell. Baseline font size ~12, this does the numbers in the boxes. 
                value = pc2_array[i]
                ax.text(x, y, fmt % value, ha="center", va="center", color=color, fontsize=17,**kw)

    # # IF ERRORS ARISE, MAINLY WITH RANDOM DATA then uncomment this section
    # if VAL_MIN > enrichedLogPVals.max():
    # 	print("VAL_MIN (",VAL_MIN,") is > than maximum logPval (",enrichedLogPVals.max(),") Will set it lower!",sep="")
    # 	VAL_MIN = enrichedLogPVals.max() - 0.5

    # 	if VAL_MIN < 0:
    # 		VAL_MIN = 0

    #  Add text
    text = ax.pcolor(
        overlaps,
        edgecolors="k",
        linewidths=LINE_WIDTH,
        vmin=VAL_MIN,
        cmap=myCmap,
        vmax=MAX_LOG,
    )
    main = ax.pcolor(
        enrichedLogPVals,
        edgecolors="k",
        linewidths=LINE_WIDTH,
        vmin=VAL_MIN,
        cmap=myCmap,
        vmax=MAX_LOG,
    )

    show_values(text, main)
    plt.colorbar(main)

    ax.set_xticks(np.arange(enrichedLogPVals.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(enrichedLogPVals.shape[0]) + 0.5, minor=False)

    #  Axis labels
    ax.set_xticklabels(yHeader, minor=False, rotation=90, fontsize=20)
    ax.set_yticklabels(xHeader, minor=False, fontsize=20)

    if ALTER_AXIS_LABEL_SIZES:
        #  Add axis labels
        plt.xlabel("Normal B module", size=AXIS_LABEL_SIZE) # replace label "" with f2 for cluster label
        plt.ylabel("CLL module", size=AXIS_LABEL_SIZE) # replace label "" with f1 for cluster label
    else:
        #  Add axis labels, default size
        plt.xlabel(f2)  # original is f2
        plt.ylabel(f1)  # original is f1

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.set_adjustable("box")

    #  Tidy up plot
    plt.axis("equal")
    plt.axis("tight")
    plt.tight_layout()

    figName = os.path.join(outF, f1 + "_vs_" + f2 + "_Log" + str(LOG_BASE) + "values.png")
    plt.savefig(figName)
    # Clear figure!  Stop memory leak.
    plt.clf()
    plt.cla()
    plt.close()

    ########################################################################################


def processGeneSets(wF, gLF, outF, outFGenes="OLgenes", clustRE=None, splitGene=False, splitGeneBy="#"):
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
            clustN = os.path.basename(path)[:-4]
            print("\t\t...", clustN)

            with open(path) as f:
                for line in f:
                    line = line.rstrip().upper()

                    if splitGene:
                        bits = line.split(splitGeneBy)
                        line = bits[0]
                        if len(bits) > 1:
                            if bits[1].isdigit():
                                line = line + "#" + bits[1]

                    genesPerCluster[foldN][clustN][line] = 1

    print("\nCalculate all overlaps")

    seen = {}
    for fold1 in sorted(genesPerCluster, key=natSort):
        for fold2 in sorted(genesPerCluster, key=natSort):
            if fold1 == fold2:
                continue

            if fold1 + fold2 in seen:
                continue

            print("\t...", fold1, "-vs-", fold2)

            if not OVERWRITE_EXISTING:
                zScoreFile = os.path.join(outF, fold1 + "_vs_" + fold2 + "_Zscores.txt")
                figName = os.path.join(outF, fold1 + "_vs_" + fold2 + "_Log" + str(LOG_BASE) + "values.pdf")
                overlapFile = os.path.join(outF, fold1 + "_vs_" + fold2 + "_Overlaps.txt")
                logPvalFile = os.path.join(outF, fold1 + "_vs_" + fold2 + "_Log" + str(LOG_BASE) + "values.txt")

                if (
                    os.path.exists(zScoreFile)
                    and os.path.exists(figName)
                    and os.path.exists(overlapFile)
                    and os.path.exists(logPvalFile)
                ):
                    print("\t\t\tAll files already exists skipping")
                    continue
                else:
                    overlapSig(fold1, fold2, genesPerCluster, natSort, outF, outFGenes, clustRE)
            else:
                overlapSig(fold1, fold2, genesPerCluster, natSort, outF, outFGenes, clustRE)

            #  Store seen
            seen[fold1 + fold2] = 1
            seen[fold2 + fold1] = 1


########################################################################################
##----------------------------------------MAIN----------------------------------------##
########################################################################################


def main(finishT="Finished!"):
    processGeneSets(
        workFolder,
        geneListsFolder,
        outFolder,
        outFGenes=OUT_GENE_PREFIX,
        clustRE=CLUST_RE,
        splitGene=SPLIT_GENE,
        splitGeneBy=SPLIT_GENE_BY,
    )

    print("\n", finishT, sep="")


##----------------------------------------------------------------------------------------##
##----------------------------------------__main__----------------------------------------##
if __name__ == "__main__":
    main()