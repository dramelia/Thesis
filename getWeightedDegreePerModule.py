#!/usr/bin/env python
"""
Uses igraph to find the degree/weighted-degree per gene within each module using module subgraphs

Input:
Edged file in gzipped tab separated format, as output by PGCNA2
Gene module list in csv formate, as output by PGCNA2

Option to include information from a single correlation file, this allows the addition of
MedianPercentile, Median_QCODexpression(VarWithin) & MADpercentile(VarAcross)

to generate the following:
ModCon = degree^2 * VarWithin  * (100 - VarAcross)/100

Note: Expression is NOT used in ModCon calculation!

Outputs can then be processed by placed in a folder and processed by:
PGCNA-Private\AuxScripts\gatherWeightedDegreePerModuleAcrossNetworksAutoRank.py
"""

import gzip
import os
import re
import sys
from collections import defaultdict

import igraph as ig

#################################################################
# -----------------------Youli D13 MYC MB-Domain  CPM-----------------#
workFolder = r"PATH/PGCNA_work/PostPGCNA_normB/MEV_Vis" # change path
edgesFP = r"PATH/PGCNA_work/PostPGCNA_normB/results/PGCNA/EPG3/GEPHI/ForPGCNA_RetainF0.3_EPG3_Edges.tsv.gz" # change path
geneModuleFile = r"PATH/PGCNA_work/PostPGCNA_normB/results/PGCNA/Selected/608.csv" # change path
outFileSuffix = "_normB_0.3_GEO"

# ----------------------------OPTIONAL----------------------------#
singleCorrFile = r"path/SINGLE_CORR/I/IRF4#ENSG00000137265.16_corr_RetainF0.3.tsv.gz" # change path
#################################################################


#################################################################
# If True will output the degree weighted by edge value (correlation), else will be integer degree
WEIGHTED_DEGREE = True
# If True will also include the information from a single corr file (generate using PGCNA2 --singleCorr/--singleCorrL) and rank by ModCon
USE_SINGLE_CORR_INFO = True

ROUND_VALUES_TO = 3  # Number of decimal places to round values to

# Current pgcna2 files
SCF_GENE_COL = 0  # 0
SCF_VAR_WITHIN_COL = 4  # 4
SCF_MEDIAN_EXP_COL = 5  # 5
SCF_VAR_ACROSS_COL = 6  # 6

#################################################################

if USE_SINGLE_CORR_INFO:
    outFileSuffix = outFileSuffix + "_ModCon"

if WEIGHTED_DEGREE:
    OUT_FILE = "geneModuleWeightedDegree" + outFileSuffix + ".txt"
else:
    OUT_FILE = "geneModuleDegree" + outFileSuffix + ".txt"
#################################################################


# ----------------------------Methods----------------------------#


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


def convertGephiToIgraph(gephiEdgeFile, splitBy="\t", header=1, outN="igraphG.temp"):
    """
    Convert edge file into format for import to igraph
    """
    baseFold, fileName = os.path.split(gephiEdgeFile)
    outFileP = os.path.join(baseFold, outN)
    outFile = open(outFileP, "w")

    for line in gzip.open(gephiEdgeFile):
        line = line.decode()
        if header:
            header -= 1
            continue

        cols = line.split(splitBy)
        print(" ".join(cols[0:3]), file=outFile)

    outFile.close()
    return outFileP


def loadGenesPerModule(gmf, splitBy=",", header=1):
    print("\nLoading info on genes per module:")
    genesPerModule = defaultdict(list)

    for line in open(gmf):
        cols = line.rstrip().split(splitBy)

        if header:
            header -= 1
            continue

        try:
            genesPerModule[int(cols[1])].append(cols[0])
        except IndexError:
            print("Issue with line:", line)
            sys.exit()

    print("\tDone")
    return genesPerModule


def loadSingleCorrInfo(
    scf,
    splitBy="\t",
    header=1,
    geneCol=0,
    varWithinCol=4,
    medianPercExpCol=5,
    varAcrossCol=6,
):
    infoPerGene = {}

    for line in gzip.open(scf):
        cols = line.decode().split(splitBy)

        if header:
            header -= 1
            continue

        gene = cols[geneCol]
        varWithin = float(cols[varWithinCol])
        varAcross = float(cols[varAcrossCol])
        medianPercentile = float(cols[medianPercExpCol])

        infoPerGene[gene] = [medianPercentile, varWithin, varAcross]

    print("\tInfo for:", len(infoPerGene), "genes loaded")
    return infoPerGene


def getDegreePerModule(
    wF,
    outFile,
    edgesFP,
    geneModuleFile,
    singleCorrFile,
    weightedDegree=True,
    roundTo=3,
    useSingleCorrInfo=False,
):
    # Convert edge file to format for igraph
    tempP = convertGephiToIgraph(edgesFP)

    #  Input graph
    print("\nConverting gephi edge file --> igraph edge file")
    g = ig.Graph().Read_Ncol(tempP, directed=False)
    os.remove(tempP)

    print("\tNodes:", len(g.vs), "Edges:", len(g.es))

    if useSingleCorrInfo:
        print("\nLoading information from single corr file")
        singleCorrInfo = loadSingleCorrInfo(
            singleCorrFile,
            geneCol=SCF_GENE_COL,
            varWithinCol=SCF_VAR_WITHIN_COL,
            varAcrossCol=SCF_VAR_ACROSS_COL,
            medianPercExpCol=SCF_MEDIAN_EXP_COL,
        )

    # Get genes per module
    genesPerModule = loadGenesPerModule(geneModuleFile)

    #################################################################
    # Get output file ready
    oFile = open(os.path.join(wF, outFile), "w")
    headerStr = "Module\tModRank\tGene"

    if weightedDegree:
        headerStr = headerStr + "\tWeightedDegree"
    else:
        headerStr = headerStr + "\tDegree"

    if useSingleCorrInfo:
        headerStr = (
            headerStr
            + "\tModCon (degree^2 * VarWithin  * (100 - VarAcross)/100)\tPercentileExp\tVarWithin\t(100-VarAcross)/100"
        )

    print(headerStr, file=oFile)
    #################################################################

    print("\nGetting per module gene degree in subnetworks")
    perModGeneDegree = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    # Get subgraphs per module
    for mod in sorted(genesPerModule):
        modGenes = genesPerModule[mod]
        subnet = g.subgraph(modGenes)
        print("\tModule:", mod, "Nodes:", len(subnet.vs), "Edges:", len(subnet.es))

        for gene in sorted(modGenes):
            if weightedDegree:
                degree = round(subnet.strength(gene, weights="weight"), roundTo)
            else:
                degree = subnet.strength(gene)

            if useSingleCorrInfo:
                geneInf = singleCorrInfo[gene]
                percentileExp = round(geneInf[0], roundTo)
                varWithin = geneInf[1]
                varAcross = (100 - geneInf[2]) / 100.00
                modScore = (degree**2) * varWithin * varAcross

                percentileExp = round(percentileExp, roundTo)
                varWithin = round(varWithin, roundTo)
                varAcross = round(varAcross, roundTo)
                modScore = round(modScore, roundTo)

                perModGeneDegree[mod][modScore][gene] = [
                    degree,
                    percentileExp,
                    varWithin,
                    varAcross,
                ]
            else:
                perModGeneDegree[mod][degree][gene] = 1

    #  Print out ordered data
    natKey = make_key_naturalSort()
    modRank = 1
    for mod in sorted(perModGeneDegree):
        for rankBy in sorted(perModGeneDegree[mod], reverse=True):
            for gene in sorted(perModGeneDegree[mod][rankBy], key=natKey):
                if useSingleCorrInfo:
                    degree, percentileExp, varWithin, varAcross = perModGeneDegree[mod][
                        rankBy
                    ][gene]
                    print(
                        mod,
                        modRank,
                        gene,
                        degree,
                        rankBy,
                        percentileExp,
                        varWithin,
                        varAcross,
                        sep="\t",
                        file=oFile,
                    )
                else:
                    print(mod, modRank, gene, rankBy, sep="\t", file=oFile)

            modRank += 1  # Genes with same score get same rank

        modRank = 1


def runMain(finishT="Finished!"):
    # global # variable.  Declare each global variable and pass new values/defaults
    pass


def main(finishT="Finished!"):
    getDegreePerModule(
        workFolder,
        OUT_FILE,
        edgesFP,
        geneModuleFile,
        singleCorrFile,
        weightedDegree=WEIGHTED_DEGREE,
        roundTo=ROUND_VALUES_TO,
        useSingleCorrInfo=USE_SINGLE_CORR_INFO,
    )
    print(finishT)


#################################################################

if __name__ == "__main__":
    main()
