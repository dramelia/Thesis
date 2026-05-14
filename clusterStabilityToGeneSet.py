#!/usr/bin/env python
"""
Simple script to convert the output from clusterStability.py into a geneset for gene signature enrichment analysis
"""


import sys
import os
import glob
import time
import re

workFolder = r"PATH/PGCNA_normB_0.3/results/PGCNA/EPG3/LEIDENALG/BEST/ClstStability" # change path
stabilityInfoF = r"Clust608" # chnage cluster number

###########################################################################################


bits = stabilityInfoF.split("_")

MIN_PERCENT_IN = 80  #  50 What percentage of clusters does a gene need to be in for inclusion.  80 upwards is good.
PREFIX = "MYNAME_" + "_".join(bits[0:1]) + "_PGCNA_"  # Gene set prefix, normally MYNAME_ (change name)

STOP_IF_GENES_MISSING = True  # If True will halt if genes missing at MIN_PERCENT_IN threshold else will just skip
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
        convert = lambda text: int(text) if text.isdigit() else text
        alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", key)]

        return alphanum_key(s)

    return nSort


def makeGeneSets(wF, sIF, minPercIn=80, splitBy="\t", prefix="CARE_", trimS=".txt_geneSums.txt", stopIt=False):

    natSort = make_key_naturalSort()
    trimA = len(trimS)

    outFile = open(os.path.join(wF, "geneSet_" + sIF + "_minCI" + str(minPercIn) + ".txt"), "w")

    geneStability = {}

    for path in sorted(glob.glob(os.path.join(wF, sIF, "*geneSums*")), key=natSort):
        fileN = os.path.basename(path)
        print("\t...", fileN)

        header = 1
        genes = {}
        for line in open(path):
            cols = line.rstrip().split(splitBy)

            if header:
                header -= 1
                continue

            gene = cols[0]
            count = int(cols[1])
            perc = float(cols[2])

            if gene in geneStability:
                print("Gene:", gene, "appears in more than one cluster!")
                sys.exit()
            else:
                geneStability[gene] = [count, perc, fileN[:-trimA]]

            if perc >= minPercIn:
                genes[gene] = 1

        if len(genes) == 0:
            print("\t\t\t#File:", fileN, "Has no genes at this threshold, will not write out!\n")
            if stopIt:
                sys.exit()
            time.sleep(0.5)
        else:
            print(
                prefix + fileN[:-trimA],
                prefix + fileN[:-trimA],
                "\t".join(sorted(genes, key=natSort)),
                sep="\t",
                file=outFile,
            )
            print("\t\tgenes added:", len(genes))

    outFile.close()

    outFile2 = open(os.path.join(wF, "totalStability_" + sIF + ".txt"), "w")
    print("Module\tGene\tPercentOverlap", file=outFile2)
    for gene in sorted(geneStability, key=natSort):
        count, perc, fileN = geneStability[gene]
        print(fileN, gene, perc, sep="\t", file=outFile2)

    outFile2.close()


########################################################################################
##----------------------------------------MAIN----------------------------------------##
########################################################################################


def main(finishT="Finished!"):

    makeGeneSets(workFolder, stabilityInfoF, minPercIn=MIN_PERCENT_IN, prefix=PREFIX, stopIt=STOP_IF_GENES_MISSING)

    print("\n", finishT, sep="")


##----------------------------------------------------------------------------------------##
##----------------------------------------__main__----------------------------------------##
if __name__ == "__main__":
    main()