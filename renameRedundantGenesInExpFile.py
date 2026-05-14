#!/usr/bin/env python
"""
Given a folder of PGCNA data files along with the meta info file will make sure that all redundancy in 
gene identifiers is removed e.g. TRA --> TRA#1 , TRA#2 etc

Importantly it retains the original order of the input data
"""
import sys
import os
import re
from collections import defaultdict
import time


workFolder = r"PATH/Post_R_Lopes" # change path
folderToProcess = r"ForPGCNA" # change folder
metaInfoName = r"#FileInfo.txt"

outFolderSuffix = "_collapseR"

#################################################################
GENE_SPLIT_CHAR = "#"
COL_SPLIT = "\t"
GENE_COL = 0

#################################################################
OUT_FOLDER = os.path.join(workFolder, folderToProcess + outFolderSuffix)
if not os.path.exists(OUT_FOLDER):
    os.makedirs(OUT_FOLDER)

#################################################################

# ----------------------------Methods----------------------------#
def loadMeta(metaFile, dataF, splitBy="\t", headerL=1):

    print("\n\nLoad meta-file (", os.path.basename(metaFile), ")", sep="")
    if not os.path.exists(metaFile):
        print("\t\t# Meta file (", metaFile, ") does not exist!", sep="")
        sys.exit()

    header = headerL
    fileInfo = {}
    for line in open(metaFile):
        cols = line.rstrip().split(splitBy)

        if header:
            header -= 1
            continue

        if line.rstrip() == "":
            continue

        totalPath = os.path.join(dataF, cols[0])
        if not os.path.exists(totalPath):
            print("\t\t# File (", totalPath, ") does not exist!, won't add to fileInfo", sep="")
        else:
            try:
                fileInfo[totalPath] = int(cols[1])
            except:
                print("Meta file line (", line.rstrip(), ") is not formed properly, skipping")

    print("\tLoaded information on:", len(fileInfo), "files")
    return fileInfo


def collapseRed(path, outFold, headerL=1, geneSplitChar="#", splitBy="\t", geneCol=0):

    fileName = os.path.basename(path)

    print("\t\t", fileName, " exists, processsing")

    genesSeen = defaultdict(dict)

    headerStr = ""
    origData = []
    headerLines = []

    for line in open(path):
        cols = line.rstrip().split(splitBy)

        if headerL:
            headerL -= 1
            headerLines.append(splitBy.join(cols))
            continue

        #  Warn about unwanted spaces in the identifier
        if " " in cols[geneCol]:
            print("Gene (", cols[geneCol], ") has spaces in it! will remove them but is this what you want?", sep="")
            cols[geneCol] = cols[geneCol].replace(" ", "")
            time.sleep(5)

        genesSeen[cols[geneCol]][str(cols[geneCol + 1 :])] = 1
        origData.append(cols)

    print("\t\t\tTotal non-redundant genes processed:", len(genesSeen))

    outFile = open(os.path.join(outFold, fileName), "w")
    print("\n".join(headerLines), file=outFile)
    totalDuplicates = 0
    removedAllSame = 0
    perGenePos = {}
    allowedRep = {}
    alreadySeen = {}

    for cols in origData:

        # Deal with duplicates that are exactly the same
        if cols[geneCol] in alreadySeen:
            if not cols[geneCol] in allowedRep:
                removedAllSame += 1
                continue

        if len(genesSeen[cols[geneCol]]) > 1:
            allowedRep[cols[geneCol]] = 1

            totalDuplicates += 1
            if cols[geneCol] not in perGenePos:
                #  We've not seen this yet
                perGenePos[cols[geneCol]] = 1
                cols[geneCol] = cols[geneCol] + geneSplitChar + "1"

            else:
                perGenePos[cols[geneCol]] += 1
                cols[geneCol] = cols[geneCol] + geneSplitChar + str(perGenePos[cols[geneCol]])

        alreadySeen[cols[geneCol]] = 1
        print(splitBy.join(cols), file=outFile)

    outFile.close()
    print("\t\t\t\tTotal duplicates:", totalDuplicates)
    print("\t\t\t\tTotal removed as duplicates exactly same:", removedAllSame)
    print("\n\n")


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


def runMain(finishT="Finished!"):
    # global # variable.  Declare each global variable and pass new values/defaults
    pass


def main(finishT="Finished!"):

    naturalSort = make_key_naturalSort()
    processFolder = os.path.join(workFolder, folderToProcess)
    metaInf = loadMeta(os.path.join(processFolder, metaInfoName), processFolder)

    print("Collapsing redundancy for files:")
    for path in sorted(metaInf, key=naturalSort):
        headerL = metaInf[path]

        print("\t", path)
        if os.path.exists(path):
            collapseRed(
                path, OUT_FOLDER, headerL=headerL, geneSplitChar=GENE_SPLIT_CHAR, splitBy=COL_SPLIT, geneCol=GENE_COL
            )

    print(finishT)


#################################################################

if __name__ == "__main__":
    main()
