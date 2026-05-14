#!/usr/bin/env python
"""
Given a folder that contains subfolders with 1 (or more) geneModuleWeightedDegree_* files (created
by getWeightedDegreePerModule.py)

will gather all the information across folders and create a final matrix of gene usage across modules
at the same time will try to autoselect the best subset of genes that covers all the modules.

remModuleFile should have one line of header and then contain the data-set (after PREFIX_REM/SUFFIX_REM applied) and modules that should be
removed (separated by tabs). e.g. :

Dataset	Module
MM	28
Bcell_APRIL_STC	12
Bcell_APRIL_STC	14
Bcell_APRIL_STC	16

#  Additional required files
> A GENCODE file to add gene type annotations.  Expects this to be a gzipped file.
> A list of TFs

##################################################################
#                  AUTOSELECT GENES METHOD						 #
Will also try to automatically select the 'best' genes by ranking them using:
autoRank = AUTOSELECT_COUNT_SCALE_FACTOR * datasetsGeneIn + sum(normalised ranks in datasetsGeneIn)

Total method:
> For each data-set/module sort by base rank (set using RF_RANK_COL)
        > Create pool of genes == AUTOSELECT_POOL_SIZE, for autoselection
        > Store in new dictionary with the key set to autoRank (previously calculated)
        > Sort the pool genes by autoRank
        > Select the AUTOSELECT_GENE_NUMBER genes
                > If gene has already been used skip to next ranked gene
                > Store used gene in dictionary to avoid duplications
> If AUTOSELECT_PRINT_NOT_CHOSEN == True then output genes that weren't selected at bottom of file
##################################################################
"""

import glob
import gzip
import os
import re
import sys
import time
from collections import defaultdict

import numpy as np

workFolder = r"PATH/PGCNA_work/PostPGCNA_normB/MEV_Vis" # change path
folderToProcess = r"CLL_GEO_gmwd"
remModuleFile = (
    None  
)

# GENCODE Annotation file
gtfFile = r"PATH/PGCNA_work/genomes/gencode.v48.EBV.primary_assembly.annotation.gtf.gz"  # HUMAN  -- latest # change path

# TF DB file
tfDbFile = r"PATH/PGCNA_work/PostPGCNA_normB/MEV_Vis/Merged_TF.txt" # change path

##################################################################
#Change the suffix blow to "_ModCon_myGeneAnnotated" if its rerquired the myGene API annotation (this may be done at a different stage for different analyses). API done post R in thi pipeline.
FIND_SUFFIXES = "_ModCon_out.txt"

RF_MOD_COL = 0  # Normally 0
RF_RANK_COL = (
    1  # Normally 1   # Actual rank i.e. from 1-->n, not expecting a value to rank by!
)
RF_EXP_COL = 5  # Normally 5
RF_GENE_COL = 8  # Normally 2, or 9 if re-annotated

# PREFIX_REM = "geneModuleWeightedDegree_" # When not re-annotated
PREFIX_REM = "geneModuleWeightedDegree_"  #  For when re-annotated
SUFFIX_REM = FIND_SUFFIXES
REM_MODULES = False

##################################################################
# FOR AUTORANKING
AUTOSELECT_GENE_NUMBER = 10  # 5  Number of genes to autoselect for each module
AUTOSELECT_POOL_SIZE = 25  # 25. Max rank from which to select from
AUTOSELECT_COUNT_SCALE_FACTOR = 0.25  # 0.25 seems good.  Scale factor; used to generate autoRank = ACSF * datasetsGeneIn + sum(normalised ranks in datasetsGeneIn)
AUTOSELECT_PRINT_NOT_CHOSEN = (
    False  # If True will print out all unselected data at end of output
)

##################################################################
# FILTERS
FILTER_GENE_TYPE_UNKNOWN = (
    True  # Will remove all genes that aren't annotated in the Gencode GTF file
)

##################################################################
# GENE TYPE ANNOTATION RELATED
GENE_SPLIT_BY = "#"  # How to split multigene IDs

#  GTF FILE RELATED
GTF_HEADER_LINES = (
    5  #  Normally 5.  Will print these out so that user can see what is being discarded
)
GTF_EXPECTED_COLS = 9  # Normally 9
GTF_PRINT_EVERY = 1e5  #  Show when every GTF_PRINT_EVERY lines have been processed

# GTF regular expressions
GTF_GENE_ID_RE = re.compile(r"gene_id\s+([^\n]+)")
GTF_GENE_TYPE_RE = re.compile(r"\s?gene_type\s+([^\n]+)")
GTF_TRANSCRIPT_ID_RE = re.compile(r"\s?transcript_id\s+([^\n]+)")
GTF_GENE_NAME_RE = re.compile(r"\s?gene_name\s+([^\n]+)")

##################################################################
# TF INFORMATION
TF_ADD_INFO = True  #  Add information about TFs to genes


##################################################################
OUT_FOLDER = "AutoselectResults"
##################################################################

##################################################################
##################################################################
OUT_FOLDER = os.path.join(workFolder, OUT_FOLDER)
if not os.path.exists(OUT_FOLDER):
    os.makedirs(OUT_FOLDER)

# ----------------------------Methods----------------------------#


def runMain(finishT="Finished!"):
    # global # variable.  Declare each global variable and pass new values/defaults
    pass


def loadGTF(
    gtfFile,
    headerL=5,
    splitBy="\t",
    splitByGeneInf=";",
    expectedCols=9,
    gtfPrintEvery=1000,
    geneInfoCol=8,
    gtfGeneIdRe=None,
    gtfGeneTypeRe=None,
    gtfTranscriptIdRe=None,
    gtfGeneNameRe=None,
):
    def readNRfile(
        path,
        splitBy="\t",
        headerL=1,
        transcriptCol=0,
        geneIdCol=1,
        geneNameCol=2,
        typeCol=3,
    ):
        print("\tReading non-redundant transcript file:")
        geneFI = defaultdict(dict)

        for line in gzip.open(path):
            cols = line.decode().rstrip().split(splitBy)

            if headerL:
                headerL -= 1
                continue

            geneFI[cols[geneNameCol]][cols[typeCol]] = 1

        #  Make sure genes only have one type assigned
        for g in geneFI:
            if len(geneFI[g]) > 1:
                print("\tMultiple types for:", g)

        print("\t\tNon-redundant gene number:", len(geneFI))
        return geneFI

    ###########################################################################################

    gtfFolder, gtfFileN = os.path.split(gtfFile)
    nrFilePath = os.path.join(gtfFolder, "nr_" + gtfFileN[:-7]) + ".txt.gz"

    if not os.path.exists(nrFilePath):
        header = headerL
        printEvery = gtfPrintEvery
        tell = printEvery
        count = 0
        startTime = time.time()

        transcriptLevelInfo = defaultdict()

        print("Parsing GTF file:", gtfFileN)
        print(
            "\t##----------------------------------------GTF HEADER----------------------------------------##"
        )
        for line in gzip.open(gtfFile):
            line = line.decode().rstrip().replace('"', "")
            cols = line.split(splitBy)

            count += 1

            if header:
                header -= 1
                print("\t\t", line)

                if not header:
                    print("\nProcessing GTF main body:")
                continue

            #  Test all lines have expected number of columns
            if len(cols) != expectedCols:
                print(
                    "Number of cols (",
                    len(cols),
                    "!= expected number (",
                    expectedCols,
                    ")",
                    sep="",
                )
                sys.exit()

            geneInfoCols = cols[geneInfoCol].split(splitByGeneInf)

            gene_id = None
            gene_type = None
            transcript_id = None
            gene_name = None

            for i in geneInfoCols:
                if gtfGeneIdRe.match(i):
                    gene_id = gtfGeneIdRe.search(i).groups()[0]

                if gtfGeneTypeRe.match(i):
                    gene_type = gtfGeneTypeRe.search(i).groups()[0]

                if gtfTranscriptIdRe.match(i):
                    transcript_id = gtfTranscriptIdRe.search(i).groups()[0]

                if gtfGeneNameRe.match(i):
                    gene_name = gtfGeneNameRe.search(i).groups()[0]

            if transcript_id:
                if not transcript_id in transcriptLevelInfo:
                    transcriptLevelInfo[transcript_id] = [gene_id, gene_name, gene_type]
                else:
                    oGeneID, oGene, oGeneType = transcriptLevelInfo[transcript_id]

                    if oGene != gene_name:
                        print(
                            "For transcript_id:",
                            transcript_id,
                            "have conflicting gene names:",
                            oGene,
                            "!=",
                            gene_name,
                        )
                        sys.exit()

                    if oGeneType != gene_type:
                        print(
                            "For transcript_id:",
                            transcript_id,
                            "have conflicting gene types:",
                            oGeneType,
                            "!=",
                            gene_type,
                        )
                        sys.exit()

            if count == printEvery:
                printEvery = printEvery + tell
                print(
                    "\t\t{0:,d}".format(int(count)),
                    "in",
                    round(time.time() - startTime, 2),
                    "seconds",
                )
        print("Total GTF lines processed:{0:,d}".format(int(count)))

        ###########################################################################################
        #  Write out a smaller non-redundant version of the GTF file to save time in the future
        print("\nWriting out non-reudndant transcript information")
        natSort = make_key_naturalSort()
        nrFile = gzip.open(nrFilePath, "wb")
        nrFile.write("transcript_id\tgene_id\tgene_name\tgene_type".encode())
        for tID in sorted(transcriptLevelInfo, key=natSort):
            gene_id, gene_name, gene_type = transcriptLevelInfo[tID]
            nrFile.write(
                "\t".join([tID, gene_id, gene_name, gene_type + "\n"]).encode()
            )

        nrFile.close()
        ###########################################################################################
    else:
        print("Non-redundant file:", nrFilePath, "already exists, will use that\n")

    #  Load file back into memory
    geneFI = readNRfile(nrFilePath)

    return geneFI


def mad(arr):
    """Median Absolute Deviation: a "Robust" version of standard deviation.
    Indices variabililty of the sample.
    https://en.wikipedia.org/wiki/Median_absolute_deviation
    """
    # arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    arr = np.array(arr)  # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med)), med


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


def make_key_naturalSort_array(arrayPos=0):
    """
    A factory function: creates a key function to use in sort.
    Sort data naturally.

    This is for an array/dictionary that contains arrays.  Where we
    want to sort the arrays contained.  Can choose what element of
    the array to sort by
    """

    def nSort(s):
        convert = lambda text: int(text) if text.isdigit() else text
        alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", key)]

        return alphanum_key(s[arrayPos])

    return nSort


def getGeneRankFilesList(wF, fTP, suffix="_ModCon.txt"):
    """
    Find the paths to the geneRankPerModule files.
    """
    geneRankFiles = defaultdict(list)
    natKey = make_key_naturalSort()

    print("\nFind list of files that have suffix=", suffix, ":")
    for dirnameT, dirnamesT, filenamesT in os.walk(os.path.join(wF, fTP)):
        for fileName in filenamesT:
            if fileName[-len(suffix) :] == suffix:
                geneRankFiles[fileName].append(os.path.join(dirnameT, fileName))

    # Print out info for quick reality check
    for fN in sorted(geneRankFiles, key=natKey):
        if (len(geneRankFiles[fN])) > 1:
            print("\n\tFile ", fN, "occurs in multiple locations:")
            for path in sorted(geneRankFiles[fN], key=natKey):
                print("\t\t", path)

            print("\n\tOnly expecting one version of each file please remove one copy!")
            sys.exit()
        print("\t", fN, sep="")

    print("\n\tTotal:", len(geneRankFiles))

    return geneRankFiles


def readRemModFile(rMF, splitBy="\t", header=1):
    remModInfo = {}

    print("\nWill remove modules in", rMF, ":")
    for line in open(rMF):
        cols = line.rstrip().split(splitBy)

        if header:
            header -= 1
            continue

        dataSetMod = cols[0] + "_M" + cols[1]
        remModInfo[dataSetMod] = 1

    print("\tTotal modules to remove:", len(remModInfo))
    return remModInfo


def loadTfInfo(tfdb, splitBy="\t", header=1, geneCol=0, tfFamilyCol=1, tfSeenCol=2):
    tfInfo = {}

    print("\nLoading TF information:")
    for line in open(tfdb):
        cols = line.rstrip().split(splitBy)

        if header:
            header -= 1
            continue

        try:
            gene = cols[geneCol]
            tfFamily = cols[tfFamilyCol]
            tfSeenCount = cols[tfSeenCol]
        except:
            print("\t## Issue with line :", line)
            continue

        if gene != "":
            tfInfo[gene] = [tfFamily, tfSeenCount]

    print("\tLoaded info on ", len(tfInfo), "TFs")

    return tfInfo


def processFiles(
    oF,
    fTP,
    fileInfo,
    geneFI,
    tfDbFile,
    splitBy="\t",
    headerL=1,
    rfModCol=0,
    rfRankCol=1,
    rfGeneCol=2,
    rfExpCol=5,
    prefixRem="geneModuleWeightedDegree_",
    suffixRem="_ModCon.txt",
    remModules=False,
    remModuleFile=None,
    autoselectGeneNumber=5,
    autoselectPoolSize=100,
    autoselectPrintNotChosen=False,
    autoselectCountScaleFactor=1,
    geneSplit="#",
    geneTypeMissing="Unknown",
    filterGeneTypeUnknown=False,
    roundExpTo=1,
    tfAddInfo=True,
):
    ##################################################################
    # setup
    if remModules:
        remModInfo = readRemModFile(remModuleFile)
        outFileSuffix = "_remMods"
    else:
        outFileSuffix = ""

    outFileSuffix = (
        outFileSuffix
        + "_autoSN"
        + str(autoselectGeneNumber)
        + "_PS"
        + str(autoselectPoolSize)
        + "_autoSF"
        + str(autoselectCountScaleFactor)
    )

    if filterGeneTypeUnknown:
        outFileSuffix = outFileSuffix + "_filtUnkGT"

    if tfAddInfo:
        outFileSuffix = outFileSuffix + "_TFinfo"

    if autoselectPrintNotChosen:
        outFileSuffix = outFileSuffix + "_keepNonSel"
    ##################################################################

    ##################################################################
    # Load TF info
    if tfAddInfo:
        tfInfo = loadTfInfo(tfDbFile)
    ##################################################################

    natKey = make_key_naturalSort()

    mergedInfo = defaultdict(dict)
    dataSetMods = defaultdict(int)

    print("\nProcessing files:")
    for fN in sorted(fileInfo, key=natKey):
        filePath = fileInfo[fN][0]
        fNsimple = fN[len(prefixRem) : -len(suffixRem)]

        # First pass to gather data
        moduleGenes = defaultdict(lambda: defaultdict(dict))
        header = headerL
        for line in open(filePath):
            cols = line.rstrip().split(splitBy)

            if header:
                header -= 1
                continue
            module = cols[rfModCol]
            geneModRank = int(cols[rfRankCol])
            gene = cols[rfGeneCol]
            expression = float(cols[rfExpCol])

            moduleGenes[module][geneModRank][gene] = expression

        #  Second pass sort data by module/rank
        #  and calculate normalised ranks
        modulesSeen = {}
        for module in sorted(moduleGenes):
            dataSetMod = fNsimple + "_M" + module
            if remModules:
                if dataSetMod in remModInfo:
                    continue

            maxRank = max(moduleGenes[module])
            minRank = min(moduleGenes[module])

            dataSetMods[dataSetMod] += 1
            modulesSeen[module] = 1

            for rank in sorted(moduleGenes[module]):
                # Normalised rank, from 1 (best) --> 0 (worst)
                normRank = (maxRank - rank) / (maxRank - minRank)

                for gene in sorted(moduleGenes[module][rank], key=natKey):
                    ##################################################################
                    # ----------------Get gene type and TF information----------------#
                    ##################################################################
                    types = {}
                    tfTypes = {}
                    tfType = ""
                    tfStatus = ""

                    for subGene in gene.split(geneSplit):
                        if subGene in geneFI:
                            for gT in sorted(geneFI[subGene], key=natKey):
                                types[gT] = 1

                        if tfAddInfo:
                            if subGene in tfInfo:
                                tfFamily, tfSeenCount = tfInfo[subGene]
                                tfTypes[tfFamily] = 1
                    if types:
                        geneType = ",".join(types.keys())
                    else:
                        geneType = geneTypeMissing

                    if tfTypes:
                        tfType = ",".join(tfTypes.keys())
                        tfStatus = "+"

                    ##################################################################

                    if filterGeneTypeUnknown and (geneType == geneTypeMissing):
                        continue

                    mergedInfo[gene][dataSetMod] = [
                        module,
                        normRank,
                        rank,
                        geneType,
                        moduleGenes[module][rank][gene],
                        tfStatus,
                        tfType,
                    ]

        print("\t", fN, "-->", fNsimple, ":", len(modulesSeen))

    ##################################################################
    # Final pass gather everything together
    datasetModOrder = sorted(dataSetMods.keys(), key=natKey)
    print("\n\tTotal modules across data-sets:", len(datasetModOrder))

    dataByRank = defaultdict(lambda: defaultdict(list))
    nonRedundant = dict()

    for gene in sorted(mergedInfo, key=natKey):
        geneInf = []
        geneRanks = []
        normGeneRanks = []
        geneExpression = []

        perDatasetInf = dict()

        for dsM in datasetModOrder:
            if dsM in mergedInfo[gene]:
                module, normRank, geneRank, geneType, expression, tfStatus, tfCount = (
                    mergedInfo[gene][dsM]
                )
                geneInf.append(str(geneRank))

                geneRanks.append(geneRank)
                normGeneRanks.append(normRank)
                geneExpression.append(expression)

                perDatasetInf[dsM] = geneRank

            else:
                geneInf.append("")

        madRank, medianRank = mad(geneRanks)
        madExp, medianExp = mad(geneExpression)

        autoRank = (autoselectCountScaleFactor * len(geneRanks)) + sum(normGeneRanks)

        # Store info for final output
        totalInfo = [
            autoRank,
            gene,
            geneType,
            len(geneRanks),
            min(geneRanks),
            medianRank,
            max(geneRanks),
            madRank,
            ";".join(map(lambda s: str(s), geneRanks)),
            min(geneExpression),
            medianExp,
            max(geneExpression),
            madExp,
            ";".join(map(lambda s: str(round(s, roundExpTo)), geneExpression)),
            tfStatus,
            tfCount,
            "\t".join(geneInf),
        ]

        # Store for final output of non selected data
        nonRedundant[gene] = totalInfo

        #  Store info for easy sorting
        #  This will obviously introduce a massive amount of redundancy!
        for dsM in perDatasetInf:
            dataByRank[dsM][perDatasetInf[dsM]].append(totalInfo)

    ##################################################################
    # -------------Per datasetModule autoselect genes-----------------#
    print(
        "\nAutoselecting top:",
        autoselectGeneNumber,
        "from pool of",
        autoselectPoolSize,
        ":\n",
    )
    outFile = open(
        os.path.join(oF, "mergedDegreePerModuleAS_" + fTP + outFileSuffix + ".txt"), "w"
    )

    print(
        "sortOrder\tselected\tautoRank\tgene\tgene-type\tcountSeen\tminRank\tmedianRank\tmaxRank\tMAD_rank\tranks\tminExp\tmedianExp\tmaxExp\tMAD_Exp\texpressionPercVals\tTF\tTF-Type",
        "\t".join(datasetModOrder),
        sep="\t",
        file=outFile,
    )

    genesSelected = {}
    natArraySort = make_key_naturalSort_array(arrayPos=1)
    sortOrder = 1  #  To allow quick resorting in Excel

    for dsM in datasetModOrder:
        byAutoRank = defaultdict(list)
        added = 0

        # Grab a pool of genes as large as autoselectPoolSize from which to work
        for rank in sorted(dataByRank[dsM]):
            anySelected = False
            for info in sorted(dataByRank[dsM][rank]):
                if not info[1] in genesSelected:  # Don't duplicate, always select new
                    byAutoRank[info[0]].append(info)
                    anySelected = True

            if anySelected:
                added += 1
            if added >= autoselectPoolSize:
                break

        selected = 0

        # Rank pool of genes and select the top autoselectGeneNumber
        for aR in sorted(byAutoRank, reverse=True):
            for info in sorted(byAutoRank[aR], key=natArraySort):
                print(
                    sortOrder,
                    "Y",
                    "\t".join(map(lambda x: str(x), info)),
                    sep="\t",
                    file=outFile,
                )
                genesSelected[info[1]] = 1
                sortOrder += 1

            selected += 1
            if selected >= autoselectGeneNumber:
                break

        print("\t", dsM, ":", selected)

    if autoselectPrintNotChosen:
        # Add remaining unseen data
        for gene in sorted(nonRedundant, key=natKey):
            if not gene in genesSelected:
                info = nonRedundant[gene]
                print(
                    sortOrder,
                    "",
                    "\t".join(map(lambda x: str(x), info)),
                    sep="\t",
                    file=outFile,
                )
                sortOrder += 1

    outFile.close()
    ##################################################################
    ##################################################################


def main(finishT="Finished!"):
    # Prepare Gencode GTF file for annotation
    print("\nLoad Gencode GTF file for gene annotation information:")
    geneFI = loadGTF(
        gtfFile,
        headerL=GTF_HEADER_LINES,
        expectedCols=GTF_EXPECTED_COLS,
        gtfPrintEvery=GTF_PRINT_EVERY,
        gtfGeneIdRe=GTF_GENE_ID_RE,
        gtfGeneTypeRe=GTF_GENE_TYPE_RE,
        gtfTranscriptIdRe=GTF_TRANSCRIPT_ID_RE,
        gtfGeneNameRe=GTF_GENE_NAME_RE,
    )

    fileInfo = getGeneRankFilesList(workFolder, folderToProcess, suffix=FIND_SUFFIXES)

    processFiles(
        OUT_FOLDER,
        folderToProcess,
        fileInfo,
        geneFI,
        tfDbFile,
        rfModCol=RF_MOD_COL,
        rfRankCol=RF_RANK_COL,
        rfGeneCol=RF_GENE_COL,
        rfExpCol=RF_EXP_COL,
        prefixRem=PREFIX_REM,
        suffixRem=SUFFIX_REM,
        remModules=REM_MODULES,
        remModuleFile=remModuleFile,
        autoselectGeneNumber=AUTOSELECT_GENE_NUMBER,
        autoselectPoolSize=AUTOSELECT_POOL_SIZE,
        autoselectPrintNotChosen=AUTOSELECT_PRINT_NOT_CHOSEN,
        autoselectCountScaleFactor=AUTOSELECT_COUNT_SCALE_FACTOR,
        geneSplit=GENE_SPLIT_BY,
        filterGeneTypeUnknown=FILTER_GENE_TYPE_UNKNOWN,
        tfAddInfo=TF_ADD_INFO,
    )

    print(finishT)


#################################################################

if __name__ == "__main__":
    main()
