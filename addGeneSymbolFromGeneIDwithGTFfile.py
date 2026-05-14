#!/usr/bin/env python
"""
Load Gencode GTF file and given:
> gene_id
> transcript_id

Then will find all files (using glob) that match FILE_FILT and map to gene_name (symbol), and also add gene_type column.

NOTE: expects the GTF file to be gzipped!
"""

import glob
import gzip
import os
import re
import sys
import time
from collections import defaultdict

import natsort as ns
from natsort import natsort_keygen

workFolder = r"/PATH/Post_R_Lopes" # change path
gtfFile = r"/PATH/gencode.v48.EBV.primary_assembly.annotation.gtf.gz"  # HUMAN + EBV  ## change path

##----------------------------------------CONSTANTS----------------------------------------##
# Programs like DEXSeq output multiple genes split by + want to output info for each
SPLIT_GENE_CHAR = "+"

#  FILES TO PROCESS RELATED
# PROCESS_FILE_GLOB = "*count_WithDispInf.txt"  # Add dispersion info first with addDispersionInfoToFile.py
# PROCESS_FILE_GLOB = "*count.txt"  # Add dispersion info first with addDispersionInfoToFile.py
# Add dispersion info first with addDispersionInfoToFile.py
# PROCESS_FILE_GLOB = "*WithDispInf.txt"
# PROCESS_FILE_GLOB = "*_VST-ComBat.txt"
PROCESS_FILE_GLOB = (
    "*VST.txt"  # Add dispersion info first with addDispersionInfoToFile.py
)
# PROCESS_FILE_GLOB = "*rLog.txt"  # Add dispersion info first with addDispersionInfoToFile.py
# PROCESS_FILE_GLOB = "*DiffExp.txt"  # Add dispersion info first with addDispersionInfoToFile.py
# PROCESS_FILE_GLOB = "EnsemblIDs.FPKM.txt"  # Add dispersion info first with addDispersionInfoToFile.py
# PROCESS_FILE_GLOB = "*.txt"  # Add dispersion info first with addDispersionInfoToFile.py
# PROCESS_FILE_GLOB = "*.csv"  # Add dispersion info first with addDispersionInfoToFile.py
# PROCESS_FILE_GLOB = "*_Result.txt"  # DEXSeq output
# PROCESS_FILE_GLOB = "*_LGPM100.txt"  # DEXSeq output
# PROCESS_FILE_GLOB = "*_Cos2.txt"  # DEXSeq output
# PROCESS_FILE_GLOB = "*.tsv"
# PROCESS_FILE_GLOB = "vsd_mat_Batch2.csv"

TO_EXCLUDE_SUFFIXES = [
    "samples.txt",
    "wGeneSymbols.txt",
    "wGeneSymbols.csv",
    "TXIcountSum.txt",
]  # Should always include wGeneSymbols.txt!
PROCESS_FILES_HEADER_LINES = 1
PROCESS_FILES_ID_COL = 0  # 0 .Normally 0, for csv maybe 1
# PROCESS_FILES_SPLIT_BY = "\t"  # "\t"  # Normally \t
PROCESS_FILES_SPLIT_BY = "\t"  # "\t"  # Normally \t
PROCESS_FILES_EXTRA_INFO = (
    # If True will append all transcipts per gene.  This is probably not needed as RSEM already outputs this.
    False
)

OVERWRITE_FILES = True  # If True will recreate files

#  GTF FILE RELATED
GTF_HEADER_LINES = (
    5  # Will print these out so that user can see what is being discarded
)
GTF_EXPECTED_COLS = 9  # Normally 9
GTF_PRINT_EVERY = 1e5  # Show when every GTF_PRINT_EVERY lines have been processed

# #  GTF FILE RELATED
# GTF_HEADER_LINES = 0  # Will print these out so that user can see what is being discarded
# GTF_EXPECTED_COLS = 9  # Normally 9
# GTF_PRINT_EVERY = 1e5 #  Show when every GTF_PRINT_EVERY lines have been processed

# GTF regular expressions
GTF_GENE_ID_RE = re.compile(r"\s?gene_id\s+([^\n]+)")
GTF_GENE_TYPE_RE = re.compile(r"\s?gene_type\s+([^\n]+)")
# GTF_GENE_TYPE_RE = re.compile(r"\s?gene_biotype\s+([^\n]+)")
GTF_TRANSCRIPT_ID_RE = re.compile(r"\s?transcript_id\s+([^\n]+)")
GTF_GENE_NAME_RE = re.compile(r"\s?gene_name\s+([^\n]+)")


natSort = natsort_keygen(alg=ns.LOCALE)

###########################################################################################

##----------------------------------------Methods----------------------------------------##


def runMain(finishT="Finished!"):
    # global # variable.  Declare each global variable and pass new values/defaults
    pass


def main(finishT="Finished!"):
    print(finishT)


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
        transcriptFI = defaultdict(dict)

        for line in gzip.open(path):
            cols = line.decode().rstrip().split(splitBy)

            if headerL:
                headerL -= 1
                continue

            geneFI[cols[geneIdCol]][cols[transcriptCol]] = [
                cols[geneNameCol],
                cols[typeCol],
            ]

            transcriptFI[cols[transcriptCol]][cols[geneIdCol]] = [
                cols[geneNameCol],
                cols[typeCol],
            ]

        print("\t\tNon-redundant gene number:", len(geneFI))
        print("\t\tNon-redundant transcript number:", len(transcriptFI))

        return geneFI, transcriptFI

    ###########################################################################################

    gtfFolder, gtfFileN = os.path.split(gtfFile)
    nrFilePath = os.path.join(gtfFolder, "nr_" + gtfFileN[:-7]) + ".txt.gz"

    if not os.path.exists(nrFilePath):
        header = headerL
        printEvery = gtfPrintEvery
        tell = printEvery
        count = 0
        startTime = time.perf_counter()

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
                    round(time.perf_counter() - startTime, 2),
                    "seconds",
                )
        print("Total GTF lines processed:{0:,d}".format(int(count)))

        ###########################################################################################
        #  Write out a smaller non-redundant version of the GTF file to save time in the future
        print("\nWriting out non-redundant transcript information")

        import io

        with gzip.open(nrFilePath, "wb") as nr_file_bin:
            with io.TextIOWrapper(nr_file_bin, encoding="utf-8") as nr_file:
                print("transcript_id\tgene_id\tgene_name\tgene_type", file=nr_file)
                for transcript_id in sorted(transcriptLevelInfo, key=natSort):
                    gene_id, gene_name, gene_type = transcriptLevelInfo[transcript_id]
                    print(
                        transcript_id,
                        gene_id,
                        gene_name,
                        gene_type,
                        sep="\t",
                        file=nr_file,
                    )
        ###########################################################################################
    else:
        print("Non-redundant file:", nrFilePath, "already exists, will use that\n")

    #  Load file back into memory
    geneFI, trancriptFI = readNRfile(nrFilePath)

    return geneFI, trancriptFI


def generateFileList(wF, pfg="*.txt", te=[], overwriteFiles=True):
    fileL = glob.glob(os.path.join(wF, pfg))

    finalPaths = []
    for fPath in sorted(fileL, key=natSort):
        root, name = os.path.split(fPath)
        keep = True

        for excluded in te:
            trimName = name[-len(excluded) :]
            if excluded == trimName:
                keep = False
                break

        if keep:
            if not overwriteFiles:
                finalPath = os.path.join(fPath[:-4] + "_wGeneSymbols.txt")
                if not os.path.exists(finalPath):
                    finalPaths.append(fPath)
            else:
                finalPaths.append(fPath)
    return finalPaths


def processFiles(
    fTP,
    gFI,
    tFI,
    headerL=1,
    idCol=0,
    splitBy="\t",
    addExtraInfo=False,
    splitGeneChar="+",
):
    print("\nAdding Gene info to files:")

    for path in sorted(fTP, key=natSort):
        basePath, fileN = os.path.split(path)

        print("\t", fileN)

        header = headerL

        outFile = open(os.path.join(basePath, fileN[:-4] + "_wGeneSymbols.tsv"), "w")

        isGeneLevel = False
        firstLine = True
        for line in open(path):
            cols = line.rstrip().replace('"', "").split(splitBy)

            if header:
                header -= 1

                if not header:
                    if addExtraInfo:
                        print(
                            "OfficialSymbol\tGeneType\tAllTranscripts\t"
                            + "\t".join(cols),
                            file=outFile,
                        )
                    else:
                        print(
                            "OfficialSymbol\tGeneType\t" + "\t".join(cols), file=outFile
                        )
                continue
            elif firstLine:
                if addExtraInfo:
                    print("OfficialSymbol\tGeneType\tAllTranscripts\t", file=outFile)
                else:
                    print("OfficialSymbol\tGeneType\t", file=outFile)

            tempID = cols[idCol]
            #  Work out if gene-level or transcript-level
            if firstLine:
                firstLine = False

                if splitGeneChar in tempID:
                    if tempID.split(splitGeneChar)[0] in gFI:
                        isGeneLevel = True
                else:
                    if tempID in gFI:
                        isGeneLevel = True

            if isGeneLevel:
                tempTI = []
                geneN = dict()
                geneT = dict()

                if splitGeneChar in tempID:
                    subIDCount = len(tempID.split(splitGeneChar))
                    for subID in tempID.split(splitGeneChar):
                        for tID in sorted(gFI[subID], key=natSort):
                            gene_name, gene_type = gFI[subID][tID]
                            geneN[gene_name] = 1
                            geneT[gene_type] = 1
                            tempTI.append(tID)

                    if len(geneN) > subIDCount:
                        print("Expecting (", subIDCount, ") genes only!", sep="")
                        sys.exit()

                else:
                    for tID in sorted(gFI[tempID], key=natSort):
                        gene_name, gene_type = gFI[tempID][tID]
                        geneN[gene_name] = 1
                        geneT[gene_type] = 1
                        tempTI.append(tID)

                    if len(geneN) > 1:
                        print("Not expecting multiple genes!")
                        sys.exit()

            else:
                tempTI = []
                geneN = dict()
                geneT = dict()

                if splitGeneChar in tempID:
                    subIDCount = len(tempID.split(splitGeneChar))
                    for subID in tempID.split(splitGeneChar):
                        for tID in sorted(tFI[subID], key=natSort):
                            gene_name, gene_type = tFI[subID][tID]
                            geneN[gene_name] = 1
                            geneT[gene_type] = 1
                            tempTI.append(tID)

                    if len(geneN) > subIDCount:
                        print("Expecting (", subIDCount, ") genes only!", sep="")
                        sys.exit()

                else:
                    for tID in sorted(tFI[tempID], key=natSort):
                        gene_name, gene_type = tFI[tempID][tID]
                        geneN[gene_name] = 1
                        geneT[gene_type] = 1
                        tempTI.append(tID)

                    if len(geneN) > 1:
                        print("Not expecting multiple genes!")
                        sys.exit()

            try:
                # Write out results
                geneNF, geneTF = "", ""
                if len(geneN) > 1:
                    geneNF = splitGeneChar.join(geneN.keys())
                else:
                    geneNF = list(geneN.keys())[0]

                if len(geneT) > 1:
                    geneTF = splitGeneChar.join(geneT.keys())
                else:
                    geneTF = list(geneT.keys())[0]

                if addExtraInfo:
                    print(
                        geneNF,
                        geneTF,
                        ";".join(tempTI),
                        "\t".join(cols),
                        sep="\t",
                        file=outFile,
                    )
                else:
                    print(geneNF, geneTF, "\t".join(cols), sep="\t", file=outFile)
            except IndexError:
                print("IndexError for line", line)
                sys.exit()
        outFile.close()


if __name__ == "__main__":
    #  First load GTF file
    geneFI, trancriptFI = loadGTF(
        gtfFile,
        headerL=GTF_HEADER_LINES,
        expectedCols=GTF_EXPECTED_COLS,
        gtfPrintEvery=GTF_PRINT_EVERY,
        gtfGeneIdRe=GTF_GENE_ID_RE,
        gtfGeneTypeRe=GTF_GENE_TYPE_RE,
        gtfTranscriptIdRe=GTF_TRANSCRIPT_ID_RE,
        gtfGeneNameRe=GTF_GENE_NAME_RE,
    )

    print("\nGenerating list of files to process that aren't in TO_EXCLUDE_SUFFIXES")
    filesToProcess = generateFileList(
        workFolder,
        pfg=PROCESS_FILE_GLOB,
        te=TO_EXCLUDE_SUFFIXES,
        overwriteFiles=OVERWRITE_FILES,
    )
    print("\tTotal files to process:", len(filesToProcess))

    processFiles(
        filesToProcess,
        geneFI,
        trancriptFI,
        headerL=PROCESS_FILES_HEADER_LINES,
        idCol=PROCESS_FILES_ID_COL,
        splitBy=PROCESS_FILES_SPLIT_BY,
        addExtraInfo=PROCESS_FILES_EXTRA_INFO,
        splitGeneChar=SPLIT_GENE_CHAR,
    )

    main()
