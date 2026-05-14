#!/usr/bin/env python
"""
Use http://mygene.info/ to annotate gene symbols with the latest annotations.

See http://docs.mygene.info/projects/mygene-py/en/latest/ for readme for the python mygene package
See http://mygene.info/ for the actual API documentation

########################################################################################
# RE-ANNOTATION PIPELINE
# This also feeds into a re-annotation pipeline used for Daisy-Net/PGCNA updating:

1. Run getGeneProbePairs.py
	To find all the gene/probe pairs that were used across the data-sets used to make the nework
		NOTE: will need to manually map probes to useful ID's for data-sets where the probes aren't useful as identifiers of genes (e.g. GSE32918 or GSE22895 etc)

2. Run annotateGeneUsingMyGeneInfo.py
	a> On the non-redundant list of genes that make up the network (e.g. CLL_Genes.txt)
	b> The output from getGeneProbePairs.py (after manually updating any non-informative IDs; e.g. for GSE<number>)

3. Run mergeAnnotationsAndMakeMappingFile.py
	Input two files
	a> Re-annotated non-redundant list of genes (this is used as backbone for script - should the probe annotations fail)
	b> Re-annotated probes.  Annotations from probes are priority
"""
import sys
import os
from biothings_client import get_client  # Moved over to biothings_client now
import time
import natsort as ns
from natsort import natsort_keygen
from collections import defaultdict



workFolder = r"PATH/ForPGCNA"
# Using .txt file generated from addGeneymbolsFromGeneIDwithGTFfile_AF.py. 
fileToAnnotate = r"RNAseq_CLL_RSEM_Gene_count_VST_wGeneSymbols_PGCNA.txt"
##----------------------------------------CONSTANTS----------------------------------------##
SPECIES = "human"
# SPECIES = "mouse"
ANNOTATION_FILE_HEADER_LINES = 1
GENE_COL = 0  # 0, 2

#  Decide what to do with multiple identifiers
GENE_SPLIT_BY = "#"
KEEP_BEST = 2  #  For multiple genes split by GENE_SPLIT_BY only keep the top KEEP_BEST in output.  If set to zero will keep all!
KEEP_BEST_KEEP_TIES = (
    True  #  If True will keep all genes that tie for a score even if that pushes beyond the KEEP_BEST limit
)
KEEP_BEST_ABS_MAX = (
    4  #  Make sure this is > KEEP_BEST.  If KEEP_BEST_KEEP_TIES is True then this is the maximum limit of tied genes
)

ENTREZ_PREFIX = "LOC"  # Need to know what to remove so that we can search using ENTREZ_ID e.g.
ENSEMBL_SPLIT = True  # If true will split on '.' and only retain prefix e.g ENSG00000000003.14 --> ENSG00000000003  # Important if using ENSEMBL gene codes
ENSEMBL_DONT_WARN = (
    True  # If true won't warn that identifiers aren't of correct format.  Useful if a mixed Human/EBV for example
)

RETAIN_INPUT_ORDER = (
    True  # If True will retain the input order.  Set ALLOW_REDUNDANCY should you want to allow for redundancy
)
ALLOW_REDUNDANCY = True  # If True will allow redundancy else will cull it.
ADD_ORIGINAL_COLUMNS_BACK = True  # Add original columns back.
REMOVE_QUOTATION_MARKS = True  # If true will remove all '"' from file
SCOPE_STR = "symbol,alias,entrezgene,reporter,refseq,ensembl.gene,ensembl.transcript,ensembl.protein,accession,unigene"  # e.g. "symbol,alias,entrezgene,reporter,refseq,ensembl.gene,accession"  see http://docs.mygene.info/en/latest/doc/query_service.html?highlight=scope#available-fields
# SCOPE_STR = "accession,refseq"  # e.g. "symbol,alias,entrezgene,reporter,refseq,ensembl.gene,accession"  see http://docs.mygene.info/en/latest/doc/query_service.html?highlight=scope#available-fields
# SCOPE_STR = "entrezgene"  # e.g. "symbol,alias,entrezgene,reporter,refseq,ensembl.gene,accession"  see http://docs.mygene.info/en/latest/doc/query_service.html?highlight=scope#available-fields


natSort = natsort_keygen(alg=ns.LOCALE)
###########################################################################################
##----------------------------------------Methods----------------------------------------##
###########################################################################################


def getGeneList(
    workFolder,
    fileToAnnotate,
    geneSplitBy="#",
    entrezPrefix="LOC",
    splitBy="\t",
    headerL=0,
    geneCol=0,
    retainInputOrder=False,
    ensemblSplit=False,
    ensemblDontWarn=False,
):

    queryInfo = defaultdict(dict)
    uniqueGenes = {}
    originalData = {}

    header = headerL

    for i, line in enumerate(open(os.path.join(workFolder, fileToAnnotate))):
        if REMOVE_QUOTATION_MARKS:
            cols = line.rstrip().replace('"', "").split(splitBy)
        else:
            cols = line.rstrip().split(splitBy)

        if i == 0:
            if not header:
                originalData["header"] = [""] * (len(cols) - 1)  # Deal with when no header exists
            else:
                originalData["header"] = cols

        if header:
            header -= 1
            continue

        originalData[str(i)] = cols

        gene = cols[geneCol]
        if gene == "":
            continue

        if ":" in gene:
            print("Issue with gene:", gene, "unexpected ':' in gene name, needs to be resolved!")
            sys.exit()

        if retainInputOrder:
            geneOrderInf = str(i) + "\t" + gene
        else:
            geneOrderInf = gene + "\t" + str(i)  # Assume we want to retain redundancy until later

        ########################################################################################

        if geneSplitBy in gene:
            for g in gene.split(geneSplitBy):
                if g == "":
                    continue

                if ensemblSplit:
                    if "." in g:
                        g = g.split(".")[0]

                if entrezPrefix in g:
                    #  Trim prefix
                    tgene = g[len(entrezPrefix) :]
                    queryInfo[geneOrderInf][g + ":" + tgene] = 1
                    uniqueGenes[tgene] = 1
                else:
                    queryInfo[geneOrderInf][g + ":" + g] = 1
                    uniqueGenes[g] = 1
        else:

            if ensemblSplit:
                if "." in gene:
                    if not ensemblDontWarn:
                        if "ENS" != gene[:3]:
                            print("\nUsing ENSEMBLE_SPLIT but IDs do not look like Ensemble IDS:", gene)
                            sys.exit()
                    gene = gene.split(".")[0]

            if entrezPrefix in gene:
                #  Trim prefix
                tgene = gene[len(entrezPrefix) :]
                queryInfo[geneOrderInf][gene + ":" + tgene] = 1
                uniqueGenes[tgene] = 1
            else:
                queryInfo[geneOrderInf][gene + ":" + gene] = 1
                uniqueGenes[gene] = 1

    print("\tTotal unique genes to process:", len(uniqueGenes), "\n\n")

    return queryInfo, uniqueGenes, originalData


def reAnnotateGenes(
    wF,
    qI,
    uG,
    oD,
    species="human",
    geneSplitBy="#",
    retainInputOrder=False,
    allowRedundancy=False,
    addOriginalColumnsBack=False,
    keepBest=1,
    keepBestTies=False,
    keepBestAbsMax=10,
):

    #  Create instance of mygene class
    # mg = mygene.MyGeneInfo()  # Old mygene version
    gene_client = get_client("gene")

    # results = mg.querymany(uG.keys(),species=species,scopes=SCOPE_STR)  # OLD mygene version
    results = gene_client.querymany(uG.keys(), species=species, scopes=SCOPE_STR)

    #  Organise by Query so we can match against qI more easily
    # Retain best match only
    resultByQuery = {}
    for r in results:
        q = r["query"]

        if q in resultByQuery:
            previousScore = 0
            if "_score" in resultByQuery[q]:
                previousScore = resultByQuery[q]["_score"]

            newScore = 0
            if "_score" in r:
                newScore = r["_score"]

            if newScore > previousScore:
                resultByQuery[q] = r

        else:
            resultByQuery[q] = r

    #  Process per original
    print("\n\nProcess results from mygene.info")

    finalResults = defaultdict(dict)

    seenResults = {}

    for geneInfM in sorted(qI, key=natSort):

        if retainInputOrder:
            pos, geneInf = geneInfM.split("\t")
        else:
            geneInf, pos = geneInfM.split("\t")

        ##----------------------------------------Deal with redundancy----------------------------------------##
        if geneInf in seenResults:
            if not allowRedundancy:
                continue
        seenResults[geneInf] = 1
        ########################################################################################
        ########################################################################################

        if len(qI[geneInfM]) == 1:
            try:
                origG, q = list(qI[geneInfM].keys())[0].split(":")
            except ValueError:
                print("Issue with :", geneInfM, ">>", qI[geneInfM].keys()[0].split(":"))
                sys.exit()

            if q in resultByQuery:
                if "symbol" in resultByQuery[q]:
                    resultStr = resultByQuery[q]["symbol"] + "\t" + str(resultByQuery[q]["_score"])
                    finalResults[geneInfM][resultStr] = 1
                else:
                    finalResults[geneInfM][geneInf + "\t0"] = 1
            else:
                print("Result for (", q, ") is missing from returned results!")
                sys.exit()

        else:
            scores = []
            genes = []
            resultStr = []

            for i, qInf in enumerate(sorted(qI[geneInfM], key=natSort)):
                origG, q = qInf.split(":")
                if q in resultByQuery:
                    if "symbol" in resultByQuery[q]:
                        gene = resultByQuery[q]["symbol"]
                        score = resultByQuery[q]["_score"]

                        genes.append(gene)
                        scores.append(score)
                    else:
                        genes.append(origG)
                        scores.append(0)
                else:
                    print("Result for (", q, ") is missing from returned results!")
                    sys.exit()

            if len(genes) > 1:
                #  Only retain best keepBest number of genes
                if keepBest > 0:
                    zStr = zip(scores, genes)
                    genes = []
                    scores = []
                    seen = {}

                    #  Sort by score
                    zStrSorted = sorted(zStr, reverse=True)
                    # Store best gene/score
                    best = zStrSorted.pop(0)
                    genes.append(best[1])
                    scores.append(best[0])
                    seen[best[1]] = 1

                    gotCount = 1
                    lastScore = best[0]

                    # compare rest
                    for bI in zStrSorted:
                        cScore, cGene = bI
                        if cScore > 0:  # No point adding extra non informative genes back
                            if gotCount >= keepBest:
                                if keepBestTies & (cScore == lastScore):
                                    if gotCount >= keepBestAbsMax:
                                        break
                                    else:
                                        if not cGene in seen:
                                            genes.append(cGene)
                                            scores.append(cScore)
                                            lastScore = cScore
                                            gotCount += 1
                                            seen[cGene] = 1
                                else:
                                    break
                            else:
                                if not cGene in seen:
                                    genes.append(cGene)
                                    scores.append(cScore)
                                    gotCount += 1
                                    lastScore = cScore
                                    seen[cGene] = 1

            zStr = zip(genes, map(str, scores))

            for gI in zStr:
                resultStr.append("\t".join(gI))

            #  Remove redundancy in gene string
            seenG = {}
            finalGenes = []
            for gene in genes:
                if gene not in seenG:
                    finalGenes.append(gene)

                seenG[gene] = 1

            finalResults[geneInfM][
                geneSplitBy.join(finalGenes) + "\t" + str(sum(scores)) + "\t" + "\t".join(resultStr)
            ] = 1

    dateStamp = time.strftime("%Y%m%d%I%M%S")

    resultsFile = open(os.path.join(wF, dateStamp + "_" + fileToAnnotate[:-4] + "_mygeneAnnotated.txt"), "w")

    if addOriginalColumnsBack:
        print("\t".join(oD["header"]), "OriginalAnnotation\tOfficialSymbol\tScore", sep="\t", file=resultsFile)
    else:
        print("OriginalAnnotation\tOfficialSymbol\tScore", file=resultsFile)

    for gM in sorted(finalResults, key=natSort):
        if retainInputOrder:
            pos, g = gM.split("\t")
        else:
            g, pos = gM.split("\t")

        if len(finalResults[gM]) > 1:
            print("Multiple results for (", g, ") !", finalResults[gM], sep="")
        else:
            for resultInf in finalResults[gM]:
                if addOriginalColumnsBack:
                    print("\t".join(oD[pos]), g, resultInf, sep="\t", file=resultsFile)
                else:
                    print(g, resultInf, sep="\t", file=resultsFile)
    resultsFile.close()


def runMain(finishT="Finished!"):
    # global # variable.  Declare each global variable and pass new values/defaults
    pass
    main(finishT=finishT)


########################################################################################
##----------------------------------------MAIN----------------------------------------##
########################################################################################


def main(finishT="Finished!"):

    print("Parse list of genes to process")
    queryInfo, uniqueGenes, originalData = getGeneList(
        workFolder,
        fileToAnnotate,
        geneSplitBy=GENE_SPLIT_BY,
        entrezPrefix=ENTREZ_PREFIX,
        headerL=ANNOTATION_FILE_HEADER_LINES,
        geneCol=GENE_COL,
        retainInputOrder=RETAIN_INPUT_ORDER,
        ensemblSplit=ENSEMBL_SPLIT,
        ensemblDontWarn=ENSEMBL_DONT_WARN,
    )

    print("Using mygene.info to get information")
    reAnnotateGenes(
        workFolder,
        queryInfo,
        uniqueGenes,
        originalData,
        species=SPECIES,
        geneSplitBy=GENE_SPLIT_BY,
        retainInputOrder=RETAIN_INPUT_ORDER,
        allowRedundancy=ALLOW_REDUNDANCY,
        addOriginalColumnsBack=ADD_ORIGINAL_COLUMNS_BACK,
        keepBest=KEEP_BEST,
        keepBestTies=KEEP_BEST_KEEP_TIES,
        keepBestAbsMax=KEEP_BEST_ABS_MAX,
    )

    print("\n", finishT, sep="")


##----------------------------------------------------------------------------------------##
##----------------------------------------__main__----------------------------------------##
if __name__ == "__main__":
    main()