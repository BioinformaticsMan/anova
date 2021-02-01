#!/usr/bin/python
import sys, os
import numpy as np
from scipy import stats

inputFile = sys.argv[1]
prefix    = sys.argv[2]

## 01
# define the list for anova analysis based on tissue sorce
# the first list for index and the second list for expression data
geneExprDict = {"group1":[[],[]],
                "group2":[[],[]],
                "group3":[[],[]],
                "group4":[[],[]],
                "DEGs":[],
                "rmDEGs":[]}

## match each column with group 
with open (inputFile, "r") as f_in:
    lines = f_in.readlines()
    groupInfo = lines[-1] 
    count = 0
    group_lst = (lines[-1].strip("\n")).split("\t")
    for source in group_lst:
        if (source == "\"0\""):
            geneExprDict["group1"][0].append(count)
        elif (source == "\"III-IV\""):
            geneExprDict["group2"][0].append(count)
        elif (source == "\"I-II\""):
            geneExprDict["group3"][0].append(count)
        elif (source == "\"V-VI\""):
            geneExprDict["group4"][0].append(count)
        count = count + 1

## anova analysis for each gene       
geneExprDict["anova"] = {}
DEGs = {}
allGenes = {}


## line 0 title
## last line group 
## expression data line 1 to -2
for line in lines[1:-1]:
    exp = line.strip("\n").replace("\"","").split("\t")
    geneSymbol = exp[0]
    geneExprDict["anova"]["group1"] = []
    geneExprDict["anova"]["group2"] = []
    geneExprDict["anova"]["group3"] = []
    geneExprDict["anova"]["group4"] = []
    for index in (geneExprDict["group1"][0]):geneExprDict["anova"]["group1"].append(float(exp[index]))
    for index in (geneExprDict["group2"][0]):geneExprDict["anova"]["group2"].append(float(exp[index]))
    for index in (geneExprDict["group3"][0]):geneExprDict["anova"]["group3"].append(float(exp[index]))
    for index in (geneExprDict["group4"][0]):geneExprDict["anova"]["group4"].append(float(exp[index]))
    f,p = stats.f_oneway(geneExprDict["anova"]["group1"],
                        geneExprDict["anova"]["group2"],
                        geneExprDict["anova"]["group3"],
                        geneExprDict["anova"]["group4"])

    allGenes[geneSymbol] = [f,p]
    if p < 0.05:DEGs[geneSymbol] = [f,p]
#print (geneExprDict)
def outputFile(fileName,headline,geneExprDict):
    with open(fileName,"w") as f_out:
        f_out.write(headline)
        for key in geneExprDict:
            content = key + "\t"+"\t".join(str(i) for i in geneExprDict[key])+"\n"
            f_out.write(content)
#print (DEGs)
#print (allGenes)
### output the DEGs file
DEGsFileName = prefix+"_DEGs_"+"anova.txt"
allGeneFileName = prefix+"_allGenes_"+"anova.txt"
headline = "geneSymbol"+"\t"+"F_value" +"\t"+"P_value"+"\n"
outputFile(DEGsFileName,headline,DEGs)
outputFile(allGeneFileName,headline,allGenes)
print ("DEGS number is:     ",len(DEGs))
