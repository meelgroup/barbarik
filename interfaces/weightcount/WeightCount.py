# @author Kuldeep S Meel
# Copyright (c) 2016 Kuldeep S Meel
# Please see attached License file for licensing information
# If you use this code for experiments, please cite the following paper:
# "From Weighted to UnWeighted Model Counting", Supratik Chakraborty,
# Dror Fried, Kuldeep S. Meel, Moshe Y. Vardi
# Proc. of IJCAI 2016

import os
import math
import argparse


def pushVar(variable, cnfClauses):
    cnfLen = len(cnfClauses)
    for i in range(cnfLen):
        cnfClauses[i].append(variable)
    return cnfClauses


def getCNF(variable, binStr, sign, origTotalVars):
    cnfClauses = []
    binLen = len(binStr)
    cnfClauses.append([binLen + 1 + origTotalVars])
    for i in range(binLen):
        newVar = binLen - i + origTotalVars
        if sign == False:
            newVar = -1 * (binLen - i + origTotalVars)
        if binStr[binLen - i - 1] == "0":
            cnfClauses.append([newVar])
        else:
            cnfClauses = pushVar(newVar, cnfClauses)
    pushVar(variable, cnfClauses)
    return cnfClauses


def EncodeCNF(
    variable, kWeight, iWeight, origtotalVars, origtotalClaus, independentSet, precision
):
    writeLines = ""
    totalVars = origtotalVars + iWeight
    totalClaus = origtotalClaus
    independentSet[origtotalVars] = 1
    binStr = str(bin(int(kWeight)))[2:-1]
    binLen = len(binStr)
    for i in range(iWeight - binLen - 1):
        binStr = "0" + binStr
    for i in range(iWeight - 1):
        independentSet[origtotalVars + i + 1] = 1
    complementStr = ""
    for i in range(len(binStr)):
        if binStr[i] == "0":
            complementStr += "1"
        else:
            complementStr += "0"
    origCNFClauses = getCNF(-variable, binStr, True, origtotalVars)

    for i in range(len(origCNFClauses)):
        totalClaus += 1
        for j in range(len(origCNFClauses[i])):
            writeLines += str(origCNFClauses[i][j]) + " "
        writeLines += "0\n"
    cnfClauses = getCNF(variable, complementStr, False, origtotalVars)
    for i in range(len(cnfClauses)):
        if cnfClauses[i] in origCNFClauses:
            continue
        totalClaus += 1
        for j in range(len(cnfClauses[i])):
            writeLines += str(cnfClauses[i][j]) + " "
        writeLines += "0\n"
    return writeLines, totalVars, totalClaus, independentSet


def ParseWeights(initWeight, precision):
    if initWeight == 1:
        return 1, 0
    weight = math.ceil(initWeight * pow(2, precision))
    while weight % 2 == 0:
        weight = weight / 2
        precision -= 1
    return weight, precision


#  The code is straightforward chain formula implementation
#
def Transform(inputFile, outputFile, precision):
    f = open(inputFile, "r")
    lines = f.readlines()
    f.close()
    writeLines = ""
    weightLine = ""
    independentSet = {}
    totalVars = 0
    totalClaus = 0
    indWeight = {}
    origWeight = {}
    minWeight = 0
    indepvars = ""
    hasInd = 0
    for line in lines:
        if line.strip()[:2] == "p ":
            fields = line.strip().split()
            totalVars = int(fields[2])
            totalClaus = int(fields[3])
            continue
        if (line.strip()[:5] == "c ind"):
            hasInd = 1
            indepvars = line.strip()[5:-2].split()
            for var in indepvars:
                independentSet[int(var)] = 1
            continue
        if (line.strip()[0] == "c"):
            continue
        if (
            line.strip()[0].isdigit()
            or line.strip()[0] == "-"
        ):
            writeLines += str(line)

    if hasInd == 0:
        for i in range(1, totalVars + 1):
            independentSet[i] = 1
    equalWeightVariables = 0
    for line in lines:
        if line.strip()[:2] == "w ":
            fields = line.strip()[2:].split()
            variable = int(fields[0])
            cnfWeight = float(fields[1])
            if cnfWeight < 1:
                minWeight += math.log(min(cnfWeight, 1 - cnfWeight), 2)
            origWeight[variable] = cnfWeight
            kWeight, iWeight = ParseWeights(cnfWeight, precision)
            if not ((iWeight == 0 and kWeight == 1) or (cnfWeight == 0.5)):
                weightLine, totalVars, totalClaus, independentSet = EncodeCNF(
                    variable,
                    kWeight,
                    iWeight,
                    totalVars,
                    totalClaus,
                    independentSet,
                    precision,
                )
            else:
                if iWeight == 0:
                    if kWeight == 1:
                        totalClaus += 1
                        weightLine += str(variable) + " 0\n"
                    if kWeight == 0:
                        totalClaus += 1
                        weightLine += str(-variable) + " 0\n"
                if cnfWeight == 0.5:
                    equalWeightVariables += 1
                indWeight[variable] = 1
            writeLines += weightLine
    indWriteStr = "c ind "
    f = open(outputFile, "w")
    independentSet = list(independentSet)
    independentSet.sort()
    f.write("p cnf " + str(totalVars) + " " + str(totalClaus) + " \n")
    indWriteStr += str(' '.join([str(i) for i in independentSet])) + " 0\n"
    f.write(indWriteStr)
    f.write(writeLines)
    f.close()
    return independentSet


def ensureDirectory(path):
    d = os.path.dirname(path)
    if not os.path.exists(d):
        os.makedirs(d)
    return


####################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-precision", help="Precision (value of m)", type=int, default=4
    )
    parser.add_argument("inputFile", help="input File(in Weighted CNF format)")
    args = parser.parse_args()
    precision = args.precision
    inputFile = args.inputFile
    tmpDir = os.getenv("TMPDIR", "temp")

    # initialization
    ensureDirectory(tmpDir + "/")

    initialFileSuffix = inputFile.split("/")[-1][:-4]
    outCNFFile = tmpDir + "/" + str(initialFileSuffix) + "u.cnf"

    Transform(inputFile, outCNFFile, precision)
