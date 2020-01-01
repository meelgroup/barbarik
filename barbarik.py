#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2018 Kuldeep Meel
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; version 2
# of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.

from __future__ import print_function
import sys
import os
import math
import random
import argparse
import copy
import tempfile
from random import shuffle

SAMPLER_UNIGEN = 1
SAMPLER_QUICKSAMPLER = 2
SAMPLER_STS = 3
SAMPLER_CUSTOM = 4


# returns List of Independent Variables
def parseIndSupport(indSupportFile):
    f = open(indSupportFile, 'r')
    lines = f.readlines()
    f.close()
    indList = []
    numVars = 0
    for line in lines:
        if (line.startswith('p cnf')):
            fields = line.split()
            numVars = int(fields[2])
        if (line.startswith('c ind')):
            indList.extend(line.strip().replace('c ind', '').replace(' 0', '').strip().replace('v ', '').split())
    if (len(indList) == 0):
        indList = [int(x) for x in range(1, numVars+1)]
    else:
        indList = [int(x) for x in indList]
    return indList


def getSolutionFromUniGen(inputFile, numSolutions):
    inputFileSuffix = inputFile.split('/')[-1][:-4]
    #tempOutputFile = tempfile.gettempdir()+'/'+inputFileSuffix+"_1.txt"
    #cmd = 'python ./samplers/UniGen2.py -runIndex=1 -samples='+str(numSolutions)+' '+inputFile+' '+tempfile.gettempdir()+' > /dev/null 2>&1'
    # print(cmd)

    tempOutputFile = tempfile.gettempdir()+'/'+inputFileSuffix+".txt"
    cmd = './samplers/scalmc -s 1 -v 0 --scalmc 0 --samples '+str(numSolutions)+' --sampleFile '+str(tempOutputFile)
    cmd += ' --multisample 1 '+inputFile+' > /dev/null 2>&1'
    os.system(cmd)
    f = open(tempOutputFile, 'r')
    lines = f.readlines()
    f.close()
    solList = []
    for line in lines:
        if (line.strip().startswith('v')):
            freq = int(line.strip().split(':')[-1])
            for i in range(freq):
                solList.append(line.strip().split(':')[0].replace('v', '').strip())
                if (len(solList) == numSolutions):
                    break
            if (len(solList) == numSolutions):
                break
    solreturnList = solList
    if (len(solList) > numSolutions):
        solreturnList = random.sample(solList, numSolutions)

    cmd = 'rm '+str(tempOutputFile)
    os.system(cmd)
    return solreturnList


# @CHANGE_HERE : please make changes in the below block of code
''' this is the method where you could run your sampler for testing
Arguments : input file, number of solutions to be returned, list of independent variables
output : list of solutions '''


def getSolutionFromCustomSampler(inputFile, numSolutions, indVarList):

    solreturnList = []

    ''' write your code here '''

    return solreturnList


''' END OF BLOCK '''


def getSolutionFromSampler(inputFile, numSolutions, samplerType, indVarList):
    if (samplerType == SAMPLER_UNIGEN):
        return getSolutionFromUniGen(inputFile, numSolutions)
    if (samplerType == SAMPLER_QUICKSAMPLER):
        return getSolutionFromQuickSampler(inputFile, numSolutions, indVarList)
    if (samplerType == SAMPLER_STS):
        return getSolutionFromSTS(inputFile, numSolutions, indVarList)
    if (samplerType == SAMPLER_CUSTOM):
        return getSolutionFromCustomSampler(inputFile, numSolutions, indVarList)
    else:
        print("Error")
        return None


def getSolutionFromSTS(inputFile, numSolutions, indVarList):
    kValue = 50
    samplingRounds = numSolutions/kValue + 1
    inputFileSuffix = inputFile.split('/')[-1][:-4]
    outputFile = tempfile.gettempdir()+'/'+inputFileSuffix+".out"
    cmd = './samplers/STS -k='+str(kValue)+' -nsamples='+str(samplingRounds)+' '+str(inputFile)+' > '+str(outputFile)
    os.system(cmd)
    f = open(outputFile, 'r')
    lines = f.readlines()
    f.close()
    solList = []
    shouldStart = False
    baseList = {}
    for j in range(len(lines)):
        if(lines[j].strip() == 'Outputting samples:' or lines[j].strip() == 'start'):
            shouldStart = True
            continue
        if (lines[j].strip().startswith('Log') or lines[j].strip() == 'end'):
            shouldStart = False
        if (shouldStart):
            i = 0

            if (lines[j].strip() not in baseList):
                baseList[lines[j].strip()] = 1
            else:
                baseList[lines[j].strip()] += 1
            sol = ''

            for x in list(lines[j].strip()):
                if (x == '0'):
                    sol += ' -'+str(indVarList[i])
                else:
                    sol += ' '+str(indVarList[i])
                i += 1
            solList.append(sol)
            if (len(solList) == numSolutions):
                break
    if (len(solList) != numSolutions):
        print(len(solList))
        print("STS Did not find required number of solutions")
        exit(1)
    cmd = 'rm '+outputFile
    os.system(cmd)
    return solList


def getSolutionFromQuickSampler(inputFile, numSolutions, indVarList):
    cmd = "./samplers/quicksampler -n "+str(numSolutions*5)+' '+str(inputFile)+' > /dev/null 2>&1'
    os.system(cmd)
    cmd = "./samplers/z3 "+str(inputFile)+' > /dev/null 2>&1'
    os.system(cmd)
    if (numSolutions > 1):
        i = 0

    f = open(inputFile+'.samples', 'r')
    lines = f.readlines()
    f.close()
    f = open(inputFile+'.samples.valid', 'r')
    validLines = f.readlines()
    f.close()
    solList = []
    for j in range(len(lines)):
        if (validLines[j].strip() == '0'):
            continue
        fields = lines[j].strip().split(':')
        sol = ''
        i = 0
        for x in list(fields[1].strip()):
            if (x == '0'):
                sol += ' -'+str(indVarList[i])
            else:
                sol += ' '+str(indVarList[i])
            i += 1
        solList.append(sol)
        if (len(solList) == numSolutions):
            break
    cmd = 'rm '+inputFile+'.samples'
    os.system(cmd)
    cmd = 'rm '+inputFile+'.samples.valid'
    os.system(cmd)
    if (len(solList) != numSolutions):
        print("Did not find required number of solutions")
        exit(1)
    return solList


def getSolutionFromMUSE(inputFile, numSolutions):
    inputFileSuffix = inputFile.split('/')[-1][:-4]
    tempOutputFile = tempfile.gettempdir()+'/'+inputFileSuffix+".out"
    cmd = './spur -q -s '+str(numSolutions)+' -out '+str(tempOutputFile)+' -cnf '+str(inputFile)
    os.system(cmd)
    f = open(tempOutputFile, 'r')
    lines = f.readlines()
    f.close()
    solList = []
    startParse = False
    for line in lines:
        if (line.startswith('#START_SAMPLES')):
            startParse = True
            continue
        if (not(startParse)):
            continue
        if (line.startswith('#END_SAMPLES')):
            startParse = False
            continue
        fields = line.strip().split(',')
        solCount = int(fields[0])
        sol = ' '
        i = 1
        for x in list(fields[1]):
            if (x == '0'):
                sol += ' -'+str(i)
            else:
                sol += ' '+str(i)
            i += 1
        for i in range(solCount):
            solList.append(sol)
    cmd = 'rm '+tempOutputFile
    os.system(cmd)
    return solList


def getSolutionFromUniform(inputFile, numSolutions):
    return getSolutionFromMUSE(inputFile, numSolutions)


# returns list of solList,newVarList,oldVarList
def findWeightsForVariables(sampleSol, unifSol, numSolutions):
    countList = [5, 5, 5]
    newVarList = [4, 4, 4]
    sampleSol[0] = sampleSol[0].strip()
    unifSol[0] = unifSol[0].strip()
    if (sampleSol[0].endswith(' 0')):
        sampleSol[0] = sampleSol[0][:-2]
    if (unifSol[0].endswith(' 0')):
        unifSol[0] = unifSol[0][:-2]
    lenSol = len(sampleSol[0].split())
    logSol = math.log(5, 2)*3+1
    for i in range(min(int(math.log(numSolutions, 2))+4, lenSol-3, 5)):
        countNum = 31  # random.randint(1,64)
        countList.append(countNum)
        newVarList.append(6)
        logSol += math.log(countNum, 2)
    rExtList = []
    oldVarList = []
    sampleVarList = random.sample(sampleSol[0].split(), len(countList))
    unifVarList = []
    unifSolMap = unifSol[0].split()
    for i in sampleVarList:
        unifVarList.append(unifSolMap[abs(int(i))-1])
    oldVarList.append(sampleVarList)
    oldVarList.append(unifVarList)
    rExtList.append(countList)
    rExtList.append(newVarList)
    rExtList.append(oldVarList)
    return rExtList


def pushVar(variable, cnfClauses):
    cnfLen = len(cnfClauses)
    for i in range(cnfLen):
        cnfClauses[i].append(variable)
    return cnfClauses


def getCNF(variable, binStr, sign, origTotalVars):
    cnfClauses = []
    binLen = len(binStr)
    if (sign == False):
        cnfClauses.append([-(binLen+1+origTotalVars)])
    else:
        cnfClauses.append([binLen+1+origTotalVars])
    for i in range(binLen):
        newVar = int(binLen-i+origTotalVars)
        if (sign == False):
            newVar = -1*(binLen-i+origTotalVars)
        if (binStr[binLen-i-1] == '0'):
            cnfClauses.append([newVar])
        else:
            cnfClauses = pushVar(newVar, cnfClauses)
    pushVar(variable, cnfClauses)
    return cnfClauses


def constructChainFormula(originalVar, solCount, newVars, origTotalVars, invert):
    writeLines = ''
    binStr = str(bin(int(solCount)))[2:-1]
    binLen = len(binStr)
    for i in range(newVars-binLen-1):
        binStr = '0'+binStr

    firstCNFClauses = getCNF(-int(originalVar), binStr, invert, origTotalVars)
    addedClauseNum = 0
    for i in range(len(firstCNFClauses)):
        addedClauseNum += 1
        for j in range(len(firstCNFClauses[i])):
            writeLines += str(firstCNFClauses[i][j])+' '
        writeLines += '0\n'
    CNFClauses = []
    for i in range(len(CNFClauses)):
        if (CNFClauses[i] in firstCNFClauses):
            continue
        addedClauseNum += 1
        for j in range(len(CNFClauses[i])):
            writeLines += str(CNFClauses[i][j])+' '
        writeLines += '0\n'
    return (writeLines, addedClauseNum)


# @returns whether new file was created and the list of independent variables
def constructNewFile(inputFile, tempFile, sampleSol, unifSol, rExtList, indVarList):
    sampleMap = {}
    unifMap = {}
    diffIndex = -1
    for i in sampleSol.strip().split():
        if (not(abs(int(i)) in indVarList)):
            continue
        if (int(i) != 0):
            sampleMap[abs(int(i))] = int(int(i)/abs(int(i)))
    for j in unifSol.strip().split():
        if (int(j) != 0):
            if (not(abs(int(j)) in indVarList)):
                continue

            if (sampleMap[abs(int(j))] != int(j)/abs(int(j))):
                diffIndex = abs(int(j))
            unifMap[abs(int(j))] = int(int(j)/abs(int(j)))

    if (diffIndex == -1):
        return False, -1, -1
    solClause = ''
    f = open(inputFile, 'r')
    lines = f.readlines()
    f.close()
    countList = rExtList[0]
    newVarList = rExtList[1]
    sumNewVar = int(sum(newVarList))
    oldClauseStr = ''
    for line in lines:
        if (line.strip().startswith('p cnf')):
            numVar = int(line.strip().split()[2])
            numClaus = int(line.strip().split()[3])
        else:
            if (not(line.strip().startswith('c'))):
                fields = line.strip().split()
                for x in list(fields):
                    if (int(x) == 0):
                        continue
                    sign = int(int(x)/abs(int(x)))

                    oldClauseStr += str(sign*(abs(int(x))+sumNewVar))+' '
                oldClauseStr += ' 0\n'
    origNumClause = numClaus
    # Adding constraints to ensure only two clauses
    indStr = 'c ind '
    indLen = 0
    indIter = 1
    for i in indVarList:
        if (int(i) != diffIndex):
            numClaus += 2
            solClause += str(-(diffIndex+sumNewVar)*sampleMap[diffIndex])+' '+str(sampleMap[int(i)]*int(i+sumNewVar))+' 0\n'
            solClause += str(-(diffIndex+sumNewVar)*unifMap[diffIndex])+' '+str(unifMap[int(i)]*int(i+sumNewVar))+' 0\n'

    invert = True
    seenVars = []
    for oldVarList in rExtList[2]:
        currentNumVar = 0
        for i in range(len(oldVarList)):
            addedClause = ''
            addedClauseNum = 0
            if (not(int(oldVarList[i]) in seenVars)):
                sign = int(oldVarList[i])/abs(int(oldVarList[i]))
                (addedClause, addedClauseNum) = constructChainFormula(sign*(abs(int(oldVarList[i]))+sumNewVar),
                                                                      int(countList[i]), int(newVarList[i]), currentNumVar, invert)
            seenVars.append(int(oldVarList[i]))
            currentNumVar += int(newVarList[i])
            numClaus += addedClauseNum
            solClause += addedClause
        invert = not(invert)
    oldIndVarList = [x+sumNewVar for x in indVarList]
    tempIndVarList = copy.copy(oldIndVarList)
    for i in range(1, currentNumVar+1):
        if (indIter % 10 == 0):
            indStr += ' 0\nc ind '
        indStr += str(i)+' '
        indIter += 1
        tempIndVarList.append(i)
    for i in oldIndVarList:
        if (indIter % 10 == 0):
            indStr += ' 0\nc ind '
        indStr += str(i)+' '
        indIter += 1

    indStr += ' 0\n'
    currentNumVar += numVar

    headStr = 'p cnf '+str(currentNumVar)+' '+str(numClaus)+'\n'
    writeStr = headStr + indStr
    writeStr += solClause
    writeStr += oldClauseStr

    f = open(tempFile, 'w')
    f.write(writeStr)
    f.close()
    return True, tempIndVarList, oldIndVarList


# Returns 1 if uniform and 0 otherwise
def testUniformity(solList, indVarList, numSolutions, loThresh, hiThresh, outputFile):
    solMap = {}
    baseMap = {}
    for sol in solList:
        solution = ''
        solFields = sol.split()
        for entry in solFields:
            if ((abs(int(entry))) in indVarList):
                solution += entry+' '
        #solution = solution[:-1]
        if (solution in solMap.keys()):
            solMap[solution] += 1
        else:
            solMap[solution] = 1
        if (sol not in baseMap.keys()):
            baseMap[sol] = 1
        else:
            baseMap[sol] += 1

    if (not(bool(solMap))):
        print("No Solutions were given to the test")
        exit(1)

    key = next(iter(solMap))

    f = open(outputFile, 'a')
    f.write("baseMap:{4} numSolutions:{3} SolutionsCount:{0} loThresh:{1} hiThresh:{2}\n".format(
        solMap[key], loThresh, hiThresh, numSolutions, len(baseMap.keys())))
    f.close()
    if (solMap[key] >= loThresh and solMap[key] <= hiThresh):
        return True
    else:
        return False


def barbarik():
    parser = argparse.ArgumentParser()
    parser.add_argument('--eta', type=float, help="default = 0.9", default=0.9, dest='eta')
    parser.add_argument('--epsilon', type=float, help="default = 0.6", default=0.6, dest='epsilon')
    parser.add_argument('--delta', type=float, help="default = 0.05", default=0.05, dest='delta')
    parser.add_argument('--sampler', type=int, help=str(SAMPLER_UNIGEN)+" for UniGen;\n" +
                        str(SAMPLER_QUICKSAMPLER)+" for QuickSampler;\n"+str(SAMPLER_STS)+" for STS;\n", default=SAMPLER_STS, dest='sampler')
    parser.add_argument('--reverse', type=int, default=0, help="order to search in", dest='searchOrder')
    parser.add_argument('--minSamples', type=int, default=0, help="min samples", dest='minSamples')
    parser.add_argument('--maxSamples', type=int, default=sys.maxsize, help="max samples", dest='maxSamples')
    parser.add_argument('--seed', type=int, required=True, dest='seed')
    parser.add_argument('--exp', type=int, help="number of experiments", dest='exp', default=1)
    parser.add_argument("input", help="input file")
    parser.add_argument("output", help="output file")

    args = parser.parse_args()
    inputFile = args.input
    inputFileSuffix = inputFile.split('/')[-1][:-4]

    eta = args.eta
    epsilon = args.epsilon
    delta = args.delta
    numExperiments = args.exp
    if (numExperiments == -1):
        numExperiments = sys.maxsize
    samplerType = args.sampler
    searchOrder = args.searchOrder
    outputFile = args.output
    f = open(outputFile, 'w')
    f.close()
    seed = args.seed
    minSamples = args.minSamples
    maxSamples = args.maxSamples
    samplerString = ''
    random.seed(seed)
    if (samplerType == SAMPLER_UNIGEN):
        samplerString = 'UniGen'
    if (samplerType == SAMPLER_QUICKSAMPLER):
        samplerString = 'QuickSampler'
    if (samplerType == SAMPLER_STS):
        samplerString = 'STS'
    if (samplerType == SAMPLER_CUSTOM):
        samplerString = 'CustomSampler'

    indVarList = parseIndSupport(inputFile)
    totalLoops = int(math.ceil(math.log(2.0/(eta+epsilon), 2))+1)
    listforTraversal = range(totalLoops, 0, -1)
    if (searchOrder == 1):
        listforTraversal = range(1, totalLoops+1, 1)
    for experiment in range(numExperiments):
        breakExperiment = False
        totalSolutionsGenerated = 0
        totalUniformSamples = 0
        thresholdSolutions = 0
        for j in listforTraversal:
            tj = math.ceil(math.pow(2, j)*(epsilon+eta)/((eta-epsilon)**2)*math.log(4.0/(eta+epsilon), 2)*(4*math.e/(math.e-1)*math.log(1.0/delta)))
            beta = (math.pow(2, j-1)+1)*(eta + epsilon)*1.0/(4+(epsilon+eta)*(math.pow(2, j-1) - 1))
            gamma = (beta-epsilon)/4
            constantFactor = math.ceil(1/(9*gamma*gamma))
            boundFactor = math.log((16)*(math.e/(math.e-1))*(1/((eta-epsilon)**2))*math.log(4/(eta+epsilon), 2)*math.log(1/delta), 2)
            f = open(outputFile, 'a')
            f.write("constantFactor:{0} boundFactor:{1} logBoundFactor:{2}\ntj:{3} totalLoops:{4} beta:{5} epsilon:{6}\n".format(
                constantFactor, boundFactor, math.log(boundFactor, 2), tj, totalLoops, beta, epsilon))
            numSolutions = int(math.ceil(constantFactor*boundFactor))
            loThresh = int((numSolutions*1.0/2)*(1-(beta+epsilon)/2))
            hiThresh = int((numSolutions*1.0/2)*(1+(beta+epsilon)/2))

            tempFile = tempfile.gettempdir()+"/"+inputFileSuffix+"_t.cnf"

            f.write("tj:%d numSolutions:%d loThresh:%d hiThresh:%d\n" % (tj, numSolutions, loThresh, hiThresh))
            f.close()
            i = 0
            while(i < int(tj)):
                i += 1
                thresholdSolutions += numSolutions
                if (thresholdSolutions < minSamples):
                    continue
                sampleSol = getSolutionFromSampler(inputFile, 1, samplerType, indVarList)
                unifSol = getSolutionFromUniform(inputFile, 1)
                totalUniformSamples += 1
                totalSolutionsGenerated += 1
                rExtList = findWeightsForVariables(sampleSol, unifSol, numSolutions)
                shakuniMix, tempIndVarList, oldIndVarList = constructNewFile(inputFile, tempFile, sampleSol[0], unifSol[0], rExtList, indVarList)
                if (not(shakuniMix)):
                    i -= 1
                    continue
                solList = getSolutionFromSampler(tempFile, numSolutions, samplerType, tempIndVarList)
                isUniform = testUniformity(solList, oldIndVarList, numSolutions, loThresh, hiThresh, outputFile)
                cmd = 'rm '+tempFile
                os.system(cmd)
                totalSolutionsGenerated += numSolutions
                f = open(outputFile, 'a')
                f.write("sampler:{0} i:{1} isUniform:{2} TotalSolutionsGenerated:{3}\n".format(samplerString, i, isUniform,
                                                                                               totalSolutionsGenerated))
                f.close()
                if (isUniform == False):
                    f = open(outputFile, 'a')
                    f.write("exp:{4} RejectIteration:{0}  Loop:{1} TotalSolutionsGenerated:{2} TotalUniformSamples:{3}\n"
                            .format(i, j, totalSolutionsGenerated, totalUniformSamples, experiment))
                    f.close()
                    breakExperiment = True
                    break
                if (thresholdSolutions > maxSamples):
                    breakExperiment = True
                    break
            if (breakExperiment):
                break
        if(not(breakExperiment)):
            f = open(outputFile, 'a')
            f.write("exp:{2} Accept:1 TotalSolutionsGenerated:{0} TotalUniformSamples:{1}\n".format
                    (totalSolutionsGenerated, totalUniformSamples, experiment))
            f.close()
        breakExperiment = False


if __name__ == "__main__":
    barbarik()
