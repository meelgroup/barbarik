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
import tempfile

from interfaces.cnf import *


class ChainFormulaSetup:
    def __init__(self, countList, newVarList, indicatorLits):
        self.countList = countList
        self.newVarList = newVarList
        self.indicatorLits = indicatorLits


# returns List of Independent Variables
def parseIndSupport(indSupportFile):
    with open(indSupportFile, 'r') as f:
        lines = f.readlines()

    indList = []
    numVars = 0
    for line in lines:
        if line.startswith('p cnf'):
            fields = line.split()
            numVars = int(fields[2])

        if line.startswith('c ind'):
            line = line.strip().replace('c ind', '').replace(' 0', '').strip().replace('v ', '')
            indList.extend(line.split())

    if len(indList) == 0:
        indList = [int(x) for x in range(1, numVars+1)]
    else:
        indList = [int(x) for x in indList]
    return indList


def chainFormulaSetup(sampleSol, unifSol, numSolutions):
    # number of solutions for each: k1, k2, k3
    # TODO rename to chainSolutions
    countList = [5, 5, 5]

    # chain formula number of variables for each
    # TODO rename to chainVars
    newVarList = [4, 4, 4]

    ##########
    # clean up the solutions
    ##########
    sampleSol = sampleSol[0].strip()
    if sampleSol.endswith(' 0'):
        sampleSol = sampleSol[:-2]
    unifSol = unifSol[0].strip()
    if unifSol.endswith(' 0'):
        unifSol = unifSol[:-2]

    # adding more chain formulas (at most 8 in total: 3 + 5)
    # these chain formulas will have 31 solutions over 6 variables
    lenSol = len(sampleSol.split())
    for i in range(min(int(math.log(numSolutions, 2))+4, lenSol-3, 5)):
        countList.append(31)
        newVarList.append(6)
    assert len(countList) == len(newVarList)

    # picking selector literals, i.e. k1, k2, k3, randomly
    sampleLitList = random.sample(sampleSol.split(), len(countList))
    unifLitList = []
    unifSolMap = unifSol.split()
    for lit in sampleLitList:
        unifLitList.append(unifSolMap[abs(int(lit))-1])

    assert len(unifLitList) == len(sampleLitList)
    for a, b in zip(unifLitList, sampleLitList):
        assert abs(int(a)) == abs(int(b))

    indicatorLits = []
    indicatorLits.append(sampleLitList)
    indicatorLits.append(unifLitList)

    #print("countList:", countList)
    #print("newVarList:", newVarList)
    #print("indicatorLits:", indicatorLits)
    return ChainFormulaSetup(countList, newVarList, indicatorLits)


def pushVar(variable, cnfClauses):
    cnfLen = len(cnfClauses)
    for i in range(cnfLen):
        cnfClauses[i].append(variable)
    return cnfClauses


def getCNF(variable, binStr, sign, origTotalVars):
    cnfClauses = []
    binLen = len(binStr)
    if sign is False:
        cnfClauses.append([-(binLen+1+origTotalVars)])
    else:
        cnfClauses.append([binLen+1+origTotalVars])

    for i in range(binLen):
        newVar = int(binLen-i+origTotalVars)
        if sign is False:
            newVar = -1*(binLen-i+origTotalVars)

        if (binStr[binLen-i-1] == '0'):
            cnfClauses.append([newVar])
        else:
            cnfClauses = pushVar(newVar, cnfClauses)
    pushVar(variable, cnfClauses)
    return cnfClauses


def constructChainFormula(originalVar, solCount, newVar, origTotalVars, invert):
    assert type(solCount) == int

    binStr = str(bin(int(solCount)))[2:-1]
    binLen = len(binStr)
    for _ in range(newVar-binLen-1):
        binStr = '0'+binStr

    firstCNFClauses = getCNF(-int(originalVar), binStr, invert, origTotalVars)
    addedClauseNum = 0
    writeLines = ''
    for cl in firstCNFClauses:
        addedClauseNum += 1
        for lit in cl:
            writeLines += "%d " % lit
        writeLines += '0\n'

    return writeLines, addedClauseNum


# returns whether new file was created and the list of TMP+OLD independent variables
def constructNewCNF(inputFile, tempFile, sampleSol, unifSol, chainFormulaConf, indVarList):
    # which variables are in pos/neg value in the sample
    sampleVal = {}
    for i in sampleSol.strip().split():
        i = int(i)
        if i != 0:
            if abs(i) not in indVarList:
                continue

            sampleVal[abs(i)] = int(i/abs(i))

    # which variables are in pos/neg value in the uniform sample
    unifVal = {}
    diffIndex = -1
    for j in unifSol.strip().split():
        j = int(j)
        if j != 0:
            if abs(j) not in indVarList:
                continue

            unifVal[abs(j)] = int(j/abs(j))

            if sampleVal[abs(j)] != unifVal[abs(j)]:
                diffIndex = abs(j)

    # the two solutions are the same
    # can't do anything, let's do another experiment
    if diffIndex == -1:
        return True, None, None

    with open(inputFile, 'r') as f:
        lines = f.readlines()

    # variables must be shifted by sumNewVar
    sumNewVar = sum(chainFormulaConf.newVarList)

    # emit the original CNF, but with shifted variables
    shiftedCNFStr = ''
    for line in lines:
        line = line.strip()
        if line.startswith('p cnf'):
            numVar = int(line.split()[2])
            numCls = int(line.split()[3])
            continue

        if line.startswith('c'):
            # comment
            continue

        for x in line.split():
            x = int(x)
            if x == 0:
                continue
            sign = int(x/abs(x))
            shiftedCNFStr += "%d " % (sign*(abs(x)+sumNewVar))
        shiftedCNFStr += ' 0\n'
    del i

    # Fixing the solution based on splittingVar
    # X = sigma1 OR X = singma2
    # All variables are set except for the index where they last differ
    solClause = ''
    splittingVar = diffIndex+sumNewVar
    for var in indVarList:
        if var != diffIndex:
            numCls += 2
            solClause += "%d " % (-splittingVar*sampleVal[diffIndex])
            solClause += "%d 0\n" % (sampleVal[var]*(var+sumNewVar))

            solClause += "%d " % (-splittingVar*unifVal[diffIndex])
            solClause += "%d 0\n" % (unifVal[var]*(var+sumNewVar))

    ##########
    # We add the N number of chain formulas
    # where chainFormulaConf.indicatorLits must be of size 2
    # and len(chainFormulaConf.indicatorLits) == len(chainFormulaConf.newVarList)
    # Adding K soluitons over Z variables, where
    #    Z = chainFormulaConf.newVarList[k]
    #    K = chainFormulaConf.countList[k]
    ##########
    invert = True
    seenLits = {}
    for indicLits in chainFormulaConf.indicatorLits:   # loop runs twice
        currentNumVar = 0
        for i in range(len(indicLits)):
            newvar = chainFormulaConf.newVarList[i]
            indicLit = int(indicLits[i])
            addedClause = ''
            addedClauseNum = 0

            # not adding the chain formula twice to the same literal
            if indicLit not in seenLits:
                sign = int(indicLit/abs(indicLit))
                addedClause, addedClauseNum = constructChainFormula(
                    sign*(abs(indicLit)+sumNewVar),
                    chainFormulaConf.countList[i], newvar, currentNumVar,
                    invert)

            seenLits[indicLit] = True
            currentNumVar += newvar
            numCls += addedClauseNum
            solClause += addedClause
        invert = not invert
    del seenLits
    del invert

    # create "c ind ..." lines
    oldIndVarList = [x+sumNewVar for x in indVarList]
    tempIndVarList = []
    indIter = 1
    indStr = 'c ind '
    for i in range(1, currentNumVar+1):
        if indIter % 10 == 0:
            indStr += ' 0\nc ind '
        indStr += "%d " % i
        indIter += 1
        tempIndVarList.append(i)

    for i in oldIndVarList:
        if indIter % 10 == 0:
            indStr += ' 0\nc ind '
        indStr += "%d " % i
        indIter += 1
        tempIndVarList.append(i)
    indStr += ' 0\n'

    # dump new CNF
    with open(tempFile, 'w') as f:
        f.write('p cnf %d %d\n' % (currentNumVar+numVar, numCls))
        f.write(indStr)
        f.write(solClause)
        #f.write("c -- old CNF below -- \n")
        f.write(shiftedCNFStr)

    #print("New file: ", tempFile)
    #exit(0)

    return False, tempIndVarList, oldIndVarList


class Experiment:
    def __init__(self, inputFile, maxSamples, samplerType):
        inputFileSuffix = inputFile.split('/')[-1][:-4]
        self.tempFile = tempfile.gettempdir() + "/" + inputFileSuffix+"_t.cnf"
        self.indVarList = parseIndSupport(inputFile)
        self.inputFile = inputFile
        self.samplerType = samplerType
        self.maxSamples = maxSamples

        self.samplerString = get_sampler_string(samplerType)

    # Returns True if uniform and False otherwise
    def testUniformity(self, solList, indVarList):
        solMap = {}
        baseMap = {}
        for sol in solList:
            solution = ''
            solFields = sol.split()
            for entry in solFields:
                if abs(int(entry)) in indVarList:
                    solution += entry+' '

            if solution in solMap.keys():
                solMap[solution] += 1
            else:
                solMap[solution] = 1

            if sol not in baseMap.keys():
                baseMap[sol] = 1
            else:
                baseMap[sol] += 1

        if not bool(solMap):
            print("No Solutions were given to the test")
            exit(1)

        key = next(iter(solMap))

        print("baseMap: {:<6} numSolutions: {:<6} SolutionsCount: {:<6} loThresh: {:<6} hiThresh: {:<6}".format(
            len(baseMap.keys()), self.numSolutions, solMap[key], self.loThresh, self.hiThresh))

        if solMap[key] >= self.loThresh and solMap[key] <= self.hiThresh:
            return True
        else:
            return False

    def one_experiment(self, j, i, tj):
        self.thresholdSolutions += self.numSolutions

        # generate a new seed value for every different (i,j,experiment)
        newSeed = int(i*tj+j)
        # get sampler's solutions
        sampleSol = SolutionRetriever.getSolutionFromSampler(
            self.inputFile, 1, self.samplerType, self.indVarList, verbosity, newSeed)
        self.totalSolutionsGenerated += 1

        # get uniform sampler's solutions
        unifSol = SolutionRetriever.getSolutionFromSampler(
            self.inputFile, 1, SAMPLER_SPUR, self.indVarList, verbosity, newSeed)
        self.totalUniformSamples += 1

        chainFormulaConf = chainFormulaSetup(sampleSol, unifSol, self.numSolutions)
        identical_sol, tempIndVarList, oldIndVarList = constructNewCNF(
            self.inputFile, self.tempFile, sampleSol[0], unifSol[0],
            chainFormulaConf, self.indVarList)

        # the two solutions were the same, couldn't construct CNF
        if identical_sol:
            return False

        # seed update
        newSeed = newSeed + 1

        # get sampler's solutions
        solList = SolutionRetriever.getSolutionFromSampler(
            self.tempFile, self.numSolutions, self.samplerType, tempIndVarList, verbosity, newSeed)
        os.unlink(self.tempFile)
        self.totalSolutionsGenerated += self.numSolutions

        isUniform = self.testUniformity(solList, oldIndVarList)

        print("sampler: {:<8s} i: {:<4d} isUniform: {:<4d} TotalSolutionsGenerated: {:<6d}".format(
            self.samplerString, i, isUniform,
            self.totalSolutionsGenerated))

        if not isUniform:
            print("RejectIteration:{0}  Loop:{1} TotalSolutionsGenerated:{2} TotalUniformSamples:{3}".format(
                i, j, self.totalSolutionsGenerated, self.totalUniformSamples))

            return True

        if self.thresholdSolutions > self.maxSamples:
            return True

        return False


if __name__ == "__main__":

    samplers = str(SAMPLER_UNIGEN) + " for UniGen\n"
    samplers += str(SAMPLER_QUICKSAMPLER) + " for QuickSampler\n"
    samplers += str(SAMPLER_STS)+ " for STS\n"
    samplers += str(SAMPLER_SPUR) + " for SPUR\n"
    samplers += str(SAMPLER_CMS) + " for CMS\n"

    parser = argparse.ArgumentParser()
    parser.add_argument('--sampler', type=int, help=samplers, default=SAMPLER_STS, dest='sampler')
    parser.add_argument('--eta', type=float, help="default = 0.9", default=0.9, dest='eta')
    parser.add_argument('--epsilon', type=float, help="default = 0.3", default=0.3, dest='epsilon')
    parser.add_argument('--delta', type=float, help="default = 0.05", default=0.05, dest='delta')
    parser.add_argument('--reverse', type=int, default=0, help="order to search in", dest='searchOrder')
    parser.add_argument('--maxSamples', type=int, default=sys.maxsize, help="max samples", dest='maxSamples')
    parser.add_argument('--seed', type=int, required=True, dest='seed')
    parser.add_argument('--verb', type=int, dest='verbose')
    parser.add_argument("input", help="input file")

    args = parser.parse_args()
    inputFile = args.input

    eta = args.eta
    epsilon = args.epsilon
    if (eta < 2*epsilon):
        print("Eta needs to be at least two times epsilon")
        exit(1)
    delta = args.delta
    searchOrder = args.searchOrder
    verbosity = args.verbose

    seed = args.seed
    random.seed(seed)
    maxSamples = args.maxSamples

    totalLoops = int(math.ceil(math.log(2.0/(eta+2*epsilon), 2))+1)
    listforTraversal = range(totalLoops, 0, -1)
    if searchOrder == 1:
        listforTraversal = range(1, totalLoops+1, 1)

    exp = Experiment(
        maxSamples=maxSamples, inputFile=inputFile,
        samplerType=args.sampler)

    breakExperiment = False
    exp.totalSolutionsGenerated = 0
    exp.totalUniformSamples = 0
    exp.thresholdSolutions = 0
    for j in listforTraversal:
        tj = math.ceil(math.pow(2, j)*(2*epsilon+eta)/((eta-2*epsilon)**2)*math.log(4.0/(eta+2*epsilon), 2)*(4*math.e/(math.e-1)*math.log(1.0/delta)))
        beta = (math.pow(2, j-1)+1)*(eta + 2*epsilon)*1.0/(4+(2*epsilon+eta)*(math.pow(2, j-1) - 1))
        gamma = (beta-2*epsilon)/4
        constantFactor = math.ceil(1/(8.79*gamma*gamma))
        boundFactor = math.log((16)*(math.e/(math.e-1))*(1/(delta*(eta-2*epsilon)**2))*math.log(4/(eta+2*epsilon), 2)*math.log(1/delta), 2)
        print("constantFactor:{:<4} boundFactor: {:<20} logBoundFactor: {:<20}".format(
            constantFactor, boundFactor, math.log(boundFactor, 2)))
        print("tj: {:<6} totalLoops: {:<5} beta: {:<10} epsilon: {:<10}".format(
            tj, totalLoops, beta, epsilon))

        exp.numSolutions = int(math.ceil(constantFactor*boundFactor))
        exp.loThresh = int((exp.numSolutions*1.0/2)*(1-(beta+2*epsilon)/2))
        exp.hiThresh = int((exp.numSolutions*1.0/2)*(1+(beta+2*epsilon)/2))
        print("numSolutions: {:<5} loThresh:{:<6} hiThresh: {:<6}".format(
            exp.numSolutions, exp.loThresh, exp.hiThresh))

        breakExperiment = False

        for i in range(int(tj)):
            breakExperiment = exp.one_experiment(j, i, tj)

            if breakExperiment:
                break

        if breakExperiment:
            break

    if not breakExperiment:
        print("exp:{2} Accept:1 TotalSolutionsGenerated:{0} TotalUniformSamples:{1}".format(
            exp.totalSolutionsGenerated,
            exp.totalUniformSamples))

    breakExperiment = False
