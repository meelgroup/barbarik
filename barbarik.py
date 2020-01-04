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

SAMPLER_UNIGEN = 1
SAMPLER_APPMC3 = 5
SAMPLER_QUICKSAMPLER = 2
SAMPLER_STS = 3
SAMPLER_CUSTOM = 4


class RExtList:
    def __init__(self, countList, newVarList, oldLitLists):
        self.countList = countList
        self.newVarList = newVarList
        self.oldLitLists = oldLitLists


class SolutionRetriver:

    @staticmethod
    def getSolutionFromSampler(inputFile, numSolutions, samplerType, indVarList):
        topass = (inputFile, numSolutions, indVarList)

        if (samplerType == SAMPLER_UNIGEN):
            return SolutionRetriver.getSolutionFromUniGen(*topass)

        if (samplerType == SAMPLER_APPMC3):
            return SolutionRetriver.getSolutionFromAppMC3(*topass)

        if (samplerType == SAMPLER_QUICKSAMPLER):
            return SolutionRetriver.getSolutionFromQuickSampler(*topass)

        if (samplerType == SAMPLER_STS):
            return SolutionRetriver.getSolutionFromSTS(*topass)

        if (samplerType == SAMPLER_CUSTOM):
            return SolutionRetriver.getSolutionFromCustomSampler(*topass)

        else:
            print("Error")
            return None

    @staticmethod
    def getSolutionFromUniGen(inputFile, numSolutions, indVarList):
        # must construct ./unigen --samples=500 --verbosity=0 --threads=1  CNF-FILE SAMPLESFILE
        inputFileSuffix = inputFile.split('/')[-1][:-4]
        tempOutputFile = tempfile.gettempdir()+'/'+inputFileSuffix+".txt"

        cmd = './samplers/unigen --samples='+str(numSolutions)
        cmd += +inputFile + ' ' + str(tempOutputFile) + ' > /dev/null 2>&1'
        os.system(cmd)

        with open(tempOutputFile, 'r') as f:
            lines = f.readlines()

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

        os.unlink(str(tempOutputFile))
        return solreturnList

    @staticmethod
    def getSolutionFromAppMC3(inputFile, numSolutions, indVarList):
        # must construct: ./approxmc3 -s 1 -v2 --sampleout /dev/null --samples 500
        inputFileSuffix = inputFile.split('/')[-1][:-4]
        tempOutputFile = tempfile.gettempdir()+'/'+inputFileSuffix+".txt"

        cmd = './samplers/approxmc3 -s 1 -v 0 --samples ' + str(numSolutions)
        cmd += ' --sampleout ' + str(tempOutputFile)
        cmd += ' ' + inputFile + ' > /dev/null 2>&1'
        # print("Calling: '%s'" % cmd)
        os.system(cmd)

        with open(tempOutputFile, 'r') as f:
            lines = f.readlines()

        solList = []
        for line in lines:
            line = line.strip()
            freq = int(line.split(':')[0])
            for i in range(freq):
                solList.append(line.split(':')[1].strip())
                if len(solList) == numSolutions:
                    break
            if len(solList) == numSolutions:
                break

        solreturnList = solList
        if len(solList) > numSolutions:
            solreturnList = random.sample(solList, numSolutions)

        os.unlink(str(tempOutputFile))
        return solreturnList

    @staticmethod
    def getSolutionFromQuickSampler(inputFile, numSolutions, indVarList):
        cmd = "./samplers/quicksampler -n "+str(numSolutions*5)+' '+str(inputFile)+' > /dev/null 2>&1'
        os.system(cmd)
        cmd = "./samplers/z3 "+str(inputFile)+' > /dev/null 2>&1'
        os.system(cmd)
        if (numSolutions > 1):
            i = 0

        with open(inputFile+'.samples', 'r') as f:
            lines = f.readlines()

        with open(inputFile+'.samples.valid', 'r') as f:
            validLines = f.readlines()

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

        os.unlink(inputFile+'.samples')
        os.unlink(inputFile+'.samples.valid')

        if (len(solList) != numSolutions):
            print("Did not find required number of solutions")
            exit(1)
        return solList

    @staticmethod
    def getSolutionFromUniform(inputFile, numSolutions):
        return SolutionRetriver.getSolutionFromMUSE(inputFile, numSolutions)

    @staticmethod
    def getSolutionFromMUSE(inputFile, numSolutions):
        inputFileSuffix = inputFile.split('/')[-1][:-4]
        tempOutputFile = tempfile.gettempdir()+'/'+inputFileSuffix+".out"
        cmd = './spur -q -s '+str(numSolutions)+' -out '+str(tempOutputFile)+' -cnf '+str(inputFile)
        os.system(cmd)

        with open(tempOutputFile, 'r') as f:
            lines = f.readlines()

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

        os.unlink(tempOutputFile)
        return solList

    @staticmethod
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

                if lines[j].strip() not in baseList:
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
                if len(solList) == numSolutions:
                    break

        if len(solList) != numSolutions:
            print(len(solList))
            print("STS Did not find required number of solutions")
            exit(1)

        os.unlink(outputFile)
        return solList

    # @CHANGE_HERE : please make changes in the below block of code
    ''' this is the method where you could run your sampler for testing
    Arguments : input file, number of solutions to be returned, list of independent variables
    output : list of solutions '''
    @staticmethod
    def getSolutionFromCustomSampler(inputFile, numSolutions, indVarList):
        solreturnList = []

        ''' write your code here '''

        return solreturnList


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


def findWeightsForVariables(sampleSol, unifSol, numSolutions):
    countList = [5, 5, 5]
    newVarList = [4, 4, 4]

    ##########
    # clean up solutions
    ##########
    sampleSol[0] = sampleSol[0].strip()
    if sampleSol[0].endswith(' 0'):
        sampleSol[0] = sampleSol[0][:-2]

    unifSol[0] = unifSol[0].strip()
    if unifSol[0].endswith(' 0'):
        unifSol[0] = unifSol[0][:-2]

    lenSol = len(sampleSol[0].split())
    logSol = math.log(5, 2)*3+1     # this is 7.965784284662087

    for i in range(min(int(math.log(numSolutions, 2))+4, lenSol-3, 5)):
        countNum = 31  # random.randint(1,64)
        countList.append(countNum)
        newVarList.append(6)
        logSol += math.log(countNum, 2)

    sampleLitList = random.sample(sampleSol[0].split(), len(countList))
    unifLitList = []
    unifSolMap = unifSol[0].split()
    for lit in sampleLitList:
        unifLitList.append(unifSolMap[abs(int(lit))-1])

    oldLitLists = []
    oldLitLists.append(sampleLitList)
    oldLitLists.append(unifLitList)

    #print("countList:", countList)
    #print("newVarList:", newVarList)
    #print("oldLitLists:", oldLitLists)
    return RExtList(countList, newVarList, oldLitLists)


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
def constructNewCNF(inputFile, tempFile, sampleSol, unifSol, rExtList, indVarList):
    # which variables are in pos/neg value in the sample
    sampleMap = {}
    for i in sampleSol.strip().split():
        i = int(i)
        if i != 0:
            if abs(i) not in indVarList:
                continue

            sampleMap[abs(i)] = int(i/abs(i))

    # which variables are in pos/neg value in the uniform sample
    unifMap = {}
    diffIndex = -1
    for j in unifSol.strip().split():
        j = int(j)
        if j != 0:
            if abs(j) not in indVarList:
                continue

            unifMap[abs(j)] = int(j/abs(j))

            if sampleMap[abs(j)] != unifMap[abs(j)]:
                diffIndex = abs(j)

    # the two solutions are the same
    # can't do anything, let's do another experiment
    if diffIndex == -1:
        return False, None, None

    with open(inputFile, 'r') as f:
        lines = f.readlines()

    # shift amount is sumNewVar
    sumNewVar = int(sum(rExtList.newVarList))

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

    # Adding constraints to ensure only two solutions
    # All variables are set except for the index where they last differ
    solClause = ''
    for var in indVarList:
        if var != diffIndex:
            numCls += 2
            solClause += "%d " % (-(diffIndex+sumNewVar)*sampleMap[diffIndex])
            solClause += "%d 0\n" % (sampleMap[var]*(var+sumNewVar))

            solClause += "%d " % (-(diffIndex+sumNewVar)*unifMap[diffIndex])
            solClause += "%d 0\n" % (unifMap[var]*(var+sumNewVar))

    ##########
    #
    ##########
    invert = True
    seenLits = {}

    # there are only 2 in fact, the one from the sampler under test
    # and one that's from the uniform sampler. One is inverted one is not.
    for oldLitList in rExtList.oldLitLists:
        currentNumVar = 0
        for i in range(len(oldLitList)):
            newvar = int(rExtList.newVarList[i])
            oldlit = int(oldLitList[i])
            addedClause = ''
            addedClauseNum = 0
            if oldlit not in seenLits:
                sign = int(oldlit/abs(oldlit))
                addedClause, addedClauseNum = constructChainFormula(
                    sign*(abs(oldlit)+sumNewVar),
                    rExtList.countList[i], newvar, currentNumVar,
                    invert)

            seenLits[oldlit] = True
            currentNumVar += newvar
            numCls += addedClauseNum
            solClause += addedClause
        invert = not invert
    del seenLits
    del invert

    # create "c ind ..." lines
    oldIndVarList = [x+sumNewVar for x in indVarList]
    tempIndVarList = copy.copy(oldIndVarList)
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

    return True, tempIndVarList, oldIndVarList


class Experiment:
    def __init__(self, inputFile, maxSamples, minSamples, samplerType):
        inputFileSuffix = inputFile.split('/')[-1][:-4]
        self.tempFile = tempfile.gettempdir() + "/" + inputFileSuffix+"_t.cnf"
        self.indVarList = parseIndSupport(inputFile)
        self.inputFile = inputFile
        self.samplerType = samplerType
        self.maxSamples = maxSamples
        self.minSamples = minSamples

        self.samplerString = None
        if samplerType == SAMPLER_UNIGEN:
            self.samplerString = 'UniGen'
        if samplerType == SAMPLER_APPMC3:
            self.samplerString = 'AppMC3'
        if samplerType == SAMPLER_QUICKSAMPLER:
            self.samplerString = 'QuickSampler'
        if samplerType == SAMPLER_STS:
            self.samplerString = 'STS'
        if samplerType == SAMPLER_CUSTOM:
            self.samplerString = 'CustomSampler'

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

    def one_experiment(self, experiment, j, i):
        self.thresholdSolutions += self.numSolutions
        if self.thresholdSolutions < self.minSamples:
            return None, None

        # get sampler's solutions
        sampleSol = SolutionRetriver.getSolutionFromSampler(
            self.inputFile, 1, self.samplerType, self.indVarList)
        self.totalSolutionsGenerated += 1

        # get uniform sampler's solutions
        unifSol = SolutionRetriver.getSolutionFromUniform(
            self.inputFile, 1)
        self.totalUniformSamples += 1

        rExtList = findWeightsForVariables(sampleSol, unifSol, self.numSolutions)
        shakuniMix, tempIndVarList, oldIndVarList = constructNewCNF(
            self.inputFile, self.tempFile, sampleSol[0], unifSol[0], rExtList, self.indVarList)

        # the two solutions were the same, couldn't construct CNF
        if not shakuniMix:
            return False, None

        # get sampler's solutions
        solList = SolutionRetriver.getSolutionFromSampler(
            self.tempFile, self.numSolutions, self.samplerType, tempIndVarList)
        os.unlink(self.tempFile)
        self.totalSolutionsGenerated += self.numSolutions

        isUniform = self.testUniformity(solList, oldIndVarList)

        print("sampler: {:<8s} i: {:<4d} isUniform: {:<4d} TotalSolutionsGenerated: {:<6d}".format(
            self.samplerString, i, isUniform,
            self.totalSolutionsGenerated))

        if not isUniform:
            print("exp:{4} RejectIteration:{0}  Loop:{1} TotalSolutionsGenerated:{2} TotalUniformSamples:{3}".format(
                i, j, self.totalSolutionsGenerated, self.totalUniformSamples, experiment))

            return True, True

        if self.thresholdSolutions > self.maxSamples:
            return True, True

        return True, False


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
    parser.add_argument('--verb', type=int, dest='verbose')
    parser.add_argument('--exp', type=int, help="number of experiments", dest='exp', default=1)
    parser.add_argument("input", help="input file")

    args = parser.parse_args()
    inputFile = args.input

    eta = args.eta
    epsilon = args.epsilon
    delta = args.delta
    numExperiments = args.exp
    if numExperiments == -1:
        numExperiments = sys.maxsize
    searchOrder = args.searchOrder
    verbose = args.verbose

    seed = args.seed
    random.seed(seed)
    minSamples = args.minSamples
    maxSamples = args.maxSamples

    totalLoops = int(math.ceil(math.log(2.0/(eta+epsilon), 2))+1)
    listforTraversal = range(totalLoops, 0, -1)
    if searchOrder == 1:
        listforTraversal = range(1, totalLoops+1, 1)

    exp = Experiment(
        minSamples=minSamples, maxSamples=maxSamples, inputFile=inputFile,
        samplerType=args.sampler)

    for experiment in range(numExperiments):
        print("Experiment: {:<5} of {:>5}".format(experiment, numExperiments))
        breakExperiment = False
        exp.totalSolutionsGenerated = 0
        exp.totalUniformSamples = 0
        exp.thresholdSolutions = 0
        for j in listforTraversal:
            tj = math.ceil(math.pow(2, j)*(epsilon+eta)/((eta-epsilon)**2)*math.log(4.0/(eta+epsilon), 2)*(4*math.e/(math.e-1)*math.log(1.0/delta)))
            beta = (math.pow(2, j-1)+1)*(eta + epsilon)*1.0/(4+(epsilon+eta)*(math.pow(2, j-1) - 1))
            gamma = (beta-epsilon)/4
            constantFactor = math.ceil(1/(9*gamma*gamma))
            boundFactor = math.log((16)*(math.e/(math.e-1))*(1/((eta-epsilon)**2))*math.log(4/(eta+epsilon), 2)*math.log(1/delta), 2)
            print("constantFactor:{:<4} boundFactor: {:<20} logBoundFactor: {:<20}".format(
                constantFactor, boundFactor, math.log(boundFactor, 2)))
            print("tj: {:<6} totalLoops: {:<5} beta: {:<10} epsilon: {:<10}".format(
                tj, totalLoops, beta, epsilon))

            exp.numSolutions = int(math.ceil(constantFactor*boundFactor))
            exp.loThresh = int((exp.numSolutions*1.0/2)*(1-(beta+epsilon)/2))
            exp.hiThresh = int((exp.numSolutions*1.0/2)*(1+(beta+epsilon)/2))
            print("numSolutions: {:<5} loThresh:{:<6} hiThresh: {:<6}".format(
                exp.numSolutions, exp.loThresh, exp.hiThresh))

            i = 0
            breakExperiment = False
            while i < int(tj) and not breakExperiment:
                i += 1
                ok, breakExperiment = exp.one_experiment(experiment, j, i)

                if ok is None:
                    continue

                if not ok:
                    i -= 1
                    continue

                if breakExperiment:
                    break

            if breakExperiment:
                break

        if not breakExperiment:
            print("exp:{2} Accept:1 TotalSolutionsGenerated:{0} TotalUniformSamples:{1}".format(
                exp.totalSolutionsGenerated,
                exp.totalUniformSamples, experiment))

        breakExperiment = False


if __name__ == "__main__":
    barbarik()
