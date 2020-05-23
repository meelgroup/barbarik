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
SAMPLER_QUICKSAMPLER = 2
SAMPLER_STS = 3
SAMPLER_CMS = 4
SAMPLER_APPMC3 = 5
SAMPLER_SPUR = 6

def get_sampler_string(samplerType):
    if samplerType == SAMPLER_UNIGEN:
        return 'UniGen'
    if samplerType == SAMPLER_APPMC3:
        return 'AppMC3'
    if samplerType == SAMPLER_QUICKSAMPLER:
        return 'QuickSampler'
    if samplerType == SAMPLER_STS:
        return 'STS'
    if samplerType == SAMPLER_CMS:
        return 'CustomSampler'
    if samplerType == SAMPLER_SPUR:
        return 'SPUR'
    print("ERROR: unknown sampler type")
    exit(-1)


class ChainFormulaSetup:
    def __init__(self, countList, newVarList, indicatorLits):
        self.countList = countList
        self.newVarList = newVarList
        self.indicatorLits = indicatorLits

def check_cnf(fname):
    with open(fname, 'r') as f:
        lines = f.readlines()

    given_vars = None
    given_cls = None
    cls = 0
    max_var = 0
    for line in lines:
        line = line.strip()

        if len(line) == 0:
            print("ERROR: CNF is incorrectly formatted, empty line!")
            return False

        line = line.split()
        line = [l.strip() for l in line]

        if line[0] == "p":
            assert len(line) == 4
            assert line[1] == "cnf"
            given_vars = int(line[2])
            given_cls = int(line[3])
            continue

        if line[0] == "c":
            continue

        cls +=1
        for l in line:
            var = abs(int(l))
            max_var = max(var, max_var)

    if max_var > given_vars:
        print("ERROR: Number of variables given is LESS than the number of variables ued")
        print("ERROR: Vars in header: %d   max var: %d" % (given_vars, max_var))
        return False

    if cls != given_cls:
        print("ERROR: Number of clauses in header is DIFFERENT than the number of clauses in the CNF")
        print("ERROR: Claues in header: %d   clauses: %d" % (given_cls, cls))
        return False

    return True


class SolutionRetriver:

    @staticmethod
    def getSolutionFromSampler(inputFile, numSolutions, samplerType, indVarList, newSeed):
        topass_withseed = (inputFile, numSolutions, indVarList, newSeed)
        ok = check_cnf(inputFile)
        if not ok:
            print("ERROR: CNF is malformatted. Sampler may give wrong solutions in this case. Exiting.")
            print("File is: %s" % inputFile)
            exit(-1)


        print("Using sampler: %s" % get_sampler_string(samplerType))
        if (samplerType == SAMPLER_UNIGEN):
            sols = SolutionRetriver.getSolutionFromUniGen(*topass_withseed)

        elif (samplerType == SAMPLER_APPMC3):
            sols = SolutionRetriver.getSolutionFromAppMC3(*topass_withseed)

        elif (samplerType == SAMPLER_QUICKSAMPLER):
            sols = SolutionRetriver.getSolutionFromQuickSampler(*topass_withseed)

        elif (samplerType == SAMPLER_STS):
            sols = SolutionRetriver.getSolutionFromSTS(*topass_withseed)

        elif (samplerType == SAMPLER_CMS):
            sols = SolutionRetriver.getSolutionFromCMSsampler(*topass_withseed)

        elif (samplerType == SAMPLER_SPUR):
            sols = SolutionRetriver.getSolutionFromSpur(*topass_withseed)

        else:
            print("Error: No such sampler!")
            exit(-1)

        # clean up the solutions
        for i in range(len(sols)):
            sols[i] = sols[i].strip()
            if sols[i].endswith(' 0'):
                sols[i] = sols[i][:-2]

        print("Number of solutions returned by sampler:", len(sols))
        if args.verbose:
            print("Solutions:", sols)
        return sols

    @staticmethod
    def getSolutionFromUniGen(inputFile, numSolutions, indVarList, newSeed):
        # must construct ./unigen --samples=500 --verbosity=0 --threads=1  CNF-FILE SAMPLESFILE
        inputFileSuffix = inputFile.split('/')[-1][:-4]
        tempOutputFile = tempfile.gettempdir()+'/'+inputFileSuffix+".txt"

        cmd = './samplers/unigen --samples='+str(numSolutions)
        cmd += ' ' + inputFile + ' ' + str(tempOutputFile) + ' > /dev/null 2>&1'
        if args.verbose:
            print("cmd: ", cmd)
        os.system(cmd)

        with open(tempOutputFile, 'r') as f:
            lines = f.readlines()

        solList = []
        for line in lines:
            line = line.strip()
            if line.startswith('v'):
                freq = int(line.split(':')[-1])
                for i in range(freq):
                    solList.append(line.split(':')[0].replace('v', '').strip())
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
    def getSolutionFromAppMC3(inputFile, numSolutions, indVarList, newSeed):
        # must construct: ./approxmc3 -s 1 -v2 --sampleout /dev/null --samples 500
        inputFileSuffix = inputFile.split('/')[-1][:-4]
        tempOutputFile = tempfile.gettempdir()+'/'+inputFileSuffix+".txt"

        cmd = './samplers/approxmc3 -s ' + str(newSeed) + ' -v 0 --samples ' + str(numSolutions)
        cmd += ' --sampleout ' + str(tempOutputFile)
        cmd += ' ' + inputFile + ' > /dev/null 2>&1'
        if args.verbose:
            print("cmd: ", cmd)
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
    def getSolutionFromQuickSampler(inputFile, numSolutions, indVarList, newSeed):
        cmd = "./samplers/quicksampler -n "+str(numSolutions*5)+' '+str(inputFile)+' > /dev/null 2>&1'
        if args.verbose:
            print("cmd: ", cmd)
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
            # valutions are 0 and 1 and in the same order as c ind.
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

        if len(solList) != numSolutions:
            print("Did not find required number of solutions")
            sys.exit(1)
        return solList

    @staticmethod
    def getSolutionFromSpur(inputFile, numSolutions, indVarList, newSeed):
        inputFileSuffix = inputFile.split('/')[-1][:-4]
        tempOutputFile = tempfile.gettempdir()+'/'+inputFileSuffix+".out"
        cmd = './samplers/spur -seed %d -q -s %d -out %s -cnf %s' % (
            newSeed, numSolutions, tempOutputFile, inputFile)
        if args.verbose:
            print("cmd: ", cmd)
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
    def getSolutionFromSTS(inputFile, numSolutions, indVarList, newSeed):
        kValue = 50
        samplingRounds = numSolutions/kValue + 1
        inputFileSuffix = inputFile.split('/')[-1][:-4]
        outputFile = tempfile.gettempdir()+'/'+inputFileSuffix+".out"
        cmd = './samplers/STS -k='+str(kValue)+' -nsamples='+str(samplingRounds)+' '+str(inputFile)
        cmd += ' > '+str(outputFile)
        if args.verbose:
            print("cmd: ", cmd)
        os.system(cmd)

        with open(outputFile, 'r') as f:
            lines = f.readlines()

        solList = []
        shouldStart = False
        #baseList = {}
        for j in range(len(lines)):
            if(lines[j].strip() == 'Outputting samples:' or lines[j].strip() == 'start'):
                shouldStart = True
                continue
            if (lines[j].strip().startswith('Log') or lines[j].strip() == 'end'):
                shouldStart = False
            if (shouldStart):


                '''if lines[j].strip() not in baseList:
                    baseList[lines[j].strip()] = 1
                else:
                    baseList[lines[j].strip()] += 1'''
                sol = ''
                i = 0
                # valutions are 0 and 1 and in the same order as c ind.
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
            sys.exit(1)

        os.unlink(outputFile)
        return solList

    @staticmethod
    def getSolutionFromCMSsampler(inputFile, numSolutions, indVarList, newSeed):
        inputFileSuffix = inputFile.split('/')[-1][:-4]
        outputFile = tempfile.gettempdir()+'/'+inputFileSuffix+".out"
        cmd = "./samplers/cryptominisat5 --restart luby --maple 0 --verb 10 --nobansol"
        cmd += " --scc 1 -n1 --presimp 0 --polar rnd --freq 0.9999"
        cmd += " --random " + str(newSeed) + " --maxsol " + str(numSolutions)
        cmd += " " + inputFile
        cmd += " --dumpresult " + outputFile + " > /dev/null 2>&1"

        if args.verbose:
            print("cmd: ", cmd)
        os.system(cmd)

        with open(outputFile, 'r') as f:
            lines = f.readlines()

        solList = []
        for line in lines:
            if line.strip() == 'SAT':
                continue

            sol = ""
            lits = line.split(" ")
            for y in indVarList:
                if str(y) in lits:
                    sol += ' ' + str(y)

                if "-" + str(y) in lits:
                    sol += ' -' + str(y)
            solList.append(sol)

        solreturnList = solList
        if len(solList) > numSolutions:
            solreturnList = random.sample(solList, numSolutions)
        if len(solList) < numSolutions:
            print("cryptominisat5 Did not find required number of solutions")
            sys.exit(1)
        os.unlink(outputFile)
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
        return False, None, None

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

    def one_experiment(self, experiment, j, i, numExperiments, tj):
        self.thresholdSolutions += self.numSolutions
        if self.thresholdSolutions < self.minSamples:
            return None, None

        # generate a new seed value for every different (i,j,experiment)
        newSeed = int(numExperiments*(i*tj+j)+experiment)
        # get sampler's solutions
        sampleSol = SolutionRetriver.getSolutionFromSampler(
            self.inputFile, 1, self.samplerType, self.indVarList, newSeed)
        self.totalSolutionsGenerated += 1

        # get uniform sampler's solutions
        # get uniform sampler's solutions
        unifSol = SolutionRetriver.getSolutionFromSampler(
            self.inputFile, 1, SAMPLER_SPUR, self.indVarList, newSeed)
        self.totalUniformSamples += 1

        chainFormulaConf = chainFormulaSetup(sampleSol, unifSol, self.numSolutions)
        shakuniMix, tempIndVarList, oldIndVarList = constructNewCNF(
            self.inputFile, self.tempFile, sampleSol[0], unifSol[0],
            chainFormulaConf, self.indVarList)

        # the two solutions were the same, couldn't construct CNF
        if not shakuniMix:
            return False, None

        # seed update
        newSeed = newSeed + 1

        # get sampler's solutions
        solList = SolutionRetriver.getSolutionFromSampler(
            self.tempFile, self.numSolutions, self.samplerType, tempIndVarList, newSeed)
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


if __name__ == "__main__":
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

            i = 0
            breakExperiment = False
            while i < int(tj) and not breakExperiment:
                i += 1
                ok, breakExperiment = exp.one_experiment(experiment, j, i, numExperiments, tj)

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
