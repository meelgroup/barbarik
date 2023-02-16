import sys
import os
from math import ceil, log, sqrt, e
from copy import deepcopy
import random
import tempfile

from interfaces.WAPS.waps import sampler as samp
import interfaces.weightcount.WeightCount as chainform

SAMPLER_UNIGEN3 = 1
SAMPLER_QUICKSAMPLER = 2
SAMPLER_STS = 3
SAMPLER_CMS = 4
SAMPLER_APPMC3 = 5
SAMPLER_SPUR = 6

SAMPLER_CUSTOM = 7


def get_sampler_string(samplerType):
    if samplerType == SAMPLER_UNIGEN3:
        return 'UniGen3'
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


# @CHANGE_HERE : please make changes in the below block of code
# this is the method where you could run your sampler for testing
# Arguments : input file, number of solutions to be returned, list of independent variables
# output : list of solutions


def getSolutionFromCustomSampler(inputFile, numSolutions, indVarList):
    solreturnList = []
    # write your code here

    return solreturnList


# def getSolutionFromSampler(seed, inputFile, numSolutions, samplerType, indVarList):
#     if samplerType == SAMPLER_UNIGEN3:
#         return getSolutionFromUniGen3(inputFile, numSolutions, indVarList)
#     if samplerType == SAMPLER_QUICKSAMPLER:
#         return getSolutionFromQuickSampler(inputFile, numSolutions, indVarList)
#     if samplerType == SAMPLER_STS:
#         return getSolutionFromSTS(seed, inputFile, numSolutions, indVarList)
#     if samplerType == SAMPLER_CUSTOM:
#         return getSolutionFromCustomSampler(inputFile, numSolutions, indVarList)
#     else:
#         print("Error")
#         return None


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
            line = line.strip().replace('c ind', '').replace(
                ' 0', '').strip().replace('v ', '')
            indList.extend(line.split())

    if len(indList) == 0:
        indList = [int(x) for x in range(1, numVars+1)]
    else:
        indList = [int(x) for x in indList]
    return indList


def hoeffding(p1, p2, delta):
    assert (p1 < p2)
    return int(ceil(2*log(1/delta)/(p2-p1)**2))


def parseWeights(inputFile, indVarList):
    f = open(inputFile, "r")
    lines = f.readlines()
    f.close()
    weight_map = {}

    for line in lines:
        if line.startswith("w"):
            variable, weight = line[2:].strip().split()
            variable = int(variable)
            weight = float(weight)
            if (0.0 < weight < 1.0):
                if (variable in indVarList):
                    weight_map[variable] = weight
            else:
                print("Error: weights should only be in (0,1) ")
                exit(-1)
    return weight_map


def weightof(weight_map, sample, UserIndVarList):
    sample_w = deepcopy(sample)

    weight = 1.0

    for i in range(len(sample_w)):
        if sample_w[i] > 0:
            weight *= weight_map.get(abs(sample_w[i]), 0.5)
        else:
            weight *= 1-weight_map.get(abs(sample_w[i]), 0.5)
    return weight


def bucketof(weight_map, sample, UserIndVarList):
    return int(ceil(log(weightof(weight_map, sample, UserIndVarList), 2)))


def findWeightsForVariables(sampleSol, idealSol, numSolutions):
    """
    Finds rExtList
    """
    countList = []
    newVarList = []
    lenSol = len(sampleSol)
    tau = min(3, lenSol)

    for _ in range(tau):
        countList.append(5)
        newVarList.append(4)
    rExtList = []
    oldVarList = []

    indexes = random.sample(range(len(sampleSol)), len(countList))

    idealVarList = [idealSol[i] for i in indexes]
    sampleVarList = [sampleSol[i] for i in indexes]

    assert len(idealVarList) == len(sampleVarList)
    for a, b in zip(idealVarList, sampleVarList):
        assert abs(int(a)) == abs(int(b))

    oldVarList.append(sampleVarList)
    oldVarList.append(idealVarList)
    rExtList.append(countList)
    rExtList.append(newVarList)
    rExtList.append(oldVarList)

    return rExtList

# @returns whether new file was created and the list of independent variables


def constructNewFile(inputFile, tempFile, sampleSol, unifSol, rExtList, origIndVarList):
    sampleMap = {}
    unifMap = {}
    diffIndex = -1  # ensures that sampleSol != unifSol when projected on indVarList
    for i in sampleSol:
        if not (abs(int(i)) in origIndVarList):
            continue
        if int(i) != 0:
            sampleMap[abs(int(i))] = int(int(i) / abs(int(i)))
    for j in unifSol:
        if int(j) != 0:
            if not (abs(int(j)) in origIndVarList):
                continue

            if sampleMap[abs(int(j))] != int(j) / abs(int(j)):
                diffIndex = abs(int(j))
            unifMap[abs(int(j))] = int(int(j) / abs(int(j)))

    if diffIndex == -1:
        print("Error: both samples are the same")
        print(sampleSol, unifSol)
        exit(-1)

    solClause = ""
    f = open(inputFile, "r")
    lines = f.readlines()
    f.close()
    countList = rExtList[0]
    newVarList = rExtList[1]
    sumNewVar = int(sum(newVarList))
    oldClauseStr = ""
    numVar = 0
    for line in lines:
        if line.strip().startswith("p cnf"):
            numVar = int(line.strip().split()[2])
            numClause = int(line.strip().split()[3])
        else:
            if line.strip().startswith("w"):
                oldClauseStr += line.strip()+"\n"
            elif not (line.strip().startswith("c")):
                oldClauseStr += line.strip()+"\n"
    # Adding constraints to ensure only two clauses
    for i in origIndVarList:
        if int(i) != diffIndex:
            numClause += 2
            solClause += (
                str(-(diffIndex) * sampleMap[diffIndex])
                + " "
                + str(sampleMap[int(i)] * int(i))
                + " 0\n"
            )
            solClause += (
                str(-(diffIndex) * unifMap[diffIndex])
                + " "
                + str(unifMap[int(i)] * int(i))
                + " 0\n"
            )

    invert = True
    seenVars = []
    currentNumVar = numVar

    for oldVarList in rExtList[2]:
        for i in range(len(oldVarList)):
            addedClause = ""
            addedClauseNum = 0

            if True or not (int(oldVarList[i]) in seenVars):
                sign = int(oldVarList[i]) / abs(int(oldVarList[i]))
                addedClause, addedClauseNum = constructChainFormula(
                    sign * (abs(int(oldVarList[i]))),
                    int(countList[i]),
                    int(newVarList[i]),
                    currentNumVar,
                    invert,
                )
            seenVars.append(int(oldVarList[i]))
            currentNumVar += int(newVarList[i])
            numClause += addedClauseNum
            solClause += addedClause
        invert = True  # not (invert)

    tempIndVarList = []
    indStr = "c ind "
    indIter = 1
    for i in origIndVarList:
        if indIter % 10 == 0:
            indStr += " 0\nc ind "
        indStr += str(i) + " "
        indIter += 1
        tempIndVarList.append(i)
    for i in range(numVar+1, currentNumVar + 1):
        if indIter % 10 == 0:
            indStr += " 0\nc ind "
        indStr += str(i) + " "
        indIter += 1
        tempIndVarList.append(i)

    indStr += " 0\n"
    headStr = "p cnf " + str(currentNumVar) + " " + str(numClause) + "\n"
    writeStr = headStr + indStr
    writeStr += solClause
    writeStr += oldClauseStr

    f = open(tempFile, "w")
    f.write(writeStr)
    f.close()
    return tempIndVarList


def constructKernel(inputFile, tempFile, samplerSample, idealSample, numSolutions, origIndVarList):
    rExtList = findWeightsForVariables(
        samplerSample, idealSample, numSolutions)
    tempIndVarList = constructNewFile(
        inputFile, tempFile, samplerSample, idealSample, rExtList, origIndVarList)
    return tempIndVarList


def insideBucket(K, eps, eps2, eta, delta, UserInputFile, inputFile, samplerType, indVarList, UserIndVarList, weight_map, ideal, totalSolsGenerated, seed):

    print("eta", eta)
    print("eps2", eps2)

    if (0.99*eta - 3.25*eps2 - 2*eps/(1-eps) < 0):
        print("Error: cannot test for these params")
        exit(1)

    tempFile = inputFile[:-6] + "_t.cnf"

    eps1 = (0.99*eta - 3.25*eps2 - 2*eps/(1-eps))/1.05 + 2*eps/(1-eps)
    alpha = (eps1+2*eps/(1-eps))/2

    M = int(ceil(sqrt(K)/(0.99*eta - 3.25*eps2 - eps1)))
    print("M #of samples in each round of InBucket = ", M)

    print("denom", (0.99*eta - 3.25*eps2 - eps1))

    T = int(ceil(log(2/delta)/log(10/(10 - eps1 + alpha))))
    print("T #of iteration = " + str(T))

    assert(M > 0)
    assert(T > 0)

    print("indVarList", indVarList)

    print("Getting "+str(M*T)+" samples from Ideal")
    total_weight = ideal.weight
    AllPsamples = ideal.getSolutionFromIdeal(M*T)

    lower_threshold_probability = -int(log(total_weight, 2))

    for i in range(T):
        seed += 1

        Psamples = AllPsamples[i*M:(i+1)*M]
        Qsamples = getSolutionFromSampler(
            seed, inputFile, M, samplerType, indVarList)

        projectedQsamples = []
        for sample in Qsamples:
            projectedQsample = []
            for s in sample:
                if abs(s) in UserIndVarList:
                    projectedQsample.append(s)
            projectedQsamples.append(projectedQsample)

        BucketedQsamples = {}
        bucketsofQ = []
        for sample in projectedQsamples:
            bucketID = -bucketof(weight_map, sample,
                                 UserIndVarList) - lower_threshold_probability
            assert (bucketID >= 0)
            if bucketID <= K:
                BucketedQsamples[bucketID] = sample
                bucketsofQ.append(bucketID)

        BucketedPsamples = {}
        bucketsofP = []
        for sample in Psamples:
            bucketID = -bucketof(weight_map, sample,
                                 UserIndVarList) - lower_threshold_probability
            assert (bucketID >= 0)
            if bucketID <= K:
                BucketedPsamples[bucketID] = sample
                bucketsofP.append(bucketID)

        commonBuckets = list(set(bucketsofP).intersection(set(bucketsofQ)))

        collisions = len(commonBuckets)
        print("# of collisions ", collisions)

        if collisions == 0:
            print("No collisions in round", i)
            continue

        # R = hoeffding(L,H,delta/(2*collisions*T))

        for bucketID in commonBuckets:

            Psample = BucketedPsamples[bucketID]
            Qsample = BucketedQsamples[bucketID]

            if Psample == Qsample:
                continue

            # the weight of Psample in P
            P_p = weightof(weight_map, Psample, UserIndVarList)
            # the weightof of Qsample in P
            P_q = weightof(weight_map, Qsample, UserIndVarList)

            H = P_p/(P_p+P_q*(1+2*eps/(1-eps)))
            L = P_p/(P_p+P_q*(1+alpha))

            R = hoeffding(L, H, delta/(4*collisions*T))
            print("R #of PAIRCOND queries "+str(R))

            if (totalSolsGenerated + R > maxSamples):
                print("Looking for more than 10**7 solutions", R)
                print("too many to ask ,so quitting here")
                print("NumSol:", totalSolsGenerated)
                exit(1)

            tempIndVarList = constructKernel(UserInputFile, tempFile, Qsample,
                                             Psample, R, UserIndVarList)
            samplinglist = list(chainform.Transform(
                tempFile, tempFile, 2))  # precision set to 4
            # print("file was constructed with these", Qsample, Psample)
            # print("samplingList:", samplinglist)
            solList = getSolutionFromSampler(
                seed, tempFile, R, samplerType, samplinglist)
            totalSolsGenerated += R

            seed += 1
            c_hat = biasFind(Psample, solList, UserIndVarList)

            cmd = "rm " + tempFile
            os.system(cmd)
            print("chat", c_hat)
            print("thresh", (H+L)/2)
            if c_hat < (H+L)/2:
                breakExperiment = True
                print("Rejected at iteration: ", i)
                print("NumSol total", totalSolsGenerated)
                return 0

        print("Accept on round: ", i)

    print("All rounds passed")
    print("NumSol total", totalSolsGenerated)

    return 1


def outsideBucket(K, theta, delta, UserInputFile, inputFile, samplerType, indVarList, UserIndVarList, weight_map, ideal, seed):

    numSamp = int(ceil(max(4*(K+1), 8*log(4/delta))/(theta)**2))

    print(4*(K+1)/theta**2, 8*log(4/delta)/theta**2)
    print("Number of samples required for OutBucket ", numSamp)

    Psamples = ideal.getSolutionFromIdeal(numSamp)
    Qsamples = getSolutionFromSampler(
        seed, inputFile, numSamp, samplerType, indVarList)

    total_weight = ideal.weight
    lower_threshold_probability = int(log(total_weight, 2))

    projectedQsamples = []
    for sample in Qsamples:
        projectedQsample = []
        for s in sample:
            if abs(s) in UserIndVarList:
                projectedQsample.append(s)
        projectedQsamples.append(projectedQsample)

    empirical_bucketP = [0]*(K+2)
    for i in Psamples:
        bucketID = -bucketof(weight_map, i, UserIndVarList) + \
            lower_threshold_probability
        assert (bucketID >= 0)
        if bucketID <= K:
            empirical_bucketP[bucketID] += 1
        else:
            empirical_bucketP[K+1] += 1

    empirical_bucketQ = [0]*(K+2)
    for i in projectedQsamples:
        bucketID = -bucketof(weight_map, i, UserIndVarList) + \
            lower_threshold_probability
        assert (bucketID >= 0)
        if bucketID <= K:
            empirical_bucketQ[bucketID] += 1
        else:
            empirical_bucketQ[K+1] += 1

    print("Bucket of P", empirical_bucketP)
    print("Bucket of Q", empirical_bucketQ)

    if (len(projectedQsamples) != len(Psamples)):
        print("Error: The length of two sample sets is not equal")
        exit(0)

    dhat = 0.0
    for i in range(K+2):
        dhat += 0.5*abs(empirical_bucketP[i]-empirical_bucketQ[i])/numSamp
    print("dhat: ", dhat)
    print("NumSol Outside", numSamp)

    return dhat, numSamp


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
    for i in sampleSol:
        if i != 0:
            if abs(i) not in indVarList:
                continue

            sampleVal[abs(i)] = int(i/abs(i))

    # which variables are in pos/neg value in the uniform sample
    unifVal = {}
    diffIndex = -1
    for j in unifSol:
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
        # f.write("c -- old CNF below -- \n")
        f.write(shiftedCNFStr)

    # print("New file: ", tempFile)
    # exit(0)

    return False, tempIndVarList, oldIndVarList


class dist_P:
    def __init__(self, inputFile) -> None:
        self.sampler


class IdealSampleRetriever:
    def __init__(self, inputFile):
        self.sampler = samp(cnfFile=inputFile)
        self.sampler.compile()
        self.sampler.parse()
        self.weight = self.sampler.annotate()

    def getSolutionFromIdeal(self, numSolutions):
        return self.getSolutionFromWAPS(numSolutions)

    def getSolutionFromWAPS(self, numSolutions):
        samples = self.sampler.sample(totalSamples=numSolutions)
        solList = list(samples)
        solList = [i.strip().split() for i in solList]
        solList = [[int(x) for x in i] for i in solList]
        return solList


class ChainFormulaSetup:
    def __init__(self, countList, newVarList, indicatorLits):
        self.countList = countList
        self.newVarList = newVarList
        self.indicatorLits = indicatorLits


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

    sampleSol = sampleSol[0]
    unifSol = unifSol[0]

    # adding more chain formulas (at most 8 in total: 3 + 5)
    # these chain formulas will have 31 solutions over 6 variables
    lenSol = len(sampleSol)
    for i in range(min(int(log(numSolutions, 2))+4, lenSol-3, 5)):
        countList.append(31)
        newVarList.append(6)
    assert len(countList) == len(newVarList)

    # picking selector literals, i.e. k1, k2, k3, randomly
    sampleLitList = random.sample(sampleSol, len(countList))
    unifLitList = []
    unifSolMap = unifSol
    for lit in sampleLitList:
        unifLitList.append(unifSolMap[abs(int(lit))-1])

    assert len(unifLitList) == len(sampleLitList)
    for a, b in zip(unifLitList, sampleLitList):
        assert abs(int(a)) == abs(int(b))

    indicatorLits = []
    indicatorLits.append(sampleLitList)
    indicatorLits.append(unifLitList)

    # print("countList:", countList)
    # print("newVarList:", newVarList)
    # print("indicatorLits:", indicatorLits)
    return ChainFormulaSetup(countList, newVarList, indicatorLits)


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

        if line[0] == "w":
            assert len(line) == 3
            var = int(line[1])
            weight = float(line[2])
            assert (0 <= weight <= 1)
            max_var = max(var, max_var)
            continue

        cls += 1
        for l in line:
            var = abs(int(l))
            max_var = max(var, max_var)

    if max_var > given_vars:
        print("ERROR: Number of variables given is LESS than the number of variables ued")
        print("ERROR: Vars in header: %d   max var: %d" %
              (given_vars, max_var))
        return False

    if cls != given_cls:
        print("ERROR: Number of clauses in header is DIFFERENT than the number of clauses in the CNF")
        print("ERROR: Claues in header: %d   clauses: %d" % (given_cls, cls))
        return False

    return True

# Returns bias wrt sample


def biasFind(sample, solList, indVarList):
    solMap = {}
    numSolutions = len(solList)
    for sol in solList:
        solution = ''
        solFields = sol
        for entry in solFields:
            if abs(entry) in indVarList:
                solution += str(entry) + " "

        if solution in solMap.keys():
            solMap[solution] += 1
        else:
            solMap[solution] = 1

    if not bool(solMap):
        print("No Solutions were given to the test")
        exit(1)

    print("c Printing counts of each sample: ", list(solMap.values()))

    solution = ""

    for i in solList[0]:
        if abs(i) in indVarList:
            solution += str(i) + " "

    # print("solList[0]", solList[0])
    # print("sample", sample)

    if (len(set(solList[0]).intersection(set(sample))) == len(sample)):
        return solMap.get(solution, 0)*1.0/numSolutions
    else:
        return 1.0-(solMap.get(solution, 0)*1.0/numSolutions)


class SolutionRetriever:

    @staticmethod
    def getSolutionFromSampler(inputFile, numSolutions, samplerType, indVarList, verbosity, newSeed):
        topass_withseed = (inputFile, numSolutions,
                           indVarList, verbosity, newSeed)
        ok = check_cnf(inputFile)
        if not ok:
            print(
                "ERROR: CNF is malformatted. Sampler may give wrong solutions in this case. Exiting.")
            print("File is: %s" % inputFile)
            exit(-1)

        print("Using sampler: %s" % get_sampler_string(samplerType))
        if (samplerType == SAMPLER_UNIGEN3):
            sols = SolutionRetriever.getSolutionFromUniGen3(*topass_withseed)

        elif (samplerType == SAMPLER_APPMC3):
            sols = SolutionRetriever.getSolutionFromAppMC3(*topass_withseed)

        elif (samplerType == SAMPLER_QUICKSAMPLER):
            sols = SolutionRetriever.getSolutionFromQuickSampler(
                *topass_withseed)

        elif (samplerType == SAMPLER_STS):
            sols = SolutionRetriever.getSolutionFromSTS(*topass_withseed)

        elif (samplerType == SAMPLER_CMS):
            sols = SolutionRetriever.getSolutionFromCMSsampler(
                *topass_withseed)

        elif (samplerType == SAMPLER_SPUR):
            sols = SolutionRetriever.getSolutionFromSpur(*topass_withseed)

        else:
            print("Error: No such sampler!")
            exit(-1)

        print("Number of solutions returned by sampler:", len(sols))
        if verbosity:
            print("Solutions:", sols)

        return sols

    @staticmethod
    def getSolutionFromUniGen3(inputFile, numSolutions, indVarList, verbosity, newSeed):
        # must construct ./unigen --samples=500 --verbosity=0 --threads=1  CNF-FILE SAMPLESFILE
        inputFileSuffix = inputFile.split('/')[-1][:-4]
        tempOutputFile = tempfile.gettempdir()+'/'+inputFileSuffix+".txt"

        cmd = './samplers/unigen --samples='+str(numSolutions)
        cmd += ' ' + inputFile + ' ' + \
            str(tempOutputFile) + ' > /dev/null 2>&1'
        if verbosity:
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
                    sol = line.split(':')[0].replace(
                        'v', '').strip().split()[:-1]
                    sol = [int(i) for i in sol]
                    solList.append(sol)

        solreturnList = solList
        if (len(solList) > numSolutions):
            solreturnList = random.sample(solList, numSolutions)

        os.unlink(str(tempOutputFile))
        return solreturnList

    @staticmethod
    def getSolutionFromAppMC3(inputFile, numSolutions, indVarList, verbosity, newSeed):
        # must construct: ./approxmc3 -s 1 -v2 --sampleout /dev/null --samples 500
        inputFileSuffix = inputFile.split('/')[-1][:-4]
        tempOutputFile = tempfile.gettempdir()+'/'+inputFileSuffix+".txt"

        cmd = './samplers/approxmc3 -s ' + \
            str(newSeed) + ' -v 0 --samples ' + str(numSolutions)
        cmd += ' --sampleout ' + str(tempOutputFile)
        cmd += ' ' + inputFile + ' > /dev/null 2>&1'
        if verbosity:
            print("cmd: ", cmd)
        os.system(cmd)

        with open(tempOutputFile, 'r') as f:
            lines = f.readlines()

        solList = []
        for line in lines:
            line = line.strip()
            freq = int(line.split(':')[0])
            for i in range(freq):
                sol = line.split(':')[1].strip().split()[:-1]
                sol = [int(i) for i in sol]
                solList.append(sol)

        solreturnList = solList
        if len(solList) > numSolutions:
            solreturnList = random.sample(solList, numSolutions)

        os.unlink(str(tempOutputFile))
        return solreturnList

    @staticmethod
    def getSolutionFromQuickSampler(inputFile, numSolutions, indVarList, verbosity, newSeed):
        cmd = "./samplers/quicksampler -n " + \
            str(numSolutions*5)+' '+str(inputFile)+' > /dev/null 2>&1'
        if verbosity:
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
            sol = []
            i = 0
            # valutions are 0 and 1 and in the same order as c ind.
            for x in list(fields[1].strip()):
                if (x == '0'):
                    sol.append(-1*indVarList[i])
                else:
                    sol.append(indVarList[i])
                i += 1
            solList.append(sol)

        solreturnList = solList
        if len(solList) > numSolutions:
            solreturnList = random.sample(solList, numSolutions)
        elif len(solreturnList) < numSolutions:
            print("Error: Quicksampler did not find required number of solutions")
            exit(1)

        os.unlink(inputFile+'.samples')
        os.unlink(inputFile+'.samples.valid')

        return solreturnList

    @staticmethod
    def getSolutionFromSpur(inputFile, numSolutions, indVarList, verbosity, newSeed):
        inputFileSuffix = inputFile.split('/')[-1][:-4]
        tempOutputFile = tempfile.gettempdir()+'/'+inputFileSuffix+".out"
        cmd = './samplers/spur -seed %d -q -s %d -out %s -cnf %s' % (
            newSeed, numSolutions, tempOutputFile, inputFile)
        if verbosity:
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
            if (not (startParse)):
                continue
            if (line.startswith('#END_SAMPLES')):
                startParse = False
                continue
            fields = line.strip().split(',')

            solCount = int(fields[0])
            sol = []

            i = 1
            for x in list(fields[1]):
                if (x == '0'):
                    sol.append(-i)
                else:
                    sol.append(i)
                i += 1

            for i in range(solCount):
                solList.append(sol)

        os.unlink(tempOutputFile)
        return solList

    @staticmethod
    def getSolutionFromSTS(inputFile, numSolutions, indVarList, verbosity, newSeed):
        kValue = 50
        samplingRounds = numSolutions/kValue + 1
        inputFileSuffix = inputFile.split('/')[-1][:-4]
        outputFile = tempfile.gettempdir()+'/'+inputFileSuffix+".out"
        cmd = './samplers/STS -k=' + \
            str(kValue)+' -nsamples='+str(samplingRounds)+' '+str(inputFile)
        cmd += ' > '+str(outputFile)
        if verbosity:
            print("cmd: ", cmd)
        os.system(cmd)

        with open(outputFile, 'r') as f:
            lines = f.readlines()

        solList = []
        shouldStart = False
        for j in range(len(lines)):
            if (lines[j].strip() == 'Outputting samples:' or lines[j].strip() == 'start'):
                shouldStart = True
                continue
            if (lines[j].strip().startswith('Log') or lines[j].strip() == 'end'):
                shouldStart = False
            if (shouldStart):

                i = 0
                sol = []
                # valutions are 0 and 1 and in the same order as c ind.
                for x in list(lines[j].strip()):
                    if (x == '0'):
                        sol.append(-indVarList[i])
                    else:
                        sol.append(indVarList[i])
                    i += 1
                solList.append(sol)

        if len(solList) < numSolutions:
            print(len(solList))
            print("STS Did not find required number of solutions")
            sys.exit(1)
        elif len(solList) > numSolutions:
            solList = random.sample(solList, numSolutions)

        os.unlink(outputFile)
        return solList

    @staticmethod
    def getSolutionFromCMSsampler(inputFile, numSolutions, indVarList, verbosity, newSeed):
        inputFileSuffix = inputFile.split('/')[-1][:-4]
        outputFile = tempfile.gettempdir()+'/'+inputFileSuffix+".out"
        cmd = "./samplers/cryptominisat5 --restart luby --maple 0 --verb 10 --nobansol"
        cmd += " --scc 1 -n1 --presimp 0 --polar rnd --freq 0.9999"
        cmd += " --random " + str(newSeed) + " --maxsol " + str(numSolutions)
        cmd += " " + inputFile
        cmd += " --dumpresult " + outputFile + " > /dev/null 2>&1"

        if verbosity:
            print("cmd: ", cmd)
        os.system(cmd)

        with open(outputFile, 'r') as f:
            lines = f.readlines()

        solList = []
        for line in lines:
            if line.strip() == 'SAT':
                continue

            sol = []
            lits = line.split(" ")
            for y in indVarList:
                if str(y) in lits:
                    sol.append(y)

                if "-" + str(y) in lits:
                    sol.append(-y)
            solList.append(sol)

        solreturnList = solList
        if len(solList) > numSolutions:
            solreturnList = random.sample(solList, numSolutions)
        if len(solList) < numSolutions:
            print("cryptominisat5 Did not find required number of solutions")
            sys.exit(1)
        os.unlink(outputFile)
        return solreturnList


class cnf_test:
    def __init__(self, samplerType, inputFile, maxSamples):
        inputFileSuffix = inputFile.split('/')[-1][:-4]
        self.tempFile = tempfile.gettempdir() + "/" + inputFileSuffix+"_t.cnf"
        self.indVarList = parseIndSupport(inputFile)
        self.inputFile = inputFile
        self.samplerType = samplerType
        self.maxSamples = maxSamples
        self.thresholdSolutions = 0
        self.totalSamplesGenerated = 0
        self.samplerString = get_sampler_string(samplerType)

    def CM_test(self, epsilon, eta, delta, verbosity, seed):

        totalLoops = int(ceil(log(2.0/(eta+2*epsilon), 2))+1)
        listforTraversal = range(totalLoops, 0, -1)

        for j in listforTraversal:
            tj = ceil(pow(2, j)*(2*epsilon+eta)/((eta-2*epsilon)**2)
                      * log(4.0/(eta+2*epsilon), 2)*(4*e/(e-1)*log(1.0/delta)))
            beta = (pow(2, j-1)+1)*(eta + 2*epsilon)*1.0 / \
                (4+(2*epsilon+eta)*(pow(2, j-1) - 1))
            gamma = (beta-2*epsilon)/4
            constantFactor = ceil(1/(8.79*gamma*gamma))
            boundFactor = log((16)*(e/(e-1))*(1/(delta*(eta-2*epsilon)**2))
                              * log(4/(eta+2*epsilon), 2)*log(1/delta), 2)

            if verbosity:
                print("constantFactor:{:<4} boundFactor: {:<20} logBoundFactor: {:<20}".format(
                    constantFactor, boundFactor, log(boundFactor, 2)))
                print("tj: {:<6} totalLoops: {:<5} beta: {:<10} epsilon: {:<10}".format(
                    tj, totalLoops, beta, epsilon))

            self.numSolutions = int(ceil(constantFactor*boundFactor))

            loThresh = (1-(beta+2*epsilon)/2)/2
            hiThresh = (1+(beta+2*epsilon)/2)/2
            print("numSolutions: {:<5} loThresh:{:<6} hiThresh: {:<6}".format(
                self.numSolutions, loThresh, hiThresh))

            breakExperiment = False

            for i in range(int(tj)):
                breakExperiment = self.CM_inner_loop(
                    j, i, tj, loThresh, hiThresh, verbosity)

                if breakExperiment:
                    break

            if breakExperiment:
                break

        if not breakExperiment:
            print("Accept:1 TotalSamplesGenerated:{0} ".format(
                self.totalSamplesGenerated))

    def CM_inner_loop(self, j, i, tj, loThresh, hiThresh, verbosity):
        self.thresholdSolutions += self.numSolutions

        # generate a new seed value for every different (i,j,experiment)
        newSeed = int(i*tj+j)
        # get sampler's solutions
        sampleSol = SolutionRetriever.getSolutionFromSampler(
            self.inputFile, 1, self.samplerType, self.indVarList, verbosity, newSeed)
        self.totalSamplesGenerated += 1

        # get uniform sampler's solutions
        unifSol = SolutionRetriever.getSolutionFromSampler(
            self.inputFile, 1, SAMPLER_SPUR, self.indVarList, verbosity, newSeed)
        self.totalSamplesGenerated += 1

        chainFormulaConf = chainFormulaSetup(
            sampleSol, unifSol, self.numSolutions)
        identical_sol, tempIndVarList, oldIndVarList = constructNewCNF(
            self.inputFile, self.tempFile, sampleSol[0], unifSol[0],
            chainFormulaConf, self.indVarList)

        # the two solutions were the same, unbiased
        if identical_sol:
            return False

        # seed update
        newSeed = newSeed + 1

        # get sampler's solutions
        solList = SolutionRetriever.getSolutionFromSampler(
            self.tempFile, self.numSolutions, self.samplerType, tempIndVarList, verbosity, newSeed)
        os.unlink(self.tempFile)
        self.totalSamplesGenerated += self.numSolutions

        isUniform = True
        bias = biasFind(unifSol[0], solList, oldIndVarList)

        if loThresh <= bias <= hiThresh:
            isUniform = True
        else:
            isUniform = False

        print("sampler: {:<8s} i: {:<4d} isUniform: {:<4d} TotalSolutionsGenerated: {:<6d}".format(
            self.samplerString, i, isUniform,
            self.totalSamplesGenerated))

        if not isUniform:
            print("RejectIteration:{0}  Loop:{1} TotalSolutionsGenerated:{2} ".format(
                i, j, self.totalSamplesGenerated))

            return True

        if self.thresholdSolutions > self.maxSamples:
            return True

        return False

    def PM_test(inputFile, samplerType, eta, epsilon, delta, maxSamples, verbosity, seed):
        UserInputFile = args.input
        print("This is the user input:--", UserInputFile)

        inputFilePrefix = UserInputFile.split("/")[-1][:-4]
        inputFile = inputFilePrefix + "."+str(samplerType)+".cnf"

        print("This is the output file after weighted to unweighted:", inputFile)

        UserIndVarList = parseIndSupport(UserInputFile)
        indVarList = list(chainform.Transform(
            UserInputFile, inputFile, 2))  # precision set to 4

        weight_map = parseWeights(UserInputFile, UserIndVarList)

        totalSolsGenerated = 0

        numVars = len(UserIndVarList)
        K = int(ceil(numVars*log(2, 2)+log(100/eta, 2)))

        theta = eta/20
        print("K", K)

        ideal = IdealSampleRetriever(inputFile=UserInputFile)

        dhat, totalSolsGenerated = outsideBucket(
            K, theta, delta/2, UserInputFile, inputFile, samplerType, indVarList, UserIndVarList, weight_map, ideal, maxSamples, seed)

        if dhat - theta > epsilon/2:
            print("Rejected as dhat("+str(dhat)+") > eps/2 ("+str(epsilon/2)+") ")
        else:
            eps2 = dhat + theta
            insideBucket(K, epsilon, eps2, eta, delta/2, UserInputFile, inputFile, samplerType,
                         indVarList, UserIndVarList, weight_map, ideal, totalSolsGenerated, maxSamples, seed)
