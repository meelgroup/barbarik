import sys
import os 
import random
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
                assert(0 <= weight <= 1)
                max_var = max(var, max_var)
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


class SolutionRetriever:

    
    @staticmethod
    def getSolutionFromSampler(inputFile, numSolutions, samplerType, indVarList, verbosity, newSeed):
        topass_withseed = (inputFile, numSolutions, indVarList, verbosity, newSeed)
        ok = check_cnf(inputFile)
        if not ok:
            print("ERROR: CNF is malformatted. Sampler may give wrong solutions in this case. Exiting.")
            print("File is: %s" % inputFile)
            exit(-1)


        print("Using sampler: %s" % get_sampler_string(samplerType))
        if (samplerType == SAMPLER_UNIGEN):
            sols = SolutionRetriever.getSolutionFromUniGen(*topass_withseed)

        elif (samplerType == SAMPLER_APPMC3):
            sols = SolutionRetriever.getSolutionFromAppMC3(*topass_withseed)

        elif (samplerType == SAMPLER_QUICKSAMPLER):
            sols = SolutionRetriever.getSolutionFromQuickSampler(*topass_withseed)

        elif (samplerType == SAMPLER_STS):
            sols = SolutionRetriever.getSolutionFromSTS(*topass_withseed)

        elif (samplerType == SAMPLER_CMS):
            sols = SolutionRetriever.getSolutionFromCMSsampler(*topass_withseed)

        elif (samplerType == SAMPLER_SPUR):
            sols = SolutionRetriever.getSolutionFromSpur(*topass_withseed)

        else:
            print("Error: No such sampler!")
            exit(-1)

        # clean up the solutions
        for i in range(len(sols)):
            sols[i] = sols[i].strip()
            if sols[i].endswith(' 0'):
                sols[i] = sols[i][:-2]

        print("Number of solutions returned by sampler:", len(sols))
        if verbosity:
            print("Solutions:", sols)
        return sols

    @staticmethod
    def getSolutionFromUniGen(inputFile, numSolutions, indVarList, verbosity, newSeed):
        # must construct ./unigen --samples=500 --verbosity=0 --threads=1  CNF-FILE SAMPLESFILE
        inputFileSuffix = inputFile.split('/')[-1][:-4]
        tempOutputFile = tempfile.gettempdir()+'/'+inputFileSuffix+".txt"

        cmd = './samplers/unigen --samples='+str(numSolutions)
        cmd += ' ' + inputFile + ' ' + str(tempOutputFile) + ' > /dev/null 2>&1'
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
                    solList.append(line.split(':')[0].replace('v', '').strip())

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

        cmd = './samplers/approxmc3 -s ' + str(newSeed) + ' -v 0 --samples ' + str(numSolutions)
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
                solList.append(line.split(':')[1].strip())

        solreturnList = solList
        if len(solList) > numSolutions:
            solreturnList = random.sample(solList, numSolutions)

        os.unlink(str(tempOutputFile))
        return solreturnList

    @staticmethod
    def getSolutionFromQuickSampler(inputFile, numSolutions, indVarList, verbosity, newSeed):
        cmd = "./samplers/quicksampler -n "+str(numSolutions*5)+' '+str(inputFile)+' > /dev/null 2>&1'
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

        solreturnList = solList
        if len(solList) > numSolutions:
            solreturnList = random.sample(solList, numSolutions)

        os.unlink(inputFile+'.samples')
        os.unlink(inputFile+'.samples.valid')

        if len(solreturnList) != numSolutions:
            print("Did not find required number of solutions")
            sys.exit(1)

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
    def getSolutionFromSTS(inputFile, numSolutions, indVarList, verbosity, newSeed):
        kValue = 50
        samplingRounds = numSolutions/kValue + 1
        inputFileSuffix = inputFile.split('/')[-1][:-4]
        outputFile = tempfile.gettempdir()+'/'+inputFileSuffix+".out"
        cmd = './samplers/STS -k='+str(kValue)+' -nsamples='+str(samplingRounds)+' '+str(inputFile)
        cmd += ' > '+str(outputFile)
        if verbosity:
            print("cmd: ", cmd)
        os.system(cmd)

        with open(outputFile, 'r') as f:
            lines = f.readlines()

        solList = []
        shouldStart = False
        for j in range(len(lines)):
            if(lines[j].strip() == 'Outputting samples:' or lines[j].strip() == 'start'):
                shouldStart = True
                continue
            if (lines[j].strip().startswith('Log') or lines[j].strip() == 'end'):
                shouldStart = False
            if (shouldStart):

                i = 0
                sol = ''
                # valutions are 0 and 1 and in the same order as c ind.
                for x in list(lines[j].strip()):
                    if (x == '0'):
                        sol += ' -'+str(indVarList[i])
                    else:
                        sol += ' '+str(indVarList[i])
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

