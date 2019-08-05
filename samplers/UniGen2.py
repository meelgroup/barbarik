from __future__ import print_function
import sys
import os
import time
import math


def ensureDirectory(d):
    if not os.path.exists(d):
        os.makedirs(d)


def usage():
    out = ["Usage: python UniGen2.py [options] <inputFile> <outputFolder>\n",
           "Generated samples (and log files if enabled) are written to",
           " <outputFolder>.\n",
           "Options are entered as '-option=value' where 'option' is from the",
           " list below.\n",
           "If an option is not specified, the default value given below is",
           " used.\n",
           "\nBASIC OPTIONS:\n\n",
           "samples: number of samples to generate (default: 1) \n",
           "kappa: computed from desired tolerance as in Algorithm 1 of the",
           " TACAS-15 paper (default: 0.638, corresponding to epsilon=16)\n",
           "threads: number of threads to use for sampling (default: OpenMP",
           " system default)\n",
           "timeout: overall timeout in seconds (default: 72000)\n",
           "satTimeout: timeout for a single BoundedWeightSAT call in seconds",
           " (default: 3000)\n",
           "runIndex: number appended to names of output & log files (default:",
           " current time) \n",
           "\nADVANCED OPTIONS (not needed for typical use):\n\n",
           "logging: 0 to turn off logging and 1 to turn it on. If logging is",
           " turned on, the log file(s) will be present in",
           " <outputFolder>/logging\n",
           "aggregateSolutions: 0/1 to disable/enable merging samples from",
           " all threads (default: 1)\n",
           "writeSamples: 0/1 to disable/enable writing samples to a file",
           " (default: 1)\n",
           "multisampling: 0/1 to disable/enable returning multiple samples",
           " from each UniGen2 call (default: 1)\n",
           "startIteration: initial number of XOR clauses to use in UniGen2,",
           " or 0 if this should be computed using sharpSAT/ApproxMC, or",
           " -1 if it should be computed using ApproxMC, only (e.g. when",
           " doing projection counting) (default: 0)\n",
           "callsPerSolver: number of UniGen2 calls to make in a single solver",
           " (without clearing learned clauses), or 0 to use a built-in",
           " heuristic (default: 0)\n",
           "pivotAC: pivot value for ApproxMC (default: 60)\n",
           "tApproxMC: number of iterations for ApproxMC (default: 1)\n"]
    print("".join(out))
    sys.exit(1)


def getInputs():
    paramMap = {}
    error = ''
    gotInputFile = False
    acceptedParams = ['timeout', 'satTimeout', 'runIndex', 'samples',
                      'callsPerSolver', 'writeSamples', 'threads',
                      'aggregateSolutions', 'logging', 'kappa', 'multisampling',
                      'startIteration', 'pivotAC', 'tApproxMC']
    for i in range(1, len(sys.argv)):
        if sys.argv[i][0] != '-':
            if not gotInputFile:
                paramMap['inputFile'] = sys.argv[i]
                gotInputFile = True
            else:
                paramMap['outputFolder'] = sys.argv[i]
                action = 3
                return action, error, paramMap
        else:
            if sys.argv[i][1] == 'h':
                action = 0
                return action, error, paramMap
            fieldValues = sys.argv[i][1:].strip().split('=')
            if fieldValues[0] not in acceptedParams:
                action = 1
                error = "Could not understand the option " \
                        + str(fieldValues[0]) + "\n"
                return action, error, paramMap
            else:
                paramMap[fieldValues[0]] = fieldValues[1]
    action = 2
    error = ["You must specify an input file and an output folder.\n",
             "Use the option -h for help.\n"]
    return action, "".join(error), paramMap


def main():
    runIndex = str(int(time.time()))
    timeout = 72000
    satTimeout = 3000
    shouldLog = False
    logFilePrefix = ''
    kappa = 0.638
    samples = 1
    writeSamples = True
    multisampling = True
    threads = 0
    aggregateSolutions = True
    pivotAC = 60
    approxMCIterations = 1
    startIteration = 0
    callsPerSolver = 0
    action, error, paramMap = getInputs()
    if action == 0:
        usage()
        sys.exit(1)
    if action == 1 or action == 2:
        print(error)
        sys.exit(1)
    if 'runIndex' in paramMap:
        runIndex = int(paramMap['runIndex'])
    if 'samples' in paramMap:
        samples = int(paramMap['samples'])
    if 'callsPerSolver' in paramMap:
        callsPerSolver = int(paramMap['callsPerSolver'])
    if 'writeSamples' in paramMap:
        writeSamples = (paramMap['writeSamples'] == '1')
    if 'timeout' in paramMap:
        timeout = float(paramMap['timeout'])
    if 'satTimeout' in paramMap:
        satTimeout = float(paramMap['satTimeout'])
    if 'kappa' in paramMap:
        kappa = float(paramMap['kappa'])
    if 'outputFile' in paramMap:
        outputFileName = paramMap['outputFile']
    if 'logging' in paramMap:
        if paramMap['logging'] == '0':
            shouldLog = False
        if paramMap['logging'] == '1':
            shouldLog = True
    if 'multisampling' in paramMap:
        if paramMap['multisampling'] == '0':
            multisampling = False
        if paramMap['multisampling'] == '1':
            multisampling = True
    if 'aggregateSolutions' in paramMap:
        if paramMap['aggregateSolutions'] == '0':
            aggregateSolutions = False
        if paramMap['aggregateSolutions'] == '1':
            aggregateSolutions = True
    if 'threads' in paramMap:
        threads = int(paramMap['threads'])
    if 'pivotAC' in paramMap:
        pivotAC = int(paramMap['pivotAC'])
    if 'tApproxMC' in paramMap:
        approxMCIterations = int(paramMap['tApproxMC'])
    if 'startIteration' in paramMap:
        startIteration = int(paramMap['startIteration'])
    initialFileName = paramMap['inputFile']
    outputFolder = paramMap['outputFolder']
    initialFileNameSuffix = initialFileName.split('/')[-1][:-4]
    ensureDirectory(outputFolder)
    if shouldLog:
        ensureDirectory(outputFolder + os.sep + "logging" + os.sep)
        logFilePrefix = outputFolder + os.sep + "logging" + os.sep \
                        + str(initialFileNameSuffix) + '_' + str(runIndex)
    outputFileName = outputFolder + os.sep + str(initialFileNameSuffix) + '_' \
                     + str(runIndex) + ".txt"

    pivotUniGen = math.ceil(4.03 * (1 + 1 / kappa) * (1 + 1 / kappa))

    if startIteration < 0:
        startIteration = 0
        print("Doing projection counting; will use ApproxMC")
    elif startIteration == 0:
        startIteration = 0
        print("Attempting to compute startIteration with sharpSAT")
        countLogFile = outputFolder + os.sep + initialFileNameSuffix \
                       + "_" + str(runIndex) + ".count"
        cmd = "." + os.sep + "doalarm -t real 300 ./sharpSAT -q " \
              + initialFileName + " > " + countLogFile
        print(cmd)
        os.system(cmd)
        f = open(countLogFile, 'r')
        lines = f.readlines()
        f.close()
        failed = True
        if len(lines) > 1:
            if lines[1].strip() == "END":
                failed = False
                if sys.version_info > (3, 0):
                    count = int(lines[0].strip())
                else:
                    # noinspection PyCompatibility
                    count = long(lines[0].strip())
                if count == 0:
                    print("Unsatisfiable formula!")
                    return
                logCount = math.log(count, 2)
                startIteration = int(round(logCount + math.log(1.8, 2)
                                           - math.log(pivotUniGen, 2))) - 2
                if startIteration < 0:
                    print("The number of solutions is only " + str(count) + "\n"
                          + "The best technique is to just enumerate all the"
                          + " solutions and choose one")
                    return
                # print("Solution count estimate is %f * 2^%d"
                #       % (count / (2.0 ** int(logCount)), int(logCount)))
        if failed:
            print("sharpSAT failed or timed out; using ApproxMC instead")
        else:
            print("startIteration computed successfully")

    # Call to the main UniGen2 binary
    cmd = "." + os.sep + "unigen --samples=" + str(samples) \
          + " --kappa=" + str(kappa) + " --pivotUniGen=" + str(pivotUniGen) \
          + " --maxTotalTime=" + str(timeout)\
          + " --startIteration=" + str(startIteration) \
          + " --maxLoopTime=" + str(satTimeout) \
          + " --tApproxMC="+str(approxMCIterations) \
          + " --pivotAC=" + str(pivotAC) + " --gaussuntil=400 --verbosity=0"
    if shouldLog:
        cmd += " --logFile=" + logFilePrefix
    if callsPerSolver > 0:
        cmd += " --callsPerSolver="+str(callsPerSolver)
    if multisampling:
        cmd += " --multisample"
    if not writeSamples:
        cmd += " --nosolprint"
    if threads > 0:
        cmd += " --threads="+str(threads)
    if not aggregateSolutions:
        cmd += " --noAggregation"
    cmd += " " + initialFileName + " " + outputFileName
    print(cmd)
    sys.stdout.flush()
    os.system(cmd)


if __name__ == "__main__":
    main()
