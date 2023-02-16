#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import sys
from math import log, ceil, e
import random
import argparse

from interfaces.cnf import *

SAMPLER_UNIGEN3 = 1
SAMPLER_QUICKSAMPLER = 2
SAMPLER_STS = 3
SAMPLER_WAPS = 4
SAMPLER_CUSTOM = 5

FILE_CNF = 1

UNIF_TEST = 0
GEN_TEST = 1

if __name__ == "__main__":

    samplers = str(SAMPLER_UNIGEN3) + " for UniGen3\n"
    samplers += str(SAMPLER_QUICKSAMPLER) + " for QuickSampler\n"
    samplers += str(SAMPLER_STS)+ " for STS\n"
    samplers += str(SAMPLER_SPUR) + " for SPUR\n"
    samplers += str(SAMPLER_CMS) + " for CMS\n"

    idealsamplers = str(SAMPLER_WAPS) + " for WAPS\n"

    filetypes = str(FILE_CNF) + " for CNF\n"

    parser = argparse.ArgumentParser()
    parser.add_argument('--sampler', type=int, help=samplers, default=SAMPLER_SPUR, dest='sampler')
    parser.add_argument('--eta', type=float, help="default = 0.9", default=0.9, dest='eta')
    parser.add_argument('--epsilon', type=float, help="default = 0.3", default=0.3, dest='epsilon')
    parser.add_argument('--delta', type=float, help="default = 0.05", default=0.05, dest='delta')
    parser.add_argument('--ftype', type=int, help=filetypes, default=1, dest='ftype')
    parser.add_argument('--testtype', type=int, default=0, help="uniform(0) vs. general(1)", dest='testtype')
    parser.add_argument('--max', type=int, default=sys.maxsize, help="max samples", dest='max')
    parser.add_argument('--seed', type=int, required=True, dest='seed')
    parser.add_argument('--verb', type=int, dest='verbose')
    parser.add_argument("input", help="input file")

    args = parser.parse_args()
    inputFile = args.input
    ftype = args.ftype

    eta = args.eta
    epsilon = args.epsilon    
    delta = args.delta
    maxSamples = args.max

    testtype = args.testtype
    verbosity = args.verbose

    seed = args.seed
    random.seed(seed)
    
    #sanity checks
    if epsilon >= eta:
        print("epsilon must be less than eta")

    elif ftype == FILE_CNF:
        
        experiment = cnf_test(args.sampler, inputFile, maxSamples)

        if testtype == UNIF_TEST:
            experiment.CM_test(epsilon, eta, delta, verbosity, seed)

        elif testtype == GEN_TEST:
            P = ideal_sampler(inputFile)
            Q = weighted_sampler_under_test(args.samplerQ)
            PM_general_test(P,Q,eta,epsilon,delta,maxSamples,verbosity,seed)


def PM_test(inputFile,samplerType, eta,epsilon,delta,maxSamples,verbosity,seed):
    UserInputFile = args.input
    print("This is the user input:--", UserInputFile)

    inputFilePrefix = UserInputFile.split("/")[-1][:-4]
    inputFile = inputFilePrefix + "."+str(samplerType)+".cnf"

    print("This is the output file after weighted to unweighted:", inputFile)

    UserIndVarList = parseIndSupport(UserInputFile)
    indVarList = list(chainform.Transform(UserInputFile, inputFile, 2))  # precision set to 4

    weight_map = parseWeights(UserInputFile, UserIndVarList)

    totalSolsGenerated = 0

    numVars = len(UserIndVarList)
    K =  int(ceil(numVars*log(2,2)+log(100/eta,2)))

    theta = eta/20
    print("K", K)

    ideal = IdealSampleRetriever(inputFile = UserInputFile)

    dhat, totalSolsGenerated = outsideBucket(K,theta,delta/2,UserInputFile,inputFile,samplerType,indVarList,UserIndVarList,weight_map,ideal,maxSamples,seed)

    if dhat - theta  > epsilon/2:
        print("Rejected as dhat("+str(dhat)+") > eps/2 ("+str(epsilon/2)+") " )
    else:
        eps2 = dhat + theta
        insideBucket(K,epsilon,eps2,eta,delta/2,UserInputFile,inputFile,samplerType,indVarList,UserIndVarList,weight_map,ideal,totalSolsGenerated,maxSamples,seed)
