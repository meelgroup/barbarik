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
from math import log, ceil, e
import random
import argparse

from interfaces.cnf import *

SAMPLER_UNIGEN3 = 1
SAMPLER_QUICKSAMPLER = 2
SAMPLER_STS = 3
SAMPLER_CUSTOM = 4

def CM_unif_test(inputFile,eta,epsilon,delta,maxSamples,verbosity,seed):
    
    if (eta < 2*epsilon):
        print("Eta needs to be at least two times epsilon")
        exit(1)
    
    totalLoops = int(ceil(log(2.0/(eta+2*epsilon), 2))+1)
    listforTraversal = range(totalLoops, 0, -1)

    exp = unif_Experiment(
    maxSamples=maxSamples, inputFile=inputFile,
    samplerType=args.sampler)

    breakExperiment = False
    
    exp.totalSolutionsGenerated = 0
    exp.totalUniformSamples = 0
    exp.thresholdSolutions = 0

    for j in listforTraversal:
        tj = ceil(pow(2, j)*(2*epsilon+eta)/((eta-2*epsilon)**2)*log(4.0/(eta+2*epsilon), 2)*(4*e/(e-1)*log(1.0/delta)))
        beta = (pow(2, j-1)+1)*(eta + 2*epsilon)*1.0/(4+(2*epsilon+eta)*(pow(2, j-1) - 1))
        gamma = (beta-2*epsilon)/4
        constantFactor = ceil(1/(8.79*gamma*gamma))
        boundFactor = log((16)*(e/(e-1))*(1/(delta*(eta-2*epsilon)**2))*log(4/(eta+2*epsilon), 2)*log(1/delta), 2)
        
        if verbosity:
            print("constantFactor:{:<4} boundFactor: {:<20} logBoundFactor: {:<20}".format(
                constantFactor, boundFactor, log(boundFactor, 2)))
            print("tj: {:<6} totalLoops: {:<5} beta: {:<10} epsilon: {:<10}".format(
                tj, totalLoops, beta, epsilon))

        exp.numSolutions = int(ceil(constantFactor*boundFactor))
        exp.loThresh = int((exp.numSolutions*1.0/2)*(1-(beta+2*epsilon)/2))
        exp.hiThresh = int((exp.numSolutions*1.0/2)*(1+(beta+2*epsilon)/2))
        print("numSolutions: {:<5} loThresh:{:<6} hiThresh: {:<6}".format(
            exp.numSolutions, exp.loThresh, exp.hiThresh))

        breakExperiment = False

        for i in range(int(tj)):
            breakExperiment = exp.one_experiment(j, i, tj, verbosity)

            if breakExperiment:
                break

        if breakExperiment:
            break

    if not breakExperiment:
        print("exp:{2} Accept:1 TotalSolutionsGenerated:{0} TotalUniformSamples:{1}".format(
            exp.totalSolutionsGenerated,
            exp.totalUniformSamples))

def PM_test(inputFile,eta,epsilon,delta,maxSamples,verbosity,seed):
    UserInputFile = args.input
    print("This is the user input:--", UserInputFile)

    samplerType=args.sampler

    inputFilePrefix = UserInputFile.split("/")[-1][:-4]
    inputFile = inputFilePrefix + "."+str(samplerType)+".cnf"

    print("This is the output file after weighted to unweighted:", inputFile)

    UserIndVarList = parseIndSupport(UserInputFile)
    indVarList = list(chainform.Transform(UserInputFile, inputFile, 2))  # precision set to 4

    samplerString = ""

    if samplerType == SAMPLER_UNIGEN3:
        samplerString = "UniGen3"
    if samplerType == SAMPLER_QUICKSAMPLER:
        samplerString = "QuickSampler"
    if samplerType == SAMPLER_STS:
        samplerString = "STS"
    if samplerType == SAMPLER_CUSTOM:
        samplerString = "CustomSampler"

    weight_map = parseWeights(UserInputFile, UserIndVarList)

    totalSolsGenerated = 0

    numVars = len(UserIndVarList)
    K =  int(ceil(numVars*log(2,2)+log(100/eta,2)))

    theta = eta/20
    print("K", K)

    ideal = IdealSampleRetriever(inputFile = UserInputFile)

    dhat, totalSolsGenerated = outsideBucket(K,theta,delta/2,UserInputFile,inputFile,samplerType,indVarList,UserIndVarList,weight_map,ideal,seed)

    if dhat - theta  > epsilon/2:
        print("Rejected as dhat("+str(dhat)+") > eps/2 ("+str(epsilon/2)+") " )
    else:
        eps2 = dhat + theta
        insideBucket(K,epsilon,eps2,eta,delta/2,UserInputFile,inputFile,samplerType,indVarList,UserIndVarList,weight_map,ideal,totalSolsGenerated,seed)

if __name__ == "__main__":

    samplers = str(SAMPLER_UNIGEN3) + " for UniGen3\n"
    samplers += str(SAMPLER_QUICKSAMPLER) + " for QuickSampler\n"
    samplers += str(SAMPLER_STS)+ " for STS\n"
    samplers += str(SAMPLER_SPUR) + " for SPUR\n"
    samplers += str(SAMPLER_CMS) + " for CMS\n"

    parser = argparse.ArgumentParser()
    parser.add_argument('--sampler', type=int, help=samplers, default=SAMPLER_STS, dest='sampler')
    parser.add_argument('--eta', type=float, help="default = 0.9", default=0.9, dest='eta')
    parser.add_argument('--epsilon', type=float, help="default = 0.3", default=0.3, dest='epsilon')
    parser.add_argument('--delta', type=float, help="default = 0.05", default=0.05, dest='delta')
    parser.add_argument('--testtype', type=int, default=0, help="uniform(0) vs. general(1)", dest='testtype')
    parser.add_argument('--maxSamples', type=int, default=sys.maxsize, help="max samples", dest='maxSamples')
    parser.add_argument('--seed', type=int, required=True, dest='seed')
    parser.add_argument('--verb', type=int, dest='verbose')
    parser.add_argument("input", help="input file")

    args = parser.parse_args()
    inputFile = args.input

    eta = args.eta
    epsilon = args.epsilon    
    delta = args.delta
    maxSamples = args.maxSamples

    testtype = args.testtype
    verbosity = args.verbose

    seed = args.seed
    random.seed(seed)
    
    #change this default to weighted
    if testtype == 0:
        CM_unif_test(inputFile,eta,epsilon,delta,maxSamples,verbosity,seed)
    else:
        PM_test(inputFile,eta,epsilon,delta,maxSamples,verbosity,seed)