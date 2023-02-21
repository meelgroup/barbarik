#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import sys
from math import e
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
    samplers += str(SAMPLER_STS) + " for STS\n"
    samplers += str(SAMPLER_SPUR) + " for SPUR\n"
    samplers += str(SAMPLER_CMS) + " for CMS\n"

    idealsamplers = str(SAMPLER_WAPS) + " for WAPS\n"

    filetypes = str(FILE_CNF) + " for CNF\n"

    parser = argparse.ArgumentParser()
    parser.add_argument('--sampler', type=int, help=samplers,
                        default=SAMPLER_SPUR, dest='sampler')
    parser.add_argument('--eta', type=float,
                        help="default = 0.9", default=0.9, dest='eta')
    parser.add_argument('--epsilon', type=float,
                        help="default = 0.3", default=0.3, dest='epsilon')
    parser.add_argument('--delta', type=float,
                        help="default = 0.05", default=0.05, dest='delta')
    parser.add_argument('--ftype', type=int, help=filetypes,
                        default=1, dest='ftype')
    parser.add_argument('--testtype', type=int, default=0,
                        help="uniform(0) vs. general(1)", dest='testtype')
    parser.add_argument('--max', type=int, default=sys.maxsize,
                        help="max samples", dest='max')
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

    # sanity checks
    if epsilon >= eta:
        print("epsilon must be less than eta")

    elif ftype == FILE_CNF:

        experiment = cnf_test(args.sampler, inputFile, eta, epsilon, delta, maxSamples,verbosity)

        if testtype == UNIF_TEST:
            experiment.CM_test(seed)

        elif testtype == GEN_TEST:
            experiment.PM_test(seed)
