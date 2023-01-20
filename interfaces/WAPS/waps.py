# /***********[waps.py]
# Copyright (c) 2018 Rahul Gupta, Shubham Sharma, Subhajit Roy, Kuldeep Meel
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# ***********/

import argparse
import os
import pickle
import random
import re
import sys
import time
from fractions import Fraction

import numpy as np
from gmpy2 import mpq,mpfr
#import pkg_resources
random.seed(100)
def random_assignment(solution, samplingSet=[], weights=None):
    """Takes total number of variables and a partial assignment
    to return a complete assignment"""
    literals = set()
    solutionstr = ""
    vartolit = {}
    for x in solution.split():
        literals.add(abs(int(x)))
        vartolit[abs(int(x))] = int(x)
    for i in samplingSet:
        if i not in literals:
            if weights and i in weights:
                litchoice = np.random.choice([1, -1], p=[weights[i], 1 - weights[i]])
                solutionstr += str(litchoice * i) + " "
            else:
                solutionstr += str(((random.randint(0, 1) * 2) - 1) * i) + " "
        else:
            solutionstr += str(vartolit[i]) + " "
    return solutionstr


def normalizeWeights(weights):
    """
    Normalizes the weights
    Assumes that either weights defines all positive literal weights between 0 and 1 or
    defines weights for both positive and negative literals.
    """
    for key in weights:
        if -1 * key in weights:
            if key > 0:
                weights[key] = mpq(Fraction(weights[key])) / (
                    mpq(Fraction(weights[key])) + mpq(Fraction(weights[key * -1]))
                )
                del weights[-1 * key]
            else:
                weights[-1 * key] = mpq(Fraction(-1 * key)) / (
                    mpq(Fraction(weights[key])) + mpq(Fraction(weights[key * -1]))
                )
                del weights[key]


def fetchWeights(weightFile):
    """either specify all positive literal weights between 0 and 1 or
    specify weights for both positive and negative literals."""
    data = open(weightFile).read()
    lines = data.strip().split("\n")
    weights = {}
    for line in lines:
        if int(line.split(",")[0]) * (-1) in weights:
            if int(line.split(",")[0]) > 0:
                weights[int(line.split(",")[0])] = mpq(Fraction(line.split(",")[1])) / (
                    weights.pop(int(line.split(",")[0]) * (-1), None)
                    + mpq(Fraction(line.split(",")[1]))
                )
            else:
                weights[int(line.split(",")[0]) * (-1)] = weights[
                    int(line.split(",")[0]) * (-1)
                ] / (
                    weights[int(line.split(",")[0]) * (-1)]
                    + mpq(Fraction(line.split(",")[1]))
                )
        else:
            weights[int(line.split(",")[0])] = mpq(Fraction(line.split(",")[1]))
    print(weights)
    return weights


def conditionWeights(lits, weights):
    """Modifies the weight of positive literal as per condition given by list lits"""
    for lit in lits:
        weights[int(lit)] = 1
        weights[-1 * int(lit)] = 0


class Node:
    def __init__(self, label=None, children=[], decision=None):
        self.label = label
        self.children = children
        self.models = 1
        self.decisionat = decision
        self.weight = mpq("1")


class sampler:
    def __init__(
        self,
        cnfFile=None,
        DIMACScnf=None,
        dDNNFfile=None,
        samplingSet=[],
        randAssign=True,
        weights={},
        conditionVars=None,
    ):
        """
        This class contains functions and handles data for sampling.

        :param cnfFile: The path of file containing DIMACS cnf
        :param DIMACScnf: cnf available in DIMACS format
        :param dDNNFfile: specifies file containing d-DNNF to bypass compilation phase
        :param samplingSet: variables on which samples are projected, sampling is done on all variables by default
        :param randAssign: extend each sample to contain all variables from samplingSet
        :param weights: dictionary mapping literals to their weights
        :param conditionVars: list specifying list of literals which should be true in samples

        """
        self._cnfFile = cnfFile
        self._dDNNFfile = dDNNFfile
        self.randAssign = randAssign
        self.weights = weights
        self.conditionVars = conditionVars
        self.samplingSet = samplingSet
        self._cnfData = DIMACScnf

        # Note that change in any variable can be directly done without an explicit set function
        self.totalVariables = None
        self.totalClauses = None
        self.treenodes = []
        self.graph = None
        self.samples = None
        self.drawnNodes = {}
        self.totalSamples = 10
        self.isSamplingSetPresent = False

        if cnfFile:
            self._cnfData = open(cnfFile).read()

    def compile(self, cnfFile=None, samplingSet=[]):
        """Compiles the cnf to dDNNF projected over sampling set.
        :param cnfFile: The path of file containing DIMACS cnf
        :param samplingSet: variables on which samples are projected, sampling is done on all variables by default
        """
        if cnfFile:
            with open(cnfFile, "r") as f:
                text = f.read()
                self._cnfData = text
                f.close()
        elif self._cnfData:
            text = self._cnfData
        else:
            raise Exception("No cnf provided to sampler for compilation")
        if not self._cnfFile:
            self._cnfFile = "default.cnf"
        # if cnf has projecting vars, compile with Dsharp_Pcompile else D4
        pvarline = re.search(
            r"(\n[\s]*|^)p([\s]+)cnf([\s]+)([0-9]+)([\s]+)([0-9]+)", text
        )
        self.totalVariables = int(pvarline.groups(1)[3])
        self.totalClauses = int(pvarline.groups(1)[5])
        if samplingSet:
            self.samplingSet = samplingSet
            self.isSamplingSetPresent = True
        else:
            self.samplingSet = list(range(1, self.totalVariables + 1))
        print("Seperating weights from Input cnf")
        weighttext = ""
        print("Extracting the Sampling Set")
        if (
            not self.isSamplingSetPresent
        ):  # sampling set provided via arguments overrides that in file
            with open("/tmp/" + self._cnfFile.split("/")[-1] + ".pvars", "w") as f:
                samplingvars = "v "
                for ind in re.findall(r"c ind.*", text):
                    self.isSamplingSetPresent = True
                    samplingvars += " ".join(ind.split(" ")[2:-1])
                    samplingvars += " "
                samplingvars += "0"
                if self.isSamplingSetPresent:
                    self.samplingSet = list(map(int, samplingvars.split()[1:-1]))
                    # for variable in samplingvars.split()[1:-1]:
                    #     self.variables.append(int(variable))
                f.write(samplingvars)
                f.close()
        with open("/tmp/" + self._cnfFile.split("/")[-1] + ".tmp", "w") as f:
            f.write(text.replace("w", "c w"))
            f.close()
            weighttext = re.findall(r"^w[^\S\n]+.*", text, re.MULTILINE)
        for line in weighttext:
            if int(line.split()[1]) * (-1) in self.weights:
                if int(line.split()[1]) > 0:
                    self.weights[int(line.split()[1])] = mpq(
                        Fraction(line.split()[2])
                    ) / (
                        self.weights.pop(int(line.split()[1]) * (-1), None)
                        + mpq(Fraction(line.split()[2]))
                    )
                else:
                    self.weights[int(line.split()[1]) * (-1)] = self.weights[
                        int(line.split()[1]) * (-1)
                    ] / (
                        self.weights[int(line.split()[1]) * (-1)]
                        + mpq(Fraction(line.split()[2]))
                    )
            else:
                self.weights[int(line.split()[1])] = mpq(Fraction(line.split()[2]))
        dDNNF = self._cnfFile.split("/")[-1] + ".nnf"
        cmd = (
            "/usr/bin/time -o "
            + "/tmp/"
            + self._cnfFile.split("/")[-1]
            + ".timeout "
            + "--verbose WAPS/bin/d4 /tmp/"
            + self._cnfFile.split("/")[-1]
            + ".tmp "
            + " -out="
            + dDNNF
        )
        if self.isSamplingSetPresent:
            cmd = (
                "/usr/bin/time -o "
                + "/tmp/"
                + self._cnfFile.split("/")[-1]
                + ".timeout "
                + "--verbose WAPS/bin/Dsharp_PCompile -cs 2000 -pvarsfile "
                + "/tmp/"
                + self._cnfFile.split("/")[-1]
                + ".pvars"
                + " -Fnnf "
                + dDNNF
                + " /tmp/"
                + self._cnfFile.split("/")[-1]
                + ".tmp"
            )
        self._dDNNFfile = dDNNF
        start = time.time()
        os.system(cmd)
        print("The time taken by D4/Dsharp_PCompile is ", time.time() - start)

    def parse(self, dDNNFfile=None):
        """Parses the d-DNNF tree to a tree like object

        :param dDNNFfile: specifies file containing d-DNNF of the formula to sample from
        """
        if dDNNFfile:
            self._dDNNFfile = dDNNFfile
        with open(self._dDNNFfile) as f:
            treetext = f.readlines()
        nodelen = 0
        for node in treetext:
            node = node.split()
            if node[0] == "c":
                continue
            elif node[0] == "nnf":
                self.totalVariables = int(node[3])
            elif node[0] == "L":
                self.treenodes.append(Node(label=int(node[1])))
                nodelen += 1
            elif node[0] == "A":
                if node[1] == "0":
                    self.treenodes.append(Node(label="T " + str(nodelen)))
                else:
                    andnode = Node(label="A " + str(nodelen))
                    andnode.children = list(
                        map(lambda x: self.treenodes[int(x)], node[2:])
                    )
                    self.treenodes.append(andnode)
                nodelen += 1
            elif node[0] == "O":
                if node[2] == "0":
                    self.treenodes.append(Node(label="F " + str(nodelen)))
                else:
                    ornode = Node(label="O " + str(nodelen), decision=int(node[1]))
                    ornode.children = list(
                        map(lambda x: self.treenodes[int(x)], node[3:])
                    )
                    self.treenodes.append(ornode)
                nodelen += 1

    def annotate(self, weights={}, conditionVars=[]):
        """Annotates d-DNNF with weights

        :param weights: dictionary mapping literals to their weights
        :param conditionVars: list specifying list of literals which should be true in samples
        """
        if weights:
            self.weights = weights
        normalizeWeights(self.weights)
        if conditionVars:
            conditionWeights(conditionVars, self.weights)
        elif self.conditionVars:
            conditionWeights(self.conditionVars, self.weights)
        self._annotate(self.treenodes[-1], weights=self.weights)

    def _annotate(self, root, weights={}):
        """Actually annotates d-DNNF with weights"""
        if str(root.label)[0] == "A":
            root.weight = mpq("1")
            for ch in root.children:  # can perform IBCP for conditioning
                root.weight *= self._annotate(ch, weights=weights)
            return root.weight
        elif str(root.label)[0] == "O":
            root.weight = self._annotate(
                root.children[0], weights=weights
            ) + self._annotate(root.children[1], weights=weights)
            return root.weight
        else:
            try:
                int(root.label)
                if weights and abs(int(root.label)) in weights:
                    if int(root.label) > 0:
                        root.weight = weights[int(root.label)]
                    else:
                        root.weight = mpq("1") - weights[abs(int(root.label))]
                else:
                    root.weight = mpq("0.5")
            except:
                if str(root.label)[0] == "F":
                    root.weight = 0
                elif str(root.label)[0] == "T":
                    root.weight = 1
            return root.weight

    def _get_samples(self, root, indices):
        """Retrieves samples from tree rooted at root"""
        if not indices.shape[0]:
            return
        if str(root.label)[0] == "O":
            z0 = root.children[0].weight
            z1 = root.children[1].weight
            p = (mpq("1.0") * z0) / (z0 + z1)
            tosses = np.random.binomial(1, p, indices.shape[0])
            self._get_samples(
                root.children[0], np.array(indices[np.where(tosses == 1)[0]])
            )
            self._get_samples(
                root.children[1], np.array(indices[np.where(tosses == 0)[0]])
            )
        elif str(root.label)[0] == "A":
            for ch in root.children:
                self._get_samples(ch, indices)
        else:
            try:
                int(root.label)
                for index in indices:
                    self.samples[index] += str(root.label) + " "
            except:
                pass

    def sample(self, totalSamples=10, randAssign=True):
        """Samples totalSamples samples and extends them to all variables if randAssign is set to True

        :param totalSamples: Number of samples to be sampled
        :param randAssign: extends each sample to contain all variables from samplingSet"""
        self.samples = []
        if totalSamples:
            self.totalSamples = totalSamples
        for _ in range(self.totalSamples):
            self.samples.append("")
        self._get_samples(self.treenodes[-1], np.arange(0, self.totalSamples))
        if randAssign:
            if not self.isSamplingSetPresent:
                self.samplingSet = list(range(1, self.totalVariables + 1))
            self.samples = map(
                lambda x: random_assignment(
                    x, samplingSet=self.samplingSet, weights=self.weights
                ),
                self.samples,
            )
        return self.samples

    def save(self, outputFile=None):
        """
        Saves the samples in outputfile

        :param outputFile: Saves the samples in outputfile. Samples are saved in samples.txt if not specified.
        """
        if outputFile:
            self.outputFile = outputFile
        else:
            self.outputFile = "samples.txt"
        f = open(self.outputFile, "w+")
        for i in range(self.totalSamples):
            f.write(str(i + 1) + ", " + self.samples[i] + "\n")
        f.close()
        print("Samples saved to", self.outputFile)

    def save_annotation_tree(self, filename=None):
        """Saves annotated d-DNNF in pickle format

        :param filename: Saves the annotated d-DNNF by filename"""
        if not filename:
            filename = "default.pkl"
            print("The tree is getting saved in current directory as: ", filename)
        fp = open(filename, "wb")
        pickle.dump((self.samplingSet, self.totalVariables, self.treenodes), fp)
        fp.close()

    def load_annotation_tree(self, filename):
        """Loads Annotation Tree saved by save_annotation_tree()

        :param filename: Loads Annotation Tree from filename"""
        fp = open(filename, "rb")
        (self.samplingSet, self.totalVariables, self.treenodes) = pickle.load(fp)
        fp.close()


class sampler2:
    """Main class for main which defines parsing, graph drawing, counting and sampling functions"""

    def __init__(self):
        self.totalVariables = None
        self.variables = []
        self.treenodes = []
        self.graph = None
        self.samples = None
        self.drawnNodes = {}
        self.isSamplingSetPresent = False

    def drawtree(self, root):
        """Recursively draws tree for the d-DNNF"""
        rootnode = pydot.Node(str(root.label) + " " + str(root.weight))
        self.graph.add_node(rootnode)
        self.drawnNodes[root.label] = rootnode
        for ch in root.children:
            if ch.label not in self.drawnNodes:
                node = self.drawtree(ch)
                self.graph.add_edge(pydot.Edge(rootnode, node))
            else:
                self.graph.add_edge(pydot.Edge(rootnode, self.drawnNodes[ch.label]))
        return rootnode

    def parse(self, inputnnffile):
        """Parses the d-DNNF tree to a tree like object"""
        with open(inputnnffile) as f:
            treetext = f.readlines()
        nodelen = 0
        for node in treetext:
            node = node.split()
            if node[0] == "c":
                continue
            elif node[0] == "nnf":
                self.totalVariables = int(node[3])
            elif node[0] == "L":
                self.treenodes.append(Node(label=int(node[1])))
                nodelen += 1
            elif node[0] == "A":
                if node[1] == "0":
                    self.treenodes.append(Node(label="T " + str(nodelen)))
                else:
                    andnode = Node(label="A " + str(nodelen))
                    andnode.children = list(
                        map(lambda x: self.treenodes[int(x)], node[2:])
                    )
                    self.treenodes.append(andnode)
                nodelen += 1
            elif node[0] == "O":
                if node[2] == "0":
                    self.treenodes.append(Node(label="F " + str(nodelen)))
                else:
                    ornode = Node(label="O " + str(nodelen), decision=int(node[1]))
                    ornode.children = list(
                        map(lambda x: self.treenodes[int(x)], node[3:])
                    )
                    self.treenodes.append(ornode)
                nodelen += 1

    def annotate(self, root, weights=None):
        """Computes Model Counts"""
        if str(root.label)[0] == "A":
            root.weight = mpq("1")
            for ch in root.children:  # can perform IBCP for conditioning
                root.weight *= self.annotate(ch, weights=weights)
            return root.weight
        elif str(root.label)[0] == "O":
            root.weight = self.annotate(
                root.children[0], weights=weights
            ) + self.annotate(root.children[1], weights=weights)
            return root.weight
        else:
            try:
                int(root.label)
                if weights and abs(int(root.label)) in weights:
                    if int(root.label) > 0:
                        root.weight = weights[int(root.label)]
                    else:
                        root.weight = mpq("1") - weights[abs(int(root.label))]
                else:
                    root.weight = mpq("0.5")
            except:
                if str(root.label)[0] == "F":
                    root.weight = 0
                elif str(root.label)[0] == "T":
                    root.weight = 1
            return root.weight

    def getsamples(self, root, indices):
        """Generates Uniform Independent Samples"""
        if not indices.shape[0]:
            return
        if str(root.label)[0] == "O":
            z0 = root.children[0].weight
            z1 = root.children[1].weight
            p = (mpq("1.0") * z0) / (z0 + z1)
            tosses = np.random.binomial(1, p, indices.shape[0])
            self.getsamples(
                root.children[0], np.array(indices[np.where(tosses == 1)[0]])
            )
            self.getsamples(
                root.children[1], np.array(indices[np.where(tosses == 0)[0]])
            )
        elif str(root.label)[0] == "A":
            for ch in root.children:
                self.getsamples(ch, indices)
        else:
            try:
                int(root.label)
                for index in indices:
                    self.samples[index] += str(root.label) + " "
            except:
                pass


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--outputfile",
        type=str,
        default="samples.txt",
        help="output file for samples",
        dest="outputfile",
    )
    parser.add_argument(
        "--drawtree",
        type=str,
        default="",
        help="draw nnf tree in png format with given filename",
        dest="figureName",
    )
    parser.add_argument(
        "--samples", type=int, default=10, help="number of samples", dest="samples"
    )
    parser.add_argument(
        "--randAssign",
        type=int,
        default=1,
        help="randomly assign unassigned variables in a model with partial assignments",
        dest="randAssign",
    )
    parser.add_argument(
        "--saveAnnotation",
        type=str,
        default=None,
        help="specify filename for saving pickle of count annotated dDNNF for incremental sampling",
        dest="pickleName",
    )
    parser.add_argument(
        "--weights",
        type=str,
        default=None,
        help="specify a csv file which contains weights for literals",
        dest="weightFile",
    )
    parser.add_argument(
        "--conditionVars",
        type=str,
        default="",
        help="specify the literals separated by space within quotes on which you want to condition",
        dest="conditionVars",
    )
    parser.add_argument(
        "--conditionFile",
        type=str,
        default="",
        help="specify the file containing the literals on which you want to condition",
        dest="conditionFile",
    )
    parser.add_argument("--dDNNF", type=str, help="specify dDNNF file", dest="dDNNF")
    parser.add_argument(
        "--loadAnnotation",
        type=str,
        help="specify filename containing pickle of count annotated dDNNF",
        dest="loadPickle",
    )
    parser.add_argument(
        "DIMACSCNF", nargs="?", type=str, default="", help="input cnf file"
    )
    args = parser.parse_args()
    if not (args.dDNNF or args.loadPickle or args.DIMACSCNF):
        parser.error(
            "Please give at least one argument out of dDNNF, countPickle and DIMACSCNF"
        )
    draw = args.figureName
    totalsamples = args.samples
    randAssignInt = args.randAssign
    dDNNF = False
    countPickle = False
    inputFile = False
    if args.loadPickle:
        countPickle = args.loadPickle
    else:
        if args.dDNNF:
            dDNNF = args.dDNNF
        if args.DIMACSCNF:
            inputFile = args.DIMACSCNF
    savePickle = args.pickleName
    randAssign = False
    if randAssignInt == 1:
        randAssign = True
    if args.weightFile:
        weights = fetchWeights(args.weightFile)
    else:
        weights = {}
    sampler = sampler2()
    if inputFile:
        print("Seperating weights from Input cnf")
        weighttext = ""
        with open(inputFile, "r") as f:
            text = f.read()
            f.close()
        print("Extracting the Sampling Set")
        with open("/tmp/" + inputFile.split("/")[-1] + ".pvars", "w") as f:
            samplingvars = "v "
            for ind in re.findall(r"c ind.*", text):
                sampler.isSamplingSetPresent = True
                samplingvars += " ".join(ind.split(" ")[2:-1])
                samplingvars += " "
            samplingvars += "0"
            if sampler.isSamplingSetPresent:
                for variable in samplingvars.split()[1:-1]:
                    sampler.variables.append(int(variable))
            f.write(samplingvars)
            f.close()
        with open("/tmp/" + inputFile.split("/")[-1] + ".tmp", "w") as f:
            f.write(text.replace("w", "c w"))
            f.close()
            weighttext = re.findall(r"^w[^\S\n]+.*", text, re.MULTILINE)
        for line in weighttext:
            if int(line.split()[1]) * (-1) in weights:
                if int(line.split()[1]) > 0:
                    weights[int(line.split()[1])] = mpq(Fraction(line.split()[2])) / (
                        weights.pop(int(line.split()[1]) * (-1), None)
                        + mpq(Fraction(line.split()[2]))
                    )
                else:
                    weights[int(line.split()[1]) * (-1)] = weights[
                        int(line.split()[1]) * (-1)
                    ] / (
                        weights[int(line.split()[1]) * (-1)]
                        + mpq(Fraction(line.split()[2]))
                    )
            else:
                weights[int(line.split()[1])] = mpq(Fraction(line.split()[2]))
        if not args.dDNNF:
            dDNNF = inputFile.split("/")[-1] + ".nnf"

            cmd = (
                "/usr/bin/time -o "
                + "/tmp/"
                + inputFile.split("/")[-1]
                + ".timeout "
                + "--verbose WAPS/bin/d4 /tmp/"
                + inputFile.split("/")[-1]
                + ".tmp "
                + " -out="
                + dDNNF
            )
            if sampler.isSamplingSetPresent:
                cmd = (
                    "/usr/bin/time -o "
                    + "/tmp/"
                    + inputFile.split("/")[-1]
                    + ".timeout "
                    + "--verbose WAPS/bin/Dsharp_PCompile -cs 2000 -pvarsfile "
                    + "/tmp/"
                    + inputFile.split("/")[-1]
                    + ".pvars"
                    + " -Fnnf "
                    + dDNNF
                    + " /tmp/"
                    + inputFile.split("/")[-1]
                    + ".tmp"
                )
            start = time.time()
            if os.system(cmd):
                raise Exception("D4/Dsharp_PCompile not found")
            print("The time taken by D4/Dsharp_PCompile is ", time.time() - start)
    if dDNNF:
        start = time.time()
        sampler.parse(dDNNF)
        if sampler.variables == []:
            for i in range(1, sampler.totalVariables + 1):
                sampler.variables.append(i)
        print("The time taken to parse the nnf text:", time.time() - start)
        # can easily adjust code in conditionWeights to give cmd/file priority
        # right now, it simply takes union of the conditioned literals
        if args.conditionFile:
            lits = open(args.conditionFile).read().strip().split()
            conditionWeights(lits, weights)
        if args.conditionVars:
            lits = args.conditionVars.split()
            conditionWeights(lits, weights)

        start = time.time()
        modelcount = sampler.annotate(sampler.treenodes[-1], weights=weights)
        
        sampler.treenodes[-1].models = modelcount
        print("The time taken for Model Counting:", time.time() - start)
        timepickle = time.time()
        if savePickle:
            fp = open(savePickle, "wb")
            pickle.dump(
                (sampler.variables, sampler.totalVariables, sampler.treenodes), fp
            )
            fp.close()
            print(
                "Time taken to save the count annotated dDNNF pickle:",
                time.time() - timepickle,
            )
    else:
        timepickle = time.time()
        fp = open(countPickle, "rb")
        (sampler.variables, sampler.totalVariables, sampler.treenodes) = pickle.load(fp)
        fp.close()
        print("The time taken to read the pickle:", time.time() - timepickle)
        if savePickle:
            fp = open(savePickle, "wb")
            pickle.dump(
                (sampler.variables, sampler.totalVariables, sampler.treenodes), fp
            )
            fp.close()
            print(
                "Time taken to save the count annotated dDNNF pickle:",
                time.time() - timepickle,
            )
    if weights:
        print(
            "Weighted Model Count as per normalised weights limited to var in dDNNF:",
            mpfr(sampler.treenodes[-1].weight),
        )
    else:
        print(
            "Model Count limited to var in dDNNF:", mpfr(sampler.treenodes[-1].models)
        )
    if draw:
        sampler.graph = pydot.Dot(graph_type="digraph")
        sampler.drawtree(sampler.treenodes[-1])
        sampler.graph.write_png(draw)
    sampler.samples = []
    for i in range(totalsamples):
        sampler.samples.append("")
    start = time.time()
    if sampler.treenodes[-1].weight == 0:
        print("The current conditional assignment has no satisfying sample")
        exit()
    sampler.getsamples(sampler.treenodes[-1], np.arange(0, totalsamples))
    print("The time taken by sampling:", time.time() - start)
    if randAssign:
        sampler.samples = list(
            map(
                lambda x: random_assignment(
                    x, samplingSet=sampler.variables, weights=weights
                ),
                sampler.samples,
            )
        )
    f = open(args.outputfile, "w+")
    for i in range(totalsamples):
        f.write(str(i + 1) + ", " + sampler.samples[i] + "\n")
    f.close()
    print("Samples saved to", args.outputfile)


if __name__ == "__main__":
    main()
