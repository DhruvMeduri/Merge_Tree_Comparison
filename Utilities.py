import MergeTreeLibrary
import sys
import os
import csv
import matplotlib.pyplot as plt
import math

def single_compute (MeasureName, FileName1, arg=-1):
    rval = float("inf")

    if MeasureName == "TotalPers":
        rval = MergeTreeLibrary.TotalPers(FileName1, arg)

    elif MeasureName == "depth":
        T = MergeTreeLibrary.generateTree(FileName1)
        rval = T.depth()

    elif MeasureName == "feat":
        T = MergeTreeLibrary.generateTree(FileName1)
        rval = T.countFeatures()

    return rval

def pair_compute (MeasureName, FileName1, FileName2, CompIDFolderName):

    rval = float("inf")

    if not os.path.exists(CompIDFolderName):
        os.mkdir(CompIDFolderName)

    if MeasureName == "bottleneck":
        rval = MergeTreeLibrary.bottleneck(FileName1, FileName2)

    elif MeasureName == "Wasserstein":
        rval = MergeTreeLibrary.Wasserstein(FileName1, FileName2)

    elif MeasureName == "COMTED":
        rval = MergeTreeLibrary.COMTED(FileName1, FileName2, CompIDFolderName)

    elif MeasureName == "CWMTED":
        rval = MergeTreeLibrary.CWMTED(FileName1, FileName2, CompIDFolderName)

    return rval

def output_normalisation (MeasureName, FileName1, FileName2):

    rval = 1

    if MeasureName == "bottleneck":
        rval = MergeTreeLibrary.TotalPers(FileName1, -1) * MergeTreeLibrary.TotalPers(FileName2, -1)

    elif MeasureName == "COMTED":
        rval = MergeTreeLibrary.TotalPers(FileName1, 1) + MergeTreeLibrary.TotalPers(FileName2, 1)

    elif MeasureName == "CWMTED":
        rval = MergeTreeLibrary.TotalPers(FileName1, 1) + MergeTreeLibrary.TotalPers(FileName1, 1)

    return rval

def single_normalisation (MeasureName, FileName1):
    if MeasureName == "TotalPers":
        rval = MergeTreeLibrary.TotalPers(FileName1, -1)
    return rval

def draw_plot (MeasureName, FileName1, FileName2, OutputFileName, MatchingFileName = None, normalised = False):

    if MeasureName == "bottleneck" or MeasureName == "wasserstein":
        T1 = MergeTreeLibrary.generateTree(FileName1)
        T2 = MergeTreeLibrary.generateTree(FileName2)

        T1.generatePersistencePairings()
        T2.generatePersistencePairings()

        NodeCoords1 = T1.persistenceDiagram()
        NodeCoords2 = T2.persistenceDiagram()

        XArrs1 = [row[0] for row in NodeCoords1]
        YArrs1 = [row[1] for row in NodeCoords1]

        XArrs2 = [row[0] for row in NodeCoords2]
        YArrs2 = [row[1] for row in NodeCoords2]

        fog = plt.figure(OutputFileName)
        axis = fog.add_subplot(111)

        axis.scatter(XArrs1, YArrs1, c = 'b', label = FileName1)
        axis.scatter(XArrs2, YArrs2, c = 'r', label = FileName2)

        plt.legend()
        plt.savefig(OutputFileName)
        plt.close()

    elif MeasureName == "pers_dist":
        T1 = MergeTreeLibrary.generateTree(FileName1)
        pers = T1.persistenceDiagram()
        listabsvals = list()
        for point in pers:
            listabsvals.append(abs(point[0] - point[1]))
        listabsvals.sort()
        if normalised == True:
            listabsvals = [row / T1.TotalPers(-1) for row in listabsvals]
        plt.figure(OutputFileName)
        plt.scatter(range(len(listabsvals)), listabsvals, c = "b")
        plt.savefig(OutputFileName)
        plt.close()

    elif MeasureName == "COMTED" or MeasureName == "CWMTED":

        TFName1 = FileName1.strip(".").replace("/", "_")
        TFName2 = FileName2.strip(".").replace("/", "_")

        s = "digraph { \n"

        s1 = ""
        with open(FileName1) as IF1:

            n = int(IF1.readline().strip())
            s1 = "subgraph cluster_1 {\n style=filled \n; color=lightgrey; \n label=" + TFName1 + ";"
            for i in range(n):
                ni = IF1.readline().strip().split(' ')
                if ni[1] != "-1":
                    s1 = s1 + TFName1 + "_" + ni[1] + " -> " + TFName1 + "_" + ni[0] + ";\n"
            s1 = s1 + "}"

        s2 = ""
        with open(FileName2) as IF2:
            n = int(IF2.readline().strip())
            s2 = "subgraph cluster_2 {\n style=filled \n; color=lightgrey; \n label=" + TFName2 +";"
            for i in range(n):
                ni = IF2.readline().strip().split(' ')
                if ni[1] != "-1":
                    s2 = s2 + TFName2 + "_" + ni[1] + " -> " + TFName2 + "_" + ni[0] + ";\n"
            s2 = s2 + "}"

        sM = ""
        with open(MatchingFileName) as MatchFile:
            for line in MatchFile.readlines():
                MapNode = line.strip().split(",")
                sM = sM + TFName1 + "_" + MapNode[0] + " -> " + TFName2 + "_" + MapNode[1] + "\n"

        s = s + s1 + s2 + sM
        s = s + "}"

        with open(OutputFileName, "w") as OF1:
            OF1.write(s)

def pairwiseCompare (InputDir, OutputDir, T1Name, T2Name, MeasuresTBC):

    OutputDict = dict()

    CompID = T1Name + T2Name
    OutputDict["ID"] = CompID

    Tree1 = MergeTreeLibrary.generateTree(InputDir + "/" + T1Name)
    Tree2 = MergeTreeLibrary.generateTree(InputDir + "/" + T2Name)

    Tree1.generatePersistencePairings()
    Tree2.generatePersistencePairings()

    for Measure in MeasuresTBC:
        MeasureName = Measure[0]
        MeasureArgs = list()
        if len(Measure) > 1:
            MeasureArgs = Measure[1].split(" ")

        raw_compute = pair_compute(MeasureName, InputDir + "/" + T1Name, InputDir + "/" + T2Name, OutputDir + "/" + CompID)
        OutputDict[MeasureName] = raw_compute

        innormval = float("inf")
        if "in_normal" in MeasureArgs:

            Tree1Normalised = Tree1.normalise()
            Tree2Normalised = Tree2.normalise()

            Tree1Normalised.generatePersistencePairings()
            Tree2Normalised.generatePersistencePairings()

            Tree1Normalised.genTreeFile(Tree1Normalised.genCEDTInput(), InputDir + "/" + T1Name + "normalised")
            Tree2Normalised.genTreeFile(Tree2Normalised.genCEDTInput(), InputDir + "/" + T2Name + "normalised")

            innormval = pair_compute (MeasureName, InputDir + "/" + T1Name + "normalised", InputDir + "/" + T2Name + "normalised", OutputDir + "/" + CompID + "/normalised")
            if innormval < float("inf"):
                OutputDict["InputNormalised" + MeasureName] = innormval

        if "out_normal" in MeasureArgs:
            OutputDict["OutputNormalised" + MeasureName] = raw_compute / output_normalisation(MeasureName, InputDir + "/" + T1Name, InputDir + "/" + T2Name)

        if "plots" in MeasureArgs:
            OutputFileName = OutputDir + "/" + CompID + "/" + MeasureName + "map"
            stringmap = ""
            if MeasureName == "CWMTED":
                stringmap = "CW"
            elif MeasureName == "COMTED":
                stringmap = "CO"
            draw_plot(MeasureName, InputDir + "/" + T1Name, InputDir + "/" + T2Name, OutputFileName, OutputDir + "/" + CompID + "/forwardmap" + stringmap + ".txt")

    return OutputDict

def soloCompare (InputDir, OutputDir, T1Name, SingleMeasuresTBC):
    OutputDict = dict()
    CompID = T1Name
    OutputDict["ID"] = CompID
    Tree1 = MergeTreeLibrary.generateTree(InputDir + "/" + T1Name)
    Tree1.generatePersistencePairings()

    for Measure in SingleMeasuresTBC:
        MeasureName = Measure[0]
        MeasureArgs = list()
        if len(Measure) > 1:
            MeasureArgs = Measure[1].split(" ")

        if MeasureName == "TotalPers":
            power = float(MeasureArgs[1])
            raw_compute = single_compute(MeasureName, InputDir + "/" + T1Name, power)
            OutputDict[MeasureName + str(power)] = raw_compute
            if "normal" in MeasureArgs:
                OutputDict["Normalised" + MeasureName + str(power)] = raw_compute / single_normalisation(MeasureName, InputDir + "/" + T1Name)

        elif MeasureName == "depth" or MeasureName == "feat":
            raw_compute = single_compute(MeasureName, InputDir + "/" + T1Name)
            OutputDict[MeasureName] = raw_compute

        elif MeasureName == "pers_dist":
            if "normal" in MeasureArgs:
                draw_plot("pers_dist", InputDir + "/" + T1Name, None, OutputDir + "/" + T1Name, None, True)
            else:
                draw_plot("pers_dist", InputDir + "/" + T1Name, None, OutputDir + "/" + T1Name, None, False)

    return OutputDict

def GenTreeFromJt (FileNameJT):
    with open(FileNameJT) as JTFile:
        with open(FileNameJT.replace(".jt", ""), "w") as TreeFile:
            list_lines = list()
            n = JTFile.readline()
            TreeFile.write(n)
            for i in range(int(n)):
                newline = list()
                nodeLine = JTFile.readline().strip().split()
                newline.append(nodeLine[0])
                newline.append(nodeLine[1])
                fVal = nodeLine[4]
                nodeChildren = nodeLine[5:]
                newline.extend(nodeChildren)
                newline.append(fVal + "\n")
                list_lines.append(newline)

            text_list_lines = list()
            for nodeData in list_lines:
                text_list_lines.append(" ".join(nodeData))

            TreeFile.writelines(text_list_lines)
