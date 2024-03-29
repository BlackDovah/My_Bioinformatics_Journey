import Bio as Bio
import numpy as np
import pandas as pd
import sys as sys
import networkx as nx
# from Bio import pairwise2
# from Bio.pairwise2 import format_alignment
from Bio import Align
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import copy
from copy import deepcopy
import random
from random import randint
from Replication import ReverseComplement


P = "(+1 -2 -4 +3)"
# P = P.split()
# P = [int(i) for i in P]
# P = "(-1 +2 -3 -4 -5 -6 -7 +8 +9 +10 +11 +12 +13 +14 +15 -16 -17 -18 +19 -20 +21 +22 -23 -24 +25 -26 -27 +28 -29 -30 +31 +32 +33 -34 -35 +36 +37 +38 +39 -40 +41 +42 -43 +44 +45 +46 -47 +48 +49 +50 -51 +52 +53 +54 -55 -56 -57 +58 -59 -60 -61 +62 -63 +64 +65 -66 -67)"
# P = "(+1 -2 -4 +3)"
P = "(+1 -2 -3 +4)"
Q = "(+1 +2 -4 -3)"
# P = "(+9 -8 +12 +7 +1 -14 +13 +3 -5 -11 +6 -2 +10 -4)"
# Q = "(-11 +8 -10 -2 +3 +4 +13 +6 +12 +9 +5 +7 -14 -1)"
# P = P.removesuffix(")")
# P = P.removeprefix("(")
# P = P.split(")(")
# Nodes = "1 2 3 4 6 5 8 7 10 9 12 11 14 13 16 15 17 18 19 20 21 22 23 24 25 26 27 28 29 30 32 31 33 34 36 35 38 37 39 40 42 41 44 43 46 45 48 47 50 49 51 52 53 54 56 55 57 58 59 60 62 61 64 63 65 66 67 68 69 70 71 72 74 73 75 76 77 78 80 79 82 81 84 83 86 85 88 87 89 90 91 92 94 93 95 96 97 98 100 99 102 101 104 103 106 105 107 108 109 110 111 112 113 114 116 115 117 118 119 120"
# Nodes = Nodes.split()

# Graph = [(2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)]
# Graph = [(2, 4), (3, 5), (6, 8), (7, 10), (9, 12), (11, 13), (14, 15), (16, 17), (18, 19), (20, 21), (22, 24), (23, 25), (26, 28), (27, 30), (29, 32), (31, 33), (34, 35), (36, 37), (38, 40), (39, 42), (41, 44), (43, 45), (46, 47), (48, 49), (50, 51), (52, 53), (54, 56), (55, 58), (57, 59), (60, 61), (62, 64), (63, 66), (65, 68), (67, 69), (70, 72), (71, 73),
#  (74, 75), (76, 77), (78, 80), (79, 81), (82, 83), (84, 86), (85, 88), (87, 90), (89, 92), (91, 94), (93, 96), (95, 98), (97, 99), (100, 102), (101, 104), (103, 106), (105, 107), (108, 109), (110, 111), (112, 113), (114, 116), (115, 117), (118, 120), (119, 121), (122, 124), (123, 126), (125, 128), (127, 129), (130, 132), (131, 134), (133, 136), (135, 1)]
# Graph = [[(2, 3), (4, 5), (6, 8), (7, 9), (10, 11), (12, 13), (14, 16), (15, 17), (18, 20), (19, 22), (21, 24), (23, 25), (26, 28), (27, 29), (30, 31), (32, 33), (34, 36), (35, 38), (37, 40), (39, 42), (41, 44), (43, 45), (46, 47), (48, 49), (50, 52), (51, 53), (54, 55), (56, 58), (57, 59), (60, 61), (62, 63), (64, 65), (66, 67), (68, 69),
#   (70, 71), (72, 73), (74, 76), (75, 77), (78, 80), (79, 81), (82, 84), (83, 85), (86, 87), (88, 89), (90, 92), (91, 93), (94, 95), (96, 97), (98, 100), (99, 102), (101, 103), (104, 105), (106, 107), (108, 110), (109, 112), (111, 114), (113, 116), (115, 118), (117, 120), (119, 121), (122, 123), (124, 126), (125, 127), (128, 130), (129, 1)]]
# Graph = [[(2, 4), (3, 1), (7, 5), (6, 8)]]
# breakSeq = [1, 6, 3, 8]
# breakSeq = [87, 90, 74, 75]
# breakSeq = [116, 118, 73, 72]
# breakSeq = [1, 6, 3, 8]
seq1 = "TGCCCCGGTGGTGAG"
seq2 = "AAGGTCGCACCTCGT"
k = 3


def Reverse(revBlock):
    reved = [int(i)*-1 for i in reversed(revBlock)]
    return reved


def GreedySorting(Perm):
    revCount = 0
    startBlock = []
    revBlock = []
    endBlock = []
    PermLis = []
    for i in range(1, len(Perm)+1):
        searchPerm = [abs(i)for i in Perm]
        for y in Perm:
            if Perm.index(y)+1 < i:
                continue
            else:
                if abs(y) == i:
                    if y == i:
                        break
                    else:
                        Perm[i-1] = y * -1
                        PermLis.append(copy.deepcopy(Perm))
                        revCount += 1
                        break
                else:
                    endPos = searchPerm.index(abs(i))
                    startPos = Perm.index(y)
                    revBlock = Perm[startPos:endPos+1]
                    if startPos != 0:
                        startBlock = Perm[0:startPos]
                    if endPos != len(Perm):
                        endBlock = Perm[endPos+1:len(Perm)]
                    revedBlock = Reverse(revBlock)
                    Perm = startBlock + revedBlock + endBlock
                    PermLis.append(copy.deepcopy(Perm))
                    revCount += 1
                    y = Perm[i-1]
                    if abs(y) == i:
                        if y == i:
                            break
                        else:
                            Perm[i-1] = y * -1
                            PermLis.append(copy.deepcopy(Perm))
                            revCount += 1
                            break
    for i in PermLis:
        for y in i:
            if y != abs(y):
                i[i.index(y)] = str(y)
            else:
                i[i.index(y)] = "+"+str(y)
        i = (" ").join(i)
        print(i)

    # PermLis = np.array(PermLis)
    return revCount


def BreakPointCounter(Perm):
    breakPoints = 0
    Perm.append(len(Perm)+1)
    Perm.insert(0, 0)
    for i in range(len(Perm)-1):
        if Perm[i+1] - Perm[i] != 1:
            breakPoints += 1

    return breakPoints


def ChromosomeToCycle(Perm):
    nodesLis = []
    if type(Perm) == str:
        Perm = Perm.split()
    # index = 1
    # for i in Perm:
    #     if i[0:1] == "+":
    #         nodesLis.append(index)
    #         nodesLis.append(index + 1)
    #         index += 2
    #     elif i[0:1] == "-":
    #         nodesLis.append(index + 1)
    #         nodesLis.append(index)
    #         index += 2
    for i in Perm:
        if int(i) > 0:
            nodesLis.append(int(i)*2 - 1)
            nodesLis.append(int(i)*2)
        else:
            nodesLis.append(abs(int(i))*2)
            nodesLis.append(abs(int(i))*2 - 1)
    nodesLis = [str(i) for i in nodesLis]
    return "("+(" ").join(nodesLis)+")"


def CycleToChromosome(Nodes):
    chromLis = []
    index = 1
    for i in range(0, len(Nodes)-1, 2):
        Edge = Nodes[i:i+2]
        if int(("").join(Edge[1:2])) - int(("").join(Edge[0:1])) == 1:
            chromLis.append("+"+str(int(("").join(Edge[1:2]))//2))
            index += 1
        elif int(("").join(Edge[1:2])) - int(("").join(Edge[0:1])) == -1:
            chromLis.append("-"+str(int(("").join(Edge[0:1]))//2))
            index += 1
    chromLis = [str(i) for i in chromLis]
    return "("+(" ").join(chromLis)+")"


def ColoredEdges(Perm):
    cycleLis = []
    nodesLis = []
    temp = []
    temp2 = []
    temp3 = []
    if type(Perm) == str:
        Perm = Perm.removeprefix("(")
        Perm = Perm.removesuffix(")")
        Perm = Perm.split(")(")
        # Perm = (" ").join(Perm)
        # Perm = Perm.split()

        for i in Perm:
            Nodes = ChromosomeToCycle(i)
            cycleLis.append(Nodes)

        Nodes = (" ").join(cycleLis)
        Nodes = Nodes.replace(") (", " ")
        Nodes = Nodes.removeprefix("(")
        Nodes = Nodes.removesuffix(")")
        Nodes = Nodes.split()

        for y in cycleLis:
            y = y.removeprefix("(")
            y = y.removesuffix(")")
            y = y.split()
            temp.append(y)
    elif type(Perm) == list:
        temp = [[str(i) for i in Perm]]
    for i in temp:
        hold = copy.deepcopy(i[0])
        i.pop(0)
        i.append(hold)
        for y in range(0, len(i)-1, 2):
            temp2.append(i[y:y+2])
        temp3.append(temp2)
        temp2 = []
    for i in temp3:
        for y in i:
            temp2.append((int(y[0]), int(y[1])))
        nodesLis.append(temp2)
        temp2 = []
    if cycleLis != []:
        return nodesLis
    else:
        return nodesLis


def GraphToGenome(Graph):
    graphDic = {}
    graphDicLen = {}
    Perm = ""
    permLis = []
    permSet = set()
    for graph in Graph:
        for edge in graph:
            if edge not in graphDic:
                graphDic[edge] = []
                graphDicLen[edge] = 0
    for graph in Graph:
        for edge in graph:
            node_1 = edge[0]
            node_2 = edge[1]
            if node_1 % 2 == 0:
                for edge_2 in graph:
                    if node_1 - 1 in edge_2:
                        graphDic[edge_2].append(edge)
                        graphDicLen[edge_2] += 1
            else:
                for edge_2 in graph:
                    if node_1 + 1 in edge_2:
                        graphDic[edge_2].append(edge)
                        graphDicLen[edge_2] += 1

            if node_2 % 2 == 0:
                for edge_2 in graph:
                    if node_2 - 1 in edge_2:
                        graphDic[edge_2].append(edge)
                        graphDicLen[edge_2] += 1
            else:
                for edge_2 in graph:
                    if node_2 + 1 in edge_2:
                        graphDic[edge_2].append(edge)
                        graphDicLen[edge_2] += 1

    while 2 in graphDicLen.values():

        for key, value in graphDicLen.items():
            if value == 2:
                start = key
                node = start[1]
                break

        if start[0] % 2 == 0:
            permSet.add(start[0] // 2)
            Perm = Perm + " " + "+" + str(start[0] // 2)
        else:
            permSet.add((start[0]+1) // 2)
            Perm = Perm + " " + "-" + str((start[0]+1) // 2)

        if start[1] % 2 == 0:
            permSet.add(start[1] // 2)
            Perm = Perm + " " + "-" + str(start[1] // 2)
        elif start[1] % 2 == 1 and start[1]+1 != start[0]:
            permSet.add((start[1]+1) // 2)
            Perm = Perm + " " + "+" + str((start[1]+1) // 2)
        Perm = Perm.lstrip()

        while len(graphDic[start]) > 1:
            for i in graphDic[start]:
                if node % 2 == 0:
                    if node - 1 in i:
                        graphDic[start].remove(i)
                        graphDicLen[start] -= 1
                        start = i
                        for i in start:
                            if i != node - 1:
                                node = i
                                if node % 2 == 0:
                                    if node // 2 not in permSet:
                                        permSet.add(node // 2)
                                        Perm = Perm + " " + \
                                            "-" + str(node // 2)
                                        break
                                    else:
                                        if len(Perm) != 0:
                                            permLis.append(Perm)
                                            Perm = ""
                                        break
                                else:
                                    if (node+1) // 2 not in permSet:
                                        permSet.add((node+1) // 2)
                                        Perm = Perm + " " + "+" + \
                                            str((node+1) // 2)
                                        break
                                    else:
                                        if len(Perm) != 0:
                                            permLis.append(Perm)
                                            Perm = ""
                                        break
                else:
                    if node + 1 in i:
                        graphDic[start].remove(i)
                        graphDicLen[start] -= 1
                        start = i
                        for i in start:
                            if i != node + 1:
                                node = i
                                if node % 2 == 0:
                                    if node // 2 not in permSet:
                                        permSet.add(node // 2)
                                        Perm = Perm + " " + \
                                            "-" + str(node // 2)
                                        break
                                    else:
                                        if len(Perm) != 0:
                                            permLis.append(Perm)
                                            Perm = ""
                                        break
                                else:
                                    if (node+1) // 2 not in permSet:
                                        permSet.add((node+1) // 2)
                                        Perm = Perm + " " + "+" + \
                                            str((node+1) // 2)
                                        break
                                    else:
                                        if len(Perm) != 0:
                                            permLis.append(Perm)
                                            Perm = ""
                                        break
    for perm in permLis:
        permRep = copy.deepcopy(perm)
        permTemp = copy.deepcopy(perm[-3:].lstrip())
        perm = perm.replace(permTemp, "")
        perm = perm.rstrip()
        if perm == "":
            perm = "(" + permTemp +")"
        else:
            perm = "(" + permTemp + " " + perm + ")"
        permLis[permLis.index(permRep)] = perm
    return ("").join(permLis)


def TwoBreakDistance(P, Q):
    breakPointGraphDic = {}
    loops = []

    # P = P.removeprefix("(")
    # P = P.removesuffix(")")
    # P = P.split(")(")

    # Q = Q.removeprefix("(")
    # Q = Q.removesuffix(")")
    # Q = Q.split(")(")

    PG = ColoredEdges(P)
    if len(PG) > 1:
        PGtemp = copy.deepcopy(PG)
        PG = []
        for i in PGtemp:
            PG.extend(i)
    else:
        PG = PG[0]
    QG = ColoredEdges(Q)
    if len(QG) > 1:
        QGtemp = copy.deepcopy(QG)
        QG = []
        for i in QGtemp:
            QG.extend(i)
    else:
        QG = QG[0]
    blocksCount = copy.deepcopy(len(PG))
    breakPointGraphDic["PG"] = PG
    breakPointGraphDic["QG"] = QG

    while len(PG) > 0 and len(QG) > 0:
        loopCount = copy.deepcopy(len(loops))
        alter = "QG"
        loopStart = PG[0][1]
        initial = copy.deepcopy(PG[0][0])
        initialSet = 0
        temp = []
        temp.append(PG[0])
        breakPointGraphDic["PG"].remove(breakPointGraphDic["PG"][0])
        while breakPointGraphDic["PG"] != [] or breakPointGraphDic["QG"] != []:
            if len(loops) > loopCount:
                break
            for i in breakPointGraphDic[alter]:
                # if loopStart == initialSet:
                #     loops.append(temp)
                #     break
                if loopStart == i[0]:
                    loopStart = i[1]
                    initialSet = copy.deepcopy(initial)
                    temp.append(i)
                    breakPointGraphDic[alter].remove(
                        breakPointGraphDic[alter][breakPointGraphDic[alter].index(i)])
                    if loopStart == initialSet:
                        loops.append(temp)
                        break
                    if alter == "PG":
                        alter = "QG"
                    elif alter == "QG":
                        alter = "PG"

                elif loopStart == i[1]:
                    loopStart = i[0]
                    initialSet = copy.deepcopy(initial)
                    temp.append(i)
                    breakPointGraphDic[alter].remove(
                        breakPointGraphDic[alter][breakPointGraphDic[alter].index(i)])
                    if loopStart == initialSet:
                        loops.append(temp)
                        break
                    if alter == "PG":
                        alter = "QG"
                    elif alter == "QG":
                        alter = "PG"

    return blocksCount - len(loops)


def TwoBreakOnGenomeGraph(Graphs, breakSeq):
    edge_1 = ""
    edge_2 = ""
    if len(Graphs) > 1:
        Graphstemp = copy.deepcopy(Graphs)
        Graphs = []
        for i in Graphstemp:
            Graphs.extend(i)
        Graphs = [Graphs]
    for graph in Graphs:
        for index, edge in enumerate(graph):
            if breakSeq[0] in edge:
                edge_1 = edge
                edge_1_index = index
            elif breakSeq[2] in edge:
                edge_2 = edge
                edge_2_index = index

    for graph in Graphs:
        graph.remove(edge_1)
        graph.remove(edge_2)
        graph.append((breakSeq[0], breakSeq[2]))
        graph.append((breakSeq[1], breakSeq[3]))
    return Graphs


def TwoBreakOnGenome(Perm, breakSeq):
    Graphs = ColoredEdges(Perm)
    twoBreakGraph = TwoBreakOnGenomeGraph(Graphs, breakSeq)
    # twoBreakChrom, cycle_v_2 = GraphToGenome(twoBreakGraph, cycle)
    twoBreakChrom = GraphToGenome(twoBreakGraph)

    return twoBreakChrom


def ShortestRearrangementScenario(P,Q):
    breakSeq = []
    steps = [P]
    perm = P
    breakPointGraph = []
    redEdges = ColoredEdges(P)
    redEdges = redEdges[0]
    blueEdges = ColoredEdges(Q)
    blueEdges = blueEdges[0]
    lenBlueEdges = copy.deepcopy(len(blueEdges))
    breakPointGraph.extend(redEdges)
    breakPointGraph.extend(blueEdges)
    counter = -1
    distance = TwoBreakDistance(P,Q)
    while distance != 0:
        # count = 0
        # Add reversals where aplicable, and check every iteration if the Randomly selected blueEdge is in redEdges
        #instead of removing redEdges
        # for i in redEdges:
        #     if i in blueEdges or i[::-1] in blueEdges:
        #         count += 1
        # if count == len(blueEdges):
        #     return steps
        randomBlue = random.choice(range(0,len(blueEdges)-1))
        blueEdge = copy.deepcopy(blueEdges[randomBlue])
        if blueEdge in redEdges or blueEdge[::-1] in redEdges:
        # if blueEdge in redEdges:
            continue
        for edge2 in redEdges:
            if len(breakSeq) == 4:
                break
            if blueEdge[0] in edge2:
                breakSeq.append(edge2[0])
                breakSeq.append(edge2[1])
            if blueEdge[1] in edge2:
                breakSeq.append(edge2[0])
                breakSeq.append(edge2[1])
        # if (breakSeq[0], breakSeq[2]) != blueEdge or (breakSeq[1], breakSeq[3]) != blueEdge:

            # swap = (breakSeq[2], breakSeq[3])[::-1]
            # breakSeq[2] = swap[0]
            # breakSeq[3] = swap[1]
        twoBreak = TwoBreakOnGenome(perm, breakSeq)
        perm = twoBreak
        redEdges = ColoredEdges(twoBreak)
        if len(redEdges) > 1:
            redEdgestemp = copy.deepcopy(redEdges)
            redEdges = []
            for i in redEdgestemp:
                redEdges.extend(i)
        else:
            redEdges = redEdges[0]
        # breakPointGraph = []
        # breakPointGraph.extend(blueEdges)
        distance = TwoBreakDistance(perm,Q)
        # for i in redEdges:
        #     if i in blueEdges or i[::-1] in blueEdges:
        #         continue
        #     else:
        #         breakPointGraph.append(i)
        steps.append(twoBreak)
        breakSeq = []
    return steps


def Shared_kmers(seq1, seq2, k):
    kmerPos = []
    kmer1Dic = {}
    for position1 in range(len(seq1)-k+1):
        kmer1 = seq1[position1:position1+k]
        if kmer1 not in kmer1Dic.keys():
            kmer1Dic[kmer1] = [position1]
        elif kmer1 in kmer1Dic.keys():
            kmer1Dic[kmer1].append(position1)

    for position2 in range(len(seq2)-k+1):
        kmer2 = seq2[position2:position2+k]
        if kmer2 in kmer1Dic.keys():
            for pos in kmer1Dic[kmer2]:
                kmerPos.append((pos, position2))

        if ReverseComplement(kmer2) in kmer1Dic.keys():
            for pos in kmer1Dic[ReverseComplement(kmer2)]:
                kmerPos.append((pos, position2))

    return len(kmerPos)
        

# print(TwoBreakOnGenome(P, breakSeq))
# print(TwoBreakOnGenomeGraph(Graph, breakSeq))
# P = P.removeprefix("(")
# P = P.removesuffix(")")
# P = P.split(")(")
# P = (" ").join(P)
# P = P.split()
# print(ColoredEdges(P))
# print(GraphToGenome(Graph))
print(ShortestRearrangementScenario(P,Q))
# print(TwoBreakDistance(P,Q))
# print(Shared_kmers(seq1, seq2, k))
# T = [(+67 -1 +2 -3 +4 +5 +6 +7 -8 +9 +10 +11 +12 +13 +14 +15 -16 +17 +18 -19 -20 -21 +22 -23 +24
#        -25 -26 +27 -28 -29 +30 -31 +32 +33 -34 +35 -36 -59 +60 +61 -62 -63 -64 +65 -66)
#      (-58 +37 +38 +39 -40 -41 +42 -43 -44 +45 +46 -47 -48 +49 +50 -51 -52 +53 +54 -55 +56 +57)]

# U = [(-67 -1 +2 -3 -4 -5 -6 -7 +8 +9 +10 +11 +12 +13 +14 +15 -16 -17 -18 +19 -20 +21 +22 -23 -24
#        +25 -26 -27 +28 -29 -30 +31 +32 +33 -34 -35 +36 -59 -60 -61 +62 -63 +64 +65 -66)
#      (+58 +37 +38 +39 -40 +41 +42 -43 +44 +45 +46 -47 +48 +49 +50 -51 +52 +53 +54 -55 -56 -57)]
# print(T == U)


