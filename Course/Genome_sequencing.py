from utilities import *
import random
from random import randint
import itertools
import copy
from copy import deepcopy
# f = readTextFile(r"C:\Users\spidy\Downloads\dataset_198_3.txt")
# l = []
# for i in range(len(f)-25+1):
#     pat = f[i:i+25]
#     l.append(pat)

import sys
cycleDic = {}
cycle = open(r"C:\Users\spidy\Downloads\dataset_203_7.txt").read(
).strip().replace(":", " ").split("\n")
for string in cycle:
    lis = string.split()
    cycleDic[lis[0]] = lis[1:]

cycleDic1 = {}
cycleNum = "0: 3 \n1: 0 \n2: 1 6 \n3: 2 \n4: 2 \n5: 4 \n6: 5 8 \n7: 9 \n8: 7 \n9: 6".strip().replace(":",
                                                                                                     "").split("\n")
cycleNum2 = "0: 2 \n1: 3 \n2: 1 \n3: 0 4 \n6: 3 7 \n7: 8 \n8: 9 \n9: 6".strip(
).replace(":", "").split("\n")
for string in cycleNum2:
    lis = string.split()
    cycleDic1[lis[0]] = lis[1:]

KmersList = open(
    r"C:\Users\spidy\Downloads\dataset_203_7.txt").read().strip().split(" ")
Pair = open(
    r"C:\Users\spidy\Downloads\dataset_204_16.txt").read().strip().split(" ")

assessment = ["AAAT", "AATG", "ACCC", "ACGC", 'ATAC', 'ATCA', 'ATGC', 'CAAA', "CACC",
              "CATA", "CATC", "CCAG", "CCCA", 'CGCT', 'CTCA', 'GCAT', "GCTC", 'TACG', 'TCAC', 'TCAT', 'TGCA']
assessment2 = ['ACC|ATA', 'ACT|ATT', 'ATA|TGA', 'ATT|TGA', 'CAC|GAT', 'CCG|TAC', 'CGA|ACT', 'CTG|AGC',
               'CTG|TTC', 'GAA|CTT', 'GAT|CTG', 'GAT|CTG', 'TAC|GAT', 'TCT|AAG', 'TGA|GCT', 'TGA|TCT', 'TTC|GAA']
test = "6879654210326"
cyc = {'00': ['00', '01'], '10': ['00', '01'],
       '01': ['10', '11'], '11': ['10', '11']}
DNA = "CAATCCAAC"
DNA2 = "AAGATTCTCTAAGA"
DNA3 = "TAATGCCATGGGATGTT"
Pairs = ["GACC|GCGC", "ACCG|CGCC", "CCGA|GCCG", "CGAG|CCGG", "GAGC|CGGA"]
Pairs2 = ["GAGA|TTGA", "TCGT|GATG", "CGTG|ATGT", "TGGT|TGAG",
          "GTGA|TGTT", "GTGG|GTGA", "TGAG|GTTG", "GGTC|GAGA", "GTCG|AGAT"]
Kmers = ["AAT", "ATG", "ATG", "ATG", "CAT", "CCA", "GAT",
         "GCC", "GGA", "GGG", "GTT", "TAA", "TGC", "TGG", "TGT"]
kmers2 = ["ACCGA", "CCGAA", "CGAAG", "GAAGC", "AAGCT"]
kmers3 = ["ATGCG", "GCATG", "CATGC", "AGGCA", "GGCAT", "GGCAC"]
kmers4 = ["GAGG", "CAGG", "GGGG", "GGGA", "CAGG", "AGGG", "GGAG"]
kmers5 = ["ATG", "ATG", "TGT", "TGG", "CAT", "GGA", "GAT", "AGA"]
k = 12


def StringCompositionWithPossibleGaps(DNA, k, d=0):
    kmers = []
    if d == 0:
        kmers = [DNA[i:i+k] for i in range(len(DNA)-k+1)]
        kmers.sort()
    else:
        for i in range(len(DNA)-(2*k)+d+1):
            if len(DNA[i:i+k] + DNA[i+k+d:i+(2*k)+d]) == 2*k:
                kmers.append(DNA[i:i+k] + "|" + DNA[i+k+d:i+(2*k)+d])
    kmers.sort()
    return ' '.join(str(i) for i in kmers).split(" ")


def OverlapGraph(Kmers, k):
    kDict = {}
    KV = ""
    for kmer in Kmers:
        kDict[kmer] = ""
    for kmer in Kmers:
        for kmer2 in Kmers:
            if kmer[1:] == kmer2[:-1]:
                kDict[kmer] += kmer2 + ' '
    for key, value in kDict.items():
        if value != '':
            KV = KV + (str(key + ":"+" " + value)) + "\n"
    return KV


def DeBrujinGraph(Kmers, k):
    tupLis = []
    elements = {}
    graph = []
    if type(Kmers) == str:
        for i in range(len(Kmers)-k+1):
            pat = Kmers[i:i+k]
            tupLis.append((pat[:-1], pat[1:]))
    elif type(Kmers) == list:
        if type(Kmers[0]) == str:
            k = len(Kmers[0])
            for pat in Kmers:
                tupLis.append((pat[:-1], pat[1:]))
        elif type(Kmers[0]) == list:
            for pairs in Kmers:
                tupLis.append(
                    (pairs[0][:-1] + pairs[1][:-1], pairs[0][1:] + pairs[1][1:]))

    for prefix, suffix in tupLis:
        if prefix not in elements:
            elements[prefix] = [suffix]
        elif prefix in elements:
            elements[prefix].append(suffix)
    for prefix, suffix in elements.items():
        graph.append((prefix, suffix))
    return elements


def cycleFinder(graph):
    stack = []
    cycle = []
    location = random.choice(list(graph.keys()))
    # location = "0"
    rnd = ""
    while max(graph.values()) != [] or stack != []:
        if graph[location] != []:
            stack.append(location)
            rnd = location
            location = graph[location][random.randint(
                0, len(graph[location])-1)]
            # location = graph[location][0]
            graph[rnd].remove(location)

        else:
            cycle.append(location)
            for i in reversed(stack):
                location = i
                if graph[location] == []:
                    cycle.append(location)
                    stack.pop((len(stack))-1)
                    continue
                else:
                    stack.pop((len(stack))-1)
                    break
    cycle1 = []
    for i in reversed(cycle):
        cycle1.append(i)
    # cycle = ' '.join(str(i) for i in reversed(cycle))
    # cycle = open(r"C:\Work and Education\Python39\Bioinformatics\Eulerian_Cycle.txt", mode = "w")
    # cycle.write(' '.join(str(i) for i in reversed(cycle)))
    return cycle1


def EulerianGrader(graph, eulerType="cycle"):
    graph_for_method = copy.deepcopy(graph)
    grade = []
    if eulerType == "cycle":
        path_list = cycleFinder(graph_for_method)
        for index in range(len(path_list) - 1):
            elem = path_list[index]
            next_elem = path_list[index + 1]
            if elem not in graph.keys():
                print("invalid elem: " + str(elem) + " at posn " + str(index))
                grade.append("src")
            elif next_elem not in graph[elem]:
                print("invalid edge: " + str(elem) + " to " + str(next_elem))
                grade.append("dst")
            else:
                graph[elem].remove(next_elem)
                if not graph[elem]:
                    del graph[elem]
                grade.append("OK")
    else:
        path_list = pathFinder(graph_for_method)
        for index in range(len(path_list) - 1):
            elem = path_list[index]
            next_elem = path_list[index + 1]
            if next_elem not in graph[elem]:
                print("invalid edge: " + str(elem) + " to " + str(next_elem))
                grade.append("dst")
            else:
                graph[elem].remove(next_elem)
                if not graph[elem]:
                    del graph[elem]
                grade.append("OK")
    if not graph:
        print("emptied!")
    else:
        print(graph)
    return grade


def pathFinder(graph):
    startNode = ''
    endNode = ''
    startWith = []
    nodeDegrees = {}
    path = []
    stack = []

    for prefix, sufix in graph.items():
        if prefix not in nodeDegrees:
            nodeDegrees[prefix] = len(sufix)
    for sufix in graph.values():
        for i in sufix:
            if i not in nodeDegrees:
                nodeDegrees[i] = 1
            else:
                nodeDegrees[i] += 1
    for node in nodeDegrees:
        if nodeDegrees[node] % 2 != 0:
            startWith.append(node)
    for node in startWith:
        if node not in graph.keys():
            endNode = node
        elif (nodeDegrees[node] - len(graph[node])) < len(graph[node]):
            startNode = node
        else:
            endNode = node
    if (startNode and endNode) == '':
        return graph

    location = startNode
    rnd = ""
    while max(graph.values()) != [] or stack != []:
        try:
            if graph[location] != []:
                stack.append(location)
                rnd = location
                # location = graph[location][random.randint(0,len(graph[location])-1)]
                location = graph[location][0]
                graph[rnd].remove(location)
            else:
                path.append(location)
                for i in reversed(stack):
                    location = i
                    if location not in graph:
                        path.append(location)
                        stack.pop((len(stack))-1)
                        continue
                    elif graph[location] == []:
                        path.append(location)
                        stack.pop((len(stack))-1)
                        continue
                    else:
                        stack.pop((len(stack))-1)
                        break
        except KeyError:
            path.append(location)
            for i in reversed(stack):
                location = i
                if location not in graph:
                    path.append(location)
                    stack.pop((len(stack))-1)
                    continue
                elif graph[location] == []:
                    path.append(location)
                    stack.pop((len(stack))-1)
                    continue
                else:
                    stack.pop((len(stack))-1)
                    break

    path1 = []
    for i in reversed(path):
        path1.append(i)
    # path = ' '.join(str(i) for i in reversed(path))
    # path = open(r"C:\Work and Education\Python39\Bioinformatics\Eulerian_Cycle.txt", mode = "w")
    # path.write(' '.join(str(i) for i in reversed(path)))
    return path1


def stringReconstruction(Kmers, k, dataType=list):
    l = 0
    if type(Kmers) == list:
        k = len(Kmers[0])
    dataType = type(Kmers)
    stringCheck = False
    graphSaved = {}
    reconstructed = ""

    while stringCheck == False:
        if dataType == dict:
            graphSaved = copy.deepcopy(Kmers)
            path = pathFinder(Kmers)
        else:
            graph = DeBrujinGraph(Kmers, k)
            graphSaved = copy.deepcopy(graph)
            path = pathFinder(graph)
        if type(path) == dict:
            cycle = cycleFinder(path)
            for i in cycle:
                reconstructed = reconstructed + i[-1]
            check = EulerianGrader(graphSaved, "cycle")
            if set(check) != {"OK"}:
                # l += 1
                # print(l)
                continue
        else:
            if len(path[0]) > 1:
                reconstructed = path[0][0:-1]
            else:
                reconstructed = path[0]
                path.pop(path.index(path[0]))
            for i in path:
                reconstructed = reconstructed + i[-1]
            check = EulerianGrader(graphSaved, "path")
            if set(check) != {"OK"}:
                # l += 1
                # print(l)
                continue
        stringCheck = True
    return reconstructed

# Incomplete:


def k_universal(k):
    binary = ["0", "1"]
    kmers = []
    num = "0" * k
    kmers.append(num)
    while len(kmers) < 2**k:
        for i in kmers:
            for j in range(len(num)):
                for b in binary:
                    nu = i[:j] + b + i[j+1:]
                    if nu not in kmers:
                        kmers.append(nu)
    num1 = stringReconstruction(kmers, k, list)

    return num1


def StringsSpelledByGappedKmers(Kmers, k, d):
    dis = k+d
    n = len(Kmers)
    finalString = ""
    firstKmers = []
    secondKmers = []

    concatenatedKmers = []
    for pair in Kmers:
        cutOff = pair.index("|")
        concatenatedKmers.append([pair[:cutOff], pair[cutOff+1:]])
    graph = DeBrujinGraph(concatenatedKmers, k)
    path = pathFinder(graph)

    for extendedKmers in path:
        firstKmers.append(extendedKmers[:len(extendedKmers)//2])
        secondKmers.append(extendedKmers[len(extendedKmers)//2:])

    firstString = firstKmers[0][0:-1]
    for i in firstKmers:
        firstString = firstString + i[-1]

    secondString = secondKmers[0][0:-1]
    for i in secondKmers:
        secondString = secondString + i[-1]

    if firstString[dis:] == secondString[:len(firstString) - dis]:
        finalString = firstString + secondString[len(firstString) - dis:]

        return finalString
    else:
        return "there is no string spelled by the gapped patterns"

# rebuild contigsFinder and : use startNodes from pathFinder as starting points for contigs.
# use non-start or end nodes with cycle finder to find isolated cycle contigs with the condition that the nodes must have 1 indegree and 1 outdegree only.
# output should be AGA ATG ATG CAT GAT TGGA TGT


def contigsFinder(graph):
    startNode = []
    endNode = []
    startWith = []
    nodeDegrees = {}
    contigs = []
    rebuiltContigs = []
    location = ""
    stack = []
    for prefix, sufix in graph.items():
        if prefix not in nodeDegrees:
            nodeDegrees[prefix] = len(sufix)
    for sufix in graph.values():
        for i in sufix:
            if i not in nodeDegrees:
                nodeDegrees[i] = 1
            else:
                nodeDegrees[i] += 1
    for node in nodeDegrees:
        if nodeDegrees[node] % 2 != 0:
            startWith.append(node)
    for node in startWith:
        if node not in graph.keys():
            endNode.append(node)
        elif (nodeDegrees[node] - len(graph[node])) < len(graph[node]):
            startNode.append(node)
        else:
            endNode.append(node)
    # contigs = [list() for i in startNode]
    for start in startNode:
        location = start
        rnd = ""
        rnd2 = ""
        while max(graph.values()) != [] or stack != []:
            try:
                if graph[location] != []:
                    stack.append(location)
                    rnd = location
                    location = graph[location][random.randint(
                        0, len(graph[location])-1)]
                    rnd2 = graph[rnd][graph[rnd].index(location)]
                    graph[rnd].pop(graph[rnd].index(
                        graph[rnd][graph[rnd].index(rnd2)]))
                else:
                    contigs.append(location)
                    for i in reversed(stack):
                        location = i
                        if location not in graph:
                            contigs.append(location)
                            stack.pop((len(stack))-1)
                            continue
                        elif graph[location] == []:
                            contigs.append(location)
                            stack.pop((len(stack))-1)
                            continue
                        else:
                            stack.pop(stack.index(i))
                            break
            except KeyError:
                contigs.append(location)
                for i in reversed(stack):
                    location = i
                    if location not in graph:
                        contigs.append(location)
                        stack.pop((len(stack))-1)
                        continue
                    elif graph[location] == []:
                        contigs.append(location)
                        stack.pop((len(stack))-1)
                        continue
                    else:
                        stack.pop(stack.index(i))
                        break

    return contigs


# print(OverlapGraph(kmers5, 3))
# print(num2)
# print(cycleDic)
# print(stringReconstruction(assessment, 4))
# print(cycleFinder(cycleDic1))
# print(k_universal(9))
d = DeBrujinGraph(kmers5, 3)
# print(contigsFinder(d))
# print(StringCompositionWithPossibleGaps(DNA3, 4, 1))
# print(StringsSpelledByGappedKmers(assessment2, 3, 1))

# print(pathFinder({"0": ['1', '2'], "2": ['3'], "1": ['4'], "3": ['4']}))
# print(cycleFinder(cycleDic))
# print(EulerianGrader(cycleDic))
