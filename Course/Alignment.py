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

sys.setrecursionlimit(10000)

array = open(r"C:\Users\spidy\Downloads\alignmentTest.txt")

DAG = open(r"C:\Users\spidy\Downloads\DAG.txt")

PAM = {'A': {'A': 2, 'C': -2, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1, 'H': -1, 'K': -1, 'M': -1, 'L': -2, 'N': 0, 'Q': 0, 'P': 1, 'S': 1, 'R': -2, 'T': 1, 'W': -6, 'V': 0, 'Y': -3}, 'C': {'A': -2, 'C': 12, 'E': -5, 'D': -5, 'G': -3, 'F': -4, 'I': -2, 'H': -3, 'K': -5, 'M': -5, 'L': -6, 'N': -4, 'Q': -5, 'P': -3, 'S': 0, 'R': -4, 'T': -2, 'W': -8, 'V': -2, 'Y': 0}, 'E': {'A': 0, 'C': -5, 'E': 4, 'D': 3, 'G': 0, 'F': -5, 'I': -2, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': 2, 'P': -1, 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4}, 'D': {'A': 0, 'C': -5, 'E': 3, 'D': 4, 'G': 1, 'F': -6, 'I': -2, 'H': 1, 'K': 0, 'M': -3, 'L': -4, 'N': 2, 'Q': 2, 'P': -1, 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4}, 'G': {'A': 1, 'C': -3, 'E': 0, 'D': 1, 'G': 5, 'F': -5, 'I': -3, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -1, 'P': 0, 'S': 1, 'R': -3, 'T': 0, 'W': -7, 'V': -1, 'Y': -5}, 'F': {'A': -3, 'C': -4, 'E': -5, 'D': -6, 'G': -5, 'F': 9, 'I': 1, 'H': -2, 'K': -5, 'M': 0, 'L': 2, 'N': -3, 'Q': -5, 'P': -5, 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -1, 'Y': 7}, 'I': {'A': -1, 'C': -2, 'E': -2, 'D': -2, 'G': -3, 'F': 1, 'I': 5, 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -2, 'S': -1, 'R': -2, 'T': 0, 'W': -5, 'V': 4, 'Y': -1}, 'H': {'A': -1, 'C': -3, 'E': 1, 'D': 1, 'G': -2, 'F': -2, 'I': -2, 'H': 6, 'K': 0, 'M': -2, 'L': -2, 'N': 2, 'Q': 3, 'P': 0, 'S': -1, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': 0}, 'K': {'A': -1, 'C': -5, 'E': 0, 'D': 0, 'G': -2, 'F': -5, 'I': -2, 'H': 0, 'K': 5, 'M': 0, 'L': -3, 'N': 1, 'Q': 1, 'P': -1, 'S': 0, 'R': 3, 'T': 0, 'W': -3, 'V': -2, 'Y': -4}, 'M': {'A': -1, 'C': -5, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 2, 'H': -2, 'K': 0, 'M': 6, 'L': 4, 'N': -2, 'Q': -1, 'P': -2, 'S': -2, 'R': 0, 'T': -1, 'W': -4, 'V': 2, 'Y': -2},
       'L': {'A': -2, 'C': -6, 'E': -3, 'D': -4, 'G': -4, 'F': 2, 'I': 2, 'H': -2, 'K': -3, 'M': 4, 'L': 6, 'N': -3, 'Q': -2, 'P': -3, 'S': -3, 'R': -3, 'T': -2, 'W': -2, 'V': 2, 'Y': -1}, 'N': {'A': 0, 'C': -4, 'E': 1, 'D': 2, 'G': 0, 'F': -3, 'I': -2, 'H': 2, 'K': 1, 'M': -2, 'L': -3, 'N': 2, 'Q': 1, 'P': 0, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -2, 'Y': -2}, 'Q': {'A': 0, 'C': -5, 'E': 2, 'D': 2, 'G': -1, 'F': -5, 'I': -2, 'H': 3, 'K': 1, 'M': -1, 'L': -2, 'N': 1, 'Q': 4, 'P': 0, 'S': -1, 'R': 1, 'T': -1, 'W': -5, 'V': -2, 'Y': -4}, 'P': {'A': 1, 'C': -3, 'E': -1, 'D': -1, 'G': 0, 'F': -5, 'I': -2, 'H': 0, 'K': -1, 'M': -2, 'L': -3, 'N': 0, 'Q': 0, 'P': 6, 'S': 1, 'R': 0, 'T': 0, 'W': -6, 'V': -1, 'Y': -5}, 'S': {'A': 1, 'C': 0, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1, 'H': -1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': -1, 'P': 1, 'S': 2, 'R': 0, 'T': 1, 'W': -2, 'V': -1, 'Y': -3}, 'R': {'A': -2, 'C': -4, 'E': -1, 'D': -1, 'G': -3, 'F': -4, 'I': -2, 'H': 2, 'K': 3, 'M': 0, 'L': -3, 'N': 0, 'Q': 1, 'P': 0, 'S': 0, 'R': 6, 'T': -1, 'W': 2, 'V': -2, 'Y': -4}, 'T': {'A': 1, 'C': -2, 'E': 0, 'D': 0, 'G': 0, 'F': -3, 'I': 0, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 0, 'Q': -1, 'P': 0, 'S': 1, 'R': -1, 'T': 3, 'W': -5, 'V': 0, 'Y': -3}, 'W': {'A': -6, 'C': -8, 'E': -7, 'D': -7, 'G': -7, 'F': 0, 'I': -5, 'H': -3, 'K': -3, 'M': -4, 'L': -2, 'N': -4, 'Q': -5, 'P': -6, 'S': -2, 'R': 2, 'T': -5, 'W': 17, 'V': -6, 'Y': 0}, 'V': {'A': 0, 'C': -2, 'E': -2, 'D': -2, 'G': -1, 'F': -1, 'I': 4, 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -1, 'S': -1, 'R': -2, 'T': 0, 'W': -6, 'V': 4, 'Y': -2}, 'Y': {'A': -3, 'C': 0, 'E': -4, 'D': -4, 'G': -5, 'F': 7, 'I': -1, 'H': 0, 'K': -4, 'M': -2, 'L': -1, 'N': -2, 'Q': -4, 'P': -5, 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -2, 'Y': 10}}
BLOSUM62 = {'A': {'A': 4, 'C': 0, 'E': -1, 'D': -2, 'G': 0, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 0, 'W': -3, 'V': 0, 'Y': -2}, 'C': {'A': 0, 'C': 9, 'E': -4, 'D': -3, 'G': -3, 'F': -2, 'I': -1, 'H': -3, 'K': -3, 'M': -1, 'L': -1, 'N': -3, 'Q': -3, 'P': -3, 'S': -1, 'R': -3, 'T': -1, 'W': -2, 'V': -1, 'Y': -2}, 'E': {'A': -1, 'C': -4, 'E': 5, 'D': 2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': -2, 'L': -3, 'N': 0, 'Q': 2, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -3, 'V': -2, 'Y': -2}, 'D': {'A': -2, 'C': -3, 'E': 2, 'D': 6, 'G': -1, 'F': -3, 'I': -3, 'H': -1, 'K': -1, 'M': -3, 'L': -4, 'N': 1, 'Q': 0, 'P': -1, 'S': 0, 'R': -2, 'T': -1, 'W': -4, 'V': -3, 'Y': -3}, 'G': {'A': 0, 'C': -3, 'E': -2, 'D': -1, 'G': 6, 'F': -3, 'I': -4, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -2, 'P': -2, 'S': 0, 'R': -2, 'T': -2, 'W': -2, 'V': -3, 'Y': -3}, 'F': {'A': -2, 'C': -2, 'E': -3, 'D': -3, 'G': -3, 'F': 6, 'I': 0, 'H': -1, 'K': -3, 'M': 0, 'L': 0, 'N': -3, 'Q': -3, 'P': -4, 'S': -2, 'R': -3, 'T': -2, 'W': 1, 'V': -1, 'Y': 3}, 'I': {'A': -1, 'C': -1, 'E': -3, 'D': -3, 'G': -4, 'F': 0, 'I': 4, 'H': -3, 'K': -3, 'M': 1, 'L': 2, 'N': -3, 'Q': -3, 'P': -3, 'S': -2, 'R': -3, 'T': -1, 'W': -3, 'V': 3, 'Y': -1}, 'H': {'A': -2, 'C': -3, 'E': 0, 'D': -1, 'G': -2, 'F': -1, 'I': -3, 'H': 8, 'K': -1, 'M': -2, 'L': -3, 'N': 1, 'Q': 0, 'P': -2, 'S': -1, 'R': 0, 'T': -2, 'W': -2, 'V': -3, 'Y': 2}, 'K': {'A': -1, 'C': -3, 'E': 1, 'D': -1, 'G': -2, 'F': -3, 'I': -3, 'H': -1, 'K': 5, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -1, 'S': 0, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': -2}, 'M': {'A': -1, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 1, 'H': -2, 'K': -1, 'M': 5, 'L': 2, 'N': -2, 'Q': 0, 'P': -2, 'S': -1, 'R': -1, 'T': -1, 'W': -1, 'V': 1,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     'Y': -1}, 'L': {'A': -1, 'C': -1, 'E': -3, 'D': -4, 'G': -4, 'F': 0, 'I': 2, 'H': -3, 'K': -2, 'M': 2, 'L': 4, 'N': -3, 'Q': -2, 'P': -3, 'S': -2, 'R': -2, 'T': -1, 'W': -2, 'V': 1, 'Y': -1}, 'N': {'A': -2, 'C': -3, 'E': 0, 'D': 1, 'G': 0, 'F': -3, 'I': -3, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 6, 'Q': 0, 'P': -2, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -3, 'Y': -2}, 'Q': {'A': -1, 'C': -3, 'E': 2, 'D': 0, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': 0, 'L': -2, 'N': 0, 'Q': 5, 'P': -1, 'S': 0, 'R': 1, 'T': -1, 'W': -2, 'V': -2, 'Y': -1}, 'P': {'A': -1, 'C': -3, 'E': -1, 'D': -1, 'G': -2, 'F': -4, 'I': -3, 'H': -2, 'K': -1, 'M': -2, 'L': -3, 'N': -2, 'Q': -1, 'P': 7, 'S': -1, 'R': -2, 'T': -1, 'W': -4, 'V': -2, 'Y': -3}, 'S': {'A': 1, 'C': -1, 'E': 0, 'D': 0, 'G': 0, 'F': -2, 'I': -2, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 1, 'Q': 0, 'P': -1, 'S': 4, 'R': -1, 'T': 1, 'W': -3, 'V': -2, 'Y': -2}, 'R': {'A': -1, 'C': -3, 'E': 0, 'D': -2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 2, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -2, 'S': -1, 'R': 5, 'T': -1, 'W': -3, 'V': -3, 'Y': -2}, 'T': {'A': 0, 'C': -1, 'E': -1, 'D': -1, 'G': -2, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': 0, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 5, 'W': -2, 'V': 0, 'Y': -2}, 'W': {'A': -3, 'C': -2, 'E': -3, 'D': -4, 'G': -2, 'F': 1, 'I': -3, 'H': -2, 'K': -3, 'M': -1, 'L': -2, 'N': -4, 'Q': -2, 'P': -4, 'S': -3, 'R': -3, 'T': -2, 'W': 11, 'V': -3, 'Y': 2}, 'V': {'A': 0, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': -1, 'I': 3, 'H': -3, 'K': -2, 'M': 1, 'L': 1, 'N': -3, 'Q': -2, 'P': -2, 'S': -2, 'R': -3, 'T': 0, 'W': -3, 'V': 4, 'Y': -1}, 'Y': {'A': -2, 'C': -2, 'E': -2, 'D': -3, 'G': -3, 'F': 3, 'I': -1, 'H': 2, 'K': -2, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -3, 'S': -2, 'R': -2, 'T': -2, 'W': 2, 'V': -1, 'Y': 7}}

# coins = [50,25,20,10,5,1]
coins = [14, 6, 5, 3, 1]

# Create a DAG using networkx

first = DAG.readline().strip()
SandE = first.split()
start = SandE[0]
end = SandE[1]
graph = nx.DiGraph()
for line in DAG:
    line2 = line.split()
    graph.add_edge(int(line2[0]), int(line2[1]), weight=int(line2[2]))

source_node = start

# v = "CTT"
# w = "AGCATAAAGCATT"
v = "GAT"
w = "AT"


def coinChange(coins, amount):
    dp = [float('inf')] * (amount + 1)
    dp[0] = 0

    for i in range(1, amount + 1):
        for coin in coins:
            if coin <= i:
                dp[i] = min(dp[i], dp[i - coin] + 1)

    return dp[amount] if dp[amount] != float('inf') else -1


def inputManipulator(file):
    # Defining first variables
    list1 = []
    list2 = []
    swap = float('inf')

    # Packing the coordinates from line 0
    fnm = array.readline().strip()
    snm = fnm.split()
    snm = list(map(int, snm))

    # Separating matricies into two lists
    for linepos, line in enumerate(array):
        if line.strip() == "-":
            swap = linepos
            continue
        if linepos > swap:
            list2.append([line.strip()])
        elif linepos < swap:
            list1.append([line.strip()])

    for pos, i in enumerate(list1):
        for y in i:
            i = y.split()
            i = map(int, i)
            list1[pos] = list(i)
    for pos, i in enumerate(list2):
        for y in i:
            i = y.split()
            i = map(int, i)
            list2[pos] = list(i)

    # Packing the results
    package = (list1, list2, snm)
    return package


def ManhattanTourist(array):

    # Unpacking matricies and coordinates
    package = inputManipulator(array)
    unpack1, unpack2, nm = package
    x = np.array(unpack2)
    y = np.array(unpack1)
    n, m = nm

    # Main code
    mainMatrix = np.zeros((n+1)*(m+1))
    mainMatrix.shape = (n+1, m+1)
    i = 0
    j = 0
    s = (i, j)

    # For first row and first column
    for position in range(0, m):
        j += 1
        mainMatrix[i, j] = mainMatrix[i, j-1] + x[i, j-1]
    j = 0
    for position in range(0, n):
        i += 1
        mainMatrix[i, j] = mainMatrix[i-1, j] + y[i-1, j]
    i = 0

    # For the rest
    for position1 in range(0, n):
        i += 1
        j = 0
        for position2 in range(0, m):
            j += 1
            if mainMatrix[i-1, j] + y[i-1, j] > mainMatrix[i, j-1] + x[i, j-1]:
                mainMatrix[i, j] = mainMatrix[i-1, j] + y[i-1, j]
            else:
                mainMatrix[i, j] = mainMatrix[i, j-1] + x[i, j-1]
    return mainMatrix

# This function returns a graph with the move made based on the rule stated inside (Right move, Down move, or Diagonal move)


def LCSBackTrack(v, w):
    """This function returns a graph with the move made based on the rule stated inside (Right move, Down move, or Diagonal move)"""
    score = np.zeros((len(v)+1)*(len(w)+1))
    score.shape = ((len(v)+1), (len(w)+1))
    Backtrack = np.zeros((len(v)+1)*(len(w)+1))
    Backtrack.shape = ((len(v)+1), (len(w)+1))

    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            match = 0
            if v[i-1] == w[j-1]:
                match = 1
            score[i, j] = max((score[i-1, j]), (score[i, j-1]),
                            (score[i-1, j-1]+match))
            if score[i, j] == score[i-1, j]:
                Backtrack[i, j] = 3
            elif score[i, j] == score[i, j-1]:
                Backtrack[i, j] = 4
            elif score[i, j] == score[i-1, j-1] + match:
                Backtrack[i, j] = 5
    return Backtrack

# Recursive version of the alignment function


def OutputLCS(backtrack, v, i, j):
    """Recursive version of the alignment function"""
    if i == 0 or j == 0:
        return ""
    if backtrack[i, j] == 3:
        return OutputLCS(backtrack, v, i-1, j)
    elif backtrack[i, j] == 4:
        return OutputLCS(backtrack, v, i, j-1)
    else:
        return OutputLCS(backtrack, v, i-1, j-1) + v[i-1]

# Iterative version of the alignment function


def IterativeOutputLCS(Backtrack, v, w):
    """Iterative version of the alignment function"""
    LCS = ""
    i = len(v)
    j = len(w)
    while i > 0 and j > 0:
        if Backtrack[i, j] == 3:
            i = i-1
        elif Backtrack[i, j] == 4:
            j = j-1
        else:
            LCS = LCS + v[i-1]
            i = i-1
            j = j-1
    return LCS[::-1]


def GlobalAlignment_LocalAlignment(seq_1, seq_2, match_reward, mismatch_penalty, gap_penalty, alignment="global"):
    mismatch_penalty = -mismatch_penalty
    gap_penalty = -gap_penalty

    # Matrices
    main_matrix = np.zeros((len(seq_1)+1, len(seq_2)+1), dtype=int)
    match_checker_matrix = np.zeros((len(seq_1), len(seq_2)), dtype=int)

    # Matrix filling
    for i in range(len(seq_1)):
        for j in range(len(seq_2)):
            if seq_1[i] == seq_2[j]:
                match_checker_matrix[i][j] = match_reward
            else:
                match_checker_matrix[i][j] = mismatch_penalty

    # for i in range(len(seq_1)):
    #     for j in range(len(seq_2)):
    #         match_checker_matrix[i][j] = PAM[seq_1[i]][seq_2[j]]
    # print(match_checker_matrix)

    # STEP 1: Initialisation in global
    if alignment == "global":
        for i in range(1, len(seq_1)+1):
            main_matrix[i][0] = i * gap_penalty
        for j in range(len(seq_2)+1):
            main_matrix[0][j] = j * gap_penalty

    # Initialisation in local
    if alignment == "local":
        for i in range(len(seq_1)+1):
            main_matrix[i][0] = 0
        for j in range(len(seq_2)+1):
            main_matrix[0][j] = 0
    # print(main_matrix)

    # STEP 2: Matrix Filling in global
    if alignment == "global":
        for i in range(1, len(seq_1)+1):
            for j in range(1, len(seq_2)+1):
                # print(main_matrix)
                main_matrix[i][j] = max(main_matrix[i-1][j-1] + match_checker_matrix[i-1][j-1],
                                        main_matrix[i-1][j] + gap_penalty,
                                        main_matrix[i][j-1] + gap_penalty)

    # Matrix filling in local
    if alignment == "local":
        max_score_cell = ()
        max_score = 0
        for i in range(1, len(seq_1)+1):
            for j in range(1, len(seq_2)+1):
                print(main_matrix)
                main_matrix[i][j] = max(0,
                                        main_matrix[i-1][j-1] +
                                        match_checker_matrix[i-1][j-1],
                                        main_matrix[i-1][j] + gap_penalty,
                                        main_matrix[i][j-1] + gap_penalty)
                if main_matrix[i, j] > max_score:
                    max_score = main_matrix[i, j]
                    max_score_cell = (i, j)
    # print(main_matrix)
    # STEP 3: Global traceback
    aligned_1 = ""
    aligned_2 = ""
    if alignment == "global":
        ti = len(seq_1)
        tj = len(seq_2)

        while (ti > 0 or tj > 0):
            if (ti > 0 and tj > 0 and main_matrix[ti][tj] == main_matrix[ti-1][tj-1] + match_checker_matrix[ti-1][tj-1]):
                aligned_1 = seq_1[ti-1] + aligned_1
                aligned_2 = seq_2[tj-1] + aligned_2

                ti = ti - 1
                tj = tj - 1

            elif (ti > 0 and main_matrix[ti][tj] == main_matrix[ti-1][tj] + gap_penalty):
                aligned_1 = seq_1[ti-1] + aligned_1
                aligned_2 = "-" + aligned_2

                ti = ti - 1

            else:
                aligned_1 = "-" + aligned_1
                aligned_2 = seq_2[tj-1] + aligned_2

                tj = tj - 1
    # Slow recursive approach to store all possible solutions

    # if alignment == "global":
        # alignments = []
        # scores = []
        # def traceback(ti, tj, aligned_1, aligned_2):
        #         if ti == 0 and tj == 0:
        #             if not len(alignments) > 0:
        #                 totalScore = 0
        #                 for i in range(0, len(aligned_1)):
        #                     if aligned_1[i] == aligned_2[i]:
        #                         totalScore += match_reward
        #                         # totalScore += PAM[aligned_1[i]][aligned_2[i]]
        #                     elif aligned_1[i] == "-" or aligned_2[i] == "-":
        #                         totalScore += gap_penalty
        #                     else:
        #                         totalScore += mismatch_penalty
        #                         # totalScore += PAM[aligned_1[i]][aligned_2[i]]
        #                 scores.append(totalScore)
        #                 alignments.append((aligned_1, aligned_2, totalScore))
        #             else:
        #                 totalScore = 0
        #                 for i in range(0, len(aligned_1)):
        #                     if aligned_1[i] == aligned_2[i]:
        #                         totalScore += match_reward
        #                         # totalScore += PAM[aligned_1[i]][aligned_2[i]]
        #                     elif aligned_1[i] == "-" or aligned_2[i] == "-":
        #                         totalScore += gap_penalty
        #                     else:
        #                         totalScore += mismatch_penalty
        #                         # totalScore += PAM[aligned_1[i]][aligned_2[i]]
        #                 if totalScore > max(scores):
        #                     scores.append(totalScore)
        #                     alignments.append((aligned_1, aligned_2, totalScore))
        #             return

        #         if ti > 0 and tj > 0 and main_matrix[ti][tj] == main_matrix[ti-1][tj-1] + match_checker_matrix[ti-1][tj-1]:
        #             traceback(ti-1, tj-1, seq_1[ti-1] + aligned_1, seq_2[tj-1] + aligned_2)

        #         if ti > 0 and main_matrix[ti][tj] == main_matrix[ti-1][tj] + gap_penalty:
        #             traceback(ti-1, tj, seq_1[ti-1] + aligned_1, '-' + aligned_2)

        #         if tj > 0 and main_matrix[ti][tj] == main_matrix[ti][tj-1] + gap_penalty:
        #             traceback(ti, tj-1, '-' + aligned_1, seq_2[ti-1] + aligned_2)

        # traceback(len(seq_1), len(seq_2), '', '')

        # HighScore = max(scores)
        # FinalAlignments = set()
        # for alignment in alignments:
        #     align1,align2,score = alignment
        #     if score == HighScore:
        #         FinalAlignments.add(alignment)
        # return FinalAlignments

    # Local variation of traceback
    if alignment == "local":
        ti, tj = max_score_cell

        while (ti > 0 and tj > 0):
            if main_matrix[ti][tj] == 0:
                break
            if (ti > 0 and tj > 0 and main_matrix[ti][tj] == main_matrix[ti-1][tj-1] + match_checker_matrix[ti-1][tj-1]):
                aligned_1 = seq_1[ti-1] + aligned_1
                aligned_2 = seq_2[tj-1] + aligned_2

                ti = ti - 1
                tj = tj - 1

            elif (ti > 0 and main_matrix[ti][tj] == main_matrix[ti-1][tj] + gap_penalty):
                aligned_1 = seq_1[ti-1] + aligned_1
                aligned_2 = "-" + aligned_2

                ti = ti - 1
            else:
                aligned_1 = "-" + aligned_1
                aligned_2 = seq_2[tj-1] + aligned_2

                tj = tj - 1

    totalScore = 0

    for i in range(0, len(aligned_1)):
        if aligned_1[i] == aligned_2[i]:
            totalScore += match_reward
            # totalScore += PAM[aligned_1[i]][aligned_2[i]]
        elif aligned_1[i] == "-" or aligned_2[i] == "-":
            totalScore += gap_penalty

        else:
            totalScore += mismatch_penalty
            # totalScore += PAM[aligned_1[i]][aligned_2[i]]

    s, al1, al2 = str(totalScore), "\n" + aligned_1, "\n" + aligned_2
    return s, al1, al2

    # return main_matrix, np.where(main_matrix == main_matrix.max())


def find_longest_path(graph, source):
    # Perform topological sorting
    topological_order = list(nx.topological_sort(graph))

    # Initialize distance array
    distance = {node: float('-inf') for node in graph.nodes}
    distance[source] = 0

    # Traverse nodes in topologically sorted order
    for node in topological_order:
        for neighbor in graph.neighbors(node):
            # Update distance of neighbor node
            distance[neighbor] = max(
                distance[neighbor], distance[node] + graph[node][neighbor]['weight'])

    # Find the node with the maximum distance
    max_distance_node = max(distance, key=distance.get)
    max_distance = distance[max_distance_node]

    # Reconstruct longest path
    longest_path = []
    node = max_distance_node
    while node != source:
        longest_path.append(node)
        for neighbor in graph.predecessors(node):
            if distance[neighbor] + graph[neighbor][node]['weight'] == distance[node]:
                node = neighbor
                break
    longest_path.append(source)
    longest_path.reverse()
    longest_path = [str(i) for i in longest_path]

    return max_distance, longest_path


def EditDistance(v, w, match, mismatch, penalty):
    s, a1, a2 = GlobalAlignment_LocalAlignment(v, w, 1, 1, 1)
    a1 = a1.lstrip()
    a2 = a2.lstrip()
    transform = a1.lstrip()
    operations = 0
    for i in range(len(a1)):
        if a1[i] == a2[i]:
            continue
        elif a1[i] == "-":
            transform = transform[:i] + a2[i] + transform[i+1:]
            operations += 1
            continue
        elif a2[i] == "-":
            transform = transform[:i] + '-' + transform[i+1:]
            operations += 1
            continue
        elif a1[i] != a2[i]:
            transform = transform[:i] + a2[i] + transform[i+1:]
            operations += 1
            continue
    transform = transform.replace("-", '')
    o, t = operations, transform
    return o, t


def levenshtein_distance(str1, str2):
    m = len(str1)
    n = len(str2)

    # Create a matrix to store the edit distances
    distance = [[0] * (n + 1) for _ in range(m + 1)]

    # Initialize the first row and column of the matrix
    for i in range(m + 1):
        distance[i][0] = i
    for j in range(n + 1):
        distance[0][j] = j

    # Fill in the rest of the matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if str1[i - 1] == str2[j - 1]:
                distance[i][j] = distance[i - 1][j - 1]
            else:
                distance[i][j] = min(
                    distance[i - 1][j] + 1,  # Deletion
                    distance[i][j - 1] + 1,  # Insertion
                    distance[i - 1][j - 1] + 1  # Substitution
                )

    return distance[m][n]


def FittingAlignment(v, w):
    match = 1
    mismatch = -1
    gap_penalty = -1
    # Matrices
    main_matrix = np.zeros((len(v)+1, len(w)+1), dtype=int)
    match_checker_matrix = np.zeros((len(v), len(w)), dtype=int)

    # Matrix filling
    for i in range(len(v)):
        for j in range(len(w)):
            if v[i] == w[j]:
                match_checker_matrix[i][j] = match
            else:
                match_checker_matrix[i][j] = mismatch

    # for i in range(len(v)):
    #     for j in range(len(w)):
    #         match_checker_matrix[i][j] = BLOSUM62[v[i]][w[j]]
    # print(match_checker_matrix)

    # STEP 1: Initialisation in fitting
    for i in range(len(v)+1):
        main_matrix[i][0] = 0
    for j in range(len(w)+1):
        main_matrix[0][j] = j * gap_penalty

    # STEP 2: Matrix Filling
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            # print(main_matrix)
            main_matrix[i][j] = max(
                main_matrix[i-1][j-1] + match_checker_matrix[i-1][j-1],
                main_matrix[i-1][j] + gap_penalty,
                main_matrix[i][j-1] + gap_penalty)
    # STEP 3: Traceback
    high_score = 0
    high_cell = ()
    for i in range(len(v)):
        if main_matrix[i][len(w)] > high_score:
            high_score = main_matrix[i][len(w)]
            high_cell = (i, len(w))

    aligned_1 = ""
    aligned_2 = ""
    ti = high_cell[0]
    tj = high_cell[1]

    while (ti > 0 and tj > 0):
        if (ti > 0 and tj > 0 and main_matrix[ti][tj] == main_matrix[ti-1][tj-1] + match_checker_matrix[ti-1][tj-1]):
            aligned_1 = v[ti-1] + aligned_1
            aligned_2 = w[tj-1] + aligned_2

            ti = ti - 1
            tj = tj - 1

        elif (ti > 0 and main_matrix[ti][tj] == main_matrix[ti-1][tj] + gap_penalty):
            aligned_1 = v[ti-1] + aligned_1
            aligned_2 = "-" + aligned_2

            ti = ti - 1

        else:
            aligned_1 = "-" + aligned_1
            aligned_2 = w[tj-1] + aligned_2

            tj = tj - 1

    totalScore = 0

    for i in range(0, len(aligned_1)):
        if aligned_1[i] == aligned_2[i]:
            totalScore += match
            # totalScore += BLOSUM62[aligned_1[i]][aligned_2[i]]
        elif aligned_1[i] == "-" or aligned_2[i] == "-":
            totalScore += gap_penalty

        else:
            totalScore += mismatch
            # totalScore += BLOSUM62[aligned_1[i]][aligned_2[i]]

    s, al1, al2 = str(totalScore), "\n" + aligned_1, "\n" + aligned_2
    return s, al1, al2


def AlignmentWithAffineGapPenalties(v, w, match, mismatch, gapOpening, gapExtension):
    # Matrices creation
    Scoring = np.zeros((len(v), len(w)), dtype=float)
    Diag = np.zeros((len(v)+1, len(w)+1), dtype=float)
    Vert = np.zeros((len(v)+1, len(w)+1), dtype=float)
    Horiz = np.zeros((len(v)+1, len(w)+1), dtype=float)

    # Scoring matrix
    for i in range(len(v)):
        for j in range(len(w)):
            if v[i] == w[j]:
                Scoring[i][j] = match
            else:
                Scoring[i][j] = mismatch
    # print(Scoring)
    # Initialize matrices
    Diag[0][0] = 0
    for i in range(1, len(v)+1):
        Diag[i][0] = -gapOpening + (-gapExtension*(i-1))
        Vert[i][0] = -gapOpening + (-gapExtension*(i-1))
        Horiz[i][0] = -float("inf")
    for j in range(1, len(w)+1):
        Diag[0][j] = -gapOpening + (-gapExtension*(j-1))
        Horiz[0][j] = -gapOpening + (-gapExtension*(j-1))
        Vert[0][j] = -float("inf")

    # Matrices filling
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):

            Vert[i][j] = max(Diag[i-1][j] - (gapOpening), Vert[i-1]
                            [j] - gapExtension, Horiz[i-1][j] - gapOpening)
            # print(Vert)
            Horiz[i][j] = max(Diag[i][j-1] - (gapOpening), Horiz[i]
                            [j-1] - gapExtension, Vert[i][j-1] - gapOpening)
            # print(Horiz)
            Diag[i][j] = max(Diag[i-1][j-1] + Scoring[i-1]
                            [j-1], Vert[i][j], Horiz[i][j])
            # print(Diag)
    # print(Scoring)
    # print(Diag)
    # print(Vert)
    # print(Horiz)
    # Traceback (recursive approach to get all the possible alignments)
    alignments = []

    def traceback(i, j, align1, align2):
        if i == 0 and j == 0:
            alignments.append((align1, align2))
            return

        if i > 0 and j > 0 and Diag[i][j] == Diag[i-1][j-1] + Scoring[i-1][j-1]:
            traceback(i-1, j-1, v[i-1] + align1, w[j-1] + align2)

        if i > 0 and Diag[i][j] == Vert[i][j]:
            traceback(i-1, j, v[i-1] + align1, '-' + align2)
            # traceback(i, j-1, '-' + align1, v[i-1] + align2)
        if j > 0 and Diag[i][j] == Horiz[i][j]:
            traceback(i, j-1, '-' + align1, w[j-1] + align2)
            # traceback(i-1, j, w[j-1] + align1, '-' + align2)

    traceback(len(v), len(w), '', '')

    # Score calculation
    scores = []
    for alignment in alignments:
        align1, align2 = alignment
        score = 0
        for i in range(0, len(align1)):
            if align1[i] == align2[i]:
                score += match
            elif align1[i] == "-" or align2[i] == "-":
                if (align1[i-1] == "-" and align1[i] == "-") or (align2[i-1] == "-" and align2[i] == "-"):
                    score -= gapExtension
                else:
                    score -= gapOpening
            else:
                score += mismatch
        scores.append(score)
        HighScore = max(scores)
        newTup = (align1, align2, score)
        alignments[alignments.index(alignment)] = newTup
    FinalAlignments = set()
    for alignment in alignments:
        align1, align2, score = alignment
        if score == HighScore:
            FinalAlignments.add(alignment)
    return FinalAlignments


def findMiddle(v, w, match_reward, mismatch_penalty, gap_penalty):
    mismatch_penalty = -mismatch_penalty
    gap_penalty = -gap_penalty
    STR1 = len(v)+1
    STR2 = len(w)+1
    # Matrices
    main_matrix = np.zeros((len(v)+1, len(w)+1), dtype=int)
    # match_checker_matrix = np.zeros((len(v), len(w)), dtype=int)

    # Matrix filling
    # for i in range(len(v)):
    #     for j in range(len(w)):
    #         if v[i] == w[j]:
    #             match_checker_matrix[i][j] = match_reward
    #         else:
    #             match_checker_matrix[i][j] = mismatch_penalty

    # STEP 1: Initialisation in global
    # for i in range(1,len(v)+1):
    #     main_matrix[i][0] = i * gap_penalty
    # for j in range(len(w)+1):
    #     main_matrix[0][j] = j * gap_penalty

    # STEP 2: Matrix Filling in global
    length_i_First = [main_matrix[0][0]]
    for i in range(1, STR1):
        for j in range(1, STR2//2 + 1):
            # print(main_matrix)
            if v[i-1] == w[j-1]:
                main_matrix[i][j] = max(main_matrix[i-1][j-1] + match_reward,
                                        main_matrix[i-1][j] + gap_penalty,
                                        main_matrix[i][j-1] + gap_penalty)
            else:
                main_matrix[i][j] = max(main_matrix[i-1][j-1] + mismatch_penalty,
                                        main_matrix[i-1][j] + gap_penalty,
                                        main_matrix[i][j-1] + gap_penalty)
            if j == STR2//2:
                length_i_First.append(main_matrix[i][j])

    v = v[::-1]
    w = w[::-1]
    main_matrix = np.zeros((len(v)+1, len(w)+1), dtype=int)

    length_i_Second = [main_matrix[0][0]]
    for i in range(1, STR1):
        for j in range(1, STR2//2 + 1):
            # print(main_matrix)
            if v[i-1] == w[j-1]:
                main_matrix[i][j] = max(main_matrix[i-1][j-1] + match_reward,
                                        main_matrix[i-1][j] + gap_penalty,
                                        main_matrix[i][j-1] + gap_penalty)
            else:
                main_matrix[i][j] = max(main_matrix[i-1][j-1] + mismatch_penalty,
                                        main_matrix[i-1][j] + gap_penalty,
                                        main_matrix[i][j-1] + gap_penalty)
            if j == STR2//2:
                length_i_Second.append(main_matrix[i][j])
    length_i_Second = list(reversed(length_i_Second))
    length_i_Combined = []
    for i in range(len(length_i_Second)):
        length_i_Combined.append(length_i_First[i] + length_i_Second[i])

    MidNode = max(length_i_Combined)
    print(main_matrix)
    return length_i_Combined


w = "GAT"
v = "AA"
# print(graph)
# print(longest_distance)
# print(' '.join(longest_path))
# print(*GlobalAlignment_LocalAlignment(v, w, 1, 1, 2, "global"))
# alignments = GlobalAlignment_LocalAlignment(v, w, 1, 1, 2, "global")
# for alignment in alignments:
#     align1,align2,score= alignment
#     print(align2)
#     print(align1)
#     print(score)
#     break
# print(IterativeOutputLCS(LCSBackTrack(v, w), v, w))
# print(OutputLCS(LCSBackTrack(v, w), v, len(v), len(w)))
# print(LCSBackTrack(v, w))
# print(levenshtein_distance(v, w))


# for a in pairwise2.align.globalms(v,
#                                   w,
#                                   1,-1,-5,-5,one_alignment_only = True):
#     print(format_alignment(*a))

# aligner = Align.PairwiseAligner(mode='global', match_score=1, mismatch_score=-1, gap_score=-5)
# for alignment in aligner.align(v,w):
#     print(alignment.score)
#     print(alignment)

# v = "ATCTTGTATAATAAAATCGAGATTGGGCGTGGACAATACTCCTCGAAGAGACAGCCAATTACCCCTCTTCGCAGAACCAAGCCTATTGGG"
# w = "ATCTTGTATAATTCGAGATTGTGCCTGCCGACAATACGCCTCGAAGGACAAGCCACTACGCTATCCTCTTCCGTGTCAGATCCAAGAGTGCTTATTGGG"
# alignments = AlignmentWithAffineGapPenalties(v,w,1,-1,5,5)
# for alignment in alignments:
#     align1,align2,score= alignment
#     print(align1)
#     print(align2)
#     print(score)
#     break
# print(Diag)
# print(Scoring)
a = SeqRecord(Seq("CCAATACGAC"))
b = SeqRecord(Seq("GCCTTACGCT"))
c = SeqRecord(Seq("CCCTAGCGGC"))
print(MultipleSeqAlignment([a, b, c]))
# print(help(MultipleSeqAlignment))

# print(findMiddle(v,w,1, 1, 2))
