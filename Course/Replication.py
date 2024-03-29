from utilities import *

# f = read_FASTA(r"C:\Users\spidy\OneDrive\Desktop\Salmonella_enterica.txt")


def PatternCount(Genome, Pattern):
    """Counts all the occurences of a specific pattern (K-mer) in the genome.
    subroutines = none
    Subroutine in = 1. SymbolArray
                    2. FasterSymbolArray"""
    count = 0
    for i in range(len(Genome)-len(Pattern)+1):
        if Genome[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count


def FrequencyMap(Genome, k):
    """Finds all the K-mers of a specific length = K in a string and returns a dictionary showing their frequencies.
    subroutines = none
    Subroutine in = FrequentWords"""
    freq = {}
    n = len(Genome)
    for i in range(n-k+1):
        Pattern = Genome[i:i+k]
        if Pattern in freq:
            freq[Pattern] += 1
        else:
            freq[Pattern] = 0
            freq[Pattern] += 1
    return freq


def FrequentWords(Genome, k):
    """Finds all the K-mers in a string and returns a list of the most frequent ones.
    Subroutines = 1. FrequencyMap
    Subroutine in = none"""
    words = []
    freq = FrequencyMap(Genome, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words


def ReverseComplement(Pattern):
    """Returns the reverse complement of a string (5' -> 3').
    subroutines = 1. Reverse
                  2. Complement
    Subroutine in = none"""
    reverse = Reverse(Pattern)
    complement = Complement(reverse)
    return complement


def Reverse(Pattern):
    """Reverses the complement of a string.
    subroutines = none
    Subroutine in = 1. ReverComplement"""
    return Pattern[::-1]


def Complement(Pattern):
    """Returns the complement of a string (3' -> 5').
    subroutines = none
    Subroutine in = 1. ReverseComplement"""
    complement = ""
    for i in Pattern:
        if i == "A":
            complement += "T"
        elif i == "T":
            complement += "A"
        elif i == "C":
            complement += "G"
        elif i == "G":
            complement += "C"
    return complement


def PatternMatching(Pattern, Genome):
    """Shows the starting positions of a pattern in the genome.
    sbroutines = none
    Subroutine in = none"""
    positions = []
    for i in range(len(Genome)-len(Pattern)+1):
        if Genome[i:i+len(Pattern)] == Pattern:
            positions.append(Genome.find(Pattern, i, i+len(Pattern)))
    return positions


def SymbolArray(Genome, symbol):
    """Scans half the genom at a time to calculate the occurences of a specific nucleotide in it moving one nucleotide over every time.
    subroutines = 1. PatternCount
    Subroutine in = none"""
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(ExtendedGenome[i:i+(n//2)], symbol)
    return array


def FasterSymbolArray(Genome, symbol):
    """subroutines = 1. PatternCount
    Subroutine in = none"""
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(Genome[0:n//2], symbol)

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array


def SkewArray(Genome):
    """Calculates the difference in occurences between "C" and "G" by scanning 
    the genome one nucleotide at a time and adding +1 if a "G" is encountered or -1 if "C" is encountered.
    subroutines = none
    Subroutine in = 1. MinimumSkew"""
    Skew = [0]
    SkewDict = {}
    n = len(Genome)
    for i in range(n):
        if Genome[i:i+1] in "AT":
            Skew.append(Skew[i])
        elif Genome[i:i+1] == "G":
            Skew.append(Skew[i]+1)
        elif Genome[i:i+1] == "C":
            Skew.append(Skew[i]-1)
    for i in range(len(Skew)):
        SkewDict[i] = Skew[i]
    return SkewDict


def MinimumSkew(Genome):
    """Calls the SkewArray function and returns the starting positions in the scan that had the minimum value for G.
    subroutines = 1. SkewArray
    Subroutine in = none"""
    positions = []  # output variable
    Skew = SkewArray(Genome)
    minimum = min(Skew.values())
    for i in Skew:
        if Skew[i] == minimum:
            positions.append(i)
    return positions


def HammingDistance(p, q):
    """Compares two strings and returns the number of differences between them.
    subroutines = none
    Subroutine in = 1. ApproximatePatternMatching
                    2. ApproximatePatternCount"""
    HamDis = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            HamDis += 1
    return HamDis


def ApproximatePatternMatching(Text, Pattern, d):
    """Same as PatternMatching but accepting a "d" amount of mismatches = to the Hamming Distance.
    subroutines = 1. HammingDistance
    Subroutine in = none"""
    positions = []  # initializing list of positions
    for i in range(len(Text)-len(Pattern)+1):
        p = Text[i:i+len(Pattern)]
        q = Pattern
        if HammingDistance(p, q) <= d:
            positions.append(i)
    return positions


def ApproximatePatternCount(Text, Pattern, d):
    """Same as PatternCount but accepting a "d" amount of mismatches = to the Hamming Distance.
    subroutines = 1. HammingDistance
    Subroutine in = none"""
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        p = Text[i:i+len(Pattern)]
        q = Pattern
        if HammingDistance(p, q) <= d:
            count = count+1
    return count


def FindClumps(Text, k, L, t):
    """Finds Kmers of length K within a window of length L if
    they exist in a number equal to or higher than T and are considered as forming a clump.
    Subroutines: none
    Subroutine in: none"""
    freq = {}
    clumps = 0
    cCheck = 0
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        if Pattern in freq:
            freq[Pattern].append(i)
        else:
            freq[Pattern] = []
            freq[Pattern].append(i)
    for j in freq.values():
        for dis in range(2, len(j)):
            cCheck = clumps
            if j[dis] - j[dis - (t-1)] <= L-(k+1):
                clumps += 1
                if clumps > cCheck:
                    break

    return clumps


def ApproximateFrequencyMap(Genome, k, d):
    """Same as FrequencyMap but accepting a "d" amount of mismatches = to the Hamming Distance.
    subroutines = 1. HammingDistance
    Subroutine in = ApproximateFrequentWords"""
    freq = {}
    n = len(Genome)
    for i in range(n-k+1):
        Pattern = Genome[i:i+k]
        for j in freq:
            if HammingDistance(Pattern, j) <= d:
                freq[j] += 1
        else:
            freq[Pattern] = 1

    return freq


def ApproximateFrequentWords(Genome, k, d):
    """Same as FrequentWords but accepting a "d" amount of mismatches = to the Hamming Distance.
    subroutines = 1. FrequencyMap
                  2. HammingDistance
    Subroutine in = none"""
    words = []
    freq = ApproximateFrequencyMap(Genome, k, d)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words


def Neighbors(Pattern, d):
    """Takes a kmer and find all the possible different versions in the structure up to "d" number of differences in nucleotides.
    Subroutines: none
    Subroutine in: 1. FrequentWordsWithMismatches
                   2. FrequentWordsWithMismatchesAndReverseComplements"""
    symbol = ["A", "C", "G", "T"]
    neighbors = [Pattern]
    while d > 0:
        for i in neighbors:
            if d == 0:
                break
            d -= 1
            for j in range(len(i)):
                for s in symbol:
                    neigh = i[:j] + s + i[j+1:]
                    if neigh not in neighbors:
                        neighbors.append(neigh)
    return neighbors


def FrequentWordsWithMismatches(Genome, k, d):
    """Finds all the most frequent kmers, accounting for mismatches with a "d" Hamming Distance.
    Subroutines: 1. Neigbors
                 2. ApproximatePatternCount
    Subroutine in: none"""
    freq = {}
    mKmers = []
    Patterns = []
    n = len(Genome)
    for i in range(n-k+1):
        Pattern = Genome[i:i+k]
        neighbor = Neighbors(Pattern, d)
        Patterns.append(neighbor)
    for j in range(len(Patterns)):
        for s in Patterns:
            for y in s:
                if y not in freq:
                    freq[y] = 0
    for i in freq:
        freq[i] = ApproximatePatternCount(Genome, i, d)
    maxFreq = max(freq.values())
    for kmers, values in freq.items():
        if values == maxFreq:
            mKmers.append(kmers)
    return ' '.join(str(i) for i in mKmers)


def FrequentWordsWithMismatchesAndReverseComplements(Genome, k, d):
    """Finds all the most frequent kmers, accounting for mismatches with a 
    "d" Hamming Distance, as well as the revese complements of the kmers.
    Subroutines: 1. Neigbors
                 2. ApproximatePatternCount
    Subroutine in: none"""
    freq = {}
    mKmers = []
    Patterns = []
    n = len(Genome)
    for i in range(n-k+1):
        Pattern = Genome[i:i+k]
        neighbor = Neighbors(Pattern, d)
        Patterns.append(neighbor)
    for j in range(len(Patterns)):
        for s in Patterns:
            for y in s:
                if y not in freq:
                    freq[y] = 0
                    freq[ReverseComplement(y)] = 0
    for i in freq:
        freq[i] = ApproximatePatternCount(Genome, i, d)
        freq[i] += ApproximatePatternCount(ReverseComplement(Genome), i, d)
    maxFreq = max(freq.values())
    for kmers, values in freq.items():
        if values == maxFreq:
            mKmers.append(kmers)
    return ' '.join(str(i) for i in mKmers)


# print(PatternCount("", "TCAGGTTTC"))
# print(FrequencyMap("TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT", 3))
# print(FrequentWords("", 11))
# print(ReverseComplement(""))
# pattern = PatternMatching("CTTGATCAT", "")
# print(' '.join(str(i) for i in pattern))
# print(FasterSymbolArray(f'Course//Ori region of Thermotoga petrophila', "C"))
import matplotlib.pyplot as plt
plt.plot(range(len(SkewArray("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"))), SkewArray("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT").values())
plt.show()
# print(SkewArray("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"))
# print(MinimumSkew("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"))
# for i in f.values():
#     gen = i
# print(FrequentWordsWithMismatchesAndReverseComplements(gen[3764856:3765356], 9, 1))
