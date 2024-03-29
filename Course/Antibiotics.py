# from bio_seq import *
# from bio_structs import *
from DNAToolset.bio_seq import *
import itertools
import copy

string = open(
    r"C:\Work and Education\Bioinformatics\Genomes\Bacillus_brevis.txt").read().strip()
peptide = "VKLFPWFNQY"
# string = "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
# peptide = "MA"
spectrum = [0, 113, 128, 186, 241, 299, 314, 427]
spectrum1 = [int(i) for i in open(
    r"C:\Users\spidy\Downloads\dataset_104_4.txt").read().strip().split()]
test = [0, 71, 101, 113, 131, 184, 202, 214, 232, 285, 303, 315, 345, 416]


def protein_encoding(DNA, peptide):
    rev = ""
    last = ""
    transLis = []
    for i in range(len(DNA)-(len(peptide)*3)):
        codon = DNA[i:i+len(peptide)*3]
        test_dna = bio_seq(seq=codon, seq_type="DNA")
        rev = test_dna.reverse_complement()
        if ''.join(test_dna.translate_seq()) == peptide or ''.join(bio_seq(seq=rev).translate_seq()) == peptide:
            transLis.append(codon)
    last = ""
    for i in transLis:
        last = last + "\n" + i
    return transLis


def sub_peptide_generator_CYCLIC(Peptide):
    subPeptides = []
    extendedDNA = Peptide + Peptide[:len(Peptide)-2]

    for ind in range(len(Peptide)):
        for inc in range(1, len(Peptide)):
            subPeptides.append(extendedDNA[ind:ind+inc])
    subPeptides.sort(key=len)
    return subPeptides


def sub_peptide_generator_LINEAR(Peptide):
    subPeptides = []

    n = len(Peptide)
    for ind in range(len(Peptide)):
        if ind > 1:
            n -= 1
        for inc in range(1, n):
            subPeptides.append(Peptide[ind:ind+inc])
    subPeptides.sort(key=len)
    return subPeptides


def theoretical_spectrum_generator(Peptide, Type="cyclic"):
    c = -1
    theoreticalSpectrum = []
    tot = 0
    DNATot = 0
    if Type == "cyclic":
        subPeptides = sub_peptide_generator_CYCLIC(Peptide)
    else:
        subPeptides = sub_peptide_generator_LINEAR(Peptide)
    for sub in subPeptides:
        c += 1
        if tot != 0:
            theoreticalSpectrum.append(tot)
            tot = 0
        for AA in sub:
            tot += AA_Mass[AA]
        if c == len(subPeptides)-1:
            theoreticalSpectrum.append(tot)
    theoreticalSpectrum.append(0)
    for AA in Peptide:
        DNATot += AA_Mass[AA]
    theoreticalSpectrum.append(DNATot)
    theoreticalSpectrum.sort()
    # two different ways to change the data type of the items in the list
    theoreticalSpectrum = map(str, theoreticalSpectrum)
    # theoreticalSpectrum = [str(i) for i in theoreticalSpectrum]
    return list(theoreticalSpectrum)


def cyclo_peptide_sequencing(Spectrum, Type="cyclic"):
    compare = copy.deepcopy(Spectrum)
    Spectrum = [str(i) for i in Spectrum]
    finalPeptides = []
    candidatePeptides = []
    tempCand = []
    tempSpec = []

    for mass in Spectrum:
        for AA, value in AA_Mass18.items():
            if mass == str(value) and AA not in candidatePeptides:
                candidatePeptides.append(AA)
    compare2 = copy.deepcopy(
        theoretical_spectrum_generator(candidatePeptides[0]))
    compare2 = [int(i) for i in compare2]
    while max(compare2) != max(compare) and not max(compare2) > max(compare):
        for amino in candidatePeptides:
            for AA in AA_Mass18.keys():
                temp = theoretical_spectrum_generator(amino + AA, Type)
                o = len(amino + AA)
                p = 0
                for i in range(len(amino + AA)-1):
                    p += o
                    o -= 1
                if len(temp)-2 != p:
                    continue
                tempSpec = copy.deepcopy(Spectrum)
                for i in theoretical_spectrum_generator(amino + AA, Type):
                    if i in tempSpec:
                        temp.remove(i)
                        tempSpec.remove(i)
                        if temp == []:
                            tempCand.append(amino + AA)
                        else:
                            continue
        candidatePeptides = tempCand
        tempCand = []
        compare2 = copy.deepcopy(
            theoretical_spectrum_generator(candidatePeptides[0]))
        compare2 = [int(i) for i in compare2]

    # temp = copy.deepcopy(candidatePeptides)
    # for i in temp:
    #     pepMass = 0
    #     for j in i:
    #         pepMass += AA_Mass18[j]
    #     if str(pepMass) != ''.join(Spectrum[-1:]):
    #         candidatePeptides.remove(i)

    # finalPeptides = [list() for i in candidatePeptides]
    # for pos, peptide in enumerate(candidatePeptides):
    #     for AA in peptide:
    #         finalPeptides[pos].append(AA_Mass18[AA])
    # finalPeptides.sort(reverse=True)
    # for pos, i in enumerate(finalPeptides):
    #     finalPeptides[pos] = [str(j) for j in finalPeptides[pos]]
    #     finalPeptides[pos] = '-'.join(finalPeptides[pos])
    # return ' '.join(finalPeptides)
    return candidatePeptides


def cyclopeptide_scoring(Spectrum, peptide, type="pep"):
    score = 0
    Spectrum = [str(i) for i in Spectrum]
    tempSpec = copy.deepcopy(Spectrum)
    theoPep = ""
    scores = {}
    if type == "pep":
        #     if type(peptide) == list:
        #         for i in peptide:
        #             x = theoretical_spectrum_generator(i, "linear")
        #             for j in x:
        #                 if j in tempSpec:
        #                     score += 1
        #             scores[i] = score
        #             score = 0
        #         return scores

        # else:
        theoPep = theoretical_spectrum_generator(peptide, "linear")
        for i in theoPep:
            if i in tempSpec:
                score += 1
                tempSpec.remove(i)
    else:
        for i in peptide:
            if i in tempSpec:
                score += 1
                tempSpec.remove(i)
    return score


def leaderboard_cyclopeptide_sequencing(Spectrum, N, Convolution, Type="cyclic"):
    compare = copy.deepcopy(Spectrum)
    Spectrum = [str(i) for i in Spectrum]
    candidatePeptides = []
    tempCand = []
    tempCand2 = []
    tempSpec = []
    leaderBoard = []
    buildingBlocks = []
    scores = []
    high = 0
    highest = high
    # for mass in Spectrum:
    #     for AA, value in AA_Mass18.items():
    #         if mass == str(value) and AA not in buildingBlocks:
    #             buildingBlocks.append(AA)
    # candidatePeptides = buildingBlocks
    # compare2 = copy.deepcopy(theoretical_spectrum_generator(buildingBlocks[0]))
    # compare2 = [int(i) for i in compare2]
    buildingBlocks = Convolution
    while highest <= high and len(compare2) <= len(compare):
        if max(compare2) == max(compare):
            Type = "cyclic"
        else:
            Type = "linear"
        for amino in candidatePeptides:
            for AA in buildingBlocks:
                temp = theoretical_spectrum_generator(amino + AA, Type)
                temp2 = copy.deepcopy(temp)
                tempSpec = copy.deepcopy(Spectrum)
                for i in temp2:
                    if i in tempSpec:
                        temp.remove(i)
                        tempSpec.remove(i)
                tempCand.append(amino + AA)
                scoring = cyclopeptide_scoring(Spectrum, temp2, "spec")
                tempCand2.append((scoring, amino + AA))
                scores.append(scoring)

        tempCand2.sort(reverse=True)
        # test = copy.deepcopy(tempCand2)
        candidatePeptides = []
        highest = copy.deepcopy(high)
        high = copy.deepcopy(max(scores))
        scores = []
        while len(tempCand2) > 0:
            i = tempCand2[0]
            if len(candidatePeptides) == N:
                lastPep = copy.deepcopy(i[0])
                for j in tempCand2:
                    if j[0] == lastPep:
                        candidatePeptides.append(j[1])
            elif len(candidatePeptides) < N:
                candidatePeptides.append(i[1])
                tempCand2.remove(i)
                if tempCand2 == []:
                    break
            else:
                break

        tempCand = []
        tempCand2 = []
        compare2 = copy.deepcopy(
            theoretical_spectrum_generator(candidatePeptides[0]))
        compare2 = [int(i) for i in compare2]
    # return test
    hold = []
    for i in candidatePeptides:
        leaderBoard.append(i)
    for i, j in enumerate(leaderBoard):
        for l in j:
            hold.append(AA_Mass18[l])
        leaderBoard[i] = hold
        hold = []
    for i, l in enumerate(leaderBoard):
        leaderBoard[i] = '-'.join([str(j) for j in l])
    return '    '.join(leaderBoard)


def spectral_convolution(Spectrum, M):
    freqDic = {}
    Spectrum.sort(reverse=True)
    spec = copy.deepcopy(Spectrum)
    for i in spec:
        Spectrum.remove(i)
        for j in Spectrum:
            l = i - j
            if l not in freqDic:
                if l in range(57, 201):
                    freqDic[l] = 1
            else:
                freqDic[l] += 1
    # for i, j in freqDic.items():
    #     freqTup =
    answer = ""
    n = 0
    for i, j in freqDic.items():
        n = 0
        while n < j:
            n += 1
            answer = answer + " " + str(i)
    return freqDic


def convolution_cyclopeptide_sequencing(Spectrum, N, M):
    convolution = spectral_convolution(Spectrum, M)
    sequencing = leaderboard_cyclopeptide_sequencing(Spectrum, N, convolution)
    return sequencing


# print(len(protein_encoding(string.replace("\n", ""), peptide)))
print(theoretical_spectrum_generator(
    "LHRWEFPCVNLKTVRTWVIPRVIFHYSNFIFHRWSQWGDSANE", "linear"))
# print(cyclo_peptide_sequencing(spectrum1, "linear"))
# print(sub_peptide_generator_LINEAR("NQELVCTP"))
# print(leaderboard_cyclopeptide_sequencing([0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460], 10, "linear"))
# print(leaderboard_cyclopeptide_sequencing(spectrum1, 1000, "linear"))
# print(cyclopeptide_scoring(spectrum1, "PEPFVAWAIYDAIKCSKTHYN"))
# print(spectral_convolution([0, 137, 186, 323]))
# print(spectral_convolution([0,86,160,234,308,320,382]))
# print(cyclopeptide_scoring([0,57,71,71,71,104,131,202,202,202,256,333,333,403,404], "MAMA"))
# print(cyclopeptide_scoring([0,97,129,129,194,226,323,323,355,452], "PEEP"))
