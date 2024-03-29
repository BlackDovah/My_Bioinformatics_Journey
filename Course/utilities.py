def colored(Seq):
    bcolorS = {
        'A': '\033[92m',
        'C': '\033[94m',
        'G': '\033[93m',
        'T': '\033[91m',
        'U': '\033[91m',
        'reSet': '\033[0;0,'
    }

    tmpStr = ""

    for nuc in Seq:
        if nuc in bcolorS:
            tmpStr += bcolorS[nuc] + nuc
        else:
            tmpStr += bcolorS ['reSet'] + nuc
    return tmpStr + '\033[0;0m'

def readTextFile(filepath):
        with open(filepath, 'r') as f:
            return "".join([l.strip() for l in f.readline()])

def writeTextFile(filepath, seq, mode='w'):
    with open(filepath, mode) as f:
        f.write(seq + '\n')

def read_FASTA(filepath):
    with open(filepath, 'r') as f:
        FASTAFile = [l.strip() for l in f.readlines()]

    FASTADict = {}
    FASTALabel = ""

    """ creates a dictionary with the label (lines starting with '>') as keys and the following line as the value"""
    for line in FASTAFile:
        if '>' in line:
            FASTALabel = line
            FASTADict[FASTALabel] = ""
        else:
            FASTADict[FASTALabel] += line
    
    return FASTADict