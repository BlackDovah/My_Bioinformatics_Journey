import Biopython as bio
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq 

# Takes a random DNA string, and finds species that have highly similar versions of it

seqInFile = open(r"C:\\Work and Education\Bioinformatics\Genomes\random.txt").read()
seqInFile = "TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTAC\
AATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCAC\
CTACGGTAGAG"
result = NCBIWWW.qblast("blastn", "nt", seqInFile)
formatResult = NCBIXML.read(result)
E_Value = 0.01
for alignment in formatResult.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_Value:
            print("***Alignment***")
            print("Sequence:", alignment.title)
            print("Length:", alignment.length)
            print("E-Value:", hsp.expect)
            print(hsp.query)
            print(hsp.match)
            print(hsp.sbjct)
