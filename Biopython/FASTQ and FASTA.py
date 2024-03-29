from Bio import SeqIO
from io import StringIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator

#Reads counter
records = SeqIO.parse(r"good_quality.fastq", "fastq")
count = SeqIO.write(records, "THIS_IS_YOUR_OUTPUT_FILE.fasta", "fasta")
print("Converted %i records" % count)

#Number of reads below an average quality
phredQualLis = []
i = 0
for rec in SeqIO.parse(r"good_quality.fastq", "fastq"):
    phredQualLis.append(rec.letter_annotations['phred_quality'])
    
average = []
for i in phredQualLis:
    average.append(round(sum(i)/len(i)))

readBelAvr = 0
for i in average:
    if i < 24:
        readBelAvr += 1

print(readBelAvr)
  
