from DNAToolset.bio_structs import *
from collections import *
import random 

class bio_seq:
    """DNA sequence class. Default value: ATCG, DNA, No label"""

    def __init__(self, seq = "ATCG", seq_type = "DNA", label = "No Label"):
        """Sequence intialization, validation."""
        self.seq = seq.upper()
        self.label= label
        self.seq_type = seq_type
        self.is_valid = self.validate()
        assert self.is_valid, f"Provided data does not seem to be a correct {self.seq_type} sequence"

    #DNA Toolkit functions:

    def validate(self):
        """Checks if the DNA sequence is valide."""
        return set(Nucleotides[self.seq_type]).issuperset(self.seq)

    def get_seq_biotype(self):
        """Returns sequence type"""
        return self.seq_type

    def show_seq_info(self):
        """Returns 4 strings. Full sequence information"""
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}\n[Length]: {len(self.seq)}"

    def generate_rnd_seq(self, length = 10, seq_type = "DNA"):
        """Generate a random DNA sequence, proided the length"""
        seq = ''.join([random.choice(Nucleotides[seq_type])
                       for x in range(length)])
        self.__init__(seq, seq_type, "Randomly generated sequence")

    def countNucFrequency(self):
        """Counts the number of occurences of a nucleotide."""
        return dict(Counter(self.seq))
    
    def transcription(self):
        """DNA -> RNA Transcription. Replaces Thymine with Uracil."""
        if self.seq_type =="DNA":
            return self.seq.replace("T", "U")
        return"Not a DNA sequence"
    
    def complement(self):
        """Swapps Adenine with Thymine and Cytosine with Guanine and vice versa."""
        if self.seq_type == "DNA":
            mapping = str.maketrans("ATCG", "TAGC")
        else:
            mapping = str.maketrans("AUCG", "UAGC")
        return self.seq.translate(mapping)

    def reverse_complement(self):
        """Same as "complement" but also reverses the sequence."""
        if self.seq_type == "DNA":
            mapping = str.maketrans("ATCG", "TAGC")
        else:
            mapping = str.maketrans("AUCG", "UAGC")
        return self.seq.translate(mapping)[::-1]

    def gc_content(self):
        """GC content in a DNA/RNA sequence."""
        return str(round((self.seq.count('C') + self.seq.count('G')) / len(self.seq) * 100 )) + "%"

    def gc_content_subsec(self, k=20):

        """GC Content in a DNA/RNA sub-sequence len k, k=20 by default."""
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i:i + k]
            res.append(round((subseq.count('C') + subseq.count('G')) / len(subseq) * 100))
        return res 

    def translate_seq(self, init_pos=0):
        """Translates DNA codons to their aminoacids"""
        if self.seq_type == "DNA":
            return [DNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]
        elif self.seq_type == "RNA":
            return [RNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]

    def codon_usage(self, aminoacid):
        """Provides the frequencey of each codon encoding a given aminoacid in a DNA sequence"""
        tmpList = []
        if self.seq_type == "DNA":
            for i in range(0, len(self.seq) - 2, 3):
                if DNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])
        
        elif self.seq_type == "RNA":
            for i in range(0, len(self.seq) - 2, 3):
                if RNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])

        freqDict = dict(Counter(tmpList))
        totalWeight = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(freqDict[seq] / totalWeight, 2)
            freqDict[seq] = freqDict[seq] * 100
            freqDict[seq] = str(int(freqDict[seq])) + "%"
        return freqDict

    def gen_reading_frames(self):
        """Generate the six reading frames of a DNA sequence, including the reverse complement"""
        frames = []
        frames.append(self.translate_seq(0))
        frames.append(self.translate_seq(1))
        frames.append(self.translate_seq(2))
        tmp_seq = bio_seq(self.reverse_complement(), self.seq_type)
        frames.append(tmp_seq.translate_seq(0))
        frames.append(tmp_seq.translate_seq(1))
        frames.append(tmp_seq.translate_seq(2))
        del tmp_seq
        return frames

    def proteins_from_rf(self, aa_seq):
        """Compute all possible proteins in an aminoacid seq and return a list of possible proteins"""
        current_prot = []
        proteins = []
        for aa in aa_seq:
            if aa == "_":
                # Stop accumulating amino acids if _ . STOP was found
                if current_prot:
                    for p in current_prot:
                        proteins.append(p)
                    current_prot = []
            else:
                # START accumulating amino acids if M . START was found
                if aa == "M":
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += aa
        return proteins

    def all_proteins_from_orfs(self, startReadPos=0, endReadPos=0, ordered=False):
        """Compute all possible proteins for all open reading frames"""
        """Protine Search DB: https://www.ncbi.nlm.nih.gov/nuccore/NM_001185097.2"""
        """API can be used to pull protein info"""
        if endReadPos > startReadPos:
            tmp_seq = bio_seq(self.seq[startReadPos: endReadPos], self.seq_type)
            rfs = tmp_seq.gen_reading_frames()
            del tmp_seq
        else:
            rfs = self.gen_reading_frames()

        res = []
        for rf in rfs:
            prots = self.proteins_from_rf(rf)
            for p in prots:
                res.append(p)
                
        if ordered:
            return sorted(res, key=len, reverse=True)
        return res
