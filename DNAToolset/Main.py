from bio_seq import bio_seq
from utilities import read_FASTA, readTextFile, writeTextFile

# fasta = read_FASTA(r"C:\Work and Education\College\Bioinformatics\Clostridium thermoaceticum ATCC 39073 Spo0A (spo0A) gene, partial cds.txt")
# print(f'{fasta} \n')
# test_dna = bio_seq(seq = ''.join(fasta.values()), seq_type="DNA", label= ''.join(fasta.keys()))
# test_dna = bio_seq(seq = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA", seq_type='RNA')
# test_dna.generate_rnd_seq(40, "DNA")
# test_dna = bio_seq(seq = open(r"C:\Users\spidy\Downloads\dataset_96_4.txt").read().strip(), seq_type= "RNA")

# 

# print(f'{test_dna.show_seq_info()}\n')
# print(f'[1] Nucleotides Frequency: {test_dna.countNucFrequency()}\n')
# print(f'[2] Transcribed Sequence: {test_dna.transcription()}\n')
# print(f'[3] Complement: {test_dna.complement()}\n')
# print(f'[4] Reverse Complement: {test_dna.reverse_complement()}\n')
# print(f'[5] GC Content: {test_dna.gc_content()}\n')
# print(f'[6] GC Content of subsection(size = 20 by default): {test_dna.gc_content_subsec()}\n')
# print(f'[7] Translated Sequence: {test_dna.translate_seq()}\n')
# print(f'[8] Frequency of Specified Aminoacid: {test_dna.codon_usage("L")}\n')

# print("[9] Open Reading Frames:")
# for rf in test_dna.gen_reading_frames():
#     print(rf)

# print(test_dna.proteins_from_rf("ACG"))

# print(f'\n[10] All Proteins Translated From ORFs: {test_dna.all_proteins_from_orfs()}')

# writeTextFile("Results.html",  f'[1] Nucleotides Frequency:\n {test_dna.seq}\n'
#               + f'---Sequence info: \n {test_dna.show_seq_info()}\n' 
#               + f'[2] Transcribed Sequence: \n {test_dna.transcription()}\n'
#               + f'[3] Complement: \n {test_dna.complement()}\n'
#               + f'[4] Reverse Complement: \n {test_dna.reverse_complement()}\n'
#               + f'[5] GC Content: \n {test_dna.gc_content()}\n'
#               + f'[6] GC Content of subsection(size = 20 by default): \n {test_dna.gc_content_subsec()}\n'
#               + f'[7] Translated Sequence: \n {test_dna.translate_seq()}\n'
#               + f'[8] Frequency of Specified Aminoacid: \n {test_dna.codon_usage("L")}\n'
#               + f'[9] Open Reading Frames: \n {rf}\n'
#               + f'[10] All Proteins Translated From ORFs: \n {test_dna.all_proteins_from_orfs()}')

