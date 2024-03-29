class genomeToolkit:
    def __init__(self):
        print("[Genome Toolkit initiated]\n")

    def count_kmer(self, sequence, kmer):
        """
        Counts the number of times a spcific k-mers appears in a sequence.
         Includes overlapping k-mers
         
         Parameters:
            sequence (str): The DNA sequence to search in.
            kmer(str): The specific k-mer to search for in the sequence.
            
        Returns:
            int: The number of times the k-mer appears in the sequence.
        """
        kmer_count = 0
        for position in range(len(sequence) - (len(kmer) - 1)):
            if sequence[position:position + len(kmer)] == kmer:
                kmer_count += 1
        return kmer_count

    def find_most_frequent_kmers(self, sequence, k_len):
        """
        Finds the most frequent k-mers of a given length in a DNA string.

        Parameters:
            sequence (str): The DNA string to search.
            k_len (int): The length of the k-mers to search for.
        
        Returns:
            list: A list of the most frequent k-mers in the DNA string.
        """
        #1. A dictionary to store k-mer frequencies.
        kmer_frequencies = {}
        #2. A loop to iterate through the DNA string and extract k-mers of a given length,
        #while also incrementing the frequency of each k-mer.
        for i in range(len(sequence) - k_len + 1):
            kmer = sequence[i:i+k_len]
            if kmer in kmer_frequencies:
                kmer_frequencies[kmer] += 1
            else:
                kmer_frequencies[kmer] = 1
        #3. A variable to store the maximum value in the dictionary.
        highest_frequency = max(kmer_frequencies.values())
        #4. List comprehension/loop to scan the k-mer dictionary and return a list of k-mers
        #with the highest frequency.
        return [
            kmer for kmer, frequency in kmer_frequencies.items()
            if frequency == highest_frequency
        ]
