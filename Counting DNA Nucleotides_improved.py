import os
os.chdir('/Users/annaklimova/Downloads')

from collections import Counter

with open("rosalind_dna.txt", "r") as file:
    file = file.read().strip()

    nucleotides_set = ["A", "C", "G", "T"]
    counts = Counter(file)

    result = [str(counts[nuc]) for nuc in nucleotides_set]

    print(" ".join(result))