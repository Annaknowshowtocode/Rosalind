import os

os.chdir("/Users/annaklimova/Downloads")

with open("rosalind_prtm.txt", "r") as file:
    file = file.read().strip()
    protein = list(file)
    mass = 0
    protein = {
        "A": ["GCU", "GCC", "GCA", "GCG"],
        "C": ["UGU", "UGC"],
        "D": ["GAU", "GAC"],
        "E": ["GAA", "GAG"],
        "F": ["UUU", "UUC"],
        "G": ["GGU", "GGC", "GGA", "GGG"],
        "H": ["CAU", "CAC"],
        "I": ["AUU", "AUC", "AUA"],
        "K": ["AAA", "AAG"],
        "L": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
        "M": ["AUG"],
        "N": ["AAU", "AAC"],
        "P": ["CCU", "CCC", "CCA", "CCG"],
        "Q": ["CAA", "CAG"],
        "R": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "S": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
        "T": ["ACU", "ACC", "ACA", "ACG"],
        "V": ["GUU", "GUC", "GUA", "GUG"],
        "W": ["UGG"],
        "Y": ["UAU", "UAC"]
    }


    for amino in protein:
        mass += protein_mass[amino]

    print(f'{mass:.3f}')