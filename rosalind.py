import os

os.chdir("/Users/annaklimova/Downloads/")

with open("rosalind_splc.txt", "r") as file:
    fasta = {}
    current_label = ""
    for line in file:
        line = line.strip()
        if line.startswith(">"):
            current_label = line
            fasta[current_label] = ""
        else:
            fasta[current_label] += line

