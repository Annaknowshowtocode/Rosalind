import os
os.chdir("/Users/annaklimova/Downloads")
with open("rosalind_subs.txt", "r") as file:
    file = file.read().strip().split()
    seq = file[0]
    pattern = file[1]
    result = []
    window = len(pattern)
    for i in range(len(seq)):
        if pattern == seq[i:i+window]:
            result.append(str(i+1))
    print(" ".join(result))
