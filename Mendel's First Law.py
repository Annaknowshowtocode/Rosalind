import os
from itertools import count

os.chdir("/Users/annaklimova/Downloads")

# Читаем файл и сохраняем последовательности в список
with open("rosalind_iprb.txt", "r") as file:
    file = file.read().split()
    k = int(file[0])
    m = int(file[1])
    n = int(file[2])
    T = k + m + n
    AA_AA = k/T * (k - 1)/(T - 1)
    AA_Aa = (k/T * m/(T-1)) + (m/T * k/(T-1))
    Aa_Aa = m/(T) * (m-1)/(T-1)
    Aa_aa = (m/T * n/(T-1)) + (n/T * m/(T-1))
    AA_aa = (k/T * n/(T-1)) + (n/T * k/(T-1))
    P = AA_AA * 1 + AA_Aa * 1 + Aa_Aa * 0.75 + Aa_aa * 0.5 + AA_aa * 1
    print(P)