import os
os.chdir('/Users/annaklimova/Downloads')

with open("rosalind_fib.txt", "r") as file:
    file = file.read().strip().split()
    n = int(file[0])
    k = int(file[1])
    rabbits = [1] * n
    for i in range(2, n):
        rabbits[i] = rabbits[i - 1] + k * rabbits[i - 2]
    print(rabbits[n-1])
