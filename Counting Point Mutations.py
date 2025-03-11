import os
import re
from collections import Counter
os.chdir('/Users/annaklimova/Downloads')

# Чтение файла
with open("rosalind_hamm.txt", "r") as file: # Текущий заголовок
    file = file.read().strip().split()
    line_template = file[0]
    line_check = file[1]
    matching_letters = 0
    for i in range(len(line_template)):
        if line_template[i] != line_check[i]:
            matching_letters += 1
    print(matching_letters)


