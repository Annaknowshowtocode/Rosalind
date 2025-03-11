import os
import re
from collections import Counter
os.chdir('/Users/annaklimova/Downloads')

# Словарь для хранения данных: ключ — заголовок, значение — последовательность
fasta_data = {}
labels = [] # Список для хренения названий
GC_content = 0
GC_content_list = []

# Чтение файла
with open("rosalind_gc-2.txt", "r") as file:
    current_label = None  # Текущий заголовок
    for line in file:
        line = line.strip()  # Убираем лишние пробелы и символы новой строки
        if line.startswith(">"):  # Если строка начинается с ">", это заголовок
            current_label = line  # Устанавливаем текущий заголовок
            labels.append(line)
            fasta_data[current_label] = ""  # Создаём пустую последовательность для этого заголовка
        else:
            # Добавляем строку последовательности к текущему заголовку
            fasta_data[current_label] += line

    for label in labels:
        G, C = 0, 0
        for seq in  fasta_data[label]:
            G += 1 if seq == "G" else 0
            C += 1 if seq == "C" else 0
        current_GC = (((G + C) * 100) / len(fasta_data[label]))
        GC_content = max(current_GC, GC_content)
        GC_content_list.append(current_GC)
    max_index = GC_content_list.index(GC_content)
print(f'{labels[max_index][1::]} \n {float(GC_content_list[max_index]):.6f}')
print(GC_content_list)



