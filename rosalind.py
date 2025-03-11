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

first_key = next(iter(fasta))
first_value = fasta[first_key]
except_first = list(fasta.values())[1:]
result = ""
previous_index = -1  # Начинаем с -1, чтобы первый индекс всегда был больше
DNA_str = first_value

# Ищем индексы интронов в DNA_str
for subseq in except_first:
    current_index = DNA_str.find(subseq, previous_index + 1)  # Ищем интрон после предыдущего индекса
    if current_index == -1:
        continue  # Если интрон не найден, пропускаем
    if current_index > previous_index:
        previous_index = current_index
        result += str(current_index) + " "  # Добавляем индекс с пробелом

# Выводим результат
print(result.strip())  # Убираем лишний пробел в конце