import os

os.chdir("/Users/annaklimova/Downloads")
with open("rosalind_lcsm.txt", "r") as file:
    current_label = None  # Текущий заголовок
    labels = []
    fasta_data = {}
    for line in file:
        line = line.strip()  # Убираем лишние пробелы и символы новой строки
        if line.startswith(">"):  # Если строка начинается с ">", это заголовок
            current_label = line  # Устанавливаем текущий заголовок
            labels.append(line)
            fasta_data[current_label] = ""  # Создаём пустую последовательность для этого заголовка
        else:
            # Добавляем строку последовательности к текущему заголовку
            fasta_data[current_label] += line

# Первая последовательность
first_key = next(iter(fasta_data))
first_value = fasta_data[first_key]
# Остальные последовательности
except_first = list(fasta_data.values())[1:]

longest_common_substring = ""

for i in range(len(first_value), 0, -1):  # Длина окна от max до 1
    for start in range(len(first_value) - i + 1):  # Начальная позиция окна
        substring = first_value[start:start + i]

        # Проверяем, содержится ли эта подстрока во всех других последовательностях
        if all(substring in seq for seq in except_first):
            longest_common_substring = substring
            break  # Прерываем поиск, так как нашли самый длинный

    if longest_common_substring:
        break  # Прерываем внешний цикл, если уже нашли

print(longest_common_substring)

