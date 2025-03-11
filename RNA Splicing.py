import os

os.chdir("/Users/annaklimova/Downloads")
with open("rosalind_splc.txt", "r") as file:
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
mRNA = first_value
for introns in except_first:
    mRNA = mRNA.replace(introns, "")

mRNA = mRNA.replace("T", "U")
print(mRNA)

genetic_code = {
        "UUU": "F", "CUU": "L", "AUU": "I", "GUU": "V",
        "UUC": "F", "CUC": "L", "AUC": "I", "GUC": "V",
        "UUA": "L", "CUA": "L", "AUA": "I", "GUA": "V",
        "UUG": "L", "CUG": "L", "AUG": "M", "GUG": "V",
        "UCU": "S", "CCU": "P", "ACU": "T", "GCU": "A",
        "UCC": "S", "CCC": "P", "ACC": "T", "GCC": "A",
        "UCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
        "UCG": "S", "CCG": "P", "ACG": "T", "GCG": "A",
        "UAU": "Y", "CAU": "H", "AAU": "N", "GAU": "D",
        "UAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
        "UAA": "Stop", "CAA": "Q", "AAA": "K", "GAA": "E",
        "UAG": "Stop", "CAG": "Q", "AAG": "K", "GAG": "E",
        "UGU": "C", "CGU": "R", "AGU": "S", "GGU": "G",
        "UGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
        "UGA": "Stop", "CGA": "R", "AGA": "R", "GGA": "G",
        "UGG": "W", "CGG": "R", "AGG": "R", "GGG": "G"
    }
protein = []
for n in range(0, len(mRNA), 3):
    triplet = mRNA[n : n + 3]
    amino = genetic_code[triplet]
    protein.append(amino)


print("".join(protein[0:len(protein)-1]))