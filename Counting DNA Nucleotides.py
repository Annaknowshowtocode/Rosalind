import os
os.chdir('/Users/annaklimova/Downloads')
with open('rosalind_dna.txt', 'r') as file:  # Открываем файл в режиме чтения
    content = file.read()  # Читаем весь файл целиком
    content = content.strip()
    nucleotides = ["A", "C", "G", "T"]
    n = 0
    result = []
    for i in nucleotides:
        counts = content.count(str(i))
        result.append(str(counts))
        n = 0
    print(" ".join(result))