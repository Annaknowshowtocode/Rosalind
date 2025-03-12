import os
import time
import logging
import pandas as pd
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

# Настроим логирование ошибок
logging.basicConfig(filename="blast_errors.log", level=logging.ERROR, format="%(asctime)s - %(message)s")

# Определяем пути
script_dir = os.path.dirname(os.path.abspath(__file__))
blast_results_path = os.path.join(script_dir, "blast_results")
file_path = os.path.join(script_dir, "output_table.csv")
blast_db_path = "/home/annaklimova/Git/BLAST/"

# Создадим папку для BLAST-результатов, если её нет
os.makedirs(blast_results_path, exist_ok=True)

# Загрузка данных
df = pd.read_csv(file_path) if file_path.endswith(".csv") else pd.read_excel(file_path, engine="openpyxl")

# Получаем последовательности
df = df.dropna(subset=["oligo_sequence"]).reset_index(drop=True)
sequences = df["oligo_sequence"].iloc[:100]


# Функция для перевода последовательности в мРНК
def get_mRNA(sequence):
    complement = {"A": "U", "T": "A", "C": "G", "G": "C", "a": "u", "t": "a", "c": "g", "g": "c"}
    return "".join(complement.get(nuc.upper(), "") for nuc in sequence)


# Функция для локального запуска BLAST
def run_local_blast(sequence, index):
    blast_file = os.path.join(blast_results_path, f"blast_result_{index}.xml")
    query_file = os.path.join(blast_results_path, f"query_{index}.fasta")

    # Записываем последовательность во временный FASTA-файл
    with open(query_file, "w") as f:
        f.write(f">query_{index}\n{sequence}\n")

    # Команда для локального BLAST
    blast_cmd = [
        "blastn",
        "-query", query_file,
        "-db", blast_db_path,
        "-out", blast_file,
        "-outfmt", "5",  # XML-формат для совместимости с NCBIXML
        "-num_threads", "4",
        "-max_target_seqs", "1",
        "-evalue", "0.05"
    ]

    try:
        subprocess.run(blast_cmd, check=True)
        return blast_file
    except subprocess.CalledProcessError as e:
        logging.error(f"Ошибка BLAST {index}: {e}")
        return None


# Функция для парсинга результатов BLAST
def parse_blast_results(blast_file):
    from Bio.Blast import NCBIXML

    if not blast_file or not os.path.exists(blast_file):
        return ["NA"] * 3

    with open(blast_file) as result_handle:
        blast_record = NCBIXML.read(result_handle)

    if not blast_record.alignments:
        return ["NA"] * 3

    top_hit = blast_record.alignments[0]
    hsp = top_hit.hsps[0]

    return [
        top_hit.title,  # top1_hit
        hsp.identities / hsp.align_length,  # top1_real
        hsp.expect  # top1_predicted
    ]


# Запуск BLAST в многопоточном режиме
def process_sequence(index, seq):
    try:
        mRNA_seq = get_mRNA(seq.strip())
        if not mRNA_seq:
            return [index, "NA", "NA", "NA"]

        blast_file = run_local_blast(mRNA_seq, index)
        return [index] + parse_blast_results(blast_file)
    except Exception as e:
        logging.error(f"Ошибка при обработке последовательности {index}: {e}")
        return [index, "ERROR", "ERROR", "ERROR"]


# Выполняем BLAST для каждой последовательности
num_threads = min(4, len(sequences))  # Используем до 4 потоков
results = []
with ThreadPoolExecutor(max_workers=num_threads) as executor:
    future_to_index = {executor.submit(process_sequence, i, seq): i for i, seq in enumerate(sequences)}
    for future in as_completed(future_to_index):
        try:
            results.append(future.result())
        except Exception as e:
            logging.error(f"Ошибка выполнения потока: {e}")

# Добавляем результаты в DataFrame
for res in results:
    index, top1_hit, top1_real, top1_predicted = res
    df.at[index, "top1_hit"] = top1_hit
    df.at[index, "top1_real"] = top1_real
    df.at[index, "top1_predicted"] = top1_predicted

# Сохраняем обновленный файл
df.to_csv(file_path, index=False, na_rep="NA")
print(f"\n✅ Файл обновлен с результатами BLAST: {file_path}")
