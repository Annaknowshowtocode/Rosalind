import os
import time
import logging
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import ssl
import certifi
ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())
from Bio.Blast import NCBIWWW, NCBIXML

#  ПУТИ: ВСЕ ФАЙЛЫ ХРАНЯТСЯ В /Users/annaklimova/Git/Rosalind/BLAST
base_path = "/Users/annaklimova/Git/Rosalind/BLAST"
blast_results_path = os.path.join(base_path, "blast_results")
blast_queries_path = os.path.join(base_path, "blast_queries")
file_path = os.path.join(base_path, "output_table.csv")
log_path = os.path.join(base_path, "blast_errors.log")

# Создаём папки
os.makedirs(blast_results_path, exist_ok=True)
os.makedirs(blast_queries_path, exist_ok=True)

#  Настраиваем логирование
logging.basicConfig(filename=log_path, level=logging.ERROR, format="%(asctime)s - %(message)s")

#  Проверка соединения с NCBI
import urllib.request

def check_ncbi():
    try:
        urllib.request.urlopen("https://blast.ncbi.nlm.nih.gov/Blast.cgi", timeout=5)
        return True
    except Exception:
        return False

if not check_ncbi():
    print(" NCBI BLAST недоступен! Проверь интернет или VPN.")
    exit(1)

#  Загружаем данные
df = pd.read_csv(file_path) if file_path.endswith(".csv") else pd.read_excel(file_path, engine="openpyxl")
df = df.dropna(subset=["oligo_sequence"]).reset_index(drop=True)  # Убираем NaN
sequences = df["oligo_sequence"].iloc[39:60]  # Берём первые 15 последовательностей

#  Функция для перевода последовательности в мРНК
def get_mRNA(sequence):
    complement = {"A": "U", "T": "A", "C": "G", "G": "C",
                  "a": "u", "t": "a", "c": "g", "g": "c"}
    return "".join(complement.get(nuc.upper(), "") for nuc in sequence)[::-1]  # Реверсируем

#  Сохраняем последовательность в FASTA-файл
def save_fasta(sequence, index):
    fasta_file = os.path.join(blast_queries_path, f"query_{index}.fasta")
    with open(fasta_file, "w") as f:
        f.write(f">query_{index}\n{sequence}\n")
    return fasta_file

#  Отправляем запрос на NCBI BLAST и сохраняем результат
def run_remote_blast(sequence, index, max_retries=3):
    blast_file = os.path.join(blast_results_path, f"blast_result_{index}.xml")
    fasta_file = save_fasta(sequence, index)  # Сохраняем FASTA-запрос

    # Если результат уже есть, не отправляем повторный запрос
    if os.path.exists(blast_file):
        return parse_blast_results(blast_file, index)

    attempt = 0
    while attempt < max_retries:
        try:
            print(f"🔄 Отправка запроса {index} на NCBI BLAST...")
            result_handle = NCBIWWW.qblast("blastn", "nt", sequence,
                                           hitlist_size=1,
                                           expect=0.05,
                                           format_type="XML",
                                           entrez_query="mus.anna.study@gmail.com")

            # Сохраняем результат в XML-файл
            with open(blast_file, "w") as out_file:
                out_file.write(result_handle.read())

            return parse_blast_results(blast_file, index)

        except Exception as e:
            logging.error(f"Ошибка BLAST {index}, попытка {attempt + 1}: {e}")
            print(f"❌ Ошибка запроса {index}, попытка {attempt + 1}: {e}")
            time.sleep(3)  # Ждём перед повторной попыткой
            attempt += 1

    return [index, "ERROR", "ERROR", "ERROR"]

#  Парсим результаты BLAST
def parse_blast_results(blast_file, index):
    if not os.path.exists(blast_file):
        return [index, "NA", "NA", "NA"]

    try:
        with open(blast_file) as result_handle:
            blast_record = NCBIXML.read(result_handle)

        if not blast_record.alignments:
            return [index, "NA", "NA", "NA"]

        top_hit = blast_record.alignments[0]
        hsp = top_hit.hsps[0]

        return [index, top_hit.title, hsp.identities / hsp.align_length, hsp.expect]

    except Exception as e:
        logging.error(f"Ошибка парсинга XML {index}: {e}")
        print(f"❌ Ошибка парсинга XML {index}: {e}")
        return [index, "ERROR", "ERROR", "ERROR"]

# ️ Запуск BLAST в многопоточном режиме (до 2 потоков, чтобы не заблокировали)
num_threads = min(2, len(sequences))
results = []
with ThreadPoolExecutor(max_workers=num_threads) as executor:
    future_to_index = {executor.submit(run_remote_blast, get_mRNA(seq.strip()), i): i for i, seq in enumerate(sequences)}
    for future in as_completed(future_to_index):
        try:
            results.append(future.result())
        except Exception as e:
            logging.error(f"Ошибка выполнения потока: {e}")
            print(f"❌ Ошибка выполнения потока: {e}")

#  Добавляем результаты в DataFrame
for res in results:
    index, top1_hit, top1_real, top1_predicted = res
    df.at[index, "top1_hit"] = top1_hit
    df.at[index, "top1_real"] = top1_real
    df.at[index, "top1_predicted"] = top1_predicted

#  Сохраняем обновленный файл
df.to_csv(file_path, index=False, na_rep="NA")
print(f"\n✅ Файл обновлен с результатами BLAST: {file_path}")
