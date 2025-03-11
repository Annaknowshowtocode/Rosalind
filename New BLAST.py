import os
import logging
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio.Blast import NCBIXML

# Настроим логирование ошибок
logging.basicConfig(filename="blast_errors.log", level=logging.ERROR, format="%(asctime)s - %(message)s")

# Определяем пути
script_dir = os.path.dirname(os.path.abspath(__file__))
blast_results_path = os.path.join(script_dir, "blast_results")
os.makedirs(blast_results_path, exist_ok=True)

# Читаем входной файл
file_path = os.path.join(script_dir, "output_table.csv")
df = pd.read_csv(file_path) if file_path.endswith(".csv") else pd.read_excel(file_path, engine="openpyxl")
df = df.dropna(subset=["oligo_sequence"]).iloc[3:40].reset_index(drop=True)  # Выбираем строки 3-40


# Функция для проверки, является ли транскрипт предсказанным
def is_predicted(hit_id):
    return any(prefix in hit_id for prefix in ["XP_", "XM_"])


# Функция для обработки BLAST-результатов
def parse_blast_result(file):
    with open(file) as result_handle:
        blast_record = NCBIXML.read(result_handle)

    if not blast_record.alignments:
        return None

    top_hit = blast_record.alignments[0]
    top_hit_id = top_hit.hit_id
    is_top_predicted = is_predicted(top_hit_id)

    top_metrics = {
        "top1_hit": top_hit_id,
        "top1_score": top_hit.hsps[0].score,
        "top1_evalue": top_hit.hsps[0].expect,
        "top1_identity": top_hit.hsps[0].identities / top_hit.hsps[0].align_length,
        "top1_predicted": is_top_predicted
    }

    # Ищем первое вхождение target-гена
    match_hit = None
    match_rank = None
    for i, alignment in enumerate(blast_record.alignments, start=1):
        if "target_gene" in alignment.hit_def:  # <- Уточнить критерий
            match_hit = alignment
            match_rank = i
            break

    match_metrics = {}
    if match_hit:
        match_metrics = {
            "match_hit": match_hit.hit_id,
            "match_score": match_hit.hsps[0].score,
            "match_evalue": match_hit.hsps[0].expect,
            "match_identity": match_hit.hsps[0].identities / match_hit.hsps[0].align_length,
            "match_rank": match_rank,
            "match_predicted": is_predicted(match_hit.hit_id)
        }

    return {**top_metrics, **match_metrics}


# Список файлов с результатами BLAST
blast_files = [os.path.join(blast_results_path, f) for f in os.listdir(blast_results_path) if f.endswith(".xml")]

# Парсим результаты
results = []
with ThreadPoolExecutor() as executor:
    future_to_file = {executor.submit(parse_blast_result, file): file for file in blast_files}
    for future in as_completed(future_to_file):
        try:
            result = future.result()
            if result:
                results.append(result)
        except Exception as e:
            logging.error(f"Ошибка при обработке {future_to_file[future]}: {e}")

# Создаем DataFrame с результатами
blast_df = pd.DataFrame(results)

# Разделяем на предсказанные / реальные
blast_df_real = blast_df[blast_df["top1_predicted"] == False].add_suffix("_real")
blast_df_pred = blast_df[blast_df["top1_predicted"] == True].add_suffix("_predicted")

# Объединяем
final_df = pd.concat([df, blast_df_real, blast_df_pred], axis=1)

# Сохраняем
final_output = os.path.join(script_dir, "blast_parsed_results.csv")
final_df.to_csv(final_output, index=False)

print(f"Результаты сохранены в {final_output}")
