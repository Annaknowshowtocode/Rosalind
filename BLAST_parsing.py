import pandas as pd
import os
import subprocess
from Bio.Blast import NCBIXML
from Bio.Seq import Seq

# 1. Загрузка данных
input_file = "/Users/annaklimova/Downloads/output_table.csv"
output_file = "/Users/annaklimova/Downloads/blast_results.csv"

df = pd.read_csv(input_file)

# 2. Функция для создания комплементарной антипараллельной последовательности
def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

# 3. BLAST для каждой последовательности
blast_results = []

for index, row in df.iterrows():
    original_seq = row["oligo_sequence"]
    antisense_seq = reverse_complement(original_seq)

    # Запись последовательности в FASTA-файл
    fasta_filename = f"temp_query_{index}.fasta"
    with open(fasta_filename, "w") as f:
        f.write(f">query_{index}\n{antisense_seq}\n")

    # Выполнение BLAST
    blast_output = f"blast_result_{index}.xml"
    blast_cmd = "/opt/homebrew/bin/blastn -query {} -db refseq_rna -out {} -outfmt 5 -max_target_seqs 1 -evalue 0.01".format(
        fasta_filename, blast_output
    )
    subprocess.run(blast_cmd, shell=True)

    # Чтение результатов BLAST
    with open(blast_output) as result_handle:
        blast_record = NCBIXML.read(result_handle)

    best_predicted, best_real = None, None
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            hit_info = {
                "query_sequence": original_seq,
                "antisense_sequence": antisense_seq,
                "hit_id": alignment.accession,
                "hit_description": alignment.title,
                "evalue": hsp.expect,
                "identity": hsp.identities / hsp.align_length * 100,
                "alignment_length": hsp.align_length,
                "query_start": hsp.query_start,
                "query_end": hsp.query_end,
                "subject_start": hsp.sbjct_start,
                "subject_end": hsp.sbjct_end
            }

            # Определяем predicted и real
            if "PREDICTED" in alignment.title and (best_predicted is None or hit_info["evalue"] < best_predicted["evalue"]):
                best_predicted = hit_info
            elif "PREDICTED" not in alignment.title and (best_real is None or hit_info["evalue"] < best_real["evalue"]):
                best_real = hit_info

    # Добавляем результат
    blast_results.append({
        "original_sequence": original_seq,
        "antisense_sequence": antisense_seq,
        "best_predicted": best_predicted["hit_description"] if best_predicted else "None",
        "best_real": best_real["hit_description"] if best_real else "None",
        "evalue_predicted": best_predicted["evalue"] if best_predicted else "N/A",
        "evalue_real": best_real["evalue"] if best_real else "N/A"
    })

    # Удаляем временные файлы
    os.remove(fasta_filename)
    os.remove(blast_output)

# 4. Запись результатов
result_df = pd.DataFrame(blast_results)
result_df.to_csv(output_file, index=False)

print(f"Результаты сохранены в {output_file}")
