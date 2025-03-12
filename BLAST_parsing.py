from Bio.Blast import NCBIWWW

sequence = "AGCTGCTAGCTGATCG"  # Заменить на свою последовательность
result_handle = NCBIWWW.qblast("blastn", "nt", sequence)

with open("/Users/annaklimova/Git/Rosalind/blast_results/test.xml", "w") as save_file:
    save_file.write(result_handle.read())
result_handle.close()
