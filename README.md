# Identifikasi Kandidat Target Obat Berbasis Protein Membran melalui Analisis Sekuens dan Prediksi Topologi Transmembran

Pipeline bioinformatika berbasis Python untuk menganalisis sekuens protein, memprediksi segmen transmembran secara baseline, mengekstraksi fitur fisikokimia, serta memprioritaskan protein membran sebagai kandidat target obat awal.

---

## Latar Belakang

Protein membran memiliki peran penting dalam transport molekul, komunikasi sel, dan interaksi dengan lingkungan. Karena sifatnya yang mudah diakses dari luar sel, protein membran sering digunakan sebagai target dalam pengembangan obat dan vaksin.

Proyek ini bertujuan untuk membangun pipeline bioinformatika sederhana yang dapat membantu proses penyaringan awal kandidat target obat berbasis analisis sekuens protein dan prediksi topologi transmembran.

---

## Dataset

Dataset menggunakan protein *Escherichia coli* K-12 dari UniProtKB dengan status *reviewed*. Dataset terdiri atas:
- 30 protein membrane
- 30 protein non membrane

File dataset:

```bash
data/data_collection.csv
data/proteins.fasta
```

---

## Metode

Pipeline melakukan beberapa tahapan utama:
1. Membaca sekuens protein dari file FASTA.
2. Membersihkan sekuens menjadi 20 asam amino standar.
3. Menghitung fitur fisikokimia protein.
4. Memprediksi segmen transmembran menggunakan metode Kyte-Doolittle sliding window.
5. Menghasilkan skor prioritas kandidat protein membran.
6. Mengevaluasi hasil prediksi menggunakan confusion matrix.
7. Membuat visualisasi hasil analisis.

---

## Instalasi

```bash
pip install -r requirements.txt
```

---

## Menjalankan Pipeline

```bash
python src/main.py --fasta data/proteins.fasta --metadata-csv data/data_collection.csv --output-dir results --figure-dir figures
```

---

## Output

Pipeline menghasilkan beberapa file output berikut:

```text
results/protein_features.csv
results/candidate_ranking.csv 
results/evaluation_report.txt 
results/summary_report.txt
```

Visualisasi hasil disimpan pada folder:
```text
figures/
```

---

## Struktur Folder

```text
data/           # dataset FASTA dan metadata 
src/            # source code pipeline 
tests/          # unit testing 
results/        # hasil analisis 
figures/        # visualisasi hasil 
docs/           # laporan dan dokumentasi
```

---

## Batasan

- Prediksi topologi transmembran masih menggunakan pendekatan baseline berbasis hidrofobisitas.
- Pipeline belum menggunakan model prediksi lanjutan seperti DeepTMHMM atau SignalP.
- Hasil yang diperoleh merupakan prioritas kandidat awal dan masih memerlukan validasi biologis lebih lanjut.

---

## Teknologi yang Digunakan

- Python
- Biopython
- Pandas
- NumPy
- Matplotlib

---

## Kontributor

| Nama | NIM | 
|---|---| 
| Anella Utari Gunadi | 13523078 | 
| Nayla Zahira | 13523079 | 
| Muhammad Izzat Jundy | 13523092 | 
| Andrew Isra Saputra DB | 13523110 |