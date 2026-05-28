# Identifikasi Kandidat Target Obat Berbasis Protein Membran melalui Analisis Sekuens dan Prediksi Topologi Transmembran

Pipeline bioinformatika berbasis Python untuk menganalisis sekuens protein, memprediksi segmen transmembran secara baseline, mengekstraksi fitur fisikokimia, mengevaluasi hasil terhadap label dataset, dan memprioritaskan protein membran sebagai kandidat target obat awal.

## Mengapa proyek ini bermakna

Protein membran sering menjadi target biologis penting karena terlibat dalam transport molekul, sinyal seluler, adhesi, dan interaksi sel dengan lingkungan. Pada bakteri patogen, protein membran/permukaan lebih mudah diakses oleh molekul obat atau antibodi dibanding protein sitoplasmik. Pipeline ini bukan alat penemuan obat final, tetapi alat penyaringan awal untuk memilih protein yang layak dianalisis lebih lanjut.

## Input yang didukung

### 1. FASTA protein

Wajib. Contoh:

```text
data/proteins.fasta
```

### 2. Metadata CSV dari tim pengumpul data

Opsional tetapi sangat disarankan untuk hasil final. Template tersedia di:

```text
data/data_collection_template.csv
```

Kolom yang disarankan:

```text
accession,label,organism,protein_name,function_hint,subcellular_location,is_essential,human_homolog,source_url,notes
```

Lihat instruksi lengkap di:

```text
docs/instruksi_pengumpulan_data.md
```

### 3. Hasil DeepTMHMM/TMHMM yang dinormalisasi

Opsional. Template tersedia di:

```text
data/external_tm_predictions_template.csv
```

Format:

```text
protein_id,region_type,start,end,source
```

## Output utama

- `results/protein_features.csv`: fitur tiap protein.
- `results/candidate_ranking.csv`: ranking kandidat target obat.
- `results/evaluation_report.txt`: evaluasi baseline jika dataset memiliki label.
- `results/summary_report.txt`: ringkasan hasil biologis.
- `figures/*.png`: grafik ranking, distribusi prioritas, scatter plot, peta segmen transmembran, dan confusion matrix.

## Instalasi

```bash
pip install -r requirements.txt
```

## Jalankan demo offline

```bash
python src/main.py --fasta data/proteins.fasta --output-dir results --figure-dir figures
```

## Jalankan dengan data teman/tim

1. Minta teman mengisi:

```text
data/friend_metadata.csv
```

berdasarkan template `data/data_collection_template.csv`.

2. Ambil FASTA otomatis dari UniProt berdasarkan accession CSV:

```bash
python src/fetch_uniprot.py fetch-accessions --input-csv data/friend_metadata.csv --output-fasta data/friend_proteins.fasta
```

3. Jalankan pipeline final:

```bash
python src/main.py --fasta data/friend_proteins.fasta --metadata-csv data/friend_metadata.csv --output-dir results_friend --figure-dir figures_friend
```

4. Jika ada hasil DeepTMHMM/TMHMM:

```bash
python src/main.py --fasta data/friend_proteins.fasta --metadata-csv data/friend_metadata.csv --external-tm-csv data/external_tm_predictions.csv --output-dir results_friend --figure-dir figures_friend
```

## Cari kandidat protein dari UniProt

Contoh pencarian protein membran *E. coli* K-12 reviewed:

```bash
python src/fetch_uniprot.py search --query "(organism_id:83333) AND (cc_subcellular_location:membrane) AND reviewed:true" --output-tsv data/ecoli_membrane.tsv --size 30
```

Contoh pencarian protein sitoplasmik/non-membran *E. coli* K-12 reviewed:

```bash
python src/fetch_uniprot.py search --query "(organism_id:83333) AND (cc_subcellular_location:cytoplasm) AND reviewed:true" --output-tsv data/ecoli_non_membrane.tsv --size 30
```

## Metode ringkas

1. Input sekuens protein FASTA.
2. Bersihkan sekuens menjadi 20 asam amino standar.
3. Hitung fitur fisikokimia dengan Biopython ProtParam.
4. Prediksi segmen transmembran menggunakan sliding window Kyte-Doolittle hydropathy.
5. Prediksi kasar signal peptide N-terminal.
6. Gabungkan metadata biologis dari CSV jika tersedia.
7. Hitung skor kandidat berdasarkan bukti membran, hidrofobisitas, aksesibilitas, sinyal sekresi, fungsi, esensialitas, risiko homolog manusia, dan penalti stabilitas/panjang.
8. Evaluasi terhadap label dataset jika tersedia.
9. Hasilkan ranking, visualisasi, dan interpretasi biologis.

## Batasan ilmiah

- Prediksi transmembran lokal adalah baseline explainable, bukan pengganti DeepTMHMM/TMHMM.
- Signal peptide hanya diprediksi dengan heuristik sederhana; validasi ideal memakai SignalP.
- Skor kandidat adalah prioritisasi komputasional awal, bukan bukti bahwa protein pasti target obat.
- Untuk studi lanjutan, kandidat perlu divalidasi dengan DeepTMHMM, SignalP, InterProScan, BLASTp terhadap protein manusia, data esensialitas gen, dan literatur biologis.

## Struktur folder

```text
data/       input FASTA, accession CSV, template data
src/        kode pipeline
tests/      unit test sederhana
results/    output CSV dan ringkasan
figures/    output visualisasi
docs/       instruksi data, proposal, laporan final, flowchart, outline video
```

## Jalankan test

```bash
pytest -q
```
