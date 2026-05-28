# Identifikasi Kandidat Target Obat Berbasis Protein Membran melalui Analisis Sekuens dan Prediksi Topologi Transmembran

## 1. Pendahuluan

Protein membran merupakan kelompok makromolekul protein yang terintegrasi atau berasosiasi dengan membran biologis. Protein ini berperan dalam transport molekul, transduksi sinyal, adhesi, pengenalan sel, dan komunikasi antara sel dengan lingkungannya. Karena posisinya yang berada di membran atau permukaan sel, protein membran sering menjadi kandidat penting dalam studi target obat, terutama pada organisme patogen.

Permasalahan biologis yang diangkat pada proyek ini adalah bagaimana menyeleksi protein yang berpotensi menjadi kandidat target obat hanya dari data sekuens protein. Pendekatan bioinformatika memungkinkan ekstraksi fitur fisikokimia dan prediksi topologi transmembran secara cepat sebelum dilakukan validasi eksperimental.

**Rumusan masalah:** bagaimana pipeline bioinformatika berbasis Python dapat mengidentifikasi dan memprioritaskan protein membran sebagai kandidat target obat berdasarkan sekuens protein?

**Tujuan:** mengembangkan pipeline komputasi untuk mengekstraksi fitur protein, memprediksi segmen transmembran, mengevaluasi hasil prediksi, dan menghasilkan ranking kandidat target obat awal.

## 2. Metode

### 2.1 Dataset

Input utama berupa file FASTA protein. Setiap header FASTA dapat diberi metadata sederhana, misalnya:

```text
>ProteinA label=membrane | organism=Escherichia_coli | function=outer_membrane_channel
```

Label `membrane` dan `non_membrane` digunakan sebagai pembanding evaluasi baseline. Untuk analisis nyata, sekuens dapat diambil dari UniProt menggunakan daftar accession pada `data/uniprot_accessions_example.csv`.

### 2.2 Ekstraksi fitur fisikokimia

Pipeline menghitung beberapa fitur dari sekuens protein:

- panjang protein,
- berat molekul,
- titik isoelektrik,
- GRAVY/hidrofobisitas,
- instability index,
- aromaticity,
- persentase asam amino hidrofobik,
- persentase asam amino bermuatan positif dan negatif.

Fitur ini digunakan karena protein membran umumnya memiliki segmen hidrofobik yang berinteraksi dengan lipid bilayer.

### 2.3 Prediksi topologi transmembran baseline

Prediksi dilakukan dengan pendekatan sliding window menggunakan skala hidrofobisitas Kyte-Doolittle. Window sepanjang 19 asam amino dianggap kandidat segmen transmembran jika memenuhi:

1. rata-rata hidrofobisitas melebihi ambang tertentu,
2. fraksi residu hidrofobik cukup tinggi,
3. panjang segmen gabungan masuk akal untuk heliks transmembran.

Metode ini digunakan sebagai baseline yang mudah dijelaskan. Untuk validasi lanjutan, hasilnya dapat dibandingkan dengan DeepTMHMM atau TMHMM.

### 2.4 Skoring kandidat target obat

Skor kandidat dihitung dari beberapa komponen:

- bukti segmen transmembran,
- hidrofobisitas,
- aksesibilitas prediksi,
- indikasi signal peptide,
- petunjuk fungsi pada metadata,
- penalti untuk panjang ekstrem atau instability index tinggi,
- penalti jika tidak ditemukan segmen transmembran.

Protein dengan skor tinggi diprioritaskan sebagai kandidat untuk analisis lanjutan, bukan langsung dianggap sebagai target obat final.

### 2.5 Evaluasi

Jika dataset memiliki label, pipeline menghitung confusion matrix, accuracy, precision, recall, dan F1-score. Evaluasi ini menunjukkan apakah baseline prediksi transmembran cukup mampu membedakan protein membran dan non-membran pada dataset yang digunakan.

## 3. Hasil dan Diskusi

Output utama pipeline adalah file `candidate_ranking.csv`. File ini berisi ranking protein berdasarkan skor kandidat, jumlah segmen transmembran, nilai GRAVY, prediksi aksesibilitas, dan interpretasi biologis.

Kandidat prioritas tinggi umumnya memiliki setidaknya satu segmen transmembran, hidrofobisitas relatif tinggi, dan indikasi aksesibilitas membran. Protein tanpa segmen transmembran cenderung mendapat skor rendah karena kemungkinan besar merupakan protein sitoplasmik atau protein larut.

Visualisasi `gravy_vs_tm_segments.png` membantu melihat hubungan antara hidrofobisitas dan jumlah segmen transmembran. Visualisasi `tm_segment_map.png` memperlihatkan posisi perkiraan segmen transmembran pada protein kandidat teratas. Evaluasi pada `evaluation_report.txt` digunakan untuk menilai performa baseline terhadap label dataset.

Secara biologis, protein membran dengan fungsi transporter, kanal, reseptor, porin, atau protein permukaan memiliki alasan yang lebih kuat untuk dianalisis sebagai kandidat target obat awal. Namun, prioritas akhir harus mempertimbangkan faktor tambahan seperti esensialitas gen, konservasi pada patogen, ketiadaan homolog dekat pada manusia, struktur 3D, dan bukti literatur.

## 4. Kesimpulan

Proyek ini menghasilkan pipeline bioinformatika berbasis Python untuk mengidentifikasi dan memprioritaskan protein membran berdasarkan sekuens protein. Pipeline menggabungkan ekstraksi fitur fisikokimia, prediksi segmen transmembran baseline, skoring kandidat, evaluasi terhadap label, dan visualisasi hasil.

Kontribusi utama proyek adalah menyediakan kerangka penyaringan awal yang transparan dan dapat direproduksi untuk memilih protein membran yang layak dianalisis lebih lanjut sebagai kandidat target obat. Untuk pengembangan masa depan, pipeline dapat diintegrasikan dengan DeepTMHMM, SignalP, InterProScan, BLASTp, AlphaFold, dan data esensialitas gen.

## 5. Daftar Pustaka Awal

1. UniProt Consortium. UniProt: the Universal Protein Knowledgebase.
2. Kyte, J. and Doolittle, R. F. A simple method for displaying the hydropathic character of a protein.
3. Hallgren, J. et al. DeepTMHMM predicts alpha and beta transmembrane proteins using deep neural networks.
4. Teufel, F. et al. SignalP 6.0 predicts all five types of signal peptides using protein language models.
5. Blum, M. et al. InterPro: the protein sequence classification resource.
