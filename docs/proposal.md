# Proposal Singkat Proyek

## Judul
Identifikasi Kandidat Target Obat Berbasis Protein Membran melalui Analisis Sekuens dan Prediksi Topologi Transmembran

## Latar Belakang
Protein membran berperan penting dalam transport molekul, komunikasi sel, adhesi, dan interaksi sel dengan lingkungan. Karena banyak protein membran terletak pada permukaan sel atau membentuk kanal/transporter, protein ini sering menjadi kandidat awal yang menarik dalam pencarian target obat. Namun, jumlah protein dalam satu organisme dapat sangat banyak, sehingga diperlukan pendekatan komputasi untuk melakukan penyaringan awal.

## Rumusan Masalah
Bagaimana analisis sekuens protein dan prediksi topologi transmembran dapat digunakan untuk mengidentifikasi serta memprioritaskan protein membran sebagai kandidat target obat?

## Tujuan
1. Mengembangkan pipeline Python untuk membaca dan memproses sekuens protein.
2. Mengekstraksi fitur fisikokimia protein, seperti panjang, berat molekul, pI, GRAVY, dan instability index.
3. Memprediksi segmen transmembran menggunakan pendekatan hidrofobisitas berbasis sliding window.
4. Menghasilkan ranking kandidat protein membran berdasarkan skor prioritas.
5. Menyediakan visualisasi dan interpretasi biologis dari hasil komputasi.

## Dataset
Dataset input berupa sekuens protein dalam format FASTA. Untuk demo awal digunakan dataset kecil yang berisi contoh protein membran dan non-membran. Untuk pengujian lebih valid, dataset dapat diganti dengan sekuens protein dari UniProt.

## Metode Komputasi
1. FASTA parsing menggunakan Biopython.
2. Feature extraction menggunakan Biopython ProtParam.
3. Prediksi segmen transmembran menggunakan skala hidropati Kyte-Doolittle dengan sliding window 19 asam amino.
4. Rule-based scoring untuk menghasilkan prioritas kandidat.
5. Visualisasi menggunakan matplotlib.

## Output
1. Tabel fitur protein.
2. Tabel ranking kandidat target obat.
3. Ringkasan otomatis hasil analisis.
4. Grafik ranking kandidat, hubungan GRAVY dan segmen transmembran, distribusi kelas prediksi, dan peta posisi segmen transmembran.

## Batasan
Prediksi lokal berbasis hidrofobisitas merupakan baseline yang mudah dijelaskan, bukan pengganti model khusus seperti DeepTMHMM. Hasil perlu divalidasi lebih lanjut menggunakan DeepTMHMM, SignalP, InterProScan, BLAST, docking, atau eksperimen biologis.
