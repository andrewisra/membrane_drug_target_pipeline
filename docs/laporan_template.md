# Identifikasi Kandidat Target Obat Berbasis Protein Membran melalui Analisis Sekuens dan Prediksi Topologi Transmembran

## 1. Pendahuluan

Protein membran berperan penting dalam transport molekul, komunikasi sel, adhesi, dan interaksi sel dengan lingkungan. Karena lokasinya di membran atau permukaan sel, protein membran sering menjadi kandidat awal target obat. Penelitian ini mengembangkan pipeline bioinformatika berbasis Python untuk menganalisis sekuens protein dan memprioritaskan kandidat protein membran berdasarkan fitur fisikokimia dan prediksi segmen transmembran.

Rumusan masalah:
Bagaimana analisis sekuens protein dan prediksi topologi transmembran dapat digunakan untuk mengidentifikasi kandidat protein membran yang berpotensi menjadi target obat?

Tujuan:
1. Mengekstraksi fitur fisikokimia dari sekuens protein.
2. Memprediksi segmen transmembran berdasarkan pola hidrofobisitas.
3. Menghasilkan ranking kandidat protein membran sebagai prioritisasi awal target obat.

## 2. Metode

Input penelitian berupa sekuens protein dalam format FASTA. Setiap protein dianalisis menggunakan Biopython ProtParam untuk menghitung panjang protein, berat molekul, titik isoelektrik, GRAVY, instability index, aromaticity, dan persentase asam amino hidrofobik.

Prediksi segmen transmembran dilakukan menggunakan pendekatan sliding window berbasis skala hidropati Kyte-Doolittle. Window dengan rata-rata hidropati tinggi dan proporsi residu hidrofobik tinggi dianggap sebagai kandidat segmen transmembran.

Skor kandidat dihitung berdasarkan jumlah segmen transmembran, nilai hidrofobisitas, prediksi signal peptide sederhana, ukuran protein, dan instability index. Protein dengan skor tinggi diprioritaskan sebagai kandidat target obat berbasis protein membran.

## 3. Hasil dan Diskusi

Hasil utama berupa tabel fitur protein, tabel ranking kandidat, dan visualisasi. Protein dengan jumlah segmen transmembran lebih tinggi dan nilai GRAVY positif cenderung diklasifikasikan sebagai protein membran. Kandidat dengan skor tertinggi menunjukkan kombinasi fitur membran, hidrofobisitas, stabilitas, dan ukuran protein yang mendukung prioritas analisis lebih lanjut.

## 4. Kesimpulan

Pipeline ini berhasil mengintegrasikan analisis sekuens, ekstraksi fitur fisikokimia, prediksi segmen transmembran, dan sistem scoring untuk menghasilkan prioritisasi awal kandidat target obat berbasis protein membran. Untuk penelitian lanjutan, pipeline dapat dikembangkan dengan integrasi DeepTMHMM, SignalP, InterProScan, BLAST terhadap proteom manusia, serta analisis struktur 3D dan docking molekuler.

## 5. Daftar Pustaka

Isi dengan sumber resmi Biopython, UniProt, DeepTMHMM, dan literatur protein membran/target obat.
