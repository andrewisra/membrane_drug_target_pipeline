# Outline Video Presentasi Maksimal 10 Menit

## 1. Latar belakang, 1 menit

- Makromolekul yang dipilih: protein.
- Protein membran penting dalam transport, sinyal, dan interaksi sel.
- Karena lokasinya di membran/permukaan, protein ini relevan sebagai kandidat target obat.

## 2. Rumusan masalah dan tujuan, 1 menit

- Masalah: bagaimana memilih kandidat protein membran dari sekuens protein?
- Tujuan: membuat pipeline Python untuk analisis sekuens, prediksi topologi transmembran, scoring, dan ranking kandidat.

## 3. Dataset, 1 menit

- Input berupa FASTA.
- Demo memakai dataset kecil berlabel.
- Versi nyata dapat mengambil protein dari UniProt memakai accession.

## 4. Metode, 2 menit

- Ekstraksi fitur fisikokimia dengan Biopython.
- Prediksi segmen transmembran baseline dengan Kyte-Doolittle sliding window.
- Scoring kandidat berdasarkan TM segment, hidrofobisitas, aksesibilitas, signal peptide, dan fungsi.
- Evaluasi dengan confusion matrix jika label tersedia.

## 5. Demo program, 2 menit

Jalankan:

```bash
python src/main.py --fasta data/proteins.fasta --output-dir results --figure-dir figures
```

Tunjukkan:

- `candidate_ranking.csv`,
- `summary_report.txt`,
- grafik ranking,
- grafik GRAVY vs TM segment,
- peta segmen transmembran.

## 6. Hasil dan analisis biologis, 2 menit

- Kandidat tertinggi memiliki segmen transmembran dan hidrofobisitas tinggi.
- Protein tanpa segmen TM mendapat skor rendah.
- Jelaskan bahwa hasil adalah prioritisasi awal, bukan validasi target obat final.

## 7. Kesimpulan dan pengembangan, 1 menit

- Pipeline berhasil menghubungkan sekuens protein dengan fitur biologis protein membran.
- Pengembangan: DeepTMHMM, SignalP, InterProScan, BLASTp, AlphaFold, data esensialitas gen.
