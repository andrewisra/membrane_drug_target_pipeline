# Flowchart Pipeline

```text
Mulai
  |
  v
Input FASTA protein
  |
  v
Membersihkan sekuens protein
  |
  v
Ekstraksi fitur fisikokimia
(length, MW, pI, GRAVY, instability index, komposisi AA)
  |
  v
Prediksi segmen transmembran baseline
(sliding window Kyte-Doolittle)
  |
  v
Prediksi kasar signal peptide dan aksesibilitas
  |
  v
Hitung skor kandidat target obat
  |
  v
Ranking protein kandidat
  |
  +--> Evaluasi terhadap label dataset
  |
  +--> Visualisasi hasil
  |
  v
Interpretasi biologis dan rekomendasi kandidat
  |
  v
Selesai
```

## Komponen metode

- Input: file FASTA protein.
- Komputasi inti: ekstraksi fitur + prediksi segmen TM + scoring.
- Output: ranking kandidat, grafik, evaluasi, ringkasan biologis.
