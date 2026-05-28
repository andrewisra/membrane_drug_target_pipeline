# Instruksi Pengumpulan Data untuk Proyek Protein Membran

Judul proyek: **Identifikasi Kandidat Target Obat Berbasis Protein Membran melalui Analisis Sekuens dan Prediksi Topologi Transmembran**

Dokumen ini dapat dikirim langsung ke anggota kelompok yang bertugas mencari data.

## 1. Tujuan data

Data yang dibutuhkan bukan sekadar daftar nama protein. Kita membutuhkan dataset yang memungkinkan pipeline untuk:

1. membaca sekuens protein,
2. mengetahui label pembanding: `membrane` atau `non_membrane`,
3. mengetahui konteks biologis protein,
4. mengevaluasi apakah prediksi pipeline masuk akal,
5. memprioritaskan kandidat target obat secara lebih bertanggung jawab.

## 2. Studi kasus yang disarankan

Pilih **satu organisme target** agar laporan fokus. Rekomendasi paling aman:

- *Escherichia coli* K-12, organism_id UniProt: `83333`
- *Staphylococcus aureus*, organism_id UniProt: `1280`
- *Mycobacterium tuberculosis*, organism_id UniProt: `1773`

Untuk tugas kuliah, *E. coli* K-12 paling mudah karena data anotasinya banyak. Jika ingin konteks target obat lebih kuat, gunakan bakteri patogen seperti *Staphylococcus aureus* atau *Mycobacterium tuberculosis*.

## 3. Jumlah data minimal

Target minimal:

- 20 protein berlabel `membrane`
- 20 protein berlabel `non_membrane`

Target ideal:

- 30 protein `membrane`
- 30 protein `non_membrane`

Jangan terlalu banyak dulu. Kualitas anotasi lebih penting daripada jumlah.

## 4. Sumber data utama

Gunakan **UniProtKB**, terutama entri reviewed/Swiss-Prot jika tersedia. UniProt menyediakan REST API dan data dapat diunduh dalam format FASTA serta TSV, sehingga cocok untuk pipeline Python.

## 5. Kriteria protein membran

Protein dimasukkan ke label `membrane` jika salah satu kondisi berikut ada di UniProt:

- subcellular location mengandung `membrane`, `outer membrane`, `inner membrane`, `plasma membrane`, atau `cell membrane`,
- fitur UniProt memiliki transmembrane region,
- keyword mengandung `Transmembrane` atau `Transmembrane helix`,
- nama/fungsi protein jelas sebagai transporter, channel, receptor, porin, pump, atau membrane protein.

## 6. Kriteria protein non-membran / kontrol negatif

Protein dimasukkan ke label `non_membrane` jika:

- subcellular location adalah `cytoplasm`, `cytosol`, atau lokasi non-membran lain,
- tidak memiliki fitur transmembrane region,
- bukan protein sekretori/permukaan,
- berfungsi sebagai enzim sitosolik, chaperone sitosolik, protein ribosomal, atau protein metabolisme larut.

Hindari protein yang anotasinya ambigu.

## 7. Format file yang harus dikumpulkan

Kumpulkan dua file utama.

### A. File metadata CSV

Nama file yang disarankan:

```text
data/friend_metadata.csv
```

Kolom wajib:

```text
accession,label,organism,protein_name,function_hint,subcellular_location,is_essential,human_homolog,source_url,notes
```

Penjelasan kolom:

| Kolom | Wajib? | Isi |
|---|---:|---|
| `accession` | ya | UniProt accession, contoh `P0A910` |
| `label` | ya | `membrane` atau `non_membrane` |
| `organism` | ya | nama organisme |
| `protein_name` | ya | nama protein dari UniProt |
| `function_hint` | ya | ringkasan fungsi: transporter, porin, enzyme, receptor, dll. |
| `subcellular_location` | ya | lokasi seluler dari UniProt |
| `is_essential` | opsional | `yes`, `no`, atau `unknown` |
| `human_homolog` | opsional | `yes`, `no`, atau `unknown`; untuk bakteri, `yes` berarti berisiko off-target |
| `source_url` | ya | URL UniProt entry |
| `notes` | opsional | catatan singkat kenapa label dipilih |

Contoh:

```csv
accession,label,organism,protein_name,function_hint,subcellular_location,is_essential,human_homolog,source_url,notes
P0A910,membrane,Escherichia coli K-12,Outer membrane protein A,outer membrane porin/structural protein,outer membrane,yes,no,https://www.uniprot.org/uniprotkb/P0A910/entry,protein membran luar
P0A6F5,non_membrane,Escherichia coli K-12,60 kDa chaperonin,cytosolic chaperone,cytoplasm,yes,yes,https://www.uniprot.org/uniprotkb/P0A6F5/entry,kontrol non-membran
```

### B. File FASTA protein

Nama file yang disarankan:

```text
data/friend_proteins.fasta
```

File ini dapat dibuat otomatis dari `friend_metadata.csv` dengan script:

```bash
python src/fetch_uniprot.py fetch-accessions --input-csv data/friend_metadata.csv --output-fasta data/friend_proteins.fasta
```

Atau diunduh manual dari UniProt dalam format FASTA. Pastikan accession pada FASTA cocok dengan kolom `accession` di CSV.

## 8. Cara mencari protein di UniProt

Contoh query protein membran *E. coli* K-12:

```text
(organism_id:83333) AND (cc_subcellular_location:membrane) AND reviewed:true
```

Contoh query protein non-membran/sitoplasmik *E. coli* K-12:

```text
(organism_id:83333) AND (cc_subcellular_location:cytoplasm) AND reviewed:true
```

Contoh memakai script proyek:

```bash
python src/fetch_uniprot.py search --query "(organism_id:83333) AND (cc_subcellular_location:membrane) AND reviewed:true" --output-tsv data/search_membrane.tsv --size 30
python src/fetch_uniprot.py search --query "(organism_id:83333) AND (cc_subcellular_location:cytoplasm) AND reviewed:true" --output-tsv data/search_non_membrane.tsv --size 30
```

## 9. Opsional: hasil DeepTMHMM/TMHMM

Jika ada waktu, jalankan FASTA pada DeepTMHMM atau TMHMM, lalu normalisasi hasil ke CSV:

```text
data/external_tm_predictions.csv
```

Format:

```csv
protein_id,region_type,start,end,source
sp|P0A910|OMPA_ECOLI,tm,22,44,DeepTMHMM
sp|P0A910|OMPA_ECOLI,tm,60,82,DeepTMHMM
```

Kolom:

- `protein_id`: harus cocok dengan ID FASTA yang dibaca pipeline,
- `region_type`: isi `tm`, `tmhelix`, `transmembrane`, atau `beta_barrel`,
- `start`: posisi awal segmen,
- `end`: posisi akhir segmen,
- `source`: misalnya `DeepTMHMM`.

Jalankan pipeline dengan hasil eksternal:

```bash
python src/main.py --fasta data/friend_proteins.fasta --metadata-csv data/friend_metadata.csv --external-tm-csv data/external_tm_predictions.csv --output-dir results_friend --figure-dir figures_friend
```

## 10. Checklist kualitas data

Sebelum dikirim ke tim coding, pastikan:

- minimal 40 protein total,
- jumlah membrane dan non_membrane seimbang,
- semua accession valid di UniProt,
- semua protein berasal dari organisme yang sama,
- label tidak ambigu,
- FASTA dan CSV memiliki accession yang sama,
- data berasal dari sumber yang bisa dicantumkan di laporan.

## 11. Yang tidak perlu dikumpulkan

Tidak perlu mencari struktur 3D PDB, docking, atau data senyawa obat. Proyek ini adalah penyaringan awal target berbasis sekuens, bukan simulasi docking obat.
