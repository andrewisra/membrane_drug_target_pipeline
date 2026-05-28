from __future__ import annotations

import argparse

import os

import pandas as pd



try:

    from evaluation import evaluate_predictions, save_evaluation_report

    from feature_extraction import extract_features_from_fasta

    from external_predictions import summarize_external_tm

    from metadata import load_metadata, merge_metadata

    from scoring import score_candidates

    from visualization import (

    plot_candidate_ranking,

    plot_gravy_vs_tm,

    plot_class_distribution,

    plot_tm_segment_map,

        plot_confusion_matrix,

    )

except ImportError:

    from .evaluation import evaluate_predictions, save_evaluation_report

    from .feature_extraction import extract_features_from_fasta

    from .external_predictions import summarize_external_tm

    from .metadata import load_metadata, merge_metadata

    from .scoring import score_candidates

    from .visualization import (

        plot_candidate_ranking,

        plot_gravy_vs_tm,

        plot_class_distribution,

        plot_tm_segment_map,

        plot_confusion_matrix,

    )





def write_summary(df: pd.DataFrame, metrics: dict, output_path: str) -> None:

    total = len(df)

    membrane_like = int((df["tm_segment_count"] >= 1).sum())

    top = df.iloc[0] if total else None

    with open(output_path, "w", encoding="utf-8") as f:

        f.write("RINGKASAN HASIL ANALISIS\n")

        f.write("=======================\n")

        f.write(f"Total protein dianalisis: {total}\n")

        f.write(f"Protein dengan >=1 segmen transmembran: {membrane_like}\n")

        f.write(f"Protein tanpa segmen transmembran: {total - membrane_like}\n")

        f.write(f"Kandidat prioritas tinggi: {int((df['priority_class'] == 'prioritas_tinggi').sum())}\n")

        f.write(f"Kandidat prioritas sedang: {int((df['priority_class'] == 'prioritas_sedang').sum())}\n")

        f.write(f"Kandidat prioritas rendah: {int((df['priority_class'] == 'prioritas_rendah').sum())}\n\n")

        if metrics.get("has_ground_truth"):

            f.write("Evaluasi baseline terhadap label dataset:\n")

            f.write(f"- Accuracy : {metrics['accuracy']}\n")

            f.write(f"- Precision: {metrics['precision']}\n")

            f.write(f"- Recall   : {metrics['recall']}\n")

            f.write(f"- F1-score : {metrics['f1']}\n\n")

        if top is not None:

            f.write("Kandidat prioritas tertinggi:\n")

            f.write(f"- ID protein: {top['protein_id']}\n")

            f.write(f"- Skor kandidat: {top['candidate_score']} / 100\n")

            f.write(f"- Jumlah segmen transmembran: {top['tm_segment_count']}\n")

            f.write(f"- Posisi segmen TM: {top['tm_positions'] if top['tm_positions'] else '-'}\n")

            f.write(f"- GRAVY: {top['gravy']}\n")

            f.write(f"- Instability index: {top['instability_index']}\n")

            f.write(f"- Aksesibilitas prediksi: {top['predicted_accessibility']}\n")

            f.write(f"- Interpretasi: {top['interpretation']}\n")





def run_pipeline(fasta_path: str, output_dir: str, figure_dir: str, metadata_csv: str | None = None, external_tm_csv: str | None = None) -> None:

    os.makedirs(output_dir, exist_ok=True)

    os.makedirs(figure_dir, exist_ok=True)



    rows = extract_features_from_fasta(fasta_path)

    if not rows:

        raise ValueError("Tidak ada sekuens protein valid pada file FASTA.")



    features = pd.DataFrame(rows)

    metadata = load_metadata(metadata_csv)

    features = merge_metadata(features, metadata)




    if external_tm_csv:

        external_tm = summarize_external_tm(external_tm_csv)

        features = features.merge(external_tm, on="protein_id", how="left")


        features["external_tm_count"] = features["external_tm_count"].fillna(0).astype(int)

        features["external_tm_positions"] = features["external_tm_positions"].fillna("")



    ranked = score_candidates(features)

    metrics = evaluate_predictions(ranked)



    features_path = os.path.join(output_dir, "protein_features.csv")

    ranking_path = os.path.join(output_dir, "candidate_ranking.csv")

    summary_path = os.path.join(output_dir, "summary_report.txt")

    eval_path = os.path.join(output_dir, "evaluation_report.txt")



    features.to_csv(features_path, index=False)

    ranked.to_csv(ranking_path, index=False)

    write_summary(ranked, metrics, summary_path)

    save_evaluation_report(metrics, eval_path)



    plot_candidate_ranking(ranked, figure_dir)

    plot_gravy_vs_tm(ranked, figure_dir)

    plot_class_distribution(ranked, figure_dir)

    plot_tm_segment_map(ranked, figure_dir)

    plot_confusion_matrix(ranked, figure_dir)



    print("Pipeline selesai.")

    print(f"Fitur protein: {features_path}")

    print(f"Ranking kandidat: {ranking_path}")

    print(f"Ringkasan: {summary_path}")

    print(f"Evaluasi: {eval_path}")

    print(f"Grafik: {figure_dir}")





if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Pipeline identifikasi kandidat target obat berbasis protein membran.")

    parser.add_argument("--fasta", default="data/proteins.fasta", help="Path ke file FASTA protein.")

    parser.add_argument("--output-dir", default="results", help="Folder output CSV dan ringkasan.")

    parser.add_argument("--figure-dir", default="figures", help="Folder output visualisasi.")

    parser.add_argument("--metadata-csv", default=None, help="Opsional: CSV metadata dari tim pengumpul data.")

    parser.add_argument("--external-tm-csv", default=None, help="Opsional: CSV hasil DeepTMHMM/TMHMM yang sudah dinormalisasi.")

    args = parser.parse_args()

    run_pipeline(args.fasta, args.output_dir, args.figure_dir, args.metadata_csv, args.external_tm_csv)

