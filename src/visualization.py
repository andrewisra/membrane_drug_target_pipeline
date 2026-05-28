from __future__ import annotations



import os

import pandas as pd

import matplotlib.pyplot as plt





def _save(fig, path: str) -> None:

    fig.tight_layout()

    fig.savefig(path, dpi=200)

    plt.close(fig)





def plot_candidate_ranking(df: pd.DataFrame, figure_dir: str, top_n: int = 10) -> None:

    top = df.head(top_n).iloc[::-1]

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.barh(top["protein_id"], top["candidate_score"])

    ax.set_xlabel("Skor kandidat target obat")

    ax.set_ylabel("Protein")

    ax.set_title("Ranking Kandidat Protein Membran")

    ax.set_xlim(0, 100)

    _save(fig, os.path.join(figure_dir, "candidate_ranking.png"))





def plot_gravy_vs_tm(df: pd.DataFrame, figure_dir: str) -> None:

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(df["gravy"], df["tm_segment_count"], s=df["candidate_score"].clip(lower=10) * 3, alpha=0.7)

    for _, row in df.iterrows():

        ax.annotate(row["protein_id"], (row["gravy"], row["tm_segment_count"]), fontsize=8, alpha=0.8)

    ax.set_xlabel("GRAVY / rata-rata hidrofobisitas")

    ax.set_ylabel("Jumlah segmen transmembran")

    ax.set_title("Hubungan Hidrofobisitas dan Prediksi Segmen Transmembran")

    _save(fig, os.path.join(figure_dir, "gravy_vs_tm_segments.png"))





def plot_class_distribution(df: pd.DataFrame, figure_dir: str) -> None:

    counts = df["priority_class"].value_counts().reindex(["prioritas_tinggi", "prioritas_sedang", "prioritas_rendah"]).fillna(0)

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.bar(counts.index, counts.values)

    ax.set_xlabel("Kelas prioritas")

    ax.set_ylabel("Jumlah protein")

    ax.set_title("Distribusi Prioritas Kandidat")

    _save(fig, os.path.join(figure_dir, "priority_class_distribution.png"))





def plot_tm_segment_map(df: pd.DataFrame, figure_dir: str, top_n: int = 8) -> None:

    sub = df.head(top_n).copy()

    fig, ax = plt.subplots(figsize=(10, max(4, len(sub) * 0.6)))

    y_positions = range(len(sub))

    for y, (_, row) in zip(y_positions, sub.iterrows()):

        length = int(row["length_aa"])

        ax.hlines(y, 1, length, linewidth=2)

        positions = str(row.get("tm_positions", ""))

        if positions and positions != "nan":

            for seg in positions.split(";"):

                if "-" not in seg:

                    continue

                start, end = map(int, seg.split("-"))

                ax.hlines(y, start, end, linewidth=8)

        ax.text(length + 5, y, f"score={row['candidate_score']}", va="center", fontsize=8)

    ax.set_yticks(list(y_positions))

    ax.set_yticklabels(sub["protein_id"])

    ax.set_xlabel("Posisi asam amino")

    ax.set_title("Peta Segmen Transmembran pada Kandidat Teratas")

    _save(fig, os.path.join(figure_dir, "tm_segment_map.png"))





def plot_confusion_matrix(df: pd.DataFrame, figure_dir: str) -> None:

    labeled = df[df["known_label"].isin(["membrane", "non_membrane"])]

    if labeled.empty:

        return

    labels = ["membrane", "non_membrane"]

    matrix = pd.crosstab(labeled["known_label"], labeled["predicted_label"]).reindex(index=labels, columns=labels).fillna(0)

    fig, ax = plt.subplots(figsize=(5, 4))

    ax.imshow(matrix.values)

    ax.set_xticks(range(len(labels)))

    ax.set_xticklabels(labels, rotation=30, ha="right")

    ax.set_yticks(range(len(labels)))

    ax.set_yticklabels(labels)

    ax.set_xlabel("Prediksi")

    ax.set_ylabel("Label dataset")

    ax.set_title("Confusion Matrix Baseline")

    for i in range(len(labels)):

        for j in range(len(labels)):

            ax.text(j, i, int(matrix.values[i, j]), ha="center", va="center")

    _save(fig, os.path.join(figure_dir, "confusion_matrix.png"))

