from __future__ import annotations



import pandas as pd





def evaluate_predictions(df: pd.DataFrame) -> dict:

    labeled = df[df["known_label"].isin(["membrane", "non_membrane"])].copy()

    if labeled.empty:

        return {"has_ground_truth": False}



    tp = int(((labeled["known_label"] == "membrane") & (labeled["predicted_label"] == "membrane")).sum())

    tn = int(((labeled["known_label"] == "non_membrane") & (labeled["predicted_label"] == "non_membrane")).sum())

    fp = int(((labeled["known_label"] == "non_membrane") & (labeled["predicted_label"] == "membrane")).sum())

    fn = int(((labeled["known_label"] == "membrane") & (labeled["predicted_label"] == "non_membrane")).sum())



    precision = tp / (tp + fp) if (tp + fp) else 0.0

    recall = tp / (tp + fn) if (tp + fn) else 0.0

    accuracy = (tp + tn) / len(labeled) if len(labeled) else 0.0

    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0



    return {

        "has_ground_truth": True,

        "n_labeled": int(len(labeled)),

        "tp": tp,

        "tn": tn,

        "fp": fp,

        "fn": fn,

        "accuracy": round(accuracy, 3),

        "precision": round(precision, 3),

        "recall": round(recall, 3),

        "f1": round(f1, 3),

    }





def save_evaluation_report(metrics: dict, path: str) -> None:

    with open(path, "w", encoding="utf-8") as f:

        f.write("EVALUASI BASELINE PREDIKSI MEMBRAN\n")

        f.write("===================================\n")

        if not metrics.get("has_ground_truth"):

            f.write("Tidak ada label ground-truth pada dataset. Tambahkan label=membrane atau label=non_membrane pada header FASTA.\n")

            return

        f.write(f"Jumlah protein berlabel: {metrics['n_labeled']}\n")

        f.write("Confusion matrix:\n")

        f.write(f"- TP membrane benar: {metrics['tp']}\n")

        f.write(f"- TN non-membrane benar: {metrics['tn']}\n")

        f.write(f"- FP non-membrane diprediksi membrane: {metrics['fp']}\n")

        f.write(f"- FN membrane tidak terdeteksi: {metrics['fn']}\n\n")

        f.write(f"Accuracy : {metrics['accuracy']}\n")

        f.write(f"Precision: {metrics['precision']}\n")

        f.write(f"Recall   : {metrics['recall']}\n")

        f.write(f"F1-score : {metrics['f1']}\n")

