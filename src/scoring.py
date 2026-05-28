from __future__ import annotations



import pandas as pd



ACCESSIBILITY_POINTS = {

    "high": 20,

    "medium_high": 15,

    "medium": 10,

    "low": 0,

}



FUNCTION_POINTS_KEYWORDS = {

    "transporter": 12,

    "receptor": 12,

    "channel": 12,

    "pump": 10,

    "efflux": 10,

    "adhesin": 10,

    "porin": 10,

    "outer": 8,

    "surface": 8,

    "enzyme": 4,

}





def _clip(value: float, low: float, high: float) -> float:

    return max(low, min(high, value))





def _truthy(value: object) -> bool:

    return str(value).strip().lower() in {"true", "yes", "y", "1", "essential", "ada"}





def _negative(value: object) -> bool:

    return str(value).strip().lower() in {"true", "yes", "y", "1", "ada", "detected"}





def _function_score(text: str) -> int:

    lower = str(text).lower()

    score = 0

    for keyword, points in FUNCTION_POINTS_KEYWORDS.items():

        if keyword in lower:

            score = max(score, points)

    return score





def score_candidates(features: pd.DataFrame) -> pd.DataFrame:

    df = features.copy()




    df["tm_evidence_score"] = df["tm_segment_count"].apply(lambda n: int(_clip(n * 16, 0, 40)))




    df["hydrophobicity_score"] = df.apply(

        lambda r: int(_clip((r["hydrophobic_percent"] - 30) * 1.2 + max(r["gravy"], 0) * 10, 0, 20)), axis=1

    )




    df["accessibility_score"] = df["predicted_accessibility"].map(ACCESSIBILITY_POINTS).fillna(0).astype(int)




    df["signal_score"] = df["rough_signal_peptide"].apply(lambda v: 8 if bool(v) else 0)




    df["function_score"] = df.apply(

        lambda r: _function_score(f"{r.get('description','')} {r.get('function_hint','')}"), axis=1

    )




    if "is_essential" in df.columns:

        df["essentiality_score"] = df["is_essential"].apply(lambda v: 8 if _truthy(v) else 0)

    else:

        df["essentiality_score"] = 0




    if "human_homolog" in df.columns:

        df["human_homolog_penalty"] = df["human_homolog"].apply(lambda v: 8 if _negative(v) else 0)

    else:

        df["human_homolog_penalty"] = 0




    df["length_penalty"] = df["length_aa"].apply(lambda x: 8 if x < 80 or x > 1500 else 0)

    df["instability_penalty"] = df["instability_index"].apply(lambda x: 6 if x > 60 else (3 if x > 45 else 0))

    df["no_tm_penalty"] = df["tm_segment_count"].apply(lambda n: 20 if n == 0 else 0)



    raw = (

        df["tm_evidence_score"]

        + df["hydrophobicity_score"]

        + df["accessibility_score"]

        + df["signal_score"]

        + df["function_score"]

        + df["essentiality_score"]

        - df["human_homolog_penalty"]

        - df["length_penalty"]

        - df["instability_penalty"]

        - df["no_tm_penalty"]

    )

    df["candidate_score"] = raw.apply(lambda x: int(_clip(round(x), 0, 100)))



    def classify(score: int) -> str:

        if score >= 70:

            return "prioritas_tinggi"

        if score >= 45:

            return "prioritas_sedang"

        return "prioritas_rendah"



    df["priority_class"] = df["candidate_score"].apply(classify)

    df["predicted_label"] = df["tm_segment_count"].apply(lambda n: "membrane" if n >= 1 else "non_membrane")



    def interpret(row: pd.Series) -> str:

        reasons = []

        if row["tm_segment_count"] >= 1:

            reasons.append(f"memiliki {row['tm_segment_count']} segmen transmembran")

        else:

            reasons.append("tidak memiliki segmen transmembran yang terdeteksi")

        if row["gravy"] > 0:

            reasons.append("GRAVY positif/hidrofobik")

        if row["rough_signal_peptide"]:

            reasons.append("memiliki indikasi signal peptide")

        if "is_essential" in row and _truthy(row.get("is_essential")):

            reasons.append("memiliki indikasi esensialitas gen")

        if "human_homolog" in row and _negative(row.get("human_homolog")):

            reasons.append("memiliki homolog manusia sehingga perlu evaluasi off-target")

        if row["predicted_accessibility"] in ["high", "medium_high"]:

            reasons.append(f"aksesibilitas diprediksi {row['predicted_accessibility']}")

        if row["priority_class"] == "prioritas_tinggi":

            verdict = "Kandidat kuat untuk dianalisis lanjut."

        elif row["priority_class"] == "prioritas_sedang":

            verdict = "Kandidat moderat, memerlukan validasi tambahan."

        else:

            verdict = "Bukan prioritas utama berdasarkan fitur sekuens."

        return verdict + " Alasan: " + ", ".join(reasons) + "."



    df["interpretation"] = df.apply(interpret, axis=1)

    return df.sort_values(["candidate_score", "tm_segment_count", "hydrophobic_percent"], ascending=[False, False, False]).reset_index(drop=True)

