from __future__ import annotations



import re

from pathlib import Path

from typing import Optional



import pandas as pd





def normalize_accession(value: object) -> str:

    """Normalize common FASTA/UniProt identifiers to a comparable accession key."""

    if value is None:

        return ""

    text = str(value).strip()

    if not text or text.lower() == "nan":

        return ""


    if "|" in text:

        parts = text.split("|")

        if len(parts) >= 2 and parts[1]:

            return parts[1].strip().upper()


    token = text.split()[0]

    token = re.sub(r"\.\d+$", "", token)

    return token.upper()





def load_metadata(metadata_csv: Optional[str]) -> pd.DataFrame:

    """Load optional protein metadata collected by the data team.

    Recommended columns:
    accession,label,organism,protein_name,function_hint,subcellular_location,
    is_essential,human_homolog,source_url,notes

    `accession` is required. Unknown columns are preserved and can be used in the report.
    """

    if not metadata_csv:

        return pd.DataFrame()

    path = Path(metadata_csv)

    if not path.exists():

        raise FileNotFoundError(f"Metadata CSV tidak ditemukan: {metadata_csv}")

    df = pd.read_csv(path)

    if "accession" in df.columns:

        id_col = "accession"

    elif "protein_id" in df.columns:

        id_col = "protein_id"

        df["accession"] = df["protein_id"]

    else:

        raise ValueError("Metadata CSV wajib memiliki kolom 'accession' atau 'protein_id'.")

    df["accession_norm"] = df[id_col].apply(normalize_accession)

    df = df[df["accession_norm"] != ""].copy()

    df = df.drop_duplicates(subset=["accession_norm"], keep="first")

    return df





def merge_metadata(features: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:

    if metadata.empty:

        return features

    df = features.copy()

    if "accession_norm" not in df.columns:

        df["accession_norm"] = df["accession"].apply(normalize_accession)



    meta = metadata.copy()


    merged = df.merge(meta, on="accession_norm", how="left", suffixes=("", "_meta"))




    preference_pairs = [

        ("label", "known_label"),

        ("known_label", "known_label"),

        ("function_note", "function_hint"),

        ("organism", "organism"),

        ("function_hint", "function_hint"),

        ("protein_name", "protein_name"),

        ("subcellular_location", "subcellular_location"),

        ("is_essential", "is_essential"),

        ("human_homolog", "human_homolog"),

        ("source_url", "source_url"),

        ("notes", "notes"),

    ]

    for meta_col, target_col in preference_pairs:

        alt = f"{meta_col}_meta"

        if meta_col in merged.columns and target_col in merged.columns:


            continue

        if alt in merged.columns:

            if target_col not in merged.columns:

                merged[target_col] = "unknown"

            merged[target_col] = merged[alt].where(merged[alt].notna() & (merged[alt].astype(str).str.strip() != ""), merged[target_col])

        elif meta_col in merged.columns and target_col not in merged.columns:

            merged[target_col] = merged[meta_col]



    if "label" in merged.columns:

        merged["known_label"] = merged["label"].where(merged["label"].notna() & (merged["label"].astype(str).str.strip() != ""), merged.get("known_label", "unknown"))



    if "known_label" in merged.columns:

        merged["known_label"] = merged["known_label"].astype(str).str.strip().str.lower().replace({

            "membran": "membrane",

            "membrane protein": "membrane",

            "non membrane": "non_membrane",

            "non-membrane": "non_membrane",

            "cytosolic": "non_membrane",

            "soluble": "non_membrane",

        })

    return merged

