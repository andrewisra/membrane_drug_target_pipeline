from __future__ import annotations



import pandas as pd





def load_deeptmhmm_regions(path: str) -> pd.DataFrame:

    """Parse a simple DeepTMHMM-like TSV exported manually.

    Expected columns:
    protein_id, region_type, start, end

    Keep this parser simple because different web-service exports may vary. The
    project can merge external predictions once students normalize the output to
    this TSV format.
    """

    df = pd.read_csv(path)

    required = {"protein_id", "region_type", "start", "end"}

    missing = required - set(df.columns)

    if missing:

        raise ValueError(f"Missing columns in external prediction file: {sorted(missing)}")

    return df





def summarize_external_tm(path: str) -> pd.DataFrame:

    df = load_deeptmhmm_regions(path)

    tm = df[df["region_type"].str.lower().isin(["tmhelix", "transmembrane", "tm", "beta_barrel"])]

    summary = tm.groupby("protein_id").agg(

        external_tm_count=("region_type", "count"),

        external_tm_positions=("start", lambda s: ";".join(f"{int(a)}-{int(b)}" for a, b in zip(tm.loc[s.index, "start"], tm.loc[s.index, "end"])))

    ).reset_index()

    return summary

