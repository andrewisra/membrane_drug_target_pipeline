from __future__ import annotations



import re

from dataclasses import dataclass

from typing import Dict, List, Tuple



from Bio import SeqIO

from Bio.SeqUtils.ProtParam import ProteinAnalysis



try:

    from metadata import normalize_accession

except ImportError:

    from .metadata import normalize_accession



VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

HYDROPHOBIC_AA = set("AVILMFWY")

POSITIVE_AA = set("KR")

NEGATIVE_AA = set("DE")




KD_SCALE: Dict[str, float] = {

    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8,

    "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6,

    "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,

}



@dataclass(frozen=True)

class TMSegment:

    start: int

    end: int

    avg_hydropathy: float





def clean_protein_sequence(seq: str) -> str:

    """Return uppercase sequence containing only the 20 standard amino acids."""

    seq = str(seq).upper().replace("*", "")

    return "".join(aa for aa in seq if aa in VALID_AA)





def parse_field(description: str, field_name: str, default: str = "unknown") -> str:


    pattern = rf"(?:^|[\s|]){re.escape(field_name)}=([^|\s]+)"

    match = re.search(pattern, description)

    return match.group(1).strip() if match else default





def sliding_window_hydropathy(sequence: str, window: int = 19) -> List[float]:

    if len(sequence) < window:

        return []

    return [sum(KD_SCALE.get(aa, 0.0) for aa in sequence[i:i + window]) / window for i in range(len(sequence) - window + 1)]





def predict_tm_segments(

    sequence: str,

    window: int = 19,

    hydropathy_threshold: float = 1.6,

    min_hydrophobic_fraction: float = 0.55,

    merge_distance: int = 6,

    min_segment_length: int = 18,

    max_segment_length: int = 65,

) -> List[TMSegment]:

    """Transparent baseline predictor for alpha-helical transmembrane regions.

    This method detects hydrophobic windows that are long enough to span a lipid
    bilayer. It is intentionally explainable for a class project. It should be
    reported as a baseline and, when possible, compared with DeepTMHMM/TMHMM.
    """

    if len(sequence) < window:

        return []



    positive_windows: List[Tuple[int, int, float]] = []

    for i in range(len(sequence) - window + 1):

        chunk = sequence[i:i + window]

        avg_h = sum(KD_SCALE.get(aa, 0.0) for aa in chunk) / window

        hyd_frac = sum(aa in HYDROPHOBIC_AA for aa in chunk) / window

        if avg_h >= hydropathy_threshold and hyd_frac >= min_hydrophobic_fraction:

            positive_windows.append((i + 1, i + window, avg_h))



    if not positive_windows:

        return []



    merged: List[Tuple[int, int, float]] = []

    cur_start, cur_end = positive_windows[0][0], positive_windows[0][1]

    cur_scores = [positive_windows[0][2]]

    for start, end, score in positive_windows[1:]:

        if start <= cur_end + merge_distance:

            cur_end = max(cur_end, end)

            cur_scores.append(score)

        else:

            merged.append((cur_start, cur_end, sum(cur_scores) / len(cur_scores)))

            cur_start, cur_end, cur_scores = start, end, [score]

    merged.append((cur_start, cur_end, sum(cur_scores) / len(cur_scores)))



    segments: List[TMSegment] = []

    for start, end, score in merged:

        length = end - start + 1

        if min_segment_length <= length <= max_segment_length:

            segments.append(TMSegment(start, end, float(score)))

        elif length > max_segment_length:


            pos = start

            while pos + min_segment_length - 1 <= end:

                seg_end = min(pos + 22 - 1, end)

                chunk = sequence[pos - 1:seg_end]

                avg_h = sum(KD_SCALE.get(aa, 0.0) for aa in chunk) / len(chunk)

                segments.append(TMSegment(pos, seg_end, float(avg_h)))

                pos += 28

    return segments





def rough_signal_peptide_prediction(sequence: str) -> bool:

    """Simple N-terminal signal peptide heuristic.

    Signal peptides often contain a positively charged N-region followed by a
    hydrophobic H-region in the first ~30 aa. Use SignalP for a stronger result.
    """

    nterm = sequence[:35]

    if len(nterm) < 20:

        return False

    positive_n_region = sum(aa in POSITIVE_AA for aa in nterm[:10]) >= 1

    core = nterm[7:25]

    hydrophobic_core = sum(aa in HYDROPHOBIC_AA for aa in core) / max(len(core), 1) >= 0.55

    return bool(positive_n_region and hydrophobic_core)





def predict_accessibility(tm_count: int, signal_peptide: bool, description: str) -> str:

    desc = description.lower()

    if any(term in desc for term in ["outer", "surface", "secreted", "extracellular", "periplasmic"]):

        return "high"

    if signal_peptide and tm_count <= 1:

        return "medium_high"

    if tm_count >= 1:

        return "medium"

    return "low"





def extract_features_from_fasta(fasta_path: str) -> List[dict]:

    rows: List[dict] = []

    for record in SeqIO.parse(fasta_path, "fasta"):

        seq = clean_protein_sequence(str(record.seq))

        if not seq:

            continue



        analysis = ProteinAnalysis(seq)

        aa_percent = analysis.amino_acids_percent

        tm_segments = predict_tm_segments(seq)

        hydro_pct = sum(seq.count(aa) for aa in HYDROPHOBIC_AA) / len(seq) * 100

        positive_pct = sum(seq.count(aa) for aa in POSITIVE_AA) / len(seq) * 100

        negative_pct = sum(seq.count(aa) for aa in NEGATIVE_AA) / len(seq) * 100

        tm_positions = ";".join(f"{seg.start}-{seg.end}" for seg in tm_segments)

        signal = rough_signal_peptide_prediction(seq)



        rows.append({

            "protein_id": record.id,

            "accession": normalize_accession(record.id),

            "accession_norm": normalize_accession(record.id),

            "description": record.description,

            "protein_name": parse_field(record.description, "protein_name"),

            "subcellular_location": parse_field(record.description, "location"),

            "known_label": parse_field(record.description, "label"),

            "organism": parse_field(record.description, "organism"),

            "function_hint": parse_field(record.description, "function", default="unknown"),

            "length_aa": len(seq),

            "molecular_weight_da": round(analysis.molecular_weight(), 2),

            "isoelectric_point": round(analysis.isoelectric_point(), 2),

            "gravy": round(analysis.gravy(), 3),

            "instability_index": round(analysis.instability_index(), 2),

            "aromaticity": round(analysis.aromaticity(), 3),

            "hydrophobic_percent": round(hydro_pct, 2),

            "positive_percent": round(positive_pct, 2),

            "negative_percent": round(negative_pct, 2),

            "tm_segment_count": len(tm_segments),

            "tm_positions": tm_positions,

            "max_tm_hydropathy": round(max([seg.avg_hydropathy for seg in tm_segments], default=0.0), 3),

            "rough_signal_peptide": signal,

            "predicted_accessibility": predict_accessibility(len(tm_segments), signal, record.description),

            "percent_A": round(aa_percent.get("A", 0) * 100, 2),

            "percent_L": round(aa_percent.get("L", 0) * 100, 2),

            "percent_I": round(aa_percent.get("I", 0) * 100, 2),

            "percent_V": round(aa_percent.get("V", 0) * 100, 2),

        })

    return rows

