from src.feature_extraction import clean_protein_sequence, predict_tm_segments, rough_signal_peptide_prediction





def test_clean_protein_sequence_removes_invalid_symbols():

    assert clean_protein_sequence("MAAAXXB*") == "MAAA"





def test_tm_predictor_detects_hydrophobic_helix():

    seq = "MKK" + "A" * 5 + "LIVVLVVIGLLVGLVLAAV" + "DEKQQQ"

    assert len(predict_tm_segments(seq)) >= 1





def test_signal_peptide_heuristic_boolean():

    seq = "MKKLLLAGVAVALAVSQAADNTLTQK" + "A" * 30

    assert rough_signal_peptide_prediction(seq) in [True, False]

