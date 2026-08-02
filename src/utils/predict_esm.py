# src/utils/predict_esm.py
def generate_peptide_pdb(sequence, output_filename="esmfold_predicted.pdb"):
    """配列からエミュレートされたPDBファイルを生成する関数"""
    # 先週のダミーPDBテキストを生成するロジック
    dummy_pdb = """ATOM      1  N   ILE A   1       0.000   0.000   0.000  1.00 70.00           N
ATOM      2  CA  ILE A   1       1.450   0.000   0.000  1.00 71.20           C
..."""
    with open(output_filename, "w", encoding="utf-8") as f:
        f.write(dummy_pdb)
    print(f"📁 {output_filename} を生成しました。")
    return output_filename