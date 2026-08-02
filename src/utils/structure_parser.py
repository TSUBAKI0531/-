# 関数（def）の形にして、外部からpdb_pathなどを自由に受け取れるようにする
from Bio import PDB

def get_biopython_structure(structure_id, pdb_path):
    """PDBファイルを読み込んでBiopythonのStructureオブジェクトを返す関数"""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(structure_id, pdb_path)
    return structure  # パースした結果を外に「返却」する