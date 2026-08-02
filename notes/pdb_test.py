#from Bio.PDB import PDBList

# PDBListのインスタンスを作成
#pdbl = PDBList()

# '1FAT' の構造データを取得し、現在のフォルダに保存
#pdbl.retrieve_pdb_file('1FAT', pdir='.', file_format='pdb')

from Bio.PDB.PDBParser import PDBParser

# 1. パーサーの作成
parser = PDBParser()

# 2. PDBファイルの読み込み
# ※コードAでダウンロードされた実際のファイル名（例: "pdb1fat.ent"）に書き換えてください
structure = parser.get_structure("1FAT_sample", "pdb1fat.ent")

# 3. Model 0 の Chain B を取得（ショートカット利用）
chain_b = structure[0]["B"]

# 4. Chain B 内の全残基をループで回し、残基名を出力
print("--- Chain B の残基名一覧 ---")
for residue in chain_b:
    print(residue.get_resname())