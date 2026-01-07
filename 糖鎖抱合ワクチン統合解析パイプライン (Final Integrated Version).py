import json
import os
import pandas as pd
import freesasa
from rdkit import Chem
from Bio.PDB import MMCIFParser, NeighborSearch
import py3Dmol

class GlycoConjugateWorkflow:
    """
    糖鎖抱合ワクチンのAF3準備、一括解析、およびランキングを行う統合クラス
    """
    def __init__(self, job_name):
        self.job_name = job_name
        self.results_df = None

    # --- 1. 準備フェーズ: インデックス特定とJSON作成 ---

    def _get_af3_atom_name(self, mol, target_idx):
        """RDKitのインデックスをAF3形式(例: C12)に変換"""
        atom = mol.GetAtomWithIdx(target_idx)
        symbol = atom.GetSymbol()
        count = 0
        for i, a in enumerate(mol.GetAtoms()):
            if a.GetSymbol() == symbol:
                count += 1
            if i == target_idx:
                return f"{symbol}{count}"
        return None

    def find_terminal_atom(self, smiles, smarts_pattern):
        """SMARTSパターンを用いて結合原子を自動特定"""
        mol = Chem.MolFromSmiles(smiles)
        pattern = Chem.MolFromSmarts(smarts_pattern)
        matches = mol.GetSubstructMatches(pattern)
        if not matches:
            raise ValueError(f"パターン '{smarts_pattern}' が見つかりません。")
        
        target_idx = matches[0][0] # パターンの最初の原子を選択
        return self._get_af3_atom_name(mol, target_idx)

    def prepare_af3_input(self, protein_seq, taca_linker_smiles, bond_res_idx, terminal_smarts, bond_atom_protein="NZ"):
        """AF3実行用のJSONファイルを生成"""
        bond_atom_ligand = self.find_terminal_atom(taca_linker_smiles, terminal_smarts)
        
        data = {
            "name": self.job_name,
            "modelSeeds": [1],
            "sequences": [
                {"protein": {"id": "A", "sequence": protein_seq}},
                {"ligand": {"id": "B", "smiles": taca_linker_smiles}}
            ],
            "bondedAtomPairs": [
                {
                    "at1": {"resChainId": "A", "resIdx": bond_res_idx, "atomName": bond_atom_protein},
                    "at2": {"resChainId": "B", "resIdx": 1, "atomName": bond_atom_ligand}
                }
            ]
        }
        
        output_path = f"{self.job_name}.json"
        with open(output_path, "w") as f:
            json.dump(data, f, indent=4)
        print(f"✅ AF3用JSONを作成完了: {output_path} (結合原子: {bond_atom_ligand})")
        return output_path

    # --- 2. 解析フェーズ: SASAおよび相互作用解析 ---

    def analyze_interactions(self, cif_path, distance_cutoff=5.0):
        """指定されたCIFファイルから接触残基を特定"""
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("model", cif_path)
        model = structure[0]
        
        protein_atoms = [a for a in model.get_atoms() if a.get_parent().get_resname() not in ["HOH"]]
        # 糖鎖(Ligand)はAF3では通常ヘテロ原子として扱われる
        sugar_atoms = [a for a in model.get_atoms() if "H_" in a.get_parent().get_full_id()[3][0]]
        
        ns = NeighborSearch(protein_atoms)
        contact_residues = set()
        for s_atom in sugar_atoms:
            contacts = ns.search(s_atom.coord, distance_cutoff)
            for c_atom in contacts:
                res = c_atom.get_parent()
                contact_residues.add(f"{res.get_resname()}{res.id[1]}")
        return sorted(list(contact_residues))

    def _calculate_sasa(self, cif_path):
        """FreeSASAを用いて露出面積を計算"""
        structure = freesasa.Structure(cif_path)
        result = freesasa.calc(structure)
        
        # Chain Aをタンパク、それ以外を糖鎖と定義
        s_prot = freesasa.selectObjects(["prot, chain A"], structure, result)
        s_glyc = freesasa.selectObjects(["glyc, not chain A"], structure, result)
        
        return {
            "protein_sasa": s_prot["prot"],
            "glycan_sasa": s_glyc["glyc"]
        }

    def batch_analyze_models(self, num_models=5):
        """全モデルを一括解析し、抗原提示の質でランキング"""
        analysis_data = []
        
        for i in range(num_models):
            file_path = f"{self.job_name}_model_{i}.cif"
            if not os.path.exists(file_path):
                continue
            
            # SASA計算
            sasa = self._calculate_sasa(file_path)
            # 接触残基
            contacts = self.analyze_interactions(file_path)
            
            # スコアリング: (糖鎖露出面積) / (接触残基数 + 1)
            # 値が高いほど「タンパク質に埋もれず、外側を向いている」と評価
            exposure_score = sasa["glycan_sasa"] / (len(contacts) + 1)
            
            analysis_data.append({
                "Model_Index": i,
                "Glycan_SASA": round(sasa["glycan_sasa"], 2),
                "Contact_Res_Count": len(contacts),
                "Exposure_Score": round(exposure_score, 2),
                "File": file_path
            })
        
        self.results_df = pd.DataFrame(analysis_data).sort_values(by="Exposure_Score", ascending=False)
        print("\n🏆 モデル解析ランキング (Exposure_Score順):")
        print(self.results_df.to_string(index=False))
        return self.results_df

    # --- 3. 可視化フェーズ ---

    def visualize_model(self, model_idx):
        """指定したモデル番号の構造を表示"""
        file_path = f"{self.job_name}_model_{model_idx}.cif"
        if not os.path.exists(file_path):
            print("ファイルが見つかりません。")
            return
            
        view = py3Dmol.view(width=800, height=600)
        with open(file_path, 'r') as f:
            view.addModel(f.read(), 'mcif')
        
        view.setStyle({'cartoon': {'color': 'spectrum'}}) # タンパクは信頼度
        view.setStyle({'hetflag': True}, {'stick': {'colorscheme': 'magentaCarbon'}}) # 糖鎖はマゼンタ
        view.zoomTo({'hetflag': True})
        print(f"Visualizing Model {model_idx}...")
        return view.show()

# --- 実際の運用例 ---

# 1. 初期化
wf = GlycoConjugateWorkflow("TACA_Vaccine_Project")

# 2. AF3の準備 (SMILESから自動インデックス特定)
protein_seq = "MKTIIALSYIFCLVFA..." # 実際の配列
taca_linker = "CC(=O)N[C@@H]1[C@@H](O)O[C@H](CO)[C@H](O)[C@@H]1OC[C@H](NC=O)C(=O)NCCCC" 
smarts_end = "C(=O)N" # リンカーの末端（タンパク質との結合点）を指定

wf.prepare_af3_input(
    protein_seq=protein_seq,
    taca_linker_smiles=taca_linker,
    bond_res_idx=105, # Lys105に結合
    terminal_smarts=smarts_end
)

# --- ここでAlphaFold Serverなどで計算を実行 ---

# 3. 計算結果の一括解析とランキング (結果ファイルがある前提)
# df = wf.batch_analyze_models(num_models=5)

# 4. 最良モデルの可視化
# if df is not None:
#     best_idx = df.iloc[0]['Model_Index']
#     wf.visualize_model(best_idx)