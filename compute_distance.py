import os

from Bio import PDB

DISTANCE_THRESHOLD_ANGSTROM = 5.0

def is_within_distance_threshold(
    distance: float,
    threshold: float = DISTANCE_THRESHOLD_ANGSTROM,
) -> bool:
    """距離が閾値以内か判定する。"""
    return distance <= threshold

def get_biopython_structure(structure_id, pdb_path):
    parser = PDB.PDBParser(QUIET=True)
    return parser.get_structure(structure_id, pdb_path)

def main():
    print("\n==================================================")
    print("=== HLAアレル多型対応・親和性スコアリング初期型 ===")
    print("==================================================")
    
    mhc_pdb_path = "pdb1a1m.ent"
    pep_pdb_path = "esmfold_predicted.pdb"
    
    # 1. アレルデータベースの定義（Wet知見のパラメータ化）
    ALLELE_DATABASE = {
        "A0201": ["LEU", "MET"],
        "B4402": ["GLU"],
        "B2705": ["ARG"]
    }
    
    # ★ ここを書き換えることで、評価モードを動的に切り替える
    CURRENT_ALLELE = "A0201"
    print(f"🧬 現在のシミュレーション対象アレル: {CURRENT_ALLELE}")
    
    if not os.path.exists(mhc_pdb_path):
        print(f"❌ エラー: {mhc_pdb_path} が見つかりません。")
        return

    structure_mhc = get_biopython_structure("MHC_I", mhc_pdb_path)
    structure_pep = get_biopython_structure("Peptide", pep_pdb_path)
    print("🎉 構造データを正常にパースしました。")

    try:
        mhc_pocket_ids = [5, 7, 59, 84, 143]
        total_affinity_score = 0

        minimum_distance = float("inf")
        minimum_pair = None
        
        print(f"\n--- {CURRENT_ALLELE} 特異的ポケットスキャン開始 ---")
        
        for mhc_res_id in mhc_pocket_ids:
            atom_mhc = structure_mhc[0]['A'][mhc_res_id]['CA']
            
            for pep_res in structure_pep[0]['A']:
                pep_res_id = pep_res.id[1]
                atom_pep = pep_res['CA']
                distance = atom_mhc - atom_pep

                if distance < minimum_distance:
                   minimum_distance = distance
                   minimum_pair = {
                       "mhc_res_id": mhc_res_id,
                       "pep_res_id": pep_res_id,
                       "pep_res_name": pep_res.get_resname().upper(),
                    }
                
                # スクリーニング（MHC–peptide原子間距離が5.0 Å以内）
                if is_within_distance_threshold(distance):
                    score_increment = 1
                    bonus_text = ""
                    
                    # 残基名の取得
                    res_name = pep_res.get_resname().upper()
                    
                    # P2ポケットの動的ボーナス判定
                    if pep_res_id == 2:
                        if res_name in ALLELE_DATABASE[CURRENT_ALLELE]:
                            score_increment += 5
                            bonus_text = f" [★アレル特異的ボーナス発生! P2が好みの残基({res_name})と一致]"
                    
                    total_affinity_score += score_increment
                    print(f"  検出: MHC {mhc_res_id:3d} <--> ペプチド P{pep_res_id} ({res_name}) | 距離: {distance:5.2f}Å | +{score_increment}点{bonus_text}")

        if minimum_pair is not None:
           print(
               "\n🔎 最小CA原子間距離: "
               f"{minimum_distance:.2f} Å "
               f"(MHC {minimum_pair['mhc_res_id']} "
               f"<--> ペプチド P{minimum_pair['pep_res_id']} "
               f"{minimum_pair['pep_res_name']})"
         )

        print("-" * 50)
        print(f"🔥 【{CURRENT_ALLELE}】での最終予測親和性スコア: {total_affinity_score} 点")
        
    except KeyError as e:
        print(f"❌ エラー: 残基または原子が見つかりません: {e}")

if __name__ == "__main__":
    main()