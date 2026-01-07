# GlycoConjugate-Vaccine-Analysis-Pipeline

糖鎖抱合ワクチン（Glycoconjugate Vaccine）の設計から、AlphaFold 3を用いた構造予測、および予測モデルの定量的評価を一貫して行う統合解析パイプラインです。

## 概要
本プロジェクトは、キャリアタンパク質と腫瘍関連糖鎖（TACA）を結合させた免疫抗原の三次元構造を最適化することを目的としています。複雑なSMILES構造からの原子特定、AlphaFold 3用の入力JSON作成、および出力された構造モデル（CIF）のSASA（溶媒露出面積）ベースのランキング機能を備えています。

## 主な機能
- **自動原子特定**: SMARTSパターンを用いて、リンカー末端の結合原子を自動的に特定。
- **AF3入力自動生成**: 共有結合（Covalent Bond）の指定を含むAlphaFold 3用JSONファイルの作成。
- **一括構造解析**: 予測された全モデルに対し、SASA（溶媒露出面積）とタンパク質-糖鎖間接触残基を算出。
- **自動スコアリング**: 糖鎖の露出度に基づき、抗体結合に最適なモデルをランキング。
- **可視化**: `py3Dmol`を用いたJupyter上でのインタラクティブな構造確認。

## セットアップ
以下のライブラリが必要です。
```bash
pip install rdkit-pypi biopython freesasa pandas py3dmol
