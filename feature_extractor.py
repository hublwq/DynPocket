import os
import sys
import json
import argparse
import logging
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda

def setup_logger(target_dir):
    logger = logging.getLogger("FeatureExtractor")
    logger.setLevel(logging.INFO)
    if logger.handlers: return logger
    ch = logging.StreamHandler(sys.stdout)
    fh = logging.FileHandler(os.path.join(target_dir, "feature_extraction.log"), mode='w')
    formatter = logging.Formatter('[%(asctime)s] [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)
    logger.addHandler(ch)
    logger.addHandler(fh)
    return logger

def parse_fpocket_info(info_file, target_pocket_id):
    features = {}
    is_target = False
    with open(info_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("Pocket"):
                current_id = line.split()[1].replace(":", "")
                is_target = (current_id == target_pocket_id)
            elif is_target and ":" in line:
                key, val = line.split(":", 1)
                key = key.strip()
                val = val.strip()
                # 提取有强生物学意义的特征
                if key in ["Score", "Volume", "Polar SASA", "Apolar SASA", "Proportion of polar atoms"]:
                    try:
                        features[key] = float(val.split()[0])
                    except ValueError:
                        pass
    return features

def parse_args():
    parser = argparse.ArgumentParser(description="DynPocket 特征提取与可视化引擎")
    parser.add_argument('-d', '--base_dir', default="./dynpocket_out", help="主输出目录")
    parser.add_argument('-p', '--protein', required=True, help="蛋白序列名称")
    parser.add_argument('--inner_shell_residues_num', type=int, default=20, help="最高频内衬残基输出数量 (-1代表全部)")
    return parser.parse_args()

def main():
    args = parse_args()
    target_dir = os.path.join(args.base_dir, args.protein)
    rep_dir = os.path.join(target_dir, "conformations_representative")
    pocket_data_dir = os.path.join(target_dir, "pocket_data")
    vis_dir = os.path.join(target_dir, "visualization_assets")
    json_file = os.path.join(target_dir, "selected_pockets.json")
    
    os.makedirs(vis_dir, exist_ok=True)
    logger = setup_logger(target_dir)
    
    logger.info("="*60)
    logger.info("📈 DynPocket 理化特征提取与可视化模块启动")
    logger.info("="*60)

    if not os.path.exists(json_file):
        logger.error(f"❌ 找不到 {json_file}，请先运行 pocket_analyzer.py")
        sys.exit(1)

    with open(json_file, 'r') as f:
        selected_data = json.load(f)
    
    selected_pockets = selected_data.get("selected_pockets", [])
    
    pdb_files = sorted([f for f in os.listdir(rep_dir) if f.endswith(".pdb")])
    ref_pdb_path = os.path.join(rep_dir, pdb_files[0])
    u_ref = mda.Universe(ref_pdb_path)
    res_dict = {res.resnum: res.resname for res in u_ref.residues}

    logger.info("🔍 正在跨帧提取 fpocket 核心理化动力学特征...")
    feature_records = []
    inner_residues_records = []
    all_lining_residues = []

    for pocket in selected_pockets:
        label = pocket["label"]
        frame_id, pocket_id = label.split("_")
        
        info_file = os.path.join(pocket_data_dir, f"rep_frame_{frame_id}_out", f"rep_frame_{frame_id}_info.txt")
        if os.path.exists(info_file):
            features = parse_fpocket_info(info_file, pocket_id)
            features["Pocket_Label"] = label
            features["Frame"] = int(frame_id)
            feature_records.append(features)
        
        res_nums = pocket["lining_residues"]
        res_names = [f"{res_dict.get(num, 'UNK')}{num}" for num in res_nums]
        
        inner_residues_records.append({
            "Pocket_Label": label,
            "Lining_Residues": ", ".join(res_names)
        })
        all_lining_residues.extend(res_names)

    df_features = pd.DataFrame(feature_records)
    cols = ['Pocket_Label', 'Frame'] + [c for c in df_features.columns if c not in ['Pocket_Label', 'Frame']]
    df_features = df_features[cols].sort_values(by="Frame")

    numeric_cols = df_features.select_dtypes(include=[np.number]).columns.drop('Frame')
    df_mean = df_features[numeric_cols].mean().to_frame().T
    df_mean.insert(0, 'Pocket_Label', 'MEAN')
    df_var = df_features[numeric_cols].var().to_frame().T
    df_var.insert(0, 'Pocket_Label', 'VARIANCE')

    pockets_csv_path = os.path.join(target_dir, "pockets.csv")
    with open(pockets_csv_path, 'w') as f:
        df_mean.to_csv(f, index=False)
        df_var.to_csv(f, index=False, header=False)
        df_features.to_csv(f, index=False, header=False)

    inner_csv_path = os.path.join(target_dir, "inner_residues.csv")
    pd.DataFrame(inner_residues_records).to_csv(inner_csv_path, index=False)
    
    res_counter = Counter(all_lining_residues)
    top_limit = len(res_counter) if args.inner_shell_residues_num == -1 else args.inner_shell_residues_num
    top_residues = res_counter.most_common(top_limit)
    
    list_txt_path = os.path.join(target_dir, "inner_residues_list.txt")
    top_res_names = []
    with open(list_txt_path, 'w') as f:
        f.write(f"# Top {top_limit} High-Frequency Lining Residues\n")
        for res, count in top_residues:
            f.write(f"{res}\t{count}\n")
            top_res_names.append(res)

    logger.info("📊 正在绘制动态理化特征波动分析图...")
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # 散点图 + 均值趋势线 (剔除无意义的 Score)
    axes[0].scatter(df_features['Frame'], df_features['Volume'], color='royalblue', alpha=0.5, s=40, edgecolors='none', label='Pockets')
    mean_vol = df_features.groupby('Frame')['Volume'].mean()
    axes[0].plot(mean_vol.index, mean_vol.values, color='crimson', linewidth=2.5, label='Mean Trend')
    axes[0].set_ylabel("Volume (Å³)")
    axes[0].set_title("Dynamic Pocket Volume Fluctuation")
    axes[0].legend()
    axes[0].grid(True, linestyle='--', alpha=0.6)

    if 'Apolar SASA' in df_features.columns:
        axes[1].scatter(df_features['Frame'], df_features['Apolar SASA'], color='mediumseagreen', alpha=0.5, s=40, edgecolors='none')
        mean_sasa = df_features.groupby('Frame')['Apolar SASA'].mean()
        axes[1].plot(mean_sasa.index, mean_sasa.values, color='crimson', linewidth=2.5)
        axes[1].set_ylabel("Apolar SASA (Å²)")
        axes[1].set_title("Hydrophobicity (Apolar SASA) Variation")
    
    axes[1].set_xlabel("Representative Frame Number")
    axes[1].grid(True, linestyle='--', alpha=0.6)
    
    plot_path = os.path.join(vis_dir, "properties_variation_plot.png")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300)

    # 恢复原版美观 PyMOL 脚本
    pml_path = os.path.join(vis_dir, "highlight_core_pocket.pml")
    res_ids = [''.join(filter(str.isdigit, r)) for r in top_res_names]
    resi_string = "+".join(res_ids)
    
    with open(pml_path, 'w') as f:
        f.write(f"load ../conformations_representative/{os.path.basename(ref_pdb_path)}, ref_struct\n")
        f.write("bg_color white\n")
        f.write("hide all\n")
        f.write("show cartoon, ref_struct\n")
        f.write("color gray80, ref_struct\n")
        f.write(f"select core_pocket, ref_struct and resi {resi_string}\n")
        f.write("show surface, core_pocket\n")
        f.write("color red, core_pocket\n")
        f.write("set transparency, 0.4, core_pocket\n")
        f.write("zoom core_pocket\n")
        
    logger.info("🌟 所有解析任务圆满收官！")

if __name__ == "__main__":
    main()
