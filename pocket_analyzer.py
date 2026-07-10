import os
import sys
import json
import shutil
import argparse
import logging
import subprocess
from multiprocessing import Pool
import numpy as np
import MDAnalysis as mda

def setup_logger(target_dir, resume=False):
    logger = logging.getLogger("PocketAnalyzer")
    logger.setLevel(logging.INFO)
    if logger.handlers:
        return logger
    console_handler = logging.StreamHandler(sys.stdout)
    file_mode = 'a' if resume else 'w'
    file_handler = logging.FileHandler(os.path.join(target_dir, "pocket_analyzer.log"), mode=file_mode)
    formatter = logging.Formatter('[%(asctime)s] [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    console_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)
    return logger

# --- 全局状态机读写 ---
def load_progress(progress_file):
    if os.path.exists(progress_file):
        with open(progress_file, 'r') as f:
            return json.load(f)
    return {}

def save_progress(progress_file, state):
    with open(progress_file, 'w') as f:
        json.dump(state, f, indent=4)

def run_single_fpocket(args_tuple):
    """单帧 fpocket 调用工作函数"""
    pdb_path, pocket_data_dir = args_tuple
    base_name = os.path.splitext(os.path.basename(pdb_path))[0]
    cmd = ["fpocket", "-f", pdb_path]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        orig_out_dir = pdb_path.replace(".pdb", "_out")
        target_out_dir = os.path.join(pocket_data_dir, f"{base_name}_out")
        if os.path.exists(orig_out_dir):
            if os.path.exists(target_out_dir):
                shutil.rmtree(target_out_dir)
            shutil.move(orig_out_dir, target_out_dir)
            return True, base_name
    except Exception as e:
        return False, f"{base_name} 失败: {str(e)}"
    return False, f"{base_name} 未生成输出"

def parse_fpocket_lining_residues(pocket_out_dir):
    """解析 fpocket 生成的 pocket*_atm.pdb 文件，精准提取内衬残基编号"""
    pockets_dict = {}
    pockets_dir = os.path.join(pocket_out_dir, "pockets")
    
    # 兼容性寻找：有的版本直接输出在 _out 下，有的在 _out/pockets/ 下
    search_dir = pockets_dir if os.path.exists(pockets_dir) else pocket_out_dir
    
    if not os.path.exists(search_dir):
        return pockets_dict

    for f in os.listdir(search_dir):
        if f.endswith("_atm.pdb") and f.startswith("pocket"):
            # 提取口袋编号，例如 "pocket1_atm.pdb" -> "1"
            pocket_id = f.split("_")[0].replace("pocket", "")
            res_set = set()
            with open(os.path.join(search_dir, f), 'r') as pdb_file:
                for line in pdb_file:
                    if line.startswith("ATOM  ") or line.startswith("HETATM"):
                        try:
                            # 按照标准 PDB 格式规范，第 23-26 列为残基序列编号
                            res_seq = int(line[22:26].strip())
                            res_set.add(res_seq)
                        except ValueError:
                            continue
            if res_set:
                pockets_dict[pocket_id] = list(res_set)
                
    return pockets_dict

def parse_args():
    parser = argparse.ArgumentParser(description="DynPocket 口袋并行探测与拓扑双模式筛选引擎")
    parser.add_argument('-d', '--base_dir', default="./dynpocket_out", help="主输出目录")
    parser.add_argument('-p', '--protein', required=True, help="蛋白序列名称 (例如: 1AMU_B)")
    parser.add_argument('-n', '--cores', type=int, default=4, help="并行计算使用的 CPU 核心数")
    parser.add_argument('--mode', choices=['A', 'B'], default='A', help="筛选模式: A(全局聚类), B(靶向引导)")
    parser.add_argument('--ratio', type=float, default=0.3, help="质心距离筛选保留的前百分比")
    parser.add_argument('--target_residues', default=None, help="模式B专用的目标残基编号列表")
    parser.add_argument('--resume', action='store_true', help="断点续传模式") # 接入续传参数
    return parser.parse_args()

def main():
    args = parse_args()
    target_dir = os.path.join(args.base_dir, args.protein)
    rep_dir = os.path.join(target_dir, "conformations_representative")
    pocket_data_dir = os.path.join(target_dir, "pocket_data")
    progress_file = os.path.join(target_dir, "progress.json")
    
    if not os.path.exists(rep_dir) or not os.listdir(rep_dir):
        print(f"❌ 致命错误: 找不到代表帧构象目录或其中无文件: {rep_dir}")
        sys.exit(1)
        
    os.makedirs(pocket_data_dir, exist_ok=True)
    logger = setup_logger(target_dir, args.resume)
    state = load_progress(progress_file) if args.resume else load_progress("dummy.json")
    
    logger.info("="*60)
    logger.info(f"⚙️ DynPocket 口袋并行探测与筛选流启动")
    logger.info(f"处理对象: {args.protein} | 筛选模式: {args.mode} | 续传: {args.resume}")
    logger.info("="*60)

    pdb_files = sorted([os.path.join(rep_dir, f) for f in os.listdir(rep_dir) if f.endswith(".pdb")])
    ref_pdb = pdb_files[0]

    # --- 阶段 1: 并行 fpocket 计算 ---
    if state.get("fpocket_computation_done"):
        logger.info("⏭️ [续传] 跳过 fpocket 探测阶段：pocket_data/ 已存在并行输出。")
    else:
        logger.info(f"🚀 阶段一: 启动 multiprocessing 进程池 (核心数: {args.cores}) 调度 fpocket...")
        tasks = [(pdb, pocket_data_dir) for pdb in pdb_files]
        with Pool(processes=args.cores) as pool:
            results = pool.map(run_single_fpocket, tasks)
            
        success_count = sum(1 for r in results if r[0])
        logger.info(f"   ✅ 口袋探测结束！成功完成: {success_count}/{len(pdb_files)} 帧。")
        state["fpocket_computation_done"] = True
        save_progress(progress_file, state)

    # --- 阶段 2: 解析、映射与双模式筛选 ---
    if state.get("pocket_filtering_done"):
        logger.info("⏭️ [续传] 跳过口袋映射与筛选阶段：已存在最优结果字典。")
        logger.info("🌟 恭喜！当前蛋白的所有流水线分析已闭环完成！")
        sys.exit(0)

    ref_u = mda.Universe(ref_pdb)
    ref_ca = ref_u.select_atoms("name CA")
    ref_ca_coords = {atom.residue.resnum: atom.position for atom in ref_ca}

    logger.info("🗺️ 阶段二: 提取口袋内衬残基并统一映射至参考构象坐标轴...")
    all_pockets_meta = []

    for pdb_path in pdb_files:
        base_name = os.path.splitext(os.path.basename(pdb_path))[0]
        frame_id = "".join(filter(str.isdigit, base_name))
        p_out_dir = os.path.join(pocket_data_dir, f"{base_name}_out")
        
        if not os.path.exists(p_out_dir): continue
            
        pockets_dict = parse_fpocket_lining_residues(p_out_dir)
        for p_id, lining_res in pockets_dict.items():
            if not lining_res: continue
            valid_coords = [ref_ca_coords[r] for r in lining_res if r in ref_ca_coords]
            if not valid_coords: continue
            
            pocket_centroid = np.mean(valid_coords, axis=0)
            all_pockets_meta.append({
                "label": f"{frame_id}_{p_id}",
                "lining_residues": lining_res,
                "centroid": pocket_centroid
            })

    total_pockets_found = len(all_pockets_meta)
    if total_pockets_found == 0:
        logger.warning("⚠️ 未探测到任何有效口袋，流水线提前结束。")
        sys.exit(0)

    logger.info(f"🎯 阶段三: 执行模式 {args.mode} 空间距离初筛过滤...")
    pocket_centroids = np.array([p["centroid"] for p in all_pockets_meta])

    if args.mode == 'A':
        global_center = np.mean(pocket_centroids, axis=0)
        distances = np.linalg.norm(pocket_centroids - global_center, axis=1)
    else:
        target_ids = [int(x.strip()) for x in args.target_residues.split(",")]
        target_coords = [ref_ca_coords[r] for r in target_ids if r in ref_ca_coords]
        target_center = np.mean(target_coords, axis=0)
        distances = np.linalg.norm(pocket_centroids - target_center, axis=1)

    for idx, dist in enumerate(distances):
        all_pockets_meta[idx]["distance_to_center"] = float(dist)
        
    all_pockets_meta.sort(key=lambda x: x["distance_to_center"])
    keep_num = max(1, int(total_pockets_found * args.ratio))
    selected_pockets = all_pockets_meta[:keep_num]

    output_summary_file = os.path.join(target_dir, "selected_pockets.json")
    with open(output_summary_file, 'w') as f:
        json.dump({"selected_pockets": [{"label": p["label"], "distance": p["distance_to_center"], "lining_residues": p["lining_residues"]} for p in selected_pockets]}, f, indent=4)

    logger.info("="*60)
    logger.info(f"📊 筛选完成！探测总数: {total_pockets_found} -> 保留数: {len(selected_pockets)}")
    logger.info(f"   核心摘要已写入: {output_summary_file}")
    
    state["pocket_filtering_done"] = True
    save_progress(progress_file, state)

if __name__ == "__main__":
    main()
