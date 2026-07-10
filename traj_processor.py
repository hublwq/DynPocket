import os
import sys
import json
import argparse
import logging
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms

def setup_logger(target_dir, resume=False):
    logger = logging.getLogger("TrajProcessor")
    logger.setLevel(logging.INFO)
    if logger.handlers:
        return logger
    
    console_handler = logging.StreamHandler(sys.stdout)
    # 如果是续传，采用 'a' (append) 模式追加日志；否则覆盖
    file_mode = 'a' if resume else 'w'
    file_handler = logging.FileHandler(os.path.join(target_dir, "traj_processor.log"), mode=file_mode)
    
    formatter = logging.Formatter('[%(asctime)s] [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    console_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)
    
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)
    return logger

def load_progress(progress_file):
    if os.path.exists(progress_file):
        with open(progress_file, 'r') as f:
            return json.load(f)
    return {
        "raw_extraction_done": False,
        "valid_frames_count": 0,
        "rmsd_filter_done": False,
        "rep_frames_count": 0
    }

def save_progress(progress_file, state):
    with open(progress_file, 'w') as f:
        json.dump(state, f, indent=4)

def parse_args():
    parser = argparse.ArgumentParser(description="轨迹解包与贪婪 RMSD 去冗余工具")
    parser.add_argument('-d', '--base_dir', default="./dynpocket_out", help="主输出目录")
    parser.add_argument('-p', '--protein', required=True, help="蛋白序列名称")
    parser.add_argument('-t', '--threshold', type=float, default=5.0, help="贪婪筛选 RMSD 阈值 (默认: 5.0 Å)")
    parser.add_argument('-ref', '--reference', default=None, help="指定的参考 PDB 结构 (如果不提供，默认使用第 0 帧)")
    parser.add_argument('--resume', action='store_true', help="断点续传模式")
    return parser.parse_args()

def main():
    args = parse_args()
    
    target_dir = os.path.join(args.base_dir, args.protein)
    bioemu_dir = os.path.join(target_dir, "bioemu_out")
    raw_dir = os.path.join(target_dir, "conformations_raw")
    rep_dir = os.path.join(target_dir, "conformations_representative")
    progress_file = os.path.join(target_dir, "progress.json")
    
    if not os.path.exists(bioemu_dir):
        print(f"❌ 致命错误: 找不到 BioEmu 输出目录 {bioemu_dir}")
        sys.exit(1)
        
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(rep_dir, exist_ok=True)
    
    logger = setup_logger(target_dir, args.resume)
    state = load_progress(progress_file) if args.resume else load_progress("dummy.json")
    
    logger.info("="*60)
    logger.info(f"轨迹处理流水线启动")
    logger.info(f"目标蛋白: {args.protein} | 续传模式: {args.resume}")
    logger.info("="*60)

    top_path = os.path.join(bioemu_dir, "topology.pdb")
    xtc_path = os.path.join(bioemu_dir, "samples.xtc")

    if not os.path.exists(top_path) or not os.path.exists(xtc_path):
        logger.error(f"❌ 缺少核心轨迹文件！请检查 topology.pdb 和 samples.xtc")
        sys.exit(1)

    try:
        # ---------------------------------------------------------
        # 阶段一：解包到 conformations_raw/
        # ---------------------------------------------------------
        if state.get("raw_extraction_done"):
            logger.info("⏭️ [续传] 跳过解包阶段：conformations_raw/ 已存在完整数据。")
            u = mda.Universe(top_path, xtc_path) # 仍需载入对象供后续计算
            actual_frames = state["valid_frames_count"]
        else:
            logger.info("📦 阶段一: 正在将轨迹解包至 conformations_raw/ ...")
            u = mda.Universe(top_path, xtc_path)
            actual_frames = len(u.trajectory)
            
            for ts in u.trajectory:
                out_pdb = os.path.join(raw_dir, f"frame_{ts.frame:04d}.pdb")
                u.atoms.write(out_pdb)
                
            logger.info(f"   ✅ 解包完成！共提取 {actual_frames} 个有效构象帧。")
            state["raw_extraction_done"] = True
            state["valid_frames_count"] = actual_frames
            save_progress(progress_file, state)

        # ---------------------------------------------------------
        # 阶段二：Cα RMSD 贪婪去冗余
        # ---------------------------------------------------------
        if state.get("rmsd_filter_done"):
            logger.info("⏭️ [续传] 跳过去冗余阶段：conformations_representative/ 已存在结果。")
            logger.info(f"   最终代表帧集合 PDB 数量: {state['rep_frames_count']}")
        else:
            logger.info("⚙️ 阶段二: 执行 Cα RMSD 贪婪去冗余清洗...")
            
            # 1. 坐标系对齐 (C-alpha)
            if args.reference and os.path.exists(args.reference):
                logger.info(f"   [校准] 使用用户提供的参考结构: {args.reference}")
                ref_u = mda.Universe(args.reference)
                align.AlignTraj(u, ref_u, select='name CA', in_memory=True).run()
            else:
                logger.info("   [校准] 未提供参考结构，默认使用 conformations_raw/ 的第 0 帧作为参考基准。")
                align.AlignTraj(u, u, select='name CA', in_memory=True).run()

            # 提取对齐后的 C-alpha 坐标用于极速计算
            ca_atoms = u.select_atoms('name CA')
            coords_all = np.array([ca_atoms.positions for _ in u.trajectory])
            
            # 2. 初筛：保留与参考构象 (帧 0) RMSD > 5.0 Å 的进入候选池
            # 默认帧 0 是必须保留的第一张“代表帧”
            representative_indices = [0] 
            candidate_pool = []
            
            logger.info(f"   [初筛] 正在筛选与参考构象差异 > {args.threshold} Å 的帧进入候选池...")
            for i in range(1, actual_frames):
                dist_to_ref = rms.rmsd(coords_all[i], coords_all[0])
                if dist_to_ref > args.threshold:
                    candidate_pool.append(i)
                    
            logger.info(f"   [初筛完成] 共有 {len(candidate_pool)} 帧进入候选池。")

            # 3. 迭代贪心筛选
            logger.info(f"   [复筛] 正在候选池中执行贪心比对 (确保相互之间 RMSD > {args.threshold} Å)...")
            for cand_idx in candidate_pool:
                is_valid = True
                # 和目前已入选的“所有代表帧”逐一比对
                for rep_idx in representative_indices:
                    dist = rms.rmsd(coords_all[cand_idx], coords_all[rep_idx])
                    if dist <= args.threshold:
                        is_valid = False
                        break # 一旦发现和某个代表帧太相似，直接抛弃，看下一个候选者
                
                if is_valid:
                    representative_indices.append(cand_idx)

            rep_count = len(representative_indices)
            logger.info(f"   [降维完成] 贪心筛选结束，最终保留了 {rep_count} 个核心代表帧。")
            
            # 4. 落盘代表帧
            logger.info("💾 正在将代表帧复制并保存至 conformations_representative/ ...")
            for idx in representative_indices:
                u.trajectory[idx]
                out_file = os.path.join(rep_dir, f"rep_frame_{idx:04d}.pdb")
                u.atoms.write(out_file)

            # 5. 日志与状态收尾
            logger.info("="*60)
            logger.info(f"📊 任务统计报告:")
            logger.info(f"   - 原始有效帧总数: {actual_frames}")
            logger.info(f"   - 最终代表帧集合 PDB 数量: {rep_count}")
            logger.info("="*60)

            state["rmsd_filter_done"] = True
            state["rep_frames_count"] = rep_count
            save_progress(progress_file, state)

    except Exception as e:
        logger.error(f"❌ 轨迹处理中断，捕获到异常: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
