import argparse
import sys
import os
import logging
import subprocess
from utils import validate_fasta_and_pdb
from Bio import SeqIO

def setup_logger(out_dir):
    logger = logging.getLogger("DynPocket-Master")
    logger.setLevel(logging.INFO)
    if logger.handlers: return logger
    ch = logging.StreamHandler(sys.stdout)
    fh = logging.FileHandler(os.path.join(out_dir, "dynpocket_master.log"), mode='a')
    formatter = logging.Formatter('[%(asctime)s] [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)
    logger.addHandler(ch)
    logger.addHandler(fh)
    return logger

def parse_args():
    parser = argparse.ArgumentParser(description="DynPocket-Analyzer: 全栈端到端动态口袋分析流水线")
    parser.add_argument('-i', '--fasta', required=True, help="输入蛋白的 FASTA 序列文件")
    parser.add_argument('-p', '--pdb', default=None, help="参考 PDB 结构文件路径")
    parser.add_argument('-o', '--out_dir', default="./dynpocket_out", help="输出主目录")
    
    # 管线与硬件控制
    parser.add_argument('--frame_num', type=int, default=1000, help="BioEmu 采样的构象数量")
    parser.add_argument('--cores', type=int, default=4, help="fpocket 调用的并发核心数")
    parser.add_argument('--force', action='store_true', help="跳过序列长度限制等安全检查")
    parser.add_argument('--check_env', action='store_true', help="环境自检")
    parser.add_argument('--resume', action='store_true', help="开启全局断点续传")
    
    # 算法控制
    parser.add_argument('--mode', choices=['A', 'B'], default='A', help="筛选模式 A(聚类) B(靶向)")
    parser.add_argument('--ratio', type=float, default=0.3, help="口袋筛选保留比例")
    parser.add_argument('--target_residues', default=None, help="模式B专用的目标残基编号列表")
    parser.add_argument('--inner_shell_residues_num', type=int, default=20, help="输出的高频内衬残基数量")
    parser.add_argument('--preset', default="default", help="fpocket 参数预设套件 (预留扩展接口)")
    
    return parser.parse_args()

def check_environment(logger):
    """环境自检模块"""
    logger.info("🔍 执行依赖环境自检...")
    tools = ['fpocket', 'python']
    all_pass = True
    for tool in tools:
        if shutil.which(tool) is None:
            logger.error(f"❌ 缺失系统级依赖: {tool}")
            all_pass = False
    if not all_pass:
        sys.exit(1)
    logger.info("✅ 环境自检通过！")

def run_step(cmd, step_name, logger):
    logger.info(f"\n{'='*60}\n▶ 启动核心调度模块: {step_name}\n{'='*60}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        logger.error(f"❌ 致命错误：子模块 {step_name} 崩溃，主流水线安全熔断。")
        sys.exit(1)

def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    logger = setup_logger(args.out_dir)
    
    logger.info("🚀 DynPocket-Analyzer 工业级流水线启动")
    
    if args.check_env:
        import shutil
        check_environment(logger)
        sys.exit(0)

    # 1. 批量 FASTA 预解析
    fasta_records = list(SeqIO.parse(args.fasta, "fasta"))
    if not fasta_records:
        logger.error("❌ 输入的 FASTA 文件无效或为空。")
        sys.exit(1)

    logger.info(f"🧬 检测到 Multi-FASTA 模式，共包含 {len(fasta_records)} 条蛋白序列。")

    for record in fasta_records:
        protein_name = record.id
        logger.info(f"\n{'#'*60}\n# 开始处理流水线: {protein_name}\n{'#'*60}")
        
        # 将当前记录写出为临时文件供校验
        temp_fasta = os.path.join(args.out_dir, f"{protein_name}_temp.fasta")
        SeqIO.write(record, temp_fasta, "fasta")
        
        # 强制一致性与安全校验
        seq_id, fasta_seq, _ = validate_fasta_and_pdb(temp_fasta, args.pdb, logger)
        if len(fasta_seq) > 1500 and not args.force:
            logger.error(f"❌ 序列 {protein_name} 超出 1500AA 限制，已跳过。使用 --force 强行解锁。")
            continue

        # ---------------------------------------------------------
        # 步骤 A: 解包与 Cα RMSD 去冗余 (调用 traj_processor.py)
        # ---------------------------------------------------------
        cmd_traj = [
            sys.executable, "traj_processor.py", 
            "-p", protein_name,
            "-d", args.out_dir
        ]
        if args.resume: cmd_traj.append("--resume")
        if args.pdb: cmd_traj.extend(["-ref", args.pdb])
        run_step(cmd_traj, f"Phase A: {protein_name} 贪婪降维", logger)

        # ---------------------------------------------------------
        # 步骤 B: 多进程 fpocket 与双模式筛选 (调用 pocket_analyzer.py)
        # ---------------------------------------------------------
        cmd_pocket = [
            sys.executable, "pocket_analyzer.py",
            "-p", protein_name,
            "-d", args.out_dir,
            "-n", str(args.cores),
            "--mode", args.mode,
            "--ratio", str(args.ratio)
        ]
        if args.resume: cmd_pocket.append("--resume")
        if args.target_residues: cmd_pocket.extend(["--target_residues", args.target_residues])
        run_step(cmd_pocket, f"Phase B: {protein_name} 多进程探测与筛选", logger)

        # ---------------------------------------------------------
        # 步骤 C: 特征提取与可视化制图 (调用 feature_extractor.py)
        # ---------------------------------------------------------
        cmd_feature = [
            sys.executable, "feature_extractor.py",
            "-p", protein_name,
            "-d", args.out_dir,
            "--inner_shell_residues_num", str(args.inner_shell_residues_num)
        ]
        run_step(cmd_feature, f"Phase C: {protein_name} 特征清洗与图表化", logger)
        
        # 清理临时切片文件
        os.remove(temp_fasta)

    logger.info("\n🎉 批量分析管线执行完毕！")

if __name__ == "__main__":
    main()
