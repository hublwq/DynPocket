# utils.py
import os
import sys
import logging
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder

def validate_fasta_and_pdb(fasta_path: str, pdb_path: str, logger: logging.Logger):
    """
    解析 FASTA 和 PDB 文件，并进行严格的一致性比对校验。
    """
    logger.info("🔍 正在解析输入文件并执行强制一致性校验...")
    
    # 1. 解析 FASTA
    try:
        fasta_records = list(SeqIO.parse(fasta_path, "fasta"))
        if not fasta_records:
            logger.error(f"❌ FASTA 文件 {fasta_path} 解析失败或为空！")
            sys.exit(1)
        # 单序列模式下，取第一条
        target_seq_id = fasta_records[0].id
        fasta_seq = str(fasta_records[0].seq).upper()
        logger.info(f"   [FASTA] 提取序列: {target_seq_id} (长度: {len(fasta_seq)} AA)")
    except Exception as e:
        logger.error(f"❌ 无法读取 FASTA 文件: {e}")
        sys.exit(1)

    # 如果用户没有提供 PDB，直接放行
    if not pdb_path:
        logger.info("   [校验] 用户未提供参考 PDB，跳过结构序列校验。默认使用第 0 帧作为后续坐标系参考。")
        return target_seq_id, fasta_seq, None

    # 2. 解析 PDB 并提取序列
    if not os.path.exists(pdb_path):
        logger.error(f"❌ 找不到参考 PDB 文件: {pdb_path}")
        sys.exit(1)
        
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("reference", pdb_path)
        ppb = PPBuilder()
        pdb_seq = ""
        for pp in ppb.build_peptides(structure[0], aa_only=True):
            pdb_seq += str(pp.get_sequence())
        
        logger.info(f"   [PDB] 提取多肽序列成功 (长度: {len(pdb_seq)} AA)")
    except Exception as e:
        logger.error(f"❌ 无法解析 PDB 文件: {e}")
        sys.exit(1)

    # 3. 严格一致性比对
    if fasta_seq != pdb_seq:
        logger.error("❌ 致命错误：一致性校验失败！")
        logger.error("   FASTA 序列与参考 PDB 中提取的氨基酸序列不一致。")
        logger.error("   为防止后续坐标映射错位，程序已终止。请检查输入文件。")
        sys.exit(1)
        
    logger.info("✅ 强制一致性校验通过：FASTA 与 PDB 序列完美对应。")
    return target_seq_id, fasta_seq, structure
