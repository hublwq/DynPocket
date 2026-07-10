import os
import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="DynPocket 独立可视化提取器")
    parser.add_argument('protein_dir', help="蛋白主输出目录路径 (例如: dynpocket_out/1AMU_B)")
    parser.add_argument('--pocket', required=True, help="指定的口袋编号 (例如: 10_1)")
    return parser.parse_args()

def main():
    args = parse_args()
    
    inner_csv = os.path.join(args.protein_dir, "inner_residues.csv")
    if not os.path.exists(inner_csv):
        print(f"❌ 找不到 {inner_csv}，请确保已运行特征提取模块。")
        return

    df = pd.read_csv(inner_csv)
    target_row = df[df['Pocket_Label'] == args.pocket]
    
    if target_row.empty:
        print(f"❌ 在记录中未找到标签为 {args.pocket} 的口袋！")
        return
        
    residues_str = target_row.iloc[0]['Lining_Residues']
    # 提取数字编号
    res_items = [x.strip() for x in residues_str.split(",")]
    res_ids = [''.join(filter(str.isdigit, r)) for r in res_items]
    resi_string = "+".join(res_ids)
    
    frame_id = args.pocket.split('_')[0]
    target_pdb = f"rep_frame_{frame_id}.pdb"
    
    pml_out = os.path.join(args.protein_dir, "visualization_assets", f"vis_pocket_{args.pocket}.pml")
    
    with open(pml_out, 'w') as f:
        f.write(f"load ../conformations_representative/{target_pdb}, frame_{frame_id}\n")
        f.write("hide all\n")
        f.write(f"show cartoon, frame_{frame_id}\n")
        f.write(f"color white, frame_{frame_id}\n")
        f.write(f"select target_pocket, frame_{frame_id} and resi {resi_string}\n")
        f.write("show sticks, target_pocket\n")
        f.write("util.cbay target_pocket\n")
        f.write("show surface, target_pocket\n")
        f.write("set transparency, 0.5, target_pocket\n")
        f.write("zoom target_pocket\n")
        
    print(f"✅ 针对 {args.pocket} 的专属 PyMOL 高亮脚本已生成: {pml_out}")

if __name__ == "__main__":
    main()
