#!/usr/bin/env python3
"""
Chemprop 水溶性预测结果完整计算工具
输入: 包含 SMILES 和预测 LogS 的 CSV 文件
输出: 包含分子量、溶解度换算的完整结果 CSV
"""

from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd
import argparse
import sys
import os

def calculate_molecular_weight(smiles):
    """
    从SMILES计算精确分子量和平均分子量
    返回: (精确分子量, 平均分子量, 分子式)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, None
    
    exact_mw = Descriptors.ExactMolWt(mol)      # 精确分子量（单一同位素）
    avg_mw = Descriptors.MolWt(mol)             # 平均分子量（考虑同位素丰度）
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)  # 分子式
    
    return exact_mw, avg_mw, formula

def convert_solubility(logS, molecular_weight):
    """
    将预测的LogS转换为不同单位的溶解度
    
    参数:
        logS: 预测的对数溶解度 (log10 mol/L)
        molecular_weight: 分子量 (g/mol)
    
    返回:
        S_mol_L: 溶解度 (mol/L)
        S_mg_mL: 溶解度 (mg/mL)
        S_ug_mL: 溶解度 (μg/mL)
        water_per_mg: 溶解1mg所需水量 (mL)
    """
    # LogS → S (mol/L)
    S_mol_L = 10 ** logS
    
    # S (mol/L) → S (mg/mL)
    S_mg_mL = S_mol_L * molecular_weight * 1000
    
    # S (mg/mL) → S (μg/mL)
    S_ug_mL = S_mg_mL * 1000
    
    # 溶解1mg所需水量 (mL)
    if S_mg_mL > 0:
        water_per_mg = 1.0 / S_mg_mL
    else:
        water_per_mg = float('inf')
    
    return S_mol_L, S_mg_mL, S_ug_mL, water_per_mg

def classify_solubility(logS):
    """
    根据LogS对水溶性进行评级
    
    参考标准:
        > -2:  极高
        -2 ~ -4: 高
        -4 ~ -6: 中等
        -6 ~ -8: 低
        < -8:  极低
    """
    if logS > -2:
        return "极高"
    elif logS > -4:
        return "高"
    elif logS > -6:
        return "中等"
    elif logS > -8:
        return "低"
    else:
        return "极低"

def process_csv(input_file, output_file=None, smiles_col='smiles', logS_col='pred_0'):
    """
    处理CSV文件，计算所有分子的分子量和溶解度
    
    参数:
        input_file: 输入CSV文件路径
        output_file: 输出CSV文件路径（默认: 输入文件名_complete.csv）
        smiles_col: SMILES列名（默认: smiles）
        logS_col: 预测LogS列名（默认: pred_0）
    """
    # 读取CSV
    print(f"正在读取输入文件: {input_file}")
    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        print(f"❌ 读取CSV失败: {e}")
        sys.exit(1)
    
    print(f"共读取 {len(df)} 个分子")
    print(f"列名: {list(df.columns)}")
    
    # 检查必要列
    if smiles_col not in df.columns:
        print(f"❌ 找不到SMILES列 '{smiles_col}'，可用列: {list(df.columns)}")
        sys.exit(1)
    
    if logS_col not in df.columns:
        print(f"❌ 找不到预测LogS列 '{logS_col}'，可用列: {list(df.columns)}")
        # 尝试自动查找
        for col in df.columns:
            if 'pred' in col.lower() or 'log' in col.lower() or 'solubility' in col.lower():
                logS_col = col
                print(f"自动选择LogS列: '{logS_col}'")
                break
        else:
            print("❌ 无法自动识别LogS列，请使用 --logs-col 参数指定")
            sys.exit(1)
    
    # 初始化结果列表
    results = []
    
    # 处理每个分子
    print(f"\n开始计算...")
    success_count = 0
    fail_count = 0
    
    for idx, row in df.iterrows():
        smiles = row[smiles_col]
        logS = row[logS_col]
        
        # 跳过空值
        if pd.isna(smiles) or str(smiles).strip() == '':
            print(f"  行 {idx+1}: SMILES为空，跳过")
            fail_count += 1
            continue
        
        if pd.isna(logS):
            print(f"  行 {idx+1}: LogS为空，跳过")
            fail_count += 1
            continue
        
        # 计算分子量
        exact_mw, avg_mw, formula = calculate_molecular_weight(str(smiles).strip())
        
        if exact_mw is None:
            print(f"  行 {idx+1}: SMILES解析失败: {smiles[:50]}...")
            fail_count += 1
            # 仍然添加行，但分子量信息为空
            results.append({
                '原始行索引': idx + 1,
                'smiles': smiles,
                '预测LogS': logS,
                '分子式': '解析失败',
                '精确分子量 (g/mol)': None,
                '平均分子量 (g/mol)': None,
                '溶解度 S (mol/L)': None,
                '溶解度 (mg/mL)': None,
                '溶解度 (μg/mL)': None,
                '溶解1mg需水量 (mL)': None,
                '水溶性评级': '未知',
                'SMILES解析状态': '失败'
            })
            continue
        
        # 换算溶解度
        S_mol_L, S_mg_mL, S_ug_mL, water_per_mg = convert_solubility(logS, avg_mw)
        
        # 评级
        rating = classify_solubility(logS)
        
        # 保存结果
        result = {
            '原始行索引': idx + 1,
            'smiles': smiles,
            '预测LogS': logS,
            '分子式': formula,
            '精确分子量 (g/mol)': round(exact_mw, 4),
            '平均分子量 (g/mol)': round(avg_mw, 4),
            '溶解度 S (mol/L)': f"{S_mol_L:.4e}",
            '溶解度 S_num (mol/L)': S_mol_L,
            '溶解度 (mg/mL)': f"{S_mg_mL:.6f}",
            '溶解度 (μg/mL)': f"{S_ug_mL:.4f}",
            '溶解1mg需水量 (mL)': f"{water_per_mg:.2f}",
            '水溶性评级': rating,
            'SMILES解析状态': '成功'
        }
        
        results.append(result)
        success_count += 1
        
        # 打印进度
        print(f"  [{success_count}/{len(df)}] {formula} | LogS={logS:.2f} | "
              f"MW={avg_mw:.2f} | 溶解度={S_ug_mL:.4f} μg/mL | {rating}")
    
    # 生成输出文件名
    if output_file is None:
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}_complete.csv"
    
    # 保存结果
    result_df = pd.DataFrame(results)
    result_df.to_csv(output_file, index=False, encoding='utf-8-sig')
    
    # 打印统计信息
    print(f"\n{'='*80}")
    print(f"计算完成!")
    print(f"  成功: {success_count} 个分子")
    print(f"  失败: {fail_count} 个分子")
    print(f"  结果已保存至: {output_file}")
    print(f"{'='*80}")
    
    # 打印汇总统计
    if success_count > 0:
        valid_results = [r for r in results if r['SMILES解析状态'] == '成功']
        if valid_results:
            logs_values = [r['预测LogS'] for r in valid_results]
            mw_values = [r['平均分子量 (g/mol)'] for r in valid_results]
            ug_values = [float(r['溶解度 (μg/mL)']) for r in valid_results]
            
            print(f"\n📊 汇总统计:")
            print(f"  LogS范围:     {min(logs_values):.2f} ~ {max(logs_values):.2f}")
            print(f"  分子量范围:   {min(mw_values):.1f} ~ {max(mw_values):.1f} g/mol")
            print(f"  溶解度范围:   {min(ug_values):.4f} ~ {max(ug_values):.4f} μg/mL")
            
            # 评级分布
            ratings = [r['水溶性评级'] for r in valid_results]
            from collections import Counter
            rating_counts = Counter(ratings)
            print(f"\n  水溶性评级分布:")
            for rating, count in rating_counts.items():
                print(f"    {rating}: {count} 个")
            
            # 找出最好和最差
            best = max(valid_results, key=lambda r: r['预测LogS'])
            worst = min(valid_results, key=lambda r: r['预测LogS'])
            print(f"\n  水溶性最好: 行{best['原始行索引']} (LogS={best['预测LogS']:.2f}, {best['溶解度 (μg/mL)']} μg/mL)")
            print(f"  水溶性最差: 行{worst['原始行索引']} (LogS={worst['预测LogS']:.2f}, {worst['溶解度 (μg/mL)']} μg/mL)")
    
    return result_df

def main():
    parser = argparse.ArgumentParser(
        description='Chemprop水溶性预测结果完整计算工具',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
使用示例:
  # 基本用法
  python solubility_calculator.py -i predictions.csv
  
  # 指定列名
  python solubility_calculator.py -i predictions.csv --smiles-col SMILES --logs-col pred_0
  
  # 指定输出文件
  python solubility_calculator.py -i predictions.csv -o results.csv
  
  # 你的CSV需要有smiles列和pred_0列（或自定义列名）
        '''
    )
    
    parser.add_argument('-i', '--input', required=True, 
                        help='输入CSV文件路径（包含SMILES和预测LogS）')
    parser.add_argument('-o', '--output', default=None,
                        help='输出CSV文件路径（默认: 输入文件名_complete.csv）')
    parser.add_argument('--smiles-col', default='smiles',
                        help='SMILES列名（默认: smiles）')
    parser.add_argument('--logs-col', default='pred_0',
                        help='预测LogS列名（默认: pred_0）')
    
    args = parser.parse_args()
    
    process_csv(
        input_file=args.input,
        output_file=args.output,
        smiles_col=args.smiles_col,
        logS_col=args.logs_col
    )

if __name__ == '__main__':
    main()
