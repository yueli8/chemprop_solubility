import pandas as pd
import numpy as np

# 读取预测结果文件
print("读取预测结果文件...")
predictions = pd.read_csv('peg_chentao_predictions.csv')

print(f"文件包含 {len(predictions)} 个化合物")
print(f"列名: {predictions.columns.tolist()}")

# 化合物信息（按顺序对应）
compounds = ['Taxol', 'PEG1-Taxol', 'PEG3-Taxol', 'PEG5-Taxol', 'PEG10-Taxol', 'PEG5*2-Taxol']

# 分子量 (g/mol)
molecular_weights = {
    'Taxol': 853.9,
    'PEG1-Taxol': 1100,
    'PEG3-Taxol': 1200,
    'PEG5-Taxol': 1300,
    'PEG10-Taxol': 1500,
    'PEG5*2-Taxol': 1500
}

# 实验数据（用于对比）
experimental_data = {
    'Taxol': '0.35-0.7 μg/mL',
    'PEG1-Taxol': '<1 mg/mL',
    'PEG3-Taxol': '<5 mg/mL',
    'PEG5-Taxol': '<10 mg/mL',
    'PEG10-Taxol': '10-20 mg/mL',
    'PEG5*2-Taxol': '10-20 mg/mL'
}

print("\n" + "="*120)
print("logS换算为溶解度（详细计算过程）")
print("="*120)

# 创建结果列表
results = []

print(f"\n{'化合物':<12} {'logS':<10} {'10^logS':<15} {'分子量':<10} {'mol/L':<12} {'g/L':<12} {'mg/mL':<12} {'μg/mL':<12}")
print("-"*120)

for i, (idx, row) in enumerate(predictions.iterrows()):
    compound = compounds[i]
    logS = row['logS']
    mw = molecular_weights[compound]
    
    # 详细计算过程
    # 步骤1: logS → 10^logS
    ten_power_logS = 10**logS
    
    # 步骤2: 10^logS × 分子量 = mol/L × g/mol = g/L
    solubility_g_L = ten_power_logS * mw
    
    # 步骤3: g/L = mg/mL（因为 1 g/L = 1 mg/mL）
    solubility_mg_ml = solubility_g_L
    
    # 步骤4: mg/mL × 1000 = μg/mL
    solubility_ug_ml = solubility_mg_ml * 1000
    
    # 选择合适的单位显示
    if solubility_ug_ml < 1:
        display = f"{solubility_ug_ml:.4f} μg/mL"
    elif solubility_mg_ml < 1:
        display = f"{solubility_mg_ml:.4f} mg/mL"
    else:
        display = f"{solubility_mg_ml:.2f} mg/mL"
    
    # 保存详细结果
    results.append({
        'Compound': compound,
        'SMILES': row['smiles'],
        'logS': logS,
        '10^logS': ten_power_logS,
        'Molecular_Weight_g_mol': mw,
        'Solubility_mol_L': ten_power_logS,
        'Solubility_g_L': solubility_g_L,
        'Solubility_mg_mL': solubility_mg_ml,
        'Solubility_ug_mL': solubility_ug_ml,
        'Solubility_Display': display,
        'Experimental': experimental_data[compound],
        'Calculation_Formula': f"10^({logS:.4f}) × {mw} = {solubility_mg_ml:.4f} mg/mL"
    })
    
    # 打印详细计算过程
    print(f"{compound:<12} {logS:<10.4f} {ten_power_logS:<15.6e} {mw:<10.1f} {ten_power_logS:<12.6e} {solubility_g_L:<12.6f} {solubility_mg_ml:<12.6f} {solubility_ug_ml:<12.4f}")

print("-"*120)

# 创建DataFrame
df_results = pd.DataFrame(results)

# 保存为CSV（包含所有详细列）
df_results.to_csv('solubility_conversion_detailed.csv', index=False)
print(f"\n✅ 详细换算结果已保存到 solubility_conversion_detailed.csv")

# 显示详细计算过程
print("\n" + "="*120)
print("详细计算过程展示")
print("="*120)

for result in results:
    print(f"\n{result['Compound']}:")
    print(f"  Step 1: logS = {result['logS']:.6f}")
    print(f"  Step 2: 10^logS = 10^({result['logS']:.6f}) = {result['10^logS']:.6e}")
    print(f"  Step 3: 分子量 = {result['Molecular_Weight_g_mol']:.1f} g/mol")
    print(f"  Step 4: 溶解度(mol/L) = 10^logS = {result['Solubility_mol_L']:.6e} mol/L")
    print(f"  Step 5: 溶解度(g/L) = {result['10^logS']:.6e} × {result['Molecular_Weight_g_mol']:.1f} = {result['Solubility_g_L']:.6f} g/L")
    print(f"  Step 6: 溶解度(mg/mL) = {result['Solubility_mg_mL']:.6f} mg/mL")
    print(f"  Step 7: 溶解度(μg/mL) = {result['Solubility_ug_mL']:.4f} μg/mL")
    print(f"  最终结果: {result['Solubility_Display']}")
    print(f"  计算公式: {result['Calculation_Formula']}")

# 创建一个简化的展示表格
print("\n" + "="*120)
print("简化展示表格")
print("="*120)

display_df = df_results[[
    'Compound', 
    'logS', 
    'Molecular_Weight_g_mol',
    'Solubility_mol_L',
    'Solubility_mg_mL',
    'Solubility_ug_mL',
    'Solubility_Display',
    'Experimental'
]].copy()

# 格式化数值列
display_df['logS'] = display_df['logS'].round(4)
display_df['Molecular_Weight_g_mol'] = display_df['Molecular_Weight_g_mol'].round(1)
display_df['Solubility_mol_L'] = display_df['Solubility_mol_L'].apply(lambda x: f"{x:.6e}")
display_df['Solubility_mg_mL'] = display_df['Solubility_mg_mL'].round(4)
display_df['Solubility_ug_mL'] = display_df['Solubility_ug_mL'].round(2)

print(display_df.to_string(index=False))

# 保存简化版
display_df.to_csv('solubility_display_table.csv', index=False)
print(f"\n✅ 简化展示表格已保存到 solubility_display_table.csv")

# 创建计算过程说明文件
with open('calculation_process.txt', 'w', encoding='utf-8') as f:
    f.write("logS换算溶解度详细计算过程\n")
    f.write("="*80 + "\n\n")
    f.write("基本公式：\n")
    f.write("1. 溶解度(mol/L) = 10^logS\n")
    f.write("2. 溶解度(g/L) = 溶解度(mol/L) × 分子量(g/mol)\n")
    f.write("3. 溶解度(mg/mL) = 溶解度(g/L) [因为 1 g/L = 1 mg/mL]\n")
    f.write("4. 溶解度(μg/mL) = 溶解度(mg/mL) × 1000\n\n")
    f.write("="*80 + "\n\n")
    
    for result in results:
        f.write(f"{result['Compound']}:\n")
        f.write(f"  logS = {result['logS']:.6f}\n")
        f.write(f"  10^logS = {result['10^logS']:.6e}\n")
        f.write(f"  分子量 = {result['Molecular_Weight_g_mol']:.1f} g/mol\n")
        f.write(f"  溶解度 = {result['10^logS']:.6e} × {result['Molecular_Weight_g_mol']:.1f} = {result['Solubility_mg_mL']:.6f} mg/mL\n")
        f.write(f"  最终结果: {result['Solubility_Display']}\n\n")

print(f"✅ 计算过程说明已保存到 calculation_process.txt")

print("\n" + "="*120)
print("换算完成！")
print("="*120)
