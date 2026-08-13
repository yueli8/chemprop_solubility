import pandas as pd
import numpy as np

# 读取预测结果
predictions = pd.read_csv('peg_chentao_predictions.csv')

# 化合物信息
compounds = ['Taxol', 'PEG1-Taxol', 'PEG3-Taxol', 'PEG5-Taxol', 'PEG10-Taxol', 'PEG5*2-Taxol']

# 实验数据
experimental = {
    'Taxol': {'solubility': '0.35-0.7 μg/mL', 'mw': 853.9},
    'PEG1-Taxol': {'solubility': '<1 mg/mL', 'mw': 1100},
    'PEG3-Taxol': {'solubility': '<5 mg/mL', 'mw': 1200},
    'PEG5-Taxol': {'solubility': '<10 mg/mL', 'mw': 1300},
    'PEG10-Taxol': {'solubility': '10-20 mg/mL', 'mw': 1500},
    'PEG5*2-Taxol': {'solubility': '10-20 mg/mL', 'mw': 1500}
}

# 创建结果表格
print("="*100)
print("Taxol衍生物水溶性预测结果对比")
print("="*100)

results = []
for i, compound in enumerate(compounds):
    pred_logs = predictions['logS'].iloc[i]
    mw = experimental[compound]['mw']
    exp_sol = experimental[compound]['solubility']
    
    # 将预测的logS转换为mg/mL
    pred_mg_ml = 10**pred_logs * mw
    
    # 转换为合适的单位显示
    if pred_mg_ml < 0.001:
        pred_display = f"{pred_mg_ml*1000:.4f} μg/mL"
    elif pred_mg_ml < 1:
        pred_display = f"{pred_mg_ml:.4f} mg/mL"
    else:
        pred_display = f"{pred_mg_ml:.2f} mg/mL"
    
    results.append({
        'Compound': compound,
        'Experimental': exp_sol,
        'Predicted_logS': pred_logs,
        'Predicted_Solubility': pred_display,
        'Predicted_mg_ml': pred_mg_ml
    })
    
    print(f"\n{compound}:")
    print(f"  实验数据: {exp_sol}")
    print(f"  预测logS: {pred_logs:.4f}")
    print(f"  预测溶解度: {pred_display}")

# 创建DataFrame
df_results = pd.DataFrame(results)

# 保存结果
df_results.to_csv('prediction_vs_experimental.csv', index=False)

print("\n" + "="*100)
print("总结")
print("="*100)

# 分析趋势
print("\n溶解度趋势分析：")
print("-" * 50)
for i in range(len(compounds)-1):
    current = results[i]['Predicted_mg_ml']
    next_comp = results[i+1]['Predicted_mg_ml']
    improvement = next_comp / current
    print(f"{compounds[i]} → {compounds[i+1]}: 溶解度提高 {improvement:.1f} 倍")

# 关键发现
print("\n关键发现：")
print("-" * 50)
print(f"1. Taxol最难溶: {results[0]['Predicted_Solubility']}")
print(f"2. PEG修饰显著提高溶解度")
print(f"3. PEG10-Taxol溶解度最好: {results[4]['Predicted_Solubility']}")
print(f"4. 从Taxol到PEG10-Taxol，溶解度提高约 {results[4]['Predicted_mg_ml']/results[0]['Predicted_mg_ml']:.0f} 倍")

# 转换为mol/L显示
print("\n\n预测的logS值（mol/L）：")
print("-" * 50)
for i, compound in enumerate(compounds):
    print(f"{compound:<15} logS = {results[i]['Predicted_logS']:.4f}")
