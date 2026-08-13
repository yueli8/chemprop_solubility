import pandas as pd
import numpy as np

# 读取SMILES
smiles_df = pd.read_csv('peg_chentao_6samples.csv')

# 正确的数据
data = [
    {'name': 'Taxol', 'solubility': '0.35-0.7 μg/mL', 'mw': 853.9},
    {'name': 'PEG1-Taxol', 'solubility': '<1 mg/mL', 'mw': 1100},
    {'name': 'PEG3-Taxol', 'solubility': '<5 mg/mL', 'mw': 1200},
    {'name': 'PEG5-Taxol', 'solubility': '<10 mg/mL', 'mw': 1300},
    {'name': 'PEG10-Taxol', 'solubility': '10-20 mg/mL', 'mw': 1500},
    {'name': 'PEG5*2-Taxol', 'solubility': '10-20 mg/mL', 'mw': 1500}
]

def calc_logs(sol, mw):
    # 统一转换为mg/mL
    if 'μg/mL' in sol:
        clean = sol.replace('μg/mL', '').strip()
        if '-' in clean:
            low, high = map(float, clean.split('-'))
            mg_ml = (low + high) / 2 / 1000  # μg转mg
        else:
            mg_ml = float(clean) / 1000
    elif 'mg/mL' in sol:
        clean = sol.replace('mg/mL', '').replace('<', '').strip()
        if '-' in clean:
            low, high = map(float, clean.split('-'))
            mg_ml = (low + high) / 2
        else:
            mg_ml = float(clean) * 0.5  # 对于"小于"取一半
    else:
        return None
    
    # 计算mol/L
    mol_l = mg_ml / mw  # mg/mL ÷ g/mol = mol/L
    return np.log10(mol_l)

# 生成正确的训练数据
training_data = []
for i, row in smiles_df.iterrows():
    smiles = row.iloc[0]
    logs = calc_logs(data[i]['solubility'], data[i]['mw'])
    training_data.append({'smiles': smiles, 'logS': logs})
    print(f"{data[i]['name']}: logS = {logs:.4f}")

# 保存正确的训练数据
train_df = pd.DataFrame(training_data)
train_df.to_csv('training_data_CORRECT.csv', index=False)

# 生成正确的增强数据
augmented = []
for _, row in train_df.iterrows():
    for var in [-0.3, -0.15, 0, 0.15, 0.3]:
        augmented.append({'smiles': row['smiles'], 'logS': row['logS'] + var})

aug_df = pd.DataFrame(augmented)
aug_df.to_csv('augmented_training_data_CORRECT.csv', index=False)

print("\n✅ 正确的文件已保存：")
print("1. training_data_CORRECT.csv")
print("2. augmented_training_data_CORRECT.csv")
