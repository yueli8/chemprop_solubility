import pandas as pd
import numpy as np

# 读取数据
predictions = pd.read_csv('predictions_training.csv')
actual = pd.read_csv('augmented_training_data_CORRECT.csv')

print("预测文件列名:", predictions.columns.tolist())
print("实际文件列名:", actual.columns.tolist())
print("预测文件前几行:")
print(predictions.head())
print("\n实际文件前几行:")
print(actual.head())

# 合并数据
results = pd.DataFrame({
    'smiles': actual['smiles'],
    'actual_logS': actual['logS'],
    'predicted_logS': predictions['logS']  # 预测结果的列名是logS
})

# 添加化合物名称（根据顺序）
compounds = ['Taxol', 'PEG1-Taxol', 'PEG3-Taxol', 'PEG5-Taxol', 'PEG10-Taxol', 'PEG5*2-Taxol']
results['compound'] = [compounds[i//5] for i in range(len(results))]  # 每5个样本对应一个化合物

print("\n" + "="*80)
print("模型预测性能总结（按化合物）")
print("="*80)

# 按化合物分组
summary = results.groupby('compound').agg({
    'actual_logS': ['mean', 'std', 'min', 'max'],
    'predicted_logS': ['mean', 'std', 'min', 'max']
}).round(4)

print("\n详细统计：")
print(summary)

# 计算每个化合物的误差
print("\n" + "="*80)
print("各化合物预测误差")
print("="*80)

for compound in compounds:
    compound_data = results[results['compound'] == compound]
    
    actual_mean = compound_data['actual_logS'].mean()
    pred_mean = compound_data['predicted_logS'].mean()
    error = actual_mean - pred_mean
    abs_error = abs(error)
    
    print(f"\n{compound}:")
    print(f"  实际平均 logS: {actual_mean:.4f}")
    print(f"  预测平均 logS: {pred_mean:.4f}")
    print(f"  平均误差: {error:.4f}")
    print(f"  平均绝对误差: {abs_error:.4f}")
    print(f"  实际范围: [{compound_data['actual_logS'].min():.4f}, {compound_data['actual_logS'].max():.4f}]")
    print(f"  预测范围: [{compound_data['predicted_logS'].min():.4f}, {compound_data['predicted_logS'].max():.4f}]")

# 计算总体统计
print("\n" + "="*80)
print("总体性能指标")
print("="*80)

results['error'] = results['actual_logS'] - results['predicted_logS']
results['abs_error'] = results['error'].abs()

mae = results['abs_error'].mean()
rmse = np.sqrt((results['error']**2).mean())
r2 = 1 - (results['error']**2).sum() / ((results['actual_logS'] - results['actual_logS'].mean())**2).sum()

print(f"样本数量: {len(results)}")
print(f"平均绝对误差 (MAE): {mae:.4f} log单位")
print(f"均方根误差 (RMSE): {rmse:.4f} log单位")
print(f"决定系数 (R²): {r2:.4f}")

# 保存详细结果
results.to_csv('detailed_analysis.csv', index=False)
print(f"\n详细结果已保存到 detailed_analysis.csv")

# 创建可视化（如果matplotlib可用）
try:
    import matplotlib.pyplot as plt
    
    # 散点图
    plt.figure(figsize=(10, 6))
    
    # 为每个化合物使用不同颜色
    colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown']
    for i, compound in enumerate(compounds):
        compound_data = results[results['compound'] == compound]
        plt.scatter(compound_data['actual_logS'], 
                   compound_data['predicted_logS'], 
                   c=colors[i], 
                   label=compound,
                   alpha=0.7,
                   s=50)
    
    # 完美预测线
    plt.plot([-7, -1], [-7, -1], 'k--', label='Perfect Prediction', alpha=0.5)
    
    plt.xlabel('Actual logS')
    plt.ylabel('Predicted logS')
    plt.title('Model Prediction vs Actual Values')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('prediction_scatter.png', dpi=300, bbox_inches='tight')
    print("散点图已保存到 prediction_scatter.png")
    
    # 误差条形图
    plt.figure(figsize=(10, 6))
    errors_by_compound = results.groupby('compound')['abs_error'].mean()
    plt.bar(range(len(errors_by_compound)), errors_by_compound.values)
    plt.xticks(range(len(errors_by_compound)), errors_by_compound.index, rotation=45)
    plt.ylabel('Mean Absolute Error (logS)')
    plt.title('Prediction Error by Compound')
    plt.tight_layout()
    plt.savefig('error_bar.png', dpi=300, bbox_inches='tight')
    print("误差条形图已保存到 error_bar.png")
    
except ImportError:
    print("\nmatplotlib未安装，跳过可视化")

print("\n分析完成！")
