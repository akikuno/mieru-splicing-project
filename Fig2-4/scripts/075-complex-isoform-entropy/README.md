
Q. 複合体を形成する遺伝子に対して多くのスプライシングイベントの変動を受けていることがわかったが、転写産物の多様性はどのように変動しているのか（増えているのか、減っているのか、はたまた変わらないのか？）

- RefSeqIDの転写産物にアラインメントする
- Controlの転写産物の発現量と比較して、転写産物の多様性を評価する

```python
import numpy as np

def bootstrap_statistic(data, func, n_bootstrap=100):
    """ブートストラップ法による統計量の分布を推定"""
    bootstrapped_stats = []
    for _ in range(n_bootstrap):
        resampled_data = {k: np.random.choice(v, size=len(v), replace=True) for k, v in data.items()}
        bootstrapped_stats.append(func(resampled_data))
    return np.array(bootstrapped_stats)

def calculate_entropy(data):
    total_sum = sum([sum(values) for values in data.values()])
    proportions = [sum(values) / total_sum for values in data.values()]
    return -np.sum([p * np.log(p) for p in proportions])

# データ定義
sample1 = {'A': [10, 10, 10], 'B': [10, 10, 10], 'C': [30]}
sample2 = {'A': [100], 'B': [100], 'C': [150, 150]}

# エントロピーの分布を推定
entropy1_dist = bootstrap_statistic(sample1, calculate_entropy)
entropy2_dist = bootstrap_statistic(sample2, calculate_entropy)

# 信頼区間とp値の計算
mean_diff = np.mean(entropy1_dist) - np.mean(entropy2_dist)
p_value = np.sum(entropy1_dist < np.mean(entropy2_dist)) / len(entropy1_dist)

print("Sample 1 Entropy CI:", np.percentile(entropy1_dist, [2.5, 97.5]))
print("Sample 2 Entropy CI:", np.percentile(entropy2_dist, [2.5, 97.5]))
print("Mean Difference:", mean_diff)
print("P-value:", p_value)
```
