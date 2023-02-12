import numpy as np
import random

def target_distribution(x):
    # Define your target distribution here
    return np.exp(-x**2 / 2)

def metropolis_hastings(x, temperature):
    # Implement the Metropolis-Hastings algorithm here
    proposal = x + np.random.normal(0, 1)
    acceptance_prob = min(1, target_distribution(proposal) / (target_distribution(x) * temperature))
    if np.random.uniform(0, 1) < acceptance_prob:
        return proposal
    else:
        return x

def parallel_tempering(num_chains, num_iterations, temperature_schedule):
    chains = [np.zeros(num_iterations) for i in range(num_chains)]
    x = [0 for i in range(num_chains)]

    #初期化(全てのchainの一番目にランダムなxを代入している)
    for i in range(num_chains):
        chains[i][0] = x[i] = np.random.normal(0, 1)

    #MCMCを回している
    for i in range(1, num_iterations):
        for j in range(num_chains):
            x[j] = metropolis_hastings(x[j], temperature_schedule[j])
            #メトロポリスヘイスティングス法によって決定されたxをチェーンに代入している
            #並列で回している訳ではない！
            chains[j][i] = x[j]
        if i % 10 == 0:
            #10の倍数のiteration時だけ、swappingを行う
            #random.sample(range(num_chains), 2)でチェーンの数以内で二つの値を含むリストを返している
            swap_indices = random.sample(range(num_chains), 2)
            i1, i2 = swap_indices[0], swap_indices[1]
            #コストの差をとってあげる
            energy_diff = target_distribution(x[i2]) / target_distribution(x[i1])
            #コストの差*温度の比の確率でスワップ
            acceptance_prob = min(1, energy_diff * temperature_schedule[i1] / temperature_schedule[i2])
            #
            if np.random.uniform(0, 1) < acceptance_prob:
                x[i1], x[i2] = x[i2], x[i1]

    return chains

num_chains = 8
num_iterations = 1000
temperature_schedule = np.linspace(0.1, 1.0, num_chains)
chains = parallel_tempering(num_chains, num_iterations, temperature_schedule)
