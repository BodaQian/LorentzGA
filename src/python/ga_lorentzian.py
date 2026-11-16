"""基于遗传算法的洛伦兹谱线拟合 —— 纯 Python（numpy）实现。

本脚本与 src/c/genetic_algorithm.c 采用完全相同的算法设置：
    * 种群大小 100，迭代 500 代
    * 适应度 = 逐点相对误差的平均值（越小越好）
    * 锦标赛选择（规模 5）、算术交叉、逐基因变异（概率 0.1）
    * 参数搜索范围与 C 版一致
区别仅在于随机数序列不同，因此两次运行的最优参数会有微小差异。

运行（仓库根目录）:
    python src/python/ga_lorentzian.py
    python src/python/ga_lorentzian.py --seed 7 --out results/ga_seed7.csv

输出:
    results/genetic_algorithm_data.csv  （每一代每个个体的参数与适应度）
"""

import argparse
import time
from pathlib import Path

import numpy as np

# 仓库根目录：本文件位于 <repo>/src/python/ 下
ROOT = Path(__file__).resolve().parents[2]

# ------------------------- 算法参数（与 C 版一致） -------------------------
POP_SIZE = 100            # 种群大小
GENERATIONS = 500         # 迭代代数
MUT_PROB = 0.1            # 变异概率
PARAMS = 4                # 基因个数：Y0, B, lambda0, delta_lambda_L
TOURNAMENT_SIZE = 5       # 锦标赛选择的规模
DEFAULT_SEED = 12345      # 默认随机种子（与 C 版相同，便于对照）

# 每个参数的取值范围（由实验数据的形态估计得到）
BOUNDS = [
    (0.0, 0.25),          # Y0
    (1.0, 100.0),         # B
    (6550.0, 6575.0),     # lambda0 (Angstrom)
    (0.0, 10.0),          # delta_lambda_L (Angstrom)
]

# 每个参数的变异幅度（均匀分布，单位为参数自身单位）
MUT_AMPLITUDE = [
    (-0.1, 0.1),
    (-5.0, 5.0),
    (-5.0, 5.0),
    (-3.0, 3.0),
]


def load_data(path):
    """读取实验数据，只使用每行的前两列（波长、强度）。

    数据文件的第三列（若存在）是噪声-free 的参考曲线，这里不使用。
    """
    data = np.loadtxt(path, usecols=(0, 1))
    wavelength = data[:, 0]
    intensity = data[:, 1]
    return wavelength, intensity


def lorentzian(wavelength, individual):
    """洛伦兹函数 J(lambda) = Y0 + B / ((lambda-lambda0)^2 + (dlambda_L/2)^2)。"""
    y0, b, lam0, dlam = individual
    return y0 + b / ((wavelength - lam0) ** 2 + (dlam / 2.0) ** 2)


def fitness_population(population, wavelength, intensity):
    """向量化计算整个种群的适应度（平均相对误差）。

    population: (POP_SIZE, PARAMS)
    返回: (POP_SIZE,) 越小越好
    """
    y0 = population[:, 0:1]
    b = population[:, 1:2]
    lam0 = population[:, 2:3]
    dlam = population[:, 3:4]
    model = y0 + b / ((wavelength[None, :] - lam0) ** 2 + (dlam / 2.0) ** 2)
    return np.mean(np.abs(model - intensity[None, :]) / np.abs(intensity[None, :]), axis=1)


def initialize_population(rng):
    """在参数范围内均匀随机生成初代种群。"""
    return np.array([[rng.uniform(lo, hi) for lo, hi in BOUNDS] for _ in range(POP_SIZE)])


def tournament_selection(fitness, rng):
    """锦标赛选择：随机选 TOURNAMENT_SIZE 个个体，返回最优者下标。"""
    candidates = rng.integers(0, POP_SIZE, size=TOURNAMENT_SIZE)
    return int(candidates[np.argmin(fitness[candidates])])


def crossover(parent1, parent2, rng):
    """算术交叉：每个基因独立取一个权重 alpha，生成两个子代。"""
    alpha = rng.random(PARAMS)
    child1 = alpha * parent1 + (1.0 - alpha) * parent2
    child2 = alpha * parent2 + (1.0 - alpha) * parent1
    return child1, child2


def mutate(individual, rng):
    """逐基因变异：每个基因以 MUT_PROB 的概率叠加一个均匀随机扰动。"""
    for j in range(PARAMS):
        if rng.random() < MUT_PROB:
            lo, hi = MUT_AMPLITUDE[j]
            individual[j] += rng.uniform(lo, hi)


def run_ga(wavelength, intensity, seed=DEFAULT_SEED, generations=GENERATIONS, verbose=True):
    """执行遗传算法。

    返回 (整轮最优个体, 最优适应度, 每代种群列表, 每代适应度列表)。
    列表共 generations+1 项：第 0~generations 代，最后一项为最终种群。
    """
    rng = np.random.default_rng(seed)

    population = initialize_population(rng)
    fitness = fitness_population(population, wavelength, intensity)

    history = []          # 每一代的种群（POP_SIZE x PARAMS）
    fitness_history = []  # 每一代对应的适应度
    for gen in range(generations):
        history.append(population.copy())
        fitness_history.append(fitness.copy())

        new_population = np.empty_like(population)
        for i in range(0, POP_SIZE, 2):
            p1 = tournament_selection(fitness, rng)
            p2 = tournament_selection(fitness, rng)
            child1, child2 = crossover(population[p1], population[p2], rng)
            mutate(child1, rng)
            mutate(child2, rng)
            new_population[i] = child1
            new_population[i + 1] = child2

        population = new_population
        fitness = fitness_population(population, wavelength, intensity)

        if verbose and (gen % 50 == 0 or gen == generations - 1):
            print(f"Generation {gen:3d}: best MRE = {fitness.min():.6f}")

    # 记录最终一代（与 C 版输出一致：CSV 覆盖第 0~generations 代）
    history.append(population.copy())
    fitness_history.append(fitness.copy())

    # 在整轮运行中找最优个体
    all_fitness = np.concatenate(fitness_history)
    flat_idx = int(np.argmin(all_fitness))
    gen_idx, ind_idx = divmod(flat_idx, POP_SIZE)
    best = history[gen_idx][ind_idx].copy()
    return best, float(all_fitness[flat_idx]), history, fitness_history


def save_history_csv(history, fitness_history, path):
    """把每一代每个个体的参数与适应度写入 CSV（与 C 版输出格式一致）。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("Generation,Individual,Y0,B,lambda0,delta_lambda_L,Fitness\n")
        for gen, (population, fitness) in enumerate(zip(history, fitness_history)):
            for idx, ind in enumerate(population):
                f.write(
                    f"{gen},{idx},{ind[0]:.6f},{ind[1]:.6f},"
                    f"{ind[2]:.6f},{ind[3]:.6f},{fitness[idx]:.8f}\n"
                )


def calculate_electron_density(delta_lambda_L):
    """由洛伦兹半高宽 ΔλL(Å) 计算电子密度 Ne(m^-3)。"""
    return 10.0 ** (21.78515 + 1.2614 * np.log10(abs(delta_lambda_L)))


def main():
    parser = argparse.ArgumentParser(description="遗传算法拟合洛伦兹谱线（Python 版）")
    parser.add_argument("--data", type=Path, default=ROOT / "data" / "Data_Genetique.txt",
                        help="实验数据文件（默认 data/Data_Genetique.txt）")
    parser.add_argument("--out", type=Path, default=ROOT / "results" / "genetic_algorithm_data.csv",
                        help="输出 CSV 路径（默认 results/genetic_algorithm_data.csv）")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="随机种子")
    parser.add_argument("--generations", type=int, default=GENERATIONS, help="迭代代数")
    args = parser.parse_args()

    t0 = time.time()
    wavelength, intensity = load_data(args.data)
    print(f"Loaded {len(wavelength)} data points from {args.data}")

    best, best_fitness, history, fitness_history = run_ga(
        wavelength, intensity, seed=args.seed, generations=args.generations
    )

    save_history_csv(history, fitness_history, args.out)
    print(f"Population history written to {args.out}")

    print(f"\nBest individual over the whole run:")
    print(f"  Y0              = {best[0]:.5f}")
    print(f"  B               = {best[1]:.5f}")
    print(f"  lambda0         = {best[2]:.5f} A")
    print(f"  delta_lambda_L  = {best[3]:.5f} A")
    print(f"  fitness (MRE)   = {best_fitness:.6f}")
    print(f"Calculated electron density (m^-3): {calculate_electron_density(best[3]):.5e}")
    print(f"Elapsed: {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
