"""分析遗传算法的运行结果：适应度演化、参数收敛与参数分布。

读取 results/genetic_algorithm_data.csv（C 版或 Python 版均可），输出图片到
figures/ 目录：
    fitness_evolution.png            最优/平均适应度随代数的变化
    parameter_convergence.png        每代最优个体的四个参数
    parameter_evolution_heatmap.png  参数在每一代中的取值分布（热图）
    parameter_distribution_3d.png    参数在每一代中的出现频率（三维散点）

运行（仓库根目录）:
    python src/python/analyse_results.py
    python src/python/analyse_results.py --csv results/ga_seed7.csv
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # 无图形界面环境下也能保存图片
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PARAMS = ["Y0", "B", "lambda0", "delta_lambda_L"]
PARAM_LABELS = [r"$Y_0$", r"$B$", r"$\lambda_0$", r"$\Delta\lambda_L$ (Å)"]


def load_csv(path):
    """读取 GA 输出 CSV，并做基本检查。"""
    df = pd.read_csv(path)
    missing = {"Generation", "Individual", "Fitness", *PARAMS} - set(df.columns)
    if missing:
        raise ValueError(f"CSV 缺少必要的列: {sorted(missing)}")
    return df


def plot_fitness_evolution(df, outdir):
    """最优与平均适应度随代数的变化（适应度越小越好）。"""
    mean_fitness = df.groupby("Generation")["Fitness"].mean()
    best_fitness = df.groupby("Generation")["Fitness"].min()

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.semilogy(best_fitness.index, best_fitness.values, color="tab:red",
                linewidth=1.6, label="Best individual")
    ax.semilogy(mean_fitness.index, mean_fitness.values, color="tab:blue",
                linewidth=1.0, alpha=0.85, label="Population mean")
    ax.set_xlabel("Generation")
    ax.set_ylabel("Fitness (mean relative error)")
    ax.set_title("Fitness evolution over generations")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "fitness_evolution.png", dpi=150)
    plt.close(fig)


def plot_parameter_convergence(df, outdir):
    """每一代最优个体的四个参数随代数的变化。"""
    best_rows = df.loc[df.groupby("Generation")["Fitness"].idxmin()]

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    for ax, name, label in zip(axes.ravel(), PARAMS, PARAM_LABELS):
        ax.plot(best_rows["Generation"], best_rows[name], color="tab:green",
                marker=".", markersize=3, linewidth=1.0)
        ax.set_xlabel("Generation")
        ax.set_ylabel(label)
        ax.set_title(f"Best {name} over generations")
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / "parameter_convergence.png", dpi=150)
    plt.close(fig)


def plot_parameter_spread(df, outdir):
    """每一代种群中参数的分布范围（最小~最大）与均值：展示收敛过程。"""
    grouped = df.groupby("Generation")
    gen = np.array(sorted(df["Generation"].unique()))

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    for ax, name, label in zip(axes.ravel(), PARAMS, PARAM_LABELS):
        low = grouped[name].min().reindex(gen).to_numpy()
        high = grouped[name].max().reindex(gen).to_numpy()
        mean = grouped[name].mean().reindex(gen).to_numpy()
        x = gen + 1  # 便于在 log 轴上显示第 0 代

        ax.fill_between(x, low, high, color="tab:blue", alpha=0.25,
                        label="Population range (min–max)")
        ax.plot(x, mean, color="tab:blue", linewidth=1.1, label="Population mean")
        ax.set_xscale("log")
        ax.set_xlabel("Generation (log scale)")
        ax.set_ylabel(label)
        ax.set_title(f"Convergence of {name}")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(outdir / "parameter_spread.png", dpi=150)
    plt.close(fig)


def plot_parameter_distribution_3d(df, outdir, nbins=40):
    """三维散点：X=参数值，Y=代数，Z=该区间内个体数（对应报告中的图 12/13）。"""
    generations = np.sort(df["Generation"].unique())
    fig = plt.figure(figsize=(13, 10))

    for subplot, (name, label) in enumerate(zip(PARAMS, PARAM_LABELS), start=1):
        ax = fig.add_subplot(2, 2, subplot, projection="3d")
        values = df[name].to_numpy()
        lo, hi = np.percentile(values, [0.5, 99.5])
        bins = np.linspace(lo, hi, nbins + 1)
        centers = 0.5 * (bins[:-1] + bins[1:])

        xs, ys, zs = [], [], []
        for gen in generations:
            vals = df.loc[df["Generation"] == gen, name].to_numpy()
            counts, _ = np.histogram(vals, bins=bins)
            for x, count in zip(centers, counts):
                if count > 2:  # 只显示出现次数大于 2 的区间，避免亮点被零散点淹没
                    xs.append(x)
                    ys.append(gen)
                    zs.append(count)

        ax.scatter(xs, ys, zs, c=zs, cmap="viridis", s=8, alpha=0.8, edgecolors="none")
        ax.set_xlabel(label)
        ax.set_ylabel("Generation")
        ax.set_zlabel("Frequency")
        ax.set_title(f"{name}: frequency across generations", fontsize=10)

    fig.tight_layout()
    fig.savefig(outdir / "parameter_distribution_3d.png", dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="分析遗传算法输出的参数演化数据")
    parser.add_argument("--csv", type=Path,
                        default=ROOT / "results" / "genetic_algorithm_data.csv",
                        help="GA 输出 CSV（默认 results/genetic_algorithm_data.csv）")
    parser.add_argument("--outdir", type=Path, default=ROOT / "figures",
                        help="图片输出目录（默认 figures/）")
    args = parser.parse_args()

    df = load_csv(args.csv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    print(f"Loaded {len(df)} rows from {args.csv}")
    print(f"  generations: {df['Generation'].nunique()}, "
          f"individuals per generation: {df.groupby('Generation').size().max()}")

    plot_fitness_evolution(df, args.outdir)
    plot_parameter_convergence(df, args.outdir)
    plot_parameter_spread(df, args.outdir)
    plot_parameter_distribution_3d(df, args.outdir)

    best = df.loc[df["Fitness"].idxmin()]
    print("\nBest individual in the whole run:")
    for name in PARAMS:
        print(f"  {name:<14} = {best[name]:.5f}")
    print(f"  {'Fitness':<14} = {best['Fitness']:.6f}")
    print(f"\nFigures written to {args.outdir}")


if __name__ == "__main__":
    main()
