"""用遗传算法得到的最优参数绘制洛伦兹拟合曲线，并评估拟合质量。

功能：
    1. 读取实验数据（前三列：波长、实测强度、理论参考曲线）；
    2. 从 GA 输出 CSV 中取出最优个体（或直接由 --params 指定）；
    3. 计算拟合曲线、残差与误差指标（MAE / RMSE / MRE / NMSE），
       并换算电子密度 Ne；
    4. 作为对照，对数据第三列（理论曲线）做最小二乘拟合，得到"真值"参数；
    5. 输出 figures/data_overview.png 与 figures/fit_lorentzian.png，
       以及文字摘要 results/run_summary.txt。

运行（仓库根目录）:
    python src/python/fit_lorentzian.py
    python src/python/fit_lorentzian.py --params 0.21992,13.90698,6562.07706,10.25711
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parents[2]
PARAM_NAMES = ["Y0", "B", "lambda0", "delta_lambda_L"]


def load_data(path):
    """读取数据文件，返回 (波长, 实测强度, 理论参考曲线或 None)。"""
    data = np.loadtxt(path)
    wavelength = data[:, 0]
    intensity = data[:, 1]
    reference = data[:, 2] if data.shape[1] >= 3 else None
    return wavelength, intensity, reference


def lorentzian(wavelength, y0, b, lam0, dlam):
    """洛伦兹函数。"""
    return y0 + b / ((wavelength - lam0) ** 2 + (dlam / 2.0) ** 2)


def metrics(model, intensity):
    """常用误差指标。"""
    residual = model - intensity
    mae = np.mean(np.abs(residual))
    rmse = np.sqrt(np.mean(residual ** 2))
    mre = np.mean(np.abs(residual) / np.abs(intensity))
    nmse = np.sum(residual ** 2) / np.sum((intensity - intensity.mean()) ** 2)
    return {"MAE": mae, "RMSE": rmse, "MRE": mre, "NMSE": nmse}


def electron_density(dlam):
    """由 ΔλL(Å) 换算电子密度 Ne(m^-3)。"""
    return 10.0 ** (21.78515 + 1.2614 * np.log10(abs(dlam)))


def best_from_csv(csv_path):
    """从 GA 输出 CSV 中取出整轮最优个体。"""
    import pandas as pd

    df = pd.read_csv(csv_path)
    best = df.loc[df["Fitness"].idxmin()]
    return np.array([best[name] for name in PARAM_NAMES], dtype=float), best["Fitness"]


def plot_data_overview(wavelength, intensity, reference, outpath):
    """原始数据总览：实测谱线与理论参考曲线。"""
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(wavelength, intensity, color="tab:blue", linewidth=1.0,
            label="Measured data")
    if reference is not None:
        ax.plot(wavelength, reference, color="black", linewidth=1.2, linestyle="--",
                label="Reference curve (column 3)")
    ax.set_xlabel(r"Wavelength $\lambda$ (Å)")
    ax.set_ylabel("Intensity (a.u.)")
    ax.set_title(r"$H_\alpha$ line profile of an air-break arc plasma")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def plot_fit(wavelength, intensity, reference, params, metrics_result, outpath):
    """拟合曲线与残差。"""
    model = lorentzian(wavelength, *params)
    residual = model - intensity

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1]})
    ax = axes[0]
    ax.plot(wavelength, intensity, color="tab:blue", linewidth=1.0,
            label="Measured data")
    ax.plot(wavelength, model, color="tab:red", linewidth=1.6,
            label="Lorentzian fit (GA)")
    if reference is not None:
        ax.plot(wavelength, reference, color="black", linewidth=1.0,
                linestyle="--", alpha=0.75, label="Reference curve")
    ax.set_ylabel("Intensity (a.u.)")
    ax.set_title("Lorentzian fit of the $H_\\alpha$ line profile")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")

    ax = axes[1]
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.plot(wavelength, residual, color="tab:orange", linewidth=0.9)
    ax.set_xlabel(r"Wavelength $\lambda$ (Å)")
    ax.set_ylabel("Residual")
    ax.set_title(f"Residuals  (RMSE = {metrics_result['RMSE']:.4f})", fontsize=10)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="绘制洛伦兹拟合曲线并评估拟合质量")
    parser.add_argument("--data", type=Path,
                        default=ROOT / "data" / "Data_Genetique.txt",
                        help="实验数据文件")
    parser.add_argument("--csv", type=Path,
                        default=ROOT / "results" / "genetic_algorithm_data.csv",
                        help="GA 输出 CSV（用于读取最优参数）")
    parser.add_argument("--params", type=str, default=None,
                        help="直接指定参数，格式 Y0,B,lambda0,delta_lambda_L；"
                             "缺省时从 CSV 中读取最优个体")
    parser.add_argument("--outdir", type=Path, default=ROOT / "figures",
                        help="图片输出目录")
    parser.add_argument("--results-dir", type=Path, default=ROOT / "results",
                        help="文字摘要输出目录")
    args = parser.parse_args()

    wavelength, intensity, reference = load_data(args.data)

    if args.params:
        params = np.array([float(v) for v in args.params.split(",")])
        best_fitness = None
    else:
        if not args.csv.exists():
            raise SystemExit(
                f"找不到 {args.csv}；请先运行 C 版或 Python 版遗传算法，"
                f"或用 --params 直接给出参数。"
            )
        params, best_fitness = best_from_csv(args.csv)

    args.outdir.mkdir(parents=True, exist_ok=True)
    args.results_dir.mkdir(parents=True, exist_ok=True)

    model = lorentzian(wavelength, *params)
    m = metrics(model, intensity)
    ne = electron_density(params[3])

    # 与数据第三列（理论曲线）做最小二乘拟合，作为"真值"参照
    reference_params = None
    if reference is not None:
        p0 = [0.22, 12.0, 6562.0, 9.5]
        reference_params, _ = curve_fit(lorentzian, wavelength, reference, p0=p0, maxfev=20000)

    plot_data_overview(wavelength, intensity, reference, args.outdir / "data_overview.png")
    plot_fit(wavelength, intensity, reference, params, m, args.outdir / "fit_lorentzian.png")

    # ------------------------------ 文字输出 ------------------------------
    lines = []
    lines.append("=" * 66)
    lines.append("Lorentzian fit of the H-alpha line profile — summary")
    lines.append("=" * 66)
    lines.append(f"data file : {args.data}")
    lines.append(f"points    : {len(wavelength)}")
    if best_fitness is not None:
        lines.append(f"fitness   : {best_fitness:.6f}  (mean relative error)")
    lines.append("")
    lines.append("Fitted parameters (genetic algorithm):")
    for name, value in zip(PARAM_NAMES, params):
        lines.append(f"  {name:<14} = {value:.5f}")
    lines.append(f"  electron density Ne = {ne:.5e} m^-3")
    lines.append("")
    lines.append("Error metrics on the measured data:")
    for key, value in m.items():
        lines.append(f"  {key:<5} = {value:.6f}")
    if reference_params is not None:
        ne_ref = electron_density(reference_params[3])
        lines.append("")
        lines.append("Reference (least-squares fit of the noise-free 3rd column):")
        for name, value in zip(PARAM_NAMES, reference_params):
            lines.append(f"  {name:<14} = {value:.5f}")
        lines.append(f"  electron density Ne = {ne_ref:.5e} m^-3")
        lines.append("")
        lines.append("Deviation of the GA fit from the reference curve:")
        for name, value, ref in zip(PARAM_NAMES, params, reference_params):
            lines.append(f"  {name:<14} = {100.0 * (value - ref) / ref:+.2f} %")
        lines.append(f"  Ne             = {100.0 * (ne - ne_ref) / ne_ref:+.2f} %")

    summary = "\n".join(lines)
    print(summary)
    (args.results_dir / "run_summary.txt").write_text(summary + "\n", encoding="utf-8")
    print(f"\nFigures written to {args.outdir}")
    print(f"Summary written to {args.results_dir / 'run_summary.txt'}")


if __name__ == "__main__":
    main()
