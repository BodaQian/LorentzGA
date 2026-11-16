# 整理前的分析 Notebook（存档）

本目录保存项目整理前的原始分析 notebook，仅用于保留历史记录和**当时的图表**，
不是当前推荐的分析流程。

| 项目 | 内容 |
| --- | --- |
| 文件 | `analyse_genetique_original.ipynb`（原文件名 `Data_Anas.ipynb`） |
| 内容 | 读取 C 版遗传算法的输出 CSV，绘制参数三维频率分布、2D 分布与多项式拟合、平均/最优适应度曲线、参数–适应度散点、每代最优参数、trisurf 曲面等 |
| 数据来源 | 当时 C 版运行输出的 `genetic_algorithm_data.csv`，其 `Fitness` 平台约为 **0.39** |
| 是否可直接运行 | 否。notebook 中的数据路径指向原作者的本地桌面目录，且所需的旧 CSV 未纳入仓库 |
| 图 | 全部为当时运行写入 notebook 的原始输出 |

## 为什么它和现在的 `analyse_genetique.ipynb` 不一样

整理仓库时发现：原 C 版用 `fscanf("%lf %lf", ...)` 连续读取数字，而数据文件每行有 3 列，
于是"第三列 + 下一行第一列"被配成一对，数据整体错位。原始 notebook 分析的就是这份错位数据
产生的输出，因此：

- 适应度曲线平台约为 **0.39**（错位数据中只有约 1/3 的点对是真实的"波长—强度"数据，其余是与参数关系很弱的错配点，把误差整体抬高了）；
- 参数收敛到 $Y_0\approx0.2247$、$B\approx12.3$、$\lambda_0\approx6561.9$ Å、$\Delta\lambda_L\approx9.48$ Å；
- 由 $\Delta\lambda_L$ 换算的电子密度偏低约 11 %。

修正数据读取后重新运行（见仓库根目录 [README](../../README.md) 的"关于原始版本与本次整理"），
适应度平台约为 **0.045**，最优参数也落在数据第三列（无噪声参考曲线）附近。

当前分析流程见 [`../analyse_genetique.ipynb`](../analyse_genetique.ipynb)。
