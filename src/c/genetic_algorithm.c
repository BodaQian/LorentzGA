/* ============================================================================
 * Projet numérique niveau 2 — Ajustement d'un profil spectral Hα par une
 * fonction lorentzienne à l'aide d'un algorithme génétique.
 *
 * 项目：基于遗传算法的 Hα 谱线洛伦兹拟合（等离子体电子密度反演）
 * 作者：QIAN Boda (22016787) ; QIN Chaofan (22016788)
 * 单位：Université Toulouse III — Paul Sabatier / Northeastern University
 *
 * 构建（Linux / macOS / MinGW）:
 *     gcc -O2 -o genetic_algorithm genetic_algorithm.c -lm
 * 构建（Windows + MSVC，在 "x64 Native Tools Command Prompt" 中）:
 *     cl /O2 genetic_algorithm.c
 *
 * 运行（在仓库根目录）:
 *     ./genetic_algorithm [数据文件] [输出CSV]
 *   默认: data/Data_Genetique.txt  ->  genetic_algorithm_data.csv
 *
 * 数据文件说明：每行 "波长 强度 [参考曲线]"（空格分隔）。
 * 本程序只使用每行的前两列作为 (λ, I) 实验数据，第三列（若有）忽略。
 * ==========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// 定义种群参数
// Définir les paramètres de population
#define POP_SIZE 100      // 种群大小  Taille des populations
#define GENERATIONS 500   // 迭代代数  L'algèbre itérative
#define MUT_PROB 0.1      // 变异概率   Probabilité de variation
#define PARAMS 4          // 参数个数（Y0, B, lambda0, delta_lambda_L）

// 定义参数范围（由实验数据的形态估计得到）
// Définir la gamme des paramètres (estimée d'après l'allure des données)
#define Y0_MIN 0
#define Y0_MAX 0.25
#define B_MIN 1
#define B_MAX 100
#define LAMBDA0_MIN 6550
#define LAMBDA0_MAX 6575
#define DELTA_LAMBDA_L_MIN 0
#define DELTA_LAMBDA_L_MAX 10

// 锦标赛选择的大小
// Choix de la taille du tournoi
#define TOURNAMENT_SIZE 5

// 固定随机种子，保证结果可复现（改变种子即可得到不同演化轨迹）
// Graine fixe pour rendre les résultats reproductibles
#define RANDOM_SEED 12345

// 默认文件路径（可在命令行覆盖）
#define DEFAULT_DATA_PATH "data/Data_Genetique.txt"
#define DEFAULT_OUTPUT_PATH "genetic_algorithm_data.csv"

/* ----------------------------------------------------------------------------
 * 读取实验数据：逐行读取，只取每行前两个数（波长、强度）
 * 注意：不能直接用 fscanf("%lf %lf") 连续读取——当数据文件每行有三列时，
 *       会把第三列与下一行第一列配成一对，数据整体错位。
 * Lecture ligne par ligne : on ne garde que les deux premiers nombres de
 * chaque ligne, ce qui évite tout décalage si le fichier contient plus de
 * deux colonnes.
 * --------------------------------------------------------------------------*/
int load_data(const char* filename, double** wavelength, double** intensity) {

    FILE* file = fopen(filename, "r");
    if (!file) return -1;  // 文件不存在  Fichier introuvable

    int capacity = 2048;   // 动态数组初始容量  Capacité initiale du tableau
    int count = 0;
    double* w = (double*)malloc(capacity * sizeof(double));
    double* i = (double*)malloc(capacity * sizeof(double));
    if (!w || !i) { fclose(file); free(w); free(i); return -1; }

    char line[256];
    while (fgets(line, sizeof(line), file)) {
        double wave, inten;
        if (sscanf(line, "%lf %lf", &wave, &inten) != 2) continue;  // 跳过空行/注释行
        if (count == capacity) {
            capacity *= 2;
            double* nw = (double*)realloc(w, capacity * sizeof(double));
            if (!nw) { fclose(file); free(w); free(i); return -1; }
            w = nw;
            double* ni = (double*)realloc(i, capacity * sizeof(double));
            if (!ni) { fclose(file); free(w); free(i); return -1; }
            i = ni;
        }
        w[count] = wave;
        i[count] = inten;
        count++;
    }

    fclose(file);
    *wavelength = w;
    *intensity = i;
    return count;
}

// 洛伦兹函数
// Fonction de lorentz
double lorentzian(double wavelength, double Y0, double B, double lambda0, double delta_lambda_L) {
    return Y0 + (B / ((wavelength - lambda0) * (wavelength - lambda0) + (delta_lambda_L / 2) * (delta_lambda_L / 2)));
}

// 误差作为适应度函数
// Erreur en fonction de l'aptitude
double fitness_function(double* individual, double* wavelength, double* intensity, int data_count) {

    double Y0 = individual[0];
    double B = individual[1];
    double lambda0 = individual[2];
    double delta_lambda_L = individual[3];

    double error = 0.0;
    for (int k = 0; k < data_count; k++) {
        double model = lorentzian(wavelength[k], Y0, B, lambda0, delta_lambda_L);
        // 以下四种误差函数曾在报告中对比，最终选用相对误差：
        // error += pow(model - intensity[k], 2);                          // MSE
        // error += pow(model - intensity[k], 2) / pow(intensity[k], 2);   // NMSE
        // error += fabs(model - intensity[k]);                            // MAE
        error += fabs(model - intensity[k]) / fabs(intensity[k]);          // MRE (相对误差)
    }
    return error / data_count;
}

// 初始化个体
// Initialiser un individu
void initialize_individual(double* individual, int i) { // population[i][0,1,2,3]

    individual[0] = Y0_MIN + (Y0_MAX - Y0_MIN) * ((double)rand() / RAND_MAX);  // 控制取值在 0 到 1 之间  dans l'intervalle [0,1]
    individual[1] = B_MIN + (B_MAX - B_MIN) * ((double)rand() / RAND_MAX);
    individual[2] = LAMBDA0_MIN + (LAMBDA0_MAX - LAMBDA0_MIN) * ((double)rand() / RAND_MAX);
    individual[3] = DELTA_LAMBDA_L_MIN + (DELTA_LAMBDA_L_MAX - DELTA_LAMBDA_L_MIN) * ((double)rand() / RAND_MAX);

    // 打印初代每一个个体的值
    // Imprimer les valeurs de chaque individu de la première génération
    printf("Initialized individual %d: Y0 = %.5f, B = %.5f, lambda0 = %.5f, delta_lambda_L = %.5f\n",
            i, individual[0], individual[1], individual[2], individual[3]);
}

// 交叉操作（算术交叉 / croisement barycentrique）
// Croiser les opérations
void crossover(double* parent1, double* parent2, double* child1, double* child2) {

    for (int j = 0; j < PARAMS; j++) {
        double alpha = ((double)rand() / RAND_MAX); // 0 - 1
        child1[j] = alpha * parent1[j] + (1 - alpha) * parent2[j];
        child2[j] = alpha * parent2[j] + (1 - alpha) * parent1[j];
    }
}

// 变异操作：以 MUT_PROB 的概率给某个基因叠加一个小区间内的随机扰动
// Opération de variation
void mutate(double* individual) {

    for (int j = 0; j < PARAMS; j++) {
        if ((double)rand() / RAND_MAX < MUT_PROB) {

            if (j == 0) {          // Y0
                double delta = ((double)rand() / RAND_MAX) * 0.2 - 0.1;
                individual[j] += delta;
            }
            else if (j == 1) {     // B
                double delta = ((double)rand() / RAND_MAX) * 10 - 5;
                individual[j] += delta;
            }
            else if (j == 2) {     // lambda0
                double delta = ((double)rand() / RAND_MAX) * 10 - 5;
                individual[j] += delta;
            }
            else {                 // delta_lambda_L
                double delta = ((double)rand() / RAND_MAX) * 6 - 3;
                individual[j] += delta;
            }
        }
    }
}

// 锦标赛选择：随机抽 TOURNAMENT_SIZE 个个体，取适应度最小（误差最小）者
// Choix du tournoi
int tournament_selection(double population[POP_SIZE][PARAMS], double fitness[POP_SIZE]) {

    int best = rand() % POP_SIZE;
    for (int t = 1; t < TOURNAMENT_SIZE; t++) {
        int challenger = rand() % POP_SIZE;
        if (fitness[challenger] < fitness[best]) {
            best = challenger;
        }
    }
    return best;
}

// 由 ΔλL 计算电子密度 Ne（单位 m^-3）
// Calculer la densité électronique
double calculate_electron_density(double delta_lambda_L) {
    // 洛伦兹轮廓的半高宽取绝对值，避免负宽度导致 log10 无定义
    return pow(10, 21.78515 + 1.2614 * log10(fabs(delta_lambda_L)));
}

// 在最优参数上计算若干误差指标，便于与报告中的数值对照
// Calculer quelques indicateurs d'erreur pour le meilleur individu
void report_metrics(double* individual, double* wavelength, double* intensity, int data_count) {

    double mae = 0.0, sse = 0.0, mre = 0.0, mean_intensity = 0.0;
    for (int k = 0; k < data_count; k++) {
        double model = lorentzian(wavelength[k], individual[0], individual[1], individual[2], individual[3]);
        double e = model - intensity[k];
        mae += fabs(e);
        sse += e * e;
        mre += fabs(e) / fabs(intensity[k]);
        mean_intensity += intensity[k];
    }
    mean_intensity /= data_count;

    // NMSE = SSE / SST
    double sst = 0.0;
    for (int k = 0; k < data_count; k++) {
        double d = intensity[k] - mean_intensity;
        sst += d * d;
    }

    printf("Metrics (on the %d data points):\n", data_count);
    printf("  MAE  = %.6f\n", mae / data_count);
    printf("  RMSE = %.6f\n", sqrt(sse / data_count));
    printf("  MRE  = %.6f  (%.2f %%)\n", mre / data_count, 100.0 * mre / data_count);
    printf("  NMSE = %.6f\n", (sst > 0.0) ? sse / sst : 0.0);
}

// 更新"历史最优个体"
// Mettre à jour le meilleur individu rencontré depuis le début
void update_best_individual(double* individual, double fit,
                            double* best, double* best_fitness) {
    if (fit < *best_fitness) {
        for (int j = 0; j < PARAMS; j++) best[j] = individual[j];
        *best_fitness = fit;
    }
}

// Le programme principal
int main(int argc, char** argv) {

    const char* data_path = (argc > 1) ? argv[1] : DEFAULT_DATA_PATH;
    const char* output_path = (argc > 2) ? argv[2] : DEFAULT_OUTPUT_PATH;

    // 固定种子保证可复现；若要每次不同，可改用 srand((unsigned)time(NULL));
    srand(RANDOM_SEED);

    // 1. 读取实验数据
    // 1. Lire les données expérimentales
    double* wavelength = NULL;
    double* intensity = NULL;
    int data_count = load_data(data_path, &wavelength, &intensity);
    if (data_count < 0) {
        printf("Unable to read file: %s\n", data_path);
        return -1;
    }
    printf("Loaded %d data points from %s\n", data_count, data_path);

    // 2. 初始化种群
    // 2. Initialiser la population
    double population[POP_SIZE][PARAMS];
    double fitness[POP_SIZE]; // Aptitude
    for (int i = 0; i < POP_SIZE; i++) {
        initialize_individual(population[i], i);
        fitness[i] = fitness_function(population[i], wavelength, intensity, data_count);
    }

    // 保留实验过程中所有个体
    // Conservation de tous les individus au cours de l'expérience
    FILE* output_file = fopen(output_path, "w");
    if (output_file == NULL) {
        printf("Unable to open file for writing: %s\n", output_path);
        free(wavelength);
        free(intensity);
        return -1;
    }
    fprintf(output_file, "Generation,Individual,Y0,B,lambda0,delta_lambda_L,Fitness\n");

    // 记录整轮运行中出现过的最优个体（不只限于最后一代）
    // Mémoriser le meilleur individu rencontré au cours de toute l'évolution
    double best_ever[PARAMS];
    double best_ever_fitness = 1e30;
    for (int i = 0; i < POP_SIZE; i++) {
        update_best_individual(population[i], fitness[i], best_ever, &best_ever_fitness);
    }

    // 3. 遗传算法主循环
    // 3. Boucle principale d'algorithme génétique
    for (int gen = 0; gen < GENERATIONS; gen++) {

        double new_population[POP_SIZE][PARAMS];  // 新一代个体  Une nouvelle génération d'individus

        // 选择、交叉和变异
        // sélection, croisement et mutation
        for (int i = 0; i < POP_SIZE; i += 2) {

            // 将每一代数据写入文件 genetic_algorithm_data.csv
            // Écrire les données de chaque génération dans le fichier CSV
            fprintf(output_file, "%d,%d,%.6f,%.6f,%.6f,%.6f,%.8f\n", gen, i,
                population[i][0], population[i][1], population[i][2], population[i][3], fitness[i]);
            fprintf(output_file, "%d,%d,%.6f,%.6f,%.6f,%.6f,%.8f\n", gen, i + 1,
                population[i + 1][0], population[i + 1][1], population[i + 1][2], population[i + 1][3], fitness[i + 1]);

            update_best_individual(population[i], fitness[i], best_ever, &best_ever_fitness);
            update_best_individual(population[i + 1], fitness[i + 1], best_ever, &best_ever_fitness);

            // 锦标赛选择，得到父代
            // Sélection de tournoi, obtenir les parents
            int p1 = tournament_selection(population, fitness);
            int p2 = tournament_selection(population, fitness);

            double child1[PARAMS], child2[PARAMS];
            crossover(population[p1], population[p2], child1, child2);

            mutate(child1);
            mutate(child2);

            for (int j = 0; j < PARAMS; j++) {
                new_population[i][j] = child1[j];
                new_population[i + 1][j] = child2[j];
            }
        }

        // 更新种群和适应度
        // Renouvellement des populations et fitness
        for (int i = 0; i < POP_SIZE; i++) {
            for (int j = 0; j < PARAMS; j++) {
                population[i][j] = new_population[i][j];
            }
            fitness[i] = fitness_function(population[i], wavelength, intensity, data_count);
        }
    }

    // 4. 记录最终一代（第 GENERATIONS 代），并更新历史最优
    // 4. Écrire la dernière génération et mettre à jour le meilleur individu
    for (int i = 0; i < POP_SIZE; i++) {
        fprintf(output_file, "%d,%d,%.6f,%.6f,%.6f,%.6f,%.8f\n", GENERATIONS, i,
            population[i][0], population[i][1], population[i][2], population[i][3], fitness[i]);
        update_best_individual(population[i], fitness[i], best_ever, &best_ever_fitness);
    }
    fclose(output_file);

    // 5. 输出整轮最优个体并计算电子密度
    // 5. Afficher le meilleur individu et calculer la densité électronique
    printf("\nBest individual over the whole run:\n");
    printf("  Y0              = %.5f\n", best_ever[0]);
    printf("  B               = %.5f\n", best_ever[1]);
    printf("  lambda0         = %.5f A\n", best_ever[2]);
    printf("  delta_lambda_L  = %.5f A\n", best_ever[3]);
    printf("  fitness (MRE)   = %.6f\n", best_ever_fitness);

    double electron_density = calculate_electron_density(best_ever[3]);
    printf("Calculated electron density (m^-3): %.5e\n", electron_density);

    report_metrics(best_ever, wavelength, intensity, data_count);

    // 释放内存
    // Libérez votre mémoire
    free(wavelength);
    free(intensity);

    return 0;
}
