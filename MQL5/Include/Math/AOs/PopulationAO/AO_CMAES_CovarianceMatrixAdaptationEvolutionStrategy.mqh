//+——————————————————————————————————————————————————————————————————+
//|                                                       C_AO_CMAES |
//|                                  Copyright 2007-2025, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/18227

#include "#C_AO.mqh"

//————————————————————————————————————————————————————————————————————
class C_AO_CMAES : public C_AO
{
  public: //----------------------------------------------------------
  ~C_AO_CMAES () { }
  C_AO_CMAES ()
  {
    ao_name = "CMAES";
    ao_desc = "Covariance Matrix Adaptation Evolution Strategy";
    ao_link = "https://www.mql5.com/ru/articles/18227";

    // Параметры по умолчанию
    popSize         = 50;          // Размер популяции по умолчанию (lambda)
    mu              = 25;          // Количество родителей (половина популяции)
    learningRateC1  = 0.01;        // Скорость обучения для ранг-1 обновления
    learningRateCMu = 0.8;         // Скорость обучения для ранг-μ обновления
    stepSizeDamping = 0.6;         // Демпфирование для размера шага

    // Создание и инициализация массива параметров
    ArrayResize (params, 5);
    params [0].name = "popSize";         params [0].val = popSize;
    params [1].name = "mu";              params [1].val = mu;
    params [2].name = "learningRateC1";  params [2].val = learningRateC1;
    params [3].name = "learningRateCMu"; params [3].val = learningRateCMu;
    params [4].name = "stepSizeDamping"; params [4].val = stepSizeDamping;
  }

  void SetParams ()
  {
    popSize         = (int)params [0].val;
    mu              = (int)params [1].val;
    learningRateC1  = params      [2].val;
    learningRateCMu = params      [3].val;
    stepSizeDamping = params      [4].val;
  }

  bool Init (const double &rangeMinP  [],  // минимальные значения
             const double &rangeMaxP  [],  // максимальные значения
             const double &rangeStepP [],  // размер шага
             const int     epochsP = 0);

  void Moving   ();
  void Revision ();

  //------------------------------------------------------------------
  private: //---------------------------------------------------------
  // Специфичные параметры CMA-ES
  int mu;                  // Количество родителей (выбранных точек)
  double sigma;            // Размер шага
  double learningRateC1;   // Скорость обучения для ранг-1 обновления
  double learningRateCMu;  // Скорость обучения для ранг-μ обновления
  double stepSizeDamping;  // Фактор демпфирования для обновления размера шага

  // Специфичные структуры данных CMA-ES
  double weights   [];      // Веса рекомбинации
  double covMatrix [];      // Ковариационная матрица (хранится как одномерный массив)
  double B  [];             // Собственные векторы C
  double D  [];             // Собственные значения (квадратные корни) C
  double pc [];             // Путь эволюции для C
  double ps [];             // Сопряженный путь эволюции для sigma
  double mu_eff;            // Эффективная масса выбора дисперсии
  int    counteval;         // Счетчик вычислений функции с последнего разложения
  int    eigeneval;         // Счетчик генераций, когда выполнялось разложение
  double chiN;              // Ожидаемая норма N(0,I)
  int    eigenInterval;     // Интервал для разложения на собственные значения

  // Переменные для оптимизации производительности
  double cs;                // Скорость обучения для пути sigma
  double cc;                // Скорость обучения для пути ранг-1
  double damps;             // Демпфирование для sigma
  double hsig_threshold;    // Порог для функции Хевисайда

  // Вспомогательные массивы
  double y_vec     [];         // Вектор мутации
  double arindex   [];         // Массив индексов для сортировки
  double arfitness [];         // Массив фитнеса для сортировки
  double temp_vec  [];         // Временный вектор для матричных операций
  double invsqrtC_times_yw []; // Временное хранилище для C^(-1/2) * y_w

  // Переменные кэширования для Бокса-Мюллера
  double cached_normal;
  bool   has_cached;

  // Вспомогательные методы
  void   InitDistribution          ();
  void   UpdateDistribution        ();
  void   ComputeEigendecomposition ();
  double GetChiN                   ();
  void   SortPopulation            ();
  void   ComputeWeights            ();
  void   UpdateMean                ();
  bool   CheckPositiveDefinite     ();
  void   EnforcePositiveDefinite   ();
};
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
bool C_AO_CMAES::Init (const double &rangeMinP [],  // минимальные значения
                       const double &rangeMaxP [],  // максимальные значения
                       const double &rangeStepP [], // размер шага
                       const int     epochsP = 0)   // количество эпох
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //------------------------------------------------------------------
  // Инициализация кэширования Бокса-Мюллера
  has_cached = false;

  // Инициализация специфичных переменных CMA-ES
  sigma          = 0.3;  // Начальный размер шага (30% от диапазона поиска)

  // Вычисление эффективной массы выбора дисперсии
  ComputeWeights ();

  // Ожидаемая норма N(0,I)
  chiN = GetChiN ();

  // Вычисление и сохранение параметров стратегии
  cs = (mu_eff + 2.0) / (coords + mu_eff + 5.0);
  cc = (4.0 + mu_eff / coords) / (coords + 4.0 + 2.0 * mu_eff / coords);
  damps = 1.0 + 2.0 * MathMax (0.0, MathSqrt ((mu_eff - 1.0) / (coords + 1.0)) - 1.0) + cs;
  hsig_threshold = 1.4 + 2.0 / (coords + 1.0);

  // Установка интервала разложения на собственные значения - настройка для производительности
  eigenInterval = (int)(coords / (10.0 * MathSqrt (learningRateC1 + learningRateCMu)));
  eigenInterval = MathMax (1, eigenInterval);

  // Выделение массивов - выделяем только один раз
  ArrayResize (covMatrix, coords * coords);
  ArrayResize (B, coords * coords);
  ArrayResize (D, coords);
  ArrayResize (pc, coords);
  ArrayResize (ps, coords);
  ArrayResize (y_vec, coords);
  ArrayResize (arindex, popSize);
  ArrayResize (arfitness, popSize);
  ArrayResize (temp_vec, coords);
  ArrayResize (invsqrtC_times_yw, coords);

  // Инициализация путей эволюции нулями
  ArrayInitialize (pc, 0);
  ArrayInitialize (ps, 0);

  // Быстрая инициализация единичной ковариационной матрицы и разложения
  ArrayInitialize (covMatrix, 0.0);
  ArrayInitialize (B, 0.0);

  for (int i = 0; i < coords; i++)
  {
    covMatrix [i * coords + i] = 1.0;
    B [i * coords + i] = 1.0;
    D [i] = 1.0;
  }

  // Инициализация начального распределения
  InitDistribution ();

  // Сброс счетчиков вычислений
  counteval = 0;
  eigeneval = 0;

  // Принудительный пересчет фитнеса
  revision = true;

  return true;
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
void C_AO_CMAES::Moving ()
{
  // Генерация lambda новых потомков
  for (int k = 0; k < popSize; k++)
  {

    // Применение преобразования y = B*D*z
    for (int i = 0; i < coords; i++)
    {
      y_vec [i] = 0.0;
      for (int j = 0; j < coords; j++)
      {
        y_vec [i] += B [i * coords + j] * D [j] * u.PowerDistribution (0.0, -8, 8, 20);
      }

      // Генерация потомка: x_k = m + σ * y
      a [k].c [i] = cB [i] + sigma * y_vec [i];
      a [k].c [i] = u.SeInDiSp (a [k].c [i], rangeMin [i], rangeMax [i], rangeStep [i]);
    }
  }

  // Отметка для пересчета фитнеса
  revision = true;
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
void C_AO_CMAES::Revision ()
{
  if (!revision) return;

  revision = false;

  // Сортировка популяции по фитнесу
  SortPopulation ();

  // Обновление параметров распределения на основе выбранных особей
  UpdateDistribution ();

  // Обновление счетчика вычислений
  counteval++;
}
//————————————————————————————————————————————————————————————————————

//+------------------------------------------------------------------+
//| Инициализация распределения поиска                               |
//+------------------------------------------------------------------+
void C_AO_CMAES::InitDistribution ()
{
  // Установка начального среднего в центр пространства поиска
  for (int i = 0; i < coords; i++)
  {
    cB [i] = u.RNDfromCI (rangeMin [i], rangeMax [i]);
  }

  for (int k = 0; k < popSize; k++)
  {
    for (int i = 0; i < coords; i++)
    {
      // Генерация равномерно распределенной точки
      a [k].c [i] = u.RNDfromCI (rangeMin [i], rangeMax [i]);
      a [k].c [i] = u.SeInDiSp (a [k].c [i], rangeMin [i], rangeMax [i], rangeStep [i]);
    }
  }
}
//————————————————————————————————————————————————————————————————————

//+------------------------------------------------------------------+
//| Вычисление ожидаемой нормы N(0,I)                                |
//+------------------------------------------------------------------+
double C_AO_CMAES::GetChiN ()
{
  double n = (double)coords;
  return MathSqrt (n) * (1.0 - 1.0 / (4.0 * n) + 1.0 / (21.0 * n * n));
}
//————————————————————————————————————————————————————————————————————

//+------------------------------------------------------------------+
//| Сортировка популяции по фитнесу                                  |
//+------------------------------------------------------------------+
void C_AO_CMAES::SortPopulation ()
{
  // Копирование значений фитнеса и индексов
  for (int i = 0; i < popSize; i++)
  {
    arindex [i] = i;
    arfitness [i] = a [i].f;
  }

  for (int i = 1; i < popSize; i++)
  {
    double tempFitness = arfitness [i];
    double tempIndex = arindex [i];
    int j = i - 1;

    // Сортировка по убыванию (для максимизации)
    while (j >= 0 && arfitness [j] < tempFitness)
    {
      arfitness [j + 1] = arfitness [j];
      arindex [j + 1] = arindex [j];
      j--;
    }

    arfitness [j + 1] = tempFitness;
    arindex [j + 1] = tempIndex;
  }

  // Обновление лучшего решения при необходимости
  if (arfitness [0] > fB)
  {
    fB = arfitness [0];
    int bestIdx = (int)arindex [0];
    ArrayCopy (cB, a [bestIdx].c, 0, 0, coords);
  }
}
//————————————————————————————————————————————————————————————————————

//+------------------------------------------------------------------+
//| Обновление среднего с использованием взвешенной рекомбинации     |
//+------------------------------------------------------------------+
void C_AO_CMAES::UpdateMean ()
{
  // Взвешенная рекомбинация: m^(g+1) = Σ w_i * x_{i:λ}^(g+1)
  for (int j = 0; j < coords; j++)
  {
    double meanSum = 0.0;
    for (int i = 0; i < mu; i++)
    {
      int idx = (int)arindex [i];
      meanSum += weights [i] * a [idx].c [j];
    }
    cB [j] = meanSum;
  }
}
//————————————————————————————————————————————————————————————————————

//+------------------------------------------------------------------+
//| Вычисление весов взвешенной рекомбинации                         |
//+------------------------------------------------------------------+
void C_AO_CMAES::ComputeWeights ()
{
  // Выделение массива весов
  ArrayResize (weights, mu);

  // Предварительное вычисление log(mu + 0.5)
  double log_mu_plus_half = MathLog (mu + 0.5);

  // Вычисление положительных весов
  double sum = 0.0;
  for (int i = 0; i < mu; i++)
  {
    weights [i] = log_mu_plus_half - MathLog (i + 1);
    sum += weights [i];
  }

  // Нормализация весов
  double sum_weights = 0.0;
  double sum_squares = 0.0;
  for (int i = 0; i < mu; i++)
  {
    weights [i] /= sum;
    sum_weights += weights [i];
    sum_squares += weights [i] * weights [i];
  }

  // Вычисление эффективной массы выбора дисперсии
  mu_eff = sum_weights * sum_weights / sum_squares;
}
//————————————————————————————————————————————————————————————————————

//+------------------------------------------------------------------+
//| Обновление параметров распределения                              |
//+------------------------------------------------------------------+
void C_AO_CMAES::UpdateDistribution ()
{
  // Проверка необходимости разложения на собственные значения
  if (counteval - eigeneval > eigenInterval)
  {
    ComputeEigendecomposition ();
    eigeneval = counteval;
  }

  // Сохранение старого среднего
  double oldMean [];
  ArrayResize (oldMean, coords);
  ArrayCopy (oldMean, cB, 0, 0, coords);

  // Обновление среднего
  UpdateMean ();

  // Вычисление смещения среднего
  double y_w [];
  ArrayResize (y_w, coords);
  for (int j = 0; j < coords; j++)
  {
    y_w [j] = (cB [j] - oldMean [j]) / sigma;
  }

  // Вычисление C^(-1/2) * y_w
  // Шаг 1: B^T * y_w
  ArrayInitialize (temp_vec, 0.0);
  for (int i = 0; i < coords; i++)
  {
    for (int j = 0; j < coords; j++)
    {
      temp_vec [i] += B [j * coords + i] * y_w [j];
    }
  }

  // Шаг 2: D^(-1) * (B^T * y_w)
  for (int i = 0; i < coords; i++)
  {
    temp_vec [i] /= D [i];
  }

  // Шаг 3: B * D^(-1) * B^T * y_w
  ArrayInitialize (invsqrtC_times_yw, 0.0);
  for (int i = 0; i < coords; i++)
  {
    for (int j = 0; j < coords; j++)
    {
      invsqrtC_times_yw [i] += B [i * coords + j] * temp_vec [j];
    }
  }

  // Обновление пути эволюции для sigma
  double norm_ps_squared = 0.0;
  for (int i = 0; i < coords; i++)
  {
    ps [i] = (1.0 - cs) * ps [i] + MathSqrt (cs * (2.0 - cs) * mu_eff) * invsqrtC_times_yw [i];
    norm_ps_squared += ps [i] * ps [i];
  }

  // Функция Хевисайда
  double norm_ps = MathSqrt (norm_ps_squared);
  double expected_length = MathSqrt (1.0 - MathPow (1.0 - cs, 2.0 * counteval)) * chiN;
  bool hsig = norm_ps / expected_length < hsig_threshold;

  // Обновление пути эволюции для C
  double delta_hsig = hsig ? 1.0 : 0.0;
  for (int i = 0; i < coords; i++)
  {
    pc [i] = (1.0 - cc) * pc [i] + delta_hsig * MathSqrt (cc * (2.0 - cc) * mu_eff) * y_w [i];
  }

  // Подготовка ранг-1 обновления
  double c1a [];
  ArrayResize (c1a, coords * coords);
  for (int i = 0; i < coords; i++)
  {
    for (int j = 0; j <= i; j++)
    {
      c1a [i * coords + j] = c1a [j * coords + i] = pc [i] * pc [j];
    }
  }

  // Подготовка ранг-μ обновления
  double cmu [];
  ArrayResize (cmu, coords * coords);
  ArrayInitialize (cmu, 0.0);

  for (int k = 0; k < mu; k++)
  {
    int idx = (int)arindex [k];

    // Вычисление y_i = (x_i - m_old) / sigma
    for (int i = 0; i < coords; i++)
    {
      temp_vec [i] = (a [idx].c [i] - oldMean [i]) / sigma;
    }

    // Добавление взвешенного внешнего произведения
    for (int i = 0; i < coords; i++)
    {
      for (int j = 0; j <= i; j++)
      {
        double update = weights [k] * temp_vec [i] * temp_vec [j];
        cmu [i * coords + j] += update;
        if (i != j) cmu [j * coords + i] += update;
      }
    }
  }

  // Обновление ковариационной матрицы C
  double c1 = learningRateC1;
  double cmu_rate = learningRateCMu;

  // Корректировка c1 если hsig false (застой прогресса)
  if (!hsig)
  {
    c1 *= (1.0 - (1.0 - delta_hsig) * cc * (2.0 - cc));
  }

  double one_minus_c1_cmu = 1.0 - c1 - cmu_rate;

  // Обновление C с ранг-1 и ранг-μ обновлениями
  for (int i = 0; i < coords; i++)
  {
    for (int j = 0; j <= i; j++)
    {
      covMatrix [i * coords + j] = one_minus_c1_cmu * covMatrix [i * coords + j] +
                                   c1 * c1a [i * coords + j] +
                                   cmu_rate * cmu [i * coords + j];

      // Сохранение симметрии
      if (i != j)
      {
        covMatrix [j * coords + i] = covMatrix [i * coords + j];
      }
    }
  }

  // Обеспечение положительной определенности
  if (counteval % (10 * eigenInterval) == 0)
  {
    EnforcePositiveDefinite ();
  }

  // Обновление размера шага sigma
  double exponent = (cs / damps) * (norm_ps / chiN - 1.0);
  sigma *= MathExp (exponent);

  // Ограничение sigma для численной стабильности
  double min_sigma = 1e-16;
  double max_eigenvalue = D [0]; // D отсортирован по убыванию
  double max_sigma = 1e4 * MathMax (1.0, MathSqrt (max_eigenvalue));

  if (sigma < min_sigma) sigma = min_sigma;
  else
    if (sigma > max_sigma) sigma = max_sigma;
}
//————————————————————————————————————————————————————————————————————

//+------------------------------------------------------------------+
//| Вычисление разложения на собственные значения методом Якоби      |
//+------------------------------------------------------------------+
void C_AO_CMAES::ComputeEigendecomposition ()
{
  // Создание копии ковариационной матрицы для разложения
  double C_copy [];
  ArrayResize (C_copy, coords * coords);
  ArrayCopy (C_copy, covMatrix);

  // Инициализация B как единичной матрицы
  for (int i = 0; i < coords; i++)
  {
    for (int j = 0; j < coords; j++)
    {
      B [i * coords + j] = (i == j) ? 1.0 : 0.0;
    }
  }

  // Улучшенное разложение Якоби на собственные значения
  int max_iterations = 10; //50 * coords;
  double tolerance   = 0.01; //1e-14 * coords * coords;

  for (int iter = 0; iter < max_iterations; iter++)
  {
    // Поиск наибольшего внедиагонального элемента
    double max_val = 0.0;
    int p = 0, q = 1;

    for (int i = 0; i < coords - 1; i++)
    {
      for (int j = i + 1; j < coords; j++)
      {
        double val = MathAbs (C_copy [i * coords + j]);
        if (val > max_val)
        {
          max_val = val;
          p = i;
          q = j;
        }
      }
    }

    // Проверка сходимости
    if (max_val < tolerance) break;

    // Вычисление угла поворота
    double app = C_copy [p * coords + p];
    double aqq = C_copy [q * coords + q];
    double apq = C_copy [p * coords + q];

    double phi = 0.5 * MathArctan (2.0 * apq / (aqq - app + 1e-14));
    double c = MathCos (phi);
    double s = MathSin (phi);

    // Обновление элементов матрицы
    double app_new = c * c * app - 2 * c * s * apq + s * s * aqq;
    double aqq_new = s * s * app + 2 * c * s * apq + c * c * aqq;

    C_copy [p * coords + p] = app_new;
    C_copy [q * coords + q] = aqq_new;
    C_copy [p * coords + q] = C_copy [q * coords + p] = 0.0;

    // Обновление других элементов в строках/столбцах p и q
    for (int i = 0; i < coords; i++)
    {
      if (i != p && i != q)
      {
        double aip = C_copy [i * coords + p];
        double aiq = C_copy [i * coords + q];
        C_copy [i * coords + p] = C_copy [p * coords + i] = c * aip - s * aiq;
        C_copy [i * coords + q] = C_copy [q * coords + i] = s * aip + c * aiq;
      }
    }

    // Обновление собственных векторов
    for (int i = 0; i < coords; i++)
    {
      double bip = B [i * coords + p];
      double biq = B [i * coords + q];
      B [i * coords + p] = c * bip - s * biq;
      B [i * coords + q] = s * bip + c * biq;
    }
  }

  // Извлечение собственных значений и обеспечение положительности
  double min_eigenvalue = 1e-14;
  for (int i = 0; i < coords; i++)
  {
    D [i] = MathSqrt (MathMax (min_eigenvalue, C_copy [i * coords + i]));
  }

  // Сортировка собственных значений и векторов по убыванию
  for (int i = 0; i < coords - 1; i++)
  {
    int max_idx = i;
    for (int j = i + 1; j < coords; j++)
    {
      if (D [j] > D [max_idx]) max_idx = j;
    }

    if (max_idx != i)
    {
      // Обмен собственных значений
      double temp = D [i];
      D [i] = D [max_idx];
      D [max_idx] = temp;

      // Обмен собственных векторов
      for (int k = 0; k < coords; k++)
      {
        temp = B [k * coords + i];
        B [k * coords + i] = B [k * coords + max_idx];
        B [k * coords + max_idx] = temp;
      }
    }
  }
}
//————————————————————————————————————————————————————————————————————

//+------------------------------------------------------------------+
//| Проверка положительной определенности ковариационной матрицы     |
//+------------------------------------------------------------------+
bool C_AO_CMAES::CheckPositiveDefinite ()
{
  // Быстрая проверка: все диагональные элементы должны быть положительными
  for (int i = 0; i < coords; i++)
  {
    if (covMatrix [i * coords + i] <= 0) return false;
  }

  // Проверка минимального собственного значения на положительность
  double min_eigenvalue = D [coords - 1]; // D отсортирован по убыванию
  return min_eigenvalue > 1e-14;
}
//————————————————————————————————————————————————————————————————————

//+------------------------------------------------------------------+
//| Обеспечение положительной определенности ковариационной матрицы  |
//+------------------------------------------------------------------+
void C_AO_CMAES::EnforcePositiveDefinite ()
{
  // Метод 1: Добавление малого значения к диагонали
  double min_diag = 1e308; // Очень большое число
  for (int i = 0; i < coords; i++)
  {
    if (covMatrix [i * coords + i] < min_diag)
    {
      min_diag = covMatrix [i * coords + i];
    }
  }

  if (min_diag < 1e-10)
  {
    double correction = 1e-10 - min_diag;
    for (int i = 0; i < coords; i++)
    {
      covMatrix [i * coords + i] += correction;
    }
  }

  // Метод 2: Обеспечение симметрии
  for (int i = 0; i < coords; i++)
  {
    for (int j = i + 1; j < coords; j++)
    {
      double avg = (covMatrix [i * coords + j] + covMatrix [j * coords + i]) * 0.5;
      covMatrix [i * coords + j] = covMatrix [j * coords + i] = avg;
    }
  }

  // Если все еще не положительно определена, выполнить разложение и исправить
  if (!CheckPositiveDefinite ())
  {
    ComputeEigendecomposition ();

    double min_eigenvalue = 1e-10;
    for (int i = 0; i < coords; i++)
    {
      if (D [i] < MathSqrt (min_eigenvalue))
      {
        D [i] = MathSqrt (min_eigenvalue);
      }
    }

    // Восстановление C = B * D^2 * B^T
    ArrayInitialize (covMatrix, 0.0);
    for (int i = 0; i < coords; i++)
    {
      for (int j = 0; j <= i; j++)
      {
        double sum = 0.0;
        for (int k = 0; k < coords; k++)
        {
          sum += B [i * coords + k] * D [k] * D [k] * B [j * coords + k];
        }
        covMatrix [i * coords + j] = covMatrix [j * coords + i] = sum;
      }
    }
  }
}
//————————————————————————————————————————————————————————————————————