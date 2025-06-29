//+——————————————————————————————————————————————————————————————————+
//|                                                          C_AO_ES |
//|                                  Copyright 2007-2025, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/18460

#include "#C_AO.mqh"

//————————————————————————————————————————————————————————————————————
class C_AO_ES : public C_AO
{
  public: //----------------------------------------------------------
  ~C_AO_ES () { }
  C_AO_ES ()
  {
    ao_name = "ES";
    ao_desc = "Eagle Strategy";
    ao_link = "https://www.mql5.com/ru/articles/18460";

    popSize         = 100;   // размер популяции
    lambda          = 1.0;  // параметр распределения Леви (1 < λ ≤ 3)
    sphereRadius    = 0.1;  // радиус гиперсферы для локального поиска
    localIterations = 20;   // количество итераций локального поиска
    alpha           = 0.1;  // параметр рандомизации для Firefly
    beta0           = 1.2;  // начальная привлекательность

    ArrayResize (params, 6);

    params [0].name = "popSize";         params [0].val = popSize;
    params [1].name = "lambda";          params [1].val = lambda;
    params [2].name = "sphereRadius";    params [2].val = sphereRadius;
    params [3].name = "localIterations"; params [3].val = localIterations;
    params [4].name = "alpha";           params [4].val = alpha;
    params [5].name = "beta0";           params [5].val = beta0;
  }

  void SetParams ()
  {
    popSize         = (int)params [0].val;
    lambda          = params      [1].val;
    sphereRadius    = params      [2].val;
    localIterations = (int)params [3].val;
    alpha           = params      [4].val;
    beta0           = params      [5].val;
  }

  bool Init (const double &rangeMinP  [],  // минимальные значения
             const double &rangeMaxP  [],  // максимальные значения
             const double &rangeStepP [],  // шаг изменения
             const int     epochsP = 0);   // количество эпох

  void Moving   ();
  void Revision ();

  //------------------------------------------------------------------
  double lambda;          // параметр распределения Леви (1 < λ ≤ 3)
  double sphereRadius;    // радиус гиперсферы для локального поиска
  int    localIterations; // количество итераций локального поиска
  double alpha;           // параметр рандомизации
  double beta0;           // начальная привлекательность

  private: //---------------------------------------------------------
  double gamma_es;           // коэффициент поглощения света
  double levyStep [];        // массив для шагов Леви

  // Отслеживание фаз
  bool   inLocalSearchPhase; // флаг локального поиска
  int    localSearchCenter;  // центр локального поиска
  int    localSearchCounter; // счетчик итераций локального поиска

  // Отслеживание сходимости
  double prevBestFitness;    // предыдущее лучшее значение
  int    stagnationCounter;  // счетчик стагнации

  // Отслеживание эпох
  int    epochCurrent;       // текущая эпоха
  int    epochMax;           // максимальное количество эпох

  // Вспомогательные методы
  void   GlobalExploration  ();
  void   LocalExploitation  ();
  void   GenerateLevyStep   ();
  double GenerateGaussian   ();
  double Gamma              (double z);
};
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
bool C_AO_ES::Init (const double &rangeMinP  [],
                    const double &rangeMaxP  [],
                    const double &rangeStepP [],
                    const int     epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //------------------------------------------------------------------
  // Инициализация массивов
  ArrayResize (levyStep, coords);

  // Инициализация отслеживания фаз
  inLocalSearchPhase = false;
  localSearchCenter  = 0;
  localSearchCounter = 0;

  // Инициализация отслеживания сходимости
  prevBestFitness   = -DBL_MAX;
  stagnationCounter = 0;

  // Инициализация отслеживания эпох
  epochMax     = epochsP;
  epochCurrent = 0;

  // Фиксированный параметр Firefly
  gamma_es = 1.0;

  // Инициализация популяции случайным образом
  for (int i = 0; i < popSize; i++)
  {
    for (int c = 0; c < coords; c++)
    {
      a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
      a [i].c [c] = u.SeInDiSp  (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }

    a [i].f  = -DBL_MAX;
    a [i].fB = -DBL_MAX;
  }

  return true;
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
void C_AO_ES::Moving ()
{
  epochCurrent++;

  // ПРИНЯТИЕ РЕШЕНИЯ О ФАЗЕ: Чередование между глобальным и локальным поиском
  if (!inLocalSearchPhase)
  {
    // ФАЗА 1: ГЛОБАЛЬНОЕ ИССЛЕДОВАНИЕ с использованием полетов Леви
    GlobalExploration ();

    // Проверка необходимости переключения на локальный поиск
    // Переключаемся, если нашли многообещающую область (улучшение лучшей приспособленности)
    if (fB > prevBestFitness)
    {
      inLocalSearchPhase = true;
      localSearchCounter = 0;
      prevBestFitness    = fB;
      stagnationCounter  = 0;

      // Поиск лучшего агента для центрирования локального поиска
      localSearchCenter = 0;
      double bestFit = -DBL_MAX;

      for (int i = 0; i < popSize; i++)
      {
        if (a [i].f > bestFit)
        {
          bestFit = a [i].f;
          localSearchCenter = i;
        }
      }
    }
    else
    {
      stagnationCounter++;

      // При стагнации увеличиваем исследование
      if (stagnationCounter > 5)
      {
        lambda = MathMax (1.0, lambda - 0.1); // Делаем полеты Леви более агрессивными
      }
    }
  }
  else
  {
    if (u.RNDprobab () < 0.8)
    {
      // ФАЗА 2: ЛОКАЛЬНАЯ ЭКСПЛУАТАЦИЯ с использованием алгоритма Firefly
      LocalExploitation ();

      localSearchCounter++;

      // Возврат к глобальному поиску после завершения локальных итераций
      if (localSearchCounter >= localIterations)
      {
        inLocalSearchPhase = false;
        lambda = params [1].val; // Сброс lambda к исходному значению
      }
    }
    else
    {
      for (int i = 0; i < popSize; i++)
      {
        for (int c = 0; c < coords; c++)
        {
          if (u.RNDprobab () < 0.5)
          {
            a [i].c [c] = cB [c];
          }
        }
      }
    }
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
void C_AO_ES::Revision ()
{
  for (int i = 0; i < popSize; i++)
  {
    // Обновление персонального лучшего
    if (a [i].f > a [i].fB)
    {
      a [i].fB = a [i].f;
      ArrayCopy (a [i].cB, a [i].c, 0, 0, WHOLE_ARRAY);
    }

    // Обновление глобального лучшего
    if (a [i].f > fB)
    {
      fB = a [i].f;
      ArrayCopy (cB, a [i].c, 0, 0, WHOLE_ARRAY);
    }
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
// ФАЗА 1: Глобальное исследование с использованием полетов Леви
void C_AO_ES::GlobalExploration ()
{
  for (int i = 0; i < popSize; i++)
  {
    // Генерация шага Леви
    GenerateLevyStep ();

    // Обновление позиции с использованием полета Леви
    for (int c = 0; c < coords; c++)
    {
      double range = rangeMax [c] - rangeMin [c];

      // Адаптивный размер шага в зависимости от прогресса поиска
      double progress = (epochMax > 0) ? (double)epochCurrent / (double)epochMax : 0.5;
      double stepScale = 0.01 + 0.2 * (1.0 - progress); // Начинаем с больших шагов

      // Применение шага Леви
      a [i].c [c] += levyStep [c] * range * stepScale;

      // Ограничения границ
      a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
// ФАЗА 2: Локальная эксплуатация с использованием алгоритма Firefly
void C_AO_ES::LocalExploitation ()
{
  // Идентификация агентов внутри гиперсферы вокруг лучшего решения
  double agents_in_sphere [];
  ArrayResize (agents_in_sphere, 0);

  for (int i = 0; i < popSize; i++)
  {
    double normalized_dist = 0.0;

    for (int c = 0; c < coords; c++)
    {
      double diff = (a [i].c [c] - a [localSearchCenter].c [c]) / (rangeMax [c] - rangeMin [c]);
      normalized_dist += diff * diff;
    }
    normalized_dist = MathSqrt (normalized_dist);

    // Включаем агентов внутри сферы или сам центр
    if (normalized_dist <= sphereRadius || i == localSearchCenter)
    {
      int size = ArraySize (agents_in_sphere);
      ArrayResize (agents_in_sphere, size + 1);
      agents_in_sphere [size] = i;
    }
  }

  // Если слишком мало агентов, расширяем до ближайших соседей
  if (ArraySize (agents_in_sphere) < 5)
  {
    ArrayResize (agents_in_sphere, 0);

    // Вычисляем расстояния для всех агентов
    double distances [];
    ArrayResize (distances, popSize);

    for (int i = 0; i < popSize; i++)
    {
      if (i == localSearchCenter)
      {
        distances [i] = 0.0;
      }
      else
      {
        double dist = 0.0;
        for (int c = 0; c < coords; c++)
        {
          double diff = (a [i].c [c] - a [localSearchCenter].c [c]) / (rangeMax [c] - rangeMin [c]);
          dist += diff * diff;
        }
        distances [i] = MathSqrt (dist);
      }
    }

    // Берем ближайших 5 агентов или 30% популяции
    int numAgents = MathMin (popSize, MathMax (5, popSize / 3));
    ArrayResize (agents_in_sphere, numAgents);

    // Простой выбор ближайших агентов
    for (int k = 0; k < numAgents; k++)
    {
      double minDist = DBL_MAX;
      int minIdx = -1;

      for (int i = 0; i < popSize; i++)
      {
        bool already_selected = false;

        for (int j = 0; j < k; j++)
        {
          if (agents_in_sphere [j] == i)
          {
            already_selected = true;
            break;
          }
        }

        if (!already_selected && distances [i] < minDist)
        {
          minDist = distances [i];
          minIdx = i;
        }
      }

      agents_in_sphere [k] = minIdx;
    }
  }

  // Выполнение алгоритма Firefly на выбранных агентах
  int numLocalAgents = ArraySize (agents_in_sphere);

  for (int i = 0; i < numLocalAgents; i++)
  {
    int idx_i = (int)agents_in_sphere [i];

    for (int j = 0; j < numLocalAgents; j++)
    {
      if (i == j) continue;

      int idx_j = (int)agents_in_sphere [j];

      // Если агент j лучше агента i, двигаем i к j
      if (a [idx_j].f > a [idx_i].f)
      {
        // Вычисление расстояния
        double r_squared = 0.0;

        for (int c = 0; c < coords; c++)
        {
          double diff = (a [idx_j].c [c] - a [idx_i].c [c]) / (rangeMax [c] - rangeMin [c]);
          r_squared += diff * diff;
        }

        // Вычисление привлекательности
        double beta = beta0 * MathExp (-gamma_es * r_squared);

        // Перемещение агента i к агенту j
        for (int c = 0; c < coords; c++)
        {
          double range = rangeMax [c] - rangeMin [c];

          // Уравнение движения Firefly
          a [idx_i].c [c] += beta * (a [idx_j].c [c] - a [idx_i].c [c]) +
                             alpha * (u.RNDfromCI (-0.5, 0.5)) * range * 0.1;

          // Применение границ
          a [idx_i].c [c] = u.SeInDiSp (a [idx_i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
        }
      }
    }
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
// Генерация шага Леви с использованием алгоритма Мантенья
void C_AO_ES::GenerateLevyStep ()
{
  // Вычисление сигмы для алгоритма Мантенья
  double numerator   = Gamma (1.0 + lambda) * MathSin (M_PI * lambda / 2.0);
  double denominator = Gamma ((1.0 + lambda) / 2.0) * lambda * MathPow (2.0, (lambda - 1.0) / 2.0);
  double sigma = MathPow (numerator / denominator, 1.0 / lambda);

  for (int c = 0; c < coords; c++)
  {
    // Генерация u и v из нормальных распределений
    double u_val = GenerateGaussian () * sigma;
    double v_val = MathAbs (GenerateGaussian ());

    // Вычисление шага Леви
    if (v_val > 1e-10)
    {
      levyStep [c] = u_val / MathPow (v_val, 1.0 / lambda);
    }
    else
    {
      levyStep [c] = 0.0;
    }

    // Ограничение экстремальных значений
    if (MathAbs (levyStep [c]) > 10.0)
    {
      levyStep [c] = 10.0 * (levyStep [c] > 0 ? 1.0 : -1.0);
    }
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
// Генерация гауссовского случайного числа с использованием преобразования Бокса-Мюллера
double C_AO_ES::GenerateGaussian ()
{
  static bool hasSpare = false;
  static double spare;

  if (hasSpare)
  {
    hasSpare = false;
    return spare;
  }

  hasSpare = true;
  double u_val = u.RNDfromCI (0.0, 1.0);
  double v_val = u.RNDfromCI (0.0, 1.0);

  double mag = MathSqrt (-2.0 * MathLog (u_val + 1e-10));
  spare = mag * MathCos (2.0 * M_PI * v_val);

  return mag * MathSin (2.0 * M_PI * v_val);
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
// Аппроксимация гамма-функции с использованием аппроксимации Ланцоша
double C_AO_ES::Gamma (double z)
{
  if (z < 0.5)
  {
    // Формула отражения для z < 0.5
    return M_PI / (MathSin (M_PI * z) * Gamma (1.0 - z));
  }

  // Коэффициенты Ланцоша
  const double g = 7.0;
  double coef [] =
  {
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7
  };

  z -= 1.0;
  double x = coef [0];

  for (int i = 1; i < 9; i++)
  {
    x += coef [i] / (z + i);
  }

  double t = z + g + 0.5;
  double sqrt2pi = MathSqrt (2.0 * M_PI);

  return sqrt2pi * MathPow (t, z + 0.5) * MathExp (-t) * x;
}
//————————————————————————————————————————————————————————————————————