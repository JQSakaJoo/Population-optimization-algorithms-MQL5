//+——————————————————————————————————————————————————————————————————+
//|                                                          C_AO_EO |
//|                                  Copyright 2007-2025, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/18755

#include "#C_AO.mqh"
/*
//————————————————————————————————————————————————————————————————————
class C_AO_EO : public C_AO
{
  public: //----------------------------------------------------------
  ~C_AO_EO () { }
  C_AO_EO ()
  {
    ao_name = "EO";
    ao_desc = "Extremal Optimization";
    ao_link = "https://www.mql5.com/ru/articles/18755";

    popSize      = 50;      // Размер популяции
    tau          = 1.4;     // Параметр степенного распределения (τ)
    greedyStart  = 0.5;     // Доля агентов с жадной инициализацией
    eliteUpdate  = 0.3;     // Доля популяции для обновления за итерацию

    ArrayResize (params, 4);

    params [0].name = "popSize";     params [0].val = popSize;
    params [1].name = "tau";         params [1].val = tau;
    params [2].name = "greedyStart"; params [2].val = greedyStart;
    params [3].name = "eliteUpdate"; params [3].val = eliteUpdate;
  }

  void SetParams ()
  {
    popSize     = (int)params [0].val;
    tau         = params      [1].val;
    greedyStart = params      [2].val;
    eliteUpdate = params      [3].val;
  }

  bool Init (const double &rangeMinP  [],
             const double &rangeMaxP  [],
             const double &rangeStepP [],
             const int     epochsP = 0);

  void Moving   ();
  void Revision ();

  //------------------------------------------------------------------
  double tau;          // Параметр степенного распределения
  double greedyStart;  // Доля жадной инициализации
  double eliteUpdate;  // Доля обновляемых агентов

  private: //---------------------------------------------------------
  // Структуры для ранжирования
  struct RankedComponent
  {
      int    agentIdx;
      int    componentIdx;
      double fitness;      // λi - fitness компонента
  };

  struct RankedAgent
  {
      int    idx;
      double fitness;      // общий fitness агента
  };

  RankedComponent compRanks [];  // ранжированные компоненты
  RankedAgent     agentRanks []; // ранжированные агенты

  void   ApplyExtremalOptimization ();
  double CalculateComponentFitness (int agentIdx, int componentIdx);
  int    SelectRankByPowerLaw (int maxRank);
};
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Инициализация
bool C_AO_EO::Init (const double &rangeMinP  [],
                    const double &rangeMaxP  [],
                    const double &rangeStepP [],
                    const int epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //------------------------------------------------------------------
  ArrayResize (compRanks, coords);
  ArrayResize (agentRanks, popSize);

  return true;
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Основной цикл алгоритма
void C_AO_EO::Moving ()
{
  // Начальная инициализация популяции
  if (!revision)
  {
    int greedyCount = (int)(popSize * greedyStart);

    for (int i = 0; i < popSize; i++)
    {
      // Случайная инициализация для остальных
      for (int c = 0; c < coords; c++)
      {
        a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }

    revision = true;
    return;
  }

  // Применяем Extremal Optimization ---------------------------------
  ApplyExtremalOptimization ();
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Применение Extremal Optimization
void C_AO_EO::ApplyExtremalOptimization ()
{
  // Количество агентов для обновления на этой итерации
  int numUpdates = MathMax (1, (int)(popSize * eliteUpdate));

  // Обновляем выбранных агентов по принципу EO
  //for (int update = 0; update < numUpdates; update++)
  for (int update = 0; update < popSize; update++)
  {
    // Шаг 1: Выбираем агента для модификации
    // Используем ранжирование по общему fitness
    int targetAgent;

    // Ранжируем агентов по fitness (от худшего к лучшему для максимизации)
    for (int i = 0; i < popSize; i++)
    {
      agentRanks [i].idx = i;
      agentRanks [i].fitness = a [i].f;
    }

    // Сортировка (худшие в начале для максимизации)
    for (int i = 0; i < popSize - 1; i++)
    {
      for (int j = i + 1; j < popSize; j++)
      {
        if (agentRanks [i].fitness > agentRanks [j].fitness)
        {
          RankedAgent temp = agentRanks [i];
          agentRanks [i] = agentRanks [j];
          agentRanks [j] = temp;
        }
      }
    }

    // Выбираем агента согласно степенному распределению
    int rank    = SelectRankByPowerLaw (popSize);
    targetAgent = agentRanks [rank].idx;

    // Шаг 2: Ранжируем компоненты выбранного агента
    for (int c = 0; c < coords; c++)
    {
      compRanks [c].agentIdx     = targetAgent;
      compRanks [c].componentIdx = c;
      compRanks [c].fitness      = CalculateComponentFitness (targetAgent, c);
    }

    // Сортировка компонентов (худшие в начале)
    for (int i = 0; i < coords - 1; i++)
    {
      for (int j = i + 1; j < coords; j++)
      {
        if (compRanks [i].fitness > compRanks [j].fitness)
        {
          RankedComponent temp = compRanks [i];
          compRanks [i] = compRanks [j];
          compRanks [j] = temp;
        }
      }
    }

    // Шаг 3: Выбираем компонент для изменения согласно P(n) ∝ n^(-τ)
    int compRank = SelectRankByPowerLaw (coords);
    int compIdx  = compRanks [compRank].componentIdx;

    // Шаг 4: Заменяем выбранный компонент новым случайным значением
    // Это ключевой принцип EO - безусловная замена на случайное
    a [targetAgent].c [compIdx] = u.RNDfromCI (rangeMin [compIdx], rangeMax [compIdx]);

    // Проверка границ
    a [targetAgent].c [compIdx] = u.SeInDiSp (a [targetAgent].c [compIdx],
                                              rangeMin [compIdx],
                                              rangeMax [compIdx],
                                              rangeStep [compIdx]);
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Расчет fitness компонента
double C_AO_EO::CalculateComponentFitness (int agentIdx, int componentIdx)
{
  // Для общей задачи оптимизации используем простую метрику
  // λi = относительный вклад компонента в общее качество

  double fitness = 0.0;

  double range = rangeMax [componentIdx] - rangeMin [componentIdx];
  if (range > 0)
  {
    // Нормализованное отклонение
    double deviation = MathAbs (a [agentIdx].c [componentIdx] - cB [componentIdx]) / range;
    fitness = 1.0 - deviation; // Инвертируем, чтобы больше = лучше
  }

  return fitness;
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Выбор ранга согласно степенному распределению
int C_AO_EO::SelectRankByPowerLaw (int maxRank)
{
  // P(n) ∝ n^(-τ), где n - ранг от 1 до maxRank
  // Используем метод обратного преобразования

  double r = u.RNDprobab ();

  if (tau != 1.0)
  {
    // Общий случай: обратное преобразование для P(n) ∝ n^(-τ)
    double norm = (1.0 - MathPow (maxRank + 1.0, 1.0 - tau)) / (1.0 - tau);
    double x = r * norm;
    int rank = (int)MathPow ((1.0 - tau) * x + 1.0, 1.0 / (1.0 - tau)) - 1;

    if (rank >= maxRank) rank = maxRank - 1;
    if (rank < 0) rank = 0;

    return rank;
  }
  else
  {
    // Специальный случай τ = 1: P(n) ∝ 1/n
    double norm = MathLog (maxRank + 1.0);
    int rank = (int)(MathExp (r * norm) - 1.0);

    if (rank >= maxRank) rank = maxRank - 1;
    if (rank < 0) rank = 0;

    return rank;
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Обновление лучших решений
void C_AO_EO::Revision ()
{
  // Сортировка популяции для МАКСИМИЗАЦИИ
  static S_AO_Agent aT [];
  ArrayResize (aT, popSize);

  // Используем встроенную функцию сортировки
  u.Sorting (a, aT, popSize);

  // Обновление глобального лучшего решения
  if (a [0].f > fB)
  {
    ArrayCopy (cB, a [0].c, 0, 0, WHOLE_ARRAY);
    fB = a [0].f;
  }
}
//————————————————————————————————————————————————————————————————————
*/

//————————————————————————————————————————————————————————————————————
class C_AO_EOm : public C_AO
{
  public: //----------------------------------------------------------
  ~C_AO_EOm () { }
  C_AO_EOm ()
  {
    ao_name = "EO";
    ao_desc = "Extremal Optimization";
    ao_link = "https://www.mql5.com/ru/articles/18755";

    popSize        = 50;      // Размер популяции
    popRaising     = 3;       // Повышение самых худших
    mutationRate   = 0.1;     // Вероятность мутации
    powCh          = 2.0;     // Степень закона распределения отбора
    powMut         = 8.0;     // Степень закона распределения мутации

    ArrayResize (params, 5);

    params [0].name = "popSize";        params [0].val = popSize;
    params [1].name = "popRaising";     params [1].val = popRaising;
    params [2].name = "mutationRate";   params [2].val = mutationRate;
    params [3].name = "powCh";          params [3].val = powCh;
    params [4].name = "powMut";         params [4].val = powMut;
  }

  void SetParams ()
  {
    popSize        = (int)params [0].val;
    popRaising     = (int)params [1].val;
    mutationRate   = params      [2].val;
    powCh          = params      [3].val;
    powMut         = params      [4].val;
  }

  bool Init (const double &rangeMinP  [],
             const double &rangeMaxP  [],
             const double &rangeStepP [],
             const int     epochsP = 0);

  void Moving   ();
  void Revision ();

  //------------------------------------------------------------------
  int    popRaising;     // Повышение самых худших
  double mutationRate;   // Вероятность мутации
  double powCh;          // Степень закона распределения отбора
  double powMut;         // Степень закона распределения мутации

  private: //---------------------------------------------------------
  int    currentEpoch;     // текущая эпоха
  int    totalEpochs;      // общее количество эпох

  void MutateComponent (int agentIdx, int componentIdx);
};
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Инициализация
bool C_AO_EOm::Init (const double &rangeMinP  [],
                     const double &rangeMaxP  [],
                     const double &rangeStepP [],
                     const int epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //------------------------------------------------------------------
  currentEpoch = 0;
  totalEpochs  = epochsP;

  return true;
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Основной цикл алгоритма
void C_AO_EOm::Moving ()
{
  currentEpoch++;

  // Начальная инициализация популяции
  if (!revision)
  {
    for (int i = 0; i < popSize; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }

    revision = true;
    return;
  }

  //Apply Extremal Optimization---------------------------------------
  static S_AO_Agent aT []; ArrayResize (aT, popSize);

  for (int i = 0; i < popSize; i++)
  {
    aT [i].Init (coords);

    for (int c = 0; c < coords; c++)
    {
      double rnd = u.RNDprobab (); rnd = pow (rnd, powCh);
      int ind = (int)u.Scale (rnd, 0.0, 1.0, 0, popSize - 1);

      // Выбор типа мутации
      double mutType = u.RNDprobab ();

      if (mutType < mutationRate)
      {
        aT [i].c [c] = u.PowerDistribution (a [ind].c [c], rangeMin [c], rangeMax [c], powMut);
      }
      else
      {
        // Направленное движение к лучшему с шумом
        aT [i].c [c] = a [ind].c [c] + u.RNDprobab () * (cB [c] - a [ind].c [c]);
      }

      // Проверка границ
      aT [i].c [c] = u.SeInDiSp (aT [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }

  for (int i = 0; i < popSize; i++) ArrayCopy (a [i].c, aT [i].c);
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Обновление лучших решений
void C_AO_EOm::Revision ()
{
  // Сортировка популяции --------------------------------------------
  static S_AO_Agent aT []; ArrayResize (aT, popSize);
  u.Sorting (a, aT, popSize);

  // Обновление глобального лучшего решения
  if (a [0].f > fB)
  {
    ArrayCopy (cB, a [0].c, 0, 0, WHOLE_ARRAY);
    fB = a [0].f;
  }

  fW = a [popSize - 1].f;

  //------------------------------------------------------------------
  for (int i = 0; i < popRaising; i++)
  {
    a [popSize - 1 - i].f = u.RNDfromCI (fW, fB);
  }

  u.Sorting (a, aT, popSize);
}
//————————————————————————————————————————————————————————————————————
