//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_COA_chaos |
//|                                            Copyright 2007-2025, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/16729

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
// Улучшенная структура агента с дополнительными полями
struct S_COA_Agent
{
    double gamma    [];       // хаотические переменные
    double velocity [];       // скорость перемещения (для усиления инерции)
    int    stagnationCounter; // счётчик стагнации

    void Init (int coords)
    {
      ArrayResize (gamma,    coords);
      ArrayResize (velocity, coords);

      // Равномерное распределение значений для gamma
      for (int i = 0; i < coords; i++)
      {
        // Используем различные начальные значения для лучшего разнообразия
        gamma [i] = 0.1 + 0.8 * (i % coords) / (double)MathMax (1, coords - 1);

        // Инициализация скорости нулями
        velocity [i] = 0.0;
      }

      stagnationCounter = 0;
    }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
class C_AO_COA_chaos : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_COA_chaos () { }
  C_AO_COA_chaos ()
  {
    ao_name = "COA(CHAOS)";
    ao_desc = "Chaos Optimization Algorithm";
    ao_link = "https://www.mql5.com/ru/articles/16729";

    // Внутренние параметры (не настраиваются извне)
    inertia      = 0.7;
    socialFactor = 1.5;
    mutationRate = 0.05;

    // Параметры по умолчанию
    popSize = 50;
    S1      = 30;
    S2      = 20;
    sigma   = 2.0;
    t3      = 1.2;
    eps     = 0.0001;

    // Инициализация массива параметров для интерфейса C_AO
    ArrayResize (params, 6);

    params [0].name = "popSize"; params [0].val = popSize;
    params [1].name = "S1";      params [1].val = S1;
    params [2].name = "S2";      params [2].val = S2;
    params [3].name = "sigma";   params [3].val = sigma;
    params [4].name = "t3";      params [4].val = t3;
    params [5].name = "eps";     params [5].val = eps;
  }

  void SetParams ()
  {
    // Обновление внутренних параметров из массива params
    popSize = (int)params [0].val;
    S1      = (int)params [1].val;
    S2      = (int)params [2].val;
    sigma   = params      [3].val;
    t3      = params      [4].val;
    eps     = params      [5].val;
  }

  bool Init (const double &rangeMinP  [], // минимальный диапазон поиска
             const double &rangeMaxP  [], // максимальный диапазон поиска
             const double &rangeStepP [], // шаг поиска
             const int     epochsP = 0);  // количество эпох

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  // Внешние параметры алгоритма
  int    S1;             // итерации первой фазы
  int    S2;             // итерации второй фазы
  double sigma;          // параметр штрафа
  double t3;             // коэффициент корректировки alpha
  double eps;            // малое число для весовых коэффициентов

  // Внутренние параметры алгоритма
  double inertia;        // параметр инерции для движения (внутренний)
  double socialFactor;   // параметр социального влияния (внутренний)
  double mutationRate;   // вероятность мутации (внутренний)

  S_COA_Agent agent [];  // массив агентов

  private: //-------------------------------------------------------------------
  int    epochNow;
  double currentSigma;           // Динамический параметр штрафа
  double alpha             [];   // параметры поиска
  double globalBestHistory [10]; // История значений глобального лучшего решения
  int    historyIndex;

  // Вспомогательные методы
  double CalculateWeightedGradient   (int agentIdx, int coordIdx);
  double CalculateConstraintValue    (int agentIdx, int coordIdx);
  double CalculateConstraintGradient (int agentIdx, int coordIdx);
  double CalculatePenaltyFunction    (int agentIdx);

  // Метод для проверки допустимости решения
  bool IsFeasible              (int agentIdx);

  // Хаотические отображения
  double LogisticMap           (double x);
  double SineMap               (double x);
  double TentMap               (double x);
  double SelectChaosMap        (double x, int type);

  void InitialPopulation       ();
  void FirstCarrierWaveSearch  ();
  void SecondCarrierWaveSearch ();
  void ApplyMutation           (int agentIdx);
  void UpdateSigma             ();
  void UpdateBestHistory       (double newBest);
  bool IsConverged             ();
  void ResetStagnatingAgents   ();
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_COA_chaos::Init (const double &rangeMinP  [],
                     const double &rangeMaxP  [],
                     const double &rangeStepP [],
                     const int     epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  epochNow     = 0;
  currentSigma = sigma;
  historyIndex = 0;

  // Инициализация истории лучших значений
  for (int i = 0; i < 10; i++) globalBestHistory [i] = -DBL_MAX;

  // Проверка и инициализация основных массивов
  int arraySize = ArraySize (rangeMinP);
  if (arraySize <= 0 || arraySize != ArraySize (rangeMaxP) || arraySize != ArraySize (rangeStepP))
  {
    return false;
  }

  ArrayResize (agent, popSize);
  ArrayResize (alpha, coords);

  // Адаптивная инициализация alpha в зависимости от диапазона поиска
  for (int c = 0; c < coords; c++)
  {
    // alpha зависит от размера пространства поиска
    double range = rangeMax [c] - rangeMin [c];
    alpha [c] = 0.1 * range / MathSqrt (MathMax (1.0, (double)coords));
  }

  // Инициализация агентов с разнообразными стратегиями
  for (int i = 0; i < popSize; i++)
  {
    agent [i].Init (coords);

    for (int c = 0; c < coords; c++)
    {
      double position;

      // Различные стратегии инициализации
      if (i < popSize / 4)
      {
        // Равномерное распределение по пространству
        position = rangeMin [c] + (i * (rangeMax [c] - rangeMin [c])) / MathMax (1, popSize / 4);
      }
      else
        if (i < popSize / 2)
        {
          // Кластеризация вокруг нескольких точек
          int cluster = (i - popSize / 4) % 3;
          double clusterCenter = rangeMin [c] + (cluster + 1) * (rangeMax [c] - rangeMin [c]) / 4.0;
          position = clusterCenter + u.RNDfromCI (-0.1, 0.1) * (rangeMax [c] - rangeMin [c]);
        }
        else
          if (i < 3 * popSize / 4)
          {
            // Случайные позиции с смещением в сторону границ
            double r = u.RNDprobab ();
            if (r < 0.5) position = rangeMin [c] + 0.2 * r * (rangeMax [c] - rangeMin [c]);
            else position = rangeMax [c] - 0.2 * (1.0 - r) * (rangeMax [c] - rangeMin [c]);
          }
          else
          {
            // Хаотические позиции с использованием разных отображений
            int mapType = i % 3;
            double chaosValue = SelectChaosMap (agent [i].gamma [c], mapType);
            position = rangeMin [c] + chaosValue * (rangeMax [c] - rangeMin [c]);
          }

      a [i].cB [c] = u.SeInDiSp (position, rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_COA_chaos::Moving ()
{
  epochNow++;

  if (!revision)
  {
    InitialPopulation ();

    revision = true;
    return;
  }

  // Динамическое обновление параметра штрафа
  UpdateSigma ();

  // Сброс агентов, находящихся в стагнации
  if (epochNow % 5 == 0)
  {
    ResetStagnatingAgents ();
  }

  // Определяем фазу поиска
  if (epochNow <= S1)
  {
    FirstCarrierWaveSearch ();
  }
  else
    if (epochNow <= S1 + S2)
    {
      SecondCarrierWaveSearch ();
    }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_COA_chaos::Revision ()
{
  int    improvementCount = 0;
  double bestImprovement  = 0.0;

  // Оценка всех решений
  for (int i = 0; i < popSize; i++)
  {
    double newValue = CalculatePenaltyFunction (i);

    // Проверка на валидность нового значения
    if (!MathIsValidNumber (newValue)) continue;

    // Сохраняем текущее значение для следующей итерации
    a [i].f = newValue;

    // Обновление личного лучшего решения (перед глобальным для эффективности)
    if (newValue > a [i].fB)
    {
      a [i].fB = newValue;
      ArrayCopy (a [i].cB, a [i].c, 0, 0, WHOLE_ARRAY);
      improvementCount++;
    }

    // Обновление глобального лучшего решения
    if (newValue > fB)
    {
      double improvement = newValue - fB;
      fB = newValue;
      ArrayCopy (cB, a [i].c, 0, 0, WHOLE_ARRAY);

      // Обновляем историю лучших значений
      UpdateBestHistory (fB);

      bestImprovement = MathMax (bestImprovement, improvement);
      improvementCount++;
    }
  }

  // Адаптация параметров поиска в зависимости от фазы и успешности
  if (epochNow > 1)
  {
    // Коэффициент успешности поиска (предотвращение деления на ноль)
    double successRate = (double)improvementCount / MathMax (1, 2 * popSize);

    // Адаптация параметра alpha
    for (int c = 0; c < coords; c++)
    {
      double oldAlpha = alpha [c];

      if (epochNow <= S1)
      {
        // В фазе глобального поиска
        if (successRate < 0.1)
        {
          // Очень мало улучшений - увеличиваем область поиска
          alpha [c] *= t3;
        }
        else
          if (successRate < 0.3)
          {
            // Мало улучшений - слегка увеличиваем область
            alpha [c] *= 1.2;
          }
          else
            if (successRate > 0.7)
            {
              // Много улучшений - сужаем область
              alpha [c] *= 0.9;
            }
      }
      else
      {
        // В фазе локального поиска
        if (successRate < 0.05)
        {
          // Очень мало улучшений - увеличиваем область поиска
          alpha [c] *= t3;
        }
        else
          if (successRate > 0.2)
          {
            // Достаточно улучшений - сужаем область
            alpha [c] *= 0.8;
          }
      }

      // Ограничения на диапазон alpha с безопасной проверкой границ
      if (c < ArraySize (rangeMax) && c < ArraySize (rangeMin))
      {
        double maxAlpha = 0.2 * (rangeMax [c] - rangeMin [c]);
        double minAlpha = 0.0001 * (rangeMax [c] - rangeMin [c]);

        if (alpha [c] > maxAlpha)
        {
          alpha [c] = maxAlpha;
        }
        else
          if (alpha [c] < minAlpha)
          {
            alpha [c] = minAlpha;
          }
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
// Улучшенные хаотические отображения
double C_AO_COA_chaos::LogisticMap (double x)
{
  // Защита от некорректных входных значений
  if (x < 0.0 || x > 1.0 || MathIsValidNumber (x) == false)
  {
    x = 0.2 + 0.6 * u.RNDprobab ();
  }

  // x(n+1) = r*x(n)*(1-x(n))
  double r = 3.9 + 0.1 * u.RNDprobab (); // Слегка рандомизированный параметр для избежания циклов
  double result = r * x * (1.0 - x);

  // Дополнительная проверка корректности
  if (result < 0.0 || result > 1.0 || MathIsValidNumber (result) == false)
  {
    result = 0.2 + 0.6 * u.RNDprobab ();
  }

  return result;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_COA_chaos::SineMap (double x)
{
  // Защита от некорректных входных значений
  if (x < 0.0 || x > 1.0 || MathIsValidNumber (x) == false)
  {
    x = 0.2 + 0.6 * u.RNDprobab ();
  }

  // x(n+1) = sin(π*x(n))
  double result = MathSin (M_PI * x);

  // Нормализация результата к диапазону [0, 1]
  result = (result + 1.0) / 2.0;

  // Дополнительная проверка корректности
  if (result < 0.0 || result > 1.0 || MathIsValidNumber (result) == false)
  {
    result = 0.2 + 0.6 * u.RNDprobab ();
  }

  return result;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_COA_chaos::TentMap (double x)
{
  // Защита от некорректных входных значений
  if (x < 0.0 || x > 1.0 || MathIsValidNumber (x) == false)
  {
    x = 0.2 + 0.6 * u.RNDprobab ();
  }

  // Tent map: x(n+1) = μ*min(x(n), 1-x(n))
  double mu = 1.99; // Параметр близкий к 2 для хаотического поведения
  double result;

  if (x <= 0.5) result = mu * x;
  else result = mu * (1.0 - x);

  // Дополнительная проверка корректности
  if (result < 0.0 || result > 1.0 || MathIsValidNumber (result) == false)
  {
    result = 0.2 + 0.6 * u.RNDprobab ();
  }

  return result;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_COA_chaos::SelectChaosMap (double x, int type)
{
  // Выбор хаотического отображения на основе типа
  switch (type % 3)
  {
    case 0:
      return LogisticMap (x);
    case 1:
      return SineMap (x);
    case 2:
      return TentMap (x);
    default:
      return LogisticMap (x);
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_COA_chaos::InitialPopulation ()
{
  // Создаем Latin Hypercube для начальной популяции
  double latinCube []; // Одномерный массив для хранения значений гиперкуба
  double tempValues []; // Временный массив для хранения и перемешивания значений

  ArrayResize (latinCube, popSize * coords);
  ArrayResize (tempValues, popSize);

  // Генерируем Латинский гиперкуб
  for (int c = 0; c < coords; c++)
  {
    // Создаем упорядоченные значения
    for (int i = 0; i < popSize; i++)
    {
      tempValues [i] = (double)i / popSize;
    }

    // Перемешиваем значения
    for (int i = popSize - 1; i > 0; i--)
    {
      int j = (int)(u.RNDprobab () * (i + 1));
      if (j < popSize)
      {
        double temp = tempValues [i];
        tempValues [i] = tempValues [j];
        tempValues [j] = temp;
      }
    }

    // Присваиваем перемешанные значения
    for (int i = 0; i < popSize; i++)
    {
      latinCube [i * coords + c] = tempValues [i];
    }
  }

  // Преобразуем значения Латинского гиперкуба в координаты
  for (int i = 0; i < popSize; i++)
  {
    for (int c = 0; c < coords; c++)
    {
      double x = rangeMin [c] + latinCube [i * coords + c] * (rangeMax [c] - rangeMin [c]);
      a [i].c [c] = u.SeInDiSp (x, rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_COA_chaos::FirstCarrierWaveSearch ()
{
  // Адаптивный баланс между исследованием и эксплуатацией
  double globalPhase = (double)epochNow / S1;
  double explorationRate = 1.0 - globalPhase * globalPhase; // Квадратичное снижение

  // Для каждого агента
  for (int i = 0; i < popSize; i++)
  {
    // Применяем мутации с некоторой вероятностью для усиления разнообразия
    if (u.RNDprobab () < mutationRate * (1.0 + explorationRate))
    {
      ApplyMutation (i);
      continue;
    }

    for (int c = 0; c < coords; c++)
    {
      // Выбор хаотического отображения с равномерным распределением
      int mapType = ((i + c + epochNow) % 3);

      // Безопасная проверка доступа к массиву gamma
      if (c < ArraySize (agent [i].gamma))
      {
        agent [i].gamma [c] = SelectChaosMap (agent [i].gamma [c], mapType);
      }
      else
      {
        continue; // Пропускаем, если индекс некорректен
      }

      // Определяем соотношение между глобальным и локальным поиском
      double strategy = u.RNDprobab ();
      double x;

      if (strategy < explorationRate)
      {
        // Глобальный поиск с хаотическим компонентом
        x = rangeMin [c] + agent [i].gamma [c] * (rangeMax [c] - rangeMin [c]);

        // Добавляем компонент скорости для сохранения направления движения
        agent [i].velocity [c] = inertia * agent [i].velocity [c] +
                                 (1.0 - inertia) * (x - a [i].c [c]);
      }
      else
      {
        // Локальный поиск вокруг лучших решений
        double personalAttraction = u.RNDprobab ();
        double globalAttraction = u.RNDprobab ();

        // Взвешенное притяжение к лучшим решениям
        double attractionTerm = //personalAttraction * (agent [i].cPrev [c] - a [i].c [c]) +
                                personalAttraction * (a [i].cB [c] - a [i].c [c]) +
                                socialFactor * globalAttraction * (cB [c] - a [i].c [c]);

        // Хаотическое возмущение для предотвращения застревания
        double chaosRange = alpha [c] * explorationRate;
        double chaosTerm = chaosRange * (2.0 * agent [i].gamma [c] - 1.0);

        // Обновление скорости с инерцией
        agent [i].velocity [c] = inertia * agent [i].velocity [c] +
                                 (1.0 - inertia) * (attractionTerm + chaosTerm);
      }

      // Ограничиваем скорость для предотвращения слишком больших шагов
      double maxVelocity = 0.1 * (rangeMax [c] - rangeMin [c]);
      if (MathAbs (agent [i].velocity [c]) > maxVelocity)
      {
        agent [i].velocity [c] = maxVelocity * (agent [i].velocity [c] > 0 ? 1.0 : -1.0);
      }

      // Применяем скорость к позиции
      x = a [i].c [c] + agent [i].velocity [c];

      // Применяем ограничения поискового пространства
      a [i].c [c] = u.SeInDiSp (x, rangeMin [c], rangeMax [c], rangeStep [c]);

      // Проверяем ограничения и применяем плавную коррекцию
      double violation = CalculateConstraintValue (i, c);
      if (violation > eps)
      {
        double gradient = CalculateWeightedGradient (i, c);
        double correction = -gradient * violation * (1.0 - globalPhase);
        a [i].c [c] = u.SeInDiSp (a [i].c [c] + correction, rangeMin [c], rangeMax [c], rangeStep [c]);

        // Сбрасываем скорость при коррекции нарушений
        agent [i].velocity [c] *= 0.5;
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_COA_chaos::SecondCarrierWaveSearch ()
{
  // Уточняющий локальный поиск с адаптивными параметрами
  double localPhase = (double)(epochNow - S1) / S2;
  double intensificationRate = localPhase * localPhase; // Квадратичное увеличение интенсификации

  // Проверка на сходимость алгоритма
  bool isConverged = IsConverged ();

  // Для каждого агента
  for (int i = 0; i < popSize; i++)
  {
    // Если обнаружена сходимость, добавляем случайную мутацию к некоторым агентам
    if (isConverged && i % 3 == 0)
    {
      ApplyMutation (i);
      continue;
    }

    for (int c = 0; c < coords; c++)
    {
      // Выбор хаотического отображения с равномерным распределением
      int mapType = ((i * c + epochNow) % 3);
      agent [i].gamma [c] = SelectChaosMap (agent [i].gamma [c], mapType);

      // Адаптивный радиус поиска с сужением к концу оптимизации
      double adaptiveAlpha = alpha [c] * (1.0 - 0.8 * intensificationRate);

      // Выбор базовой точки с приоритетом лучших решений
      double basePoint;
      if (a [i].f > a [i].fB)
      {
        basePoint = a [i].c [c];  // Текущее положение лучше
      }
      else
      {
        double r = u.RNDprobab ();

        if (r < 0.7 * (1.0 + intensificationRate)) // Увеличиваем притяжение к глобальному лучшему
        {
          basePoint = cB [c];  // Глобальное лучшее
        }
        else
        {
          basePoint = a [i].cB [c];  // Личное лучшее
        }
      }

      // Локальный поиск с хаотическим компонентом
      double chaosOffset = adaptiveAlpha * (2.0 * agent [i].gamma [c] - 1.0);

      // Добавляем шум Леви для случайных дальних прыжков (тяжелый хвост распределения)
      double levyNoise = 0.0;
      if (u.RNDprobab () < 0.1 * (1.0 - intensificationRate))
      {
        // Упрощенное приближение шума Леви
        double u1 = u.RNDprobab ();
        double u2 = u.RNDprobab ();

        if (u2 > 0.01) // Защита от деления на очень малые числа
        {
          levyNoise = 0.01 * u1 / MathPow (u2, 0.5) * adaptiveAlpha * (rangeMax [c] - rangeMin [c]);
        }
      }

      // Обновляем скорость с инерцией
      agent [i].velocity [c] = inertia * (1.0 - 0.5 * intensificationRate) * agent [i].velocity [c] +
                               (1.0 - inertia) * (chaosOffset + levyNoise);

      // Применяем скорость к позиции
      double x = basePoint + agent [i].velocity [c];

      // Ограничиваем позицию
      a [i].c [c] = u.SeInDiSp (x, rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_COA_chaos::ApplyMutation (int agentIdx)
{
  // Определяем количество координат для мутации (от 1 до 30% координат)
  int mutationCount = 1 + (int)(u.RNDprobab () * coords * 0.3);
  mutationCount = MathMin (mutationCount, coords);

  // Создаем массив индексов для мутации без повторений
  int mutationIndices [];
  ArrayResize (mutationIndices, coords);

  // Заполняем массив индексами
  for (int i = 0; i < coords; i++)
  {
    mutationIndices [i] = i;
  }

  // Перемешиваем индексы
  for (int i = coords - 1; i > 0; i--)
  {
    int j = (int)(u.RNDprobab () * (i + 1));
    if (j <= i) // Дополнительная проверка для безопасности
    {
      int temp = mutationIndices [i];
      mutationIndices [i] = mutationIndices [j];
      mutationIndices [j] = temp;
    }
  }

  // Применяем мутации к выбранным координатам
  for (int m = 0; m < mutationCount; m++)
  {
    int c = mutationIndices [m];

    // Различные типы мутаций для разнообразия
    double r = u.RNDprobab ();
    double x;

    if (r < 0.3)
    {
      // Полная случайная мутация
      x = rangeMin [c] + u.RNDprobab () * (rangeMax [c] - rangeMin [c]);
    }
    else
      if (r < 0.6)
      {
        // Мутация относительно глобального лучшего
        double offset = (u.RNDprobab () - 0.5) * (rangeMax [c] - rangeMin [c]) * 0.2;
        x = cB [c] + offset;
      }
      else
      {
        // Мутация с использованием хаотического отображения
        agent [agentIdx].gamma [c] = SelectChaosMap (agent [agentIdx].gamma [c], (epochNow + c) % 3);
        x = rangeMin [c] + agent [agentIdx].gamma [c] * (rangeMax [c] - rangeMin [c]);
      }

    // Сбрасываем скорость
    agent [agentIdx].velocity [c] = 0.0;

    // Применяем новое значение с проверкой диапазона
    a [agentIdx].c [c] = u.SeInDiSp (x, rangeMin [c], rangeMax [c], rangeStep [c]);
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_COA_chaos::UpdateSigma ()
{
  // Динамическая адаптация параметра штрафа
  // Начинаем с малого значения и увеличиваем его при необходимости

  if (epochNow == 1)
  {
    currentSigma = sigma * 0.5;
    return;
  }

  // Подсчет количества допустимых решений
  int feasibleCount = 0;
  for (int i = 0; i < popSize; i++)
  {
    if (IsFeasible (i))
    {
      feasibleCount++;
    }
  }

  double feasibleRatio = (double)feasibleCount / MathMax (1, popSize);

  // Адаптация параметра штрафа в зависимости от доли допустимых решений
  if (feasibleRatio < 0.3)
  {
    // Слишком мало допустимых решений - увеличиваем штраф
    currentSigma *= 1.2;
  }
  else
    if (feasibleRatio > 0.7)
    {
      // Слишком много допустимых решений - уменьшаем штраф
      currentSigma *= 0.9;
    }

  // Ограничиваем значение sigma
  if (currentSigma < sigma * 0.1) currentSigma = sigma * 0.1;
  else
    if (currentSigma > sigma * 5.0) currentSigma = sigma * 5.0;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_COA_chaos::IsFeasible (int agentIdx)
{
  // Проверяем, находится ли решение в допустимой области
  for (int c = 0; c < coords; c++)
  {
    double violation = CalculateConstraintValue (agentIdx, c);
    if (violation > eps)
    {
      return false;
    }
  }
  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_COA_chaos::UpdateBestHistory (double newBest)
{
  // Защита от некорректных значений
  if (!MathIsValidNumber (newBest)) return;

  // Обновляем историю лучших значений
  globalBestHistory [historyIndex] = newBest;
  historyIndex = (historyIndex + 1) % 10;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_COA_chaos::IsConverged ()
{
  // Проверка достаточного количества данных в истории
  int validValues = 0;
  double sum      = 0.0;
  double minVal   = DBL_MAX;
  double maxVal   = -DBL_MAX;

  // Находим min, max и sum значений в истории
  for (int i = 0; i < 10; i++)
  {
    if (globalBestHistory [i] == -DBL_MAX || !MathIsValidNumber (globalBestHistory [i])) continue;

    minVal = MathMin (minVal, globalBestHistory [i]);
    maxVal = MathMax (maxVal, globalBestHistory [i]);
    sum += globalBestHistory [i];
    validValues++;
  }

  // Если недостаточно данных или все значения одинаковые
  if (validValues < 5 || minVal == maxVal) return false;

  // Вычисляем среднее значение
  double average = sum / validValues;

  // Проверка случая, когда среднее близко к нулю
  if (MathAbs (average) < eps) return MathAbs (maxVal - minVal) < eps * 10.0;

  // Относительная разница для проверки сходимости
  double relDiff = MathAbs (maxVal - minVal) / MathAbs (average);

  return relDiff < 0.001; // Порог сходимости - 0.1%
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_COA_chaos::ResetStagnatingAgents ()
{
  int resetCount = 0;

  for (int i = 0; i < popSize; i++)
  {
    // Увеличиваем счетчик стагнации, если нет улучшения
    if (a [i].f <= a [i].fB)
    {
      agent [i].stagnationCounter++;
    }
    else
    {
      agent [i].stagnationCounter = 0;
    }

    // Сбрасываем агента, если он находится в стагнации слишком долго
    if (agent [i].stagnationCounter > 5)
    {
      // Сброс только части агентов с вероятностью, зависящей от стагнации
      double resetProb = 0.2 * (1.0 + (double)agent [i].stagnationCounter / 10.0);

      if (u.RNDprobab () < resetProb)
      {
        ApplyMutation (i);
        agent [i].stagnationCounter = 0;
        resetCount++;
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_COA_chaos::CalculateWeightedGradient (int agentIdx, int coordIdx)
{
  // Вычисляем максимальное значение нарушения ограничений
  double maxViolation = eps;

  for (int c = 0; c < coords; c++)
  {
    double violation = CalculateConstraintValue (agentIdx, c);
    maxViolation = MathMax (maxViolation, violation);
  }

  // Нарушение для текущей координаты
  double violation = CalculateConstraintValue (agentIdx, coordIdx);

  // Если нет значительного нарушения, возвращаем 0
  if (violation <= eps) return 0.0;

  // Вычисляем градиент
  double gradient = CalculateConstraintGradient (agentIdx, coordIdx);

  // Нормализуем влияние в зависимости от степени нарушения
  double weight = violation / maxViolation;

  return gradient * weight;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_COA_chaos::CalculateConstraintValue (int agentIdx, int coordIdx)
{
  double x = a [agentIdx].c [coordIdx];
  double min = rangeMin [coordIdx];
  double max = rangeMax [coordIdx];

  // Сглаженная функция нарушения ограничений с плавным переходом на границе
  double violation = 0.0;

  // Проверка нижней границы
  if (x < min)
  {
    violation += min - x;
  }
  // Проверка верхней границы
  else
    if (x > max)
    {
      violation += x - max;
    }

  return violation;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_COA_chaos::CalculateConstraintGradient (int agentIdx, int coordIdx)
{
  double x = a [agentIdx].c [coordIdx];
  double min = rangeMin [coordIdx];
  double max = rangeMax [coordIdx];

  // Сглаженный градиент для лучшей стабильности
  if (x < min) return -1.0; // Нужно увеличивать значение
  else
    if (x > max) return 1.0;  // Нужно уменьшать значение

  return 0.0;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_COA_chaos::CalculatePenaltyFunction (int agentIdx)
{
  // Базовое значение целевой функции
  double baseValue = a [agentIdx].f;

  // Штрафной терм
  double penaltySum = 0.0;

  for (int c = 0; c < coords; c++)
  {
    double violation = CalculateConstraintValue (agentIdx, c);

    // Квадратичный штраф для лучшего разрешения близких к границе решений
    if (violation > eps)
    {
      penaltySum += violation * violation;
    }
  }

  // Применяем динамический коэффициент штрафа
  double penalty = currentSigma * penaltySum;

  // Для задачи максимизации: f_penalty = f - penalty
  return baseValue - penalty;
}
//——————————————————————————————————————————————————————————————————————————————