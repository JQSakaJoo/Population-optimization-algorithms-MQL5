//+——————————————————————————————————————————————————————————————————+
//|                                                         C_AO_EMA |
//|                                  Copyright 2007-2025, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/18605

#include "#C_AO.mqh"

//————————————————————————————————————————————————————————————————————
class C_AO_EMA : public C_AO
{
  public: //----------------------------------------------------------
  ~C_AO_EMA () { }
  C_AO_EMA ()
  {
    ao_name = "EMA";
    ao_desc = "Exchange Market Algorithm";
    ao_link = "https://www.mql5.com/ru/articles/18605";

    popSize   = 60;     // Размер популяции
    r1        = 1.5;    // Коэффициент поглощения для группы 2
    r2        = 0.8;    // Коэффициент поглощения для группы 3
    riskAlpha = 0.3;    // Фактор риска

    ArrayResize (params, 4);

    params [0].name = "popSize";   params [0].val = popSize;
    params [1].name = "r1";        params [1].val = r1;
    params [2].name = "r2";        params [2].val = r2;
    params [3].name = "riskAlpha"; params [3].val = riskAlpha;
  }

  void SetParams ()
  {
    popSize   = (int)params [0].val;
    r1        = params      [1].val;
    r2        = params      [2].val;
    riskAlpha = params      [3].val;

    // Проверка корректности параметров
    if (popSize < 6) popSize = 6;
    if (popSize % 3 != 0) popSize = ((popSize / 3) + 1) * 3; // Кратно 3
    if (r1 < 0.0) r1 = 0.0;
    if (r2 < 0.0) r2 = 0.0;
    if (riskAlpha < 0.0) riskAlpha = 0.0;
    if (riskAlpha > 1.0) riskAlpha = 1.0;
  }

  bool Init (const double &rangeMinP  [],
             const double &rangeMaxP  [],
             const double &rangeStepP [],
             const int     epochsP = 0);

  void Moving   ();
  void Revision ();

  //------------------------------------------------------------------
  double r1;         // Коэффициент поглощения 1
  double r2;         // Коэффициент поглощения 2
  double riskAlpha;  // Фактор риска

  private: //---------------------------------------------------------
  int    groupSize;     // размер каждой группы (popSize/3)
  int    currentEpoch;  // текущая эпоха
  int    totalEpochs;   // общее количество эпох

  S_AO_Agent tempPop []; // временная популяция для хранения новых позиций

  double GetDecayRate   ();
};
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Инициализация
bool C_AO_EMA::Init (const double &rangeMinP  [],
                     const double &rangeMaxP  [],
                     const double &rangeStepP [],
                     const int epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //------------------------------------------------------------------
  groupSize    = popSize / 3;
  currentEpoch = 0;
  totalEpochs  = epochsP;

  // Инициализация временной популяции
  ArrayResize (tempPop, popSize);
  for (int i = 0; i < popSize; i++)
  {
    tempPop [i].Init (coords);
  }

  return true;
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Основной цикл алгоритма
void C_AO_EMA::Moving ()
{
  // Начальная инициализация
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

  currentEpoch++;
  double decayRate = GetDecayRate ();

  // Копирование текущей популяции во временную
  for (int i = 0; i < popSize; i++)
  {
    ArrayCopy (tempPop [i].c, a [i].c, 0, 0, WHOLE_ARRAY);
    tempPop [i].f = a [i].f;
  }

  // ФАЗА 1: Сбалансированный рынок (Поглощающие операторы)

  // Группа 1 (элита) - не изменяется
  // Индексы: 0 ... groupSize-1

  // Группа 2 - поглощающий оператор 1
  // Индексы: groupSize ... 2*groupSize-1

  double adaptiveR1 = r1 * (1.0 - decayRate * 0.5);

  for (int i = groupSize; i < 2 * groupSize; i++)
  {
    // Каждый агент группы 2 выбирает случайного лидера из группы 1
    int leaderIdx = u.RNDminusOne (groupSize);

    for (int c = 0; c < coords; c++)
    {
      tempPop [i].c [c] = a [i].c [c] + u.RNDprobab () * adaptiveR1 * (a [leaderIdx].c [c] - a [i].c [c]);
      tempPop [i].c [c] = u.SeInDiSp (tempPop [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }

  // Группа 3 - поглощающий оператор 2
  // Индексы: 2*groupSize ... popSize-1

  adaptiveR1 = r1 * (1.0 - decayRate * 0.3);

  for (int i = 2 * groupSize; i < popSize; i++)
  {
    // Выбор лидеров из групп 1 и 2
    int leader1Idx = u.RNDminusOne (groupSize);
    int leader2Idx = groupSize + u.RNDminusOne (groupSize);

    for (int c = 0; c < coords; c++)
    {
      tempPop [i].c [c] = a [i].c [c] +
                          u.RNDprobab () * adaptiveR1 * (a [leader1Idx].c [c] - a [i].c [c]) +
                          u.RNDprobab () * adaptiveR1 * (a [leader2Idx].c [c] - a [i].c [c]);

      tempPop [i].c [c] = u.SeInDiSp (tempPop [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }

  // ФАЗА 2: Колеблющийся рынок (Поисковые операторы)

  // Группа 2 - поисковый оператор 1 (умеренный риск)
  for (int i = groupSize; i < 2 * groupSize; i++)
  {
    for (int c = 0; c < coords; c++)
    {
      double range = rangeMax [c] - rangeMin [c];

      if (u.RNDprobab () < 0.5)
      {
        tempPop [i].c [c] = cB [c];// tempPop [i].c [c] + delta;
      }
      else
      {
        // Поиск вокруг центра элиты
        double eliteCenter = 0.0;
        for (int j = 0; j < groupSize; j++)
        {
          eliteCenter += a [j].c [c];
        }
        eliteCenter /= (double)groupSize;

        double noise = riskAlpha * range * u.RNDfromCI (-0.5, 0.5) * (1.0 - decayRate * 0.5);
        tempPop [i].c [c] = eliteCenter + noise;
      }

      // Проверка границ
      tempPop [i].c [c] = u.SeInDiSp (tempPop [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }

  // Группа 3 - поисковый оператор 2 (высокий риск)
  for (int i = 2 * groupSize; i < popSize; i++)
  {
    for (int c = 0; c < coords; c++)
    {
      double range = rangeMax [c] - rangeMin [c];

      if (u.RNDprobab () < riskAlpha)
      {
        // Полная реинициализация
        tempPop [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
      }
      else
        if (u.RNDprobab () < 0.5 + riskAlpha / 2.0)
        {
          // Широкий поиск
          double searchRadius = 2.0 * riskAlpha * range * (1.0 - decayRate * 0.3);
          double delta = u.RNDfromCI (-searchRadius, searchRadius);

          tempPop [i].c [c] = tempPop [i].c [c] + delta;
        }
        else
        {
          // Оппозиционное обучение
          double worstCenter = 0.0;
          int worstCount = groupSize / 2;
          for (int j = popSize - worstCount; j < popSize; j++)
          {
            worstCenter += a [j].c [c];
          }
          worstCenter /= (double)worstCount;

          // Движение в противоположном направлении от худших
          tempPop [i].c [c] = 2.0 * tempPop [i].c [c] - worstCenter;

          // Добавляем небольшой шум
          double noise = riskAlpha * range * u.RNDfromCI (-0.1, 0.1);
          tempPop [i].c [c] += noise;
        }

      // Проверка границ
      tempPop [i].c [c] = u.SeInDiSp (tempPop [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }

  // Копирование из временной популяции в основную (кроме группы 1)
  for (int i = groupSize; i < popSize; i++)
  {
    ArrayCopy (a [i].c, tempPop [i].c, 0, 0, WHOLE_ARRAY);
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Получение коэффициента затухания
double C_AO_EMA::GetDecayRate ()
{
  if (totalEpochs <= 0) return 0.0;

  // Нелинейное затухание для лучшего баланса эксплуатации/исследования
  double progress = (double)currentEpoch / (double)totalEpochs;

  // Использование сигмоидной функции для плавного перехода
  return 1.0 / (1.0 + MathExp (-10.0 * (progress - 0.5)));
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Обновление лучших решений
void C_AO_EMA::Revision ()
{
  static S_AO_Agent aT []; ArrayResize (aT, popSize);
  u.Sorting (a, aT, popSize);
  ArrayCopy (cB, a [0].c, 0, 0, WHOLE_ARRAY);
  fB = a [0].f;
}
//————————————————————————————————————————————————————————————————————
