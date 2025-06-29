//+——————————————————————————————————————————————————————————————————+
//|                                                         C_AO_DEA |
//|                                  Copyright 2007-2025, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/18495

#include "#C_AO.mqh"

//————————————————————————————————————————————————————————————————————
struct S_Alternative
{
    double value;     // значение альтернативы
    double AF;        // накопленная пригодность для этой альтернативы
};
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
struct S_Coordinate
{
    S_Alternative alt [];  // массив альтернатив для данной координаты
    int           count;   // количество альтернатив
};
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
class C_AO_DEA : public C_AO
{
  public: //----------------------------------------------------------
  ~C_AO_DEA () { }
  C_AO_DEA ()
  {
    ao_name = "DEA";
    ao_desc = "Dolphin Echolocation Algorithm";
    ao_link = "https://www.mql5.com/ru/articles/18495";

    popSize = 100;    // NL - количество локаций (дельфинов)
    Re      = 2;      // эффективный радиус поиска
    Power   = 2.0;    // степень кривой сходимости
    PP1     = 1.0;    // фактор сходимости первой итерации

    ArrayResize (params, 4);

    params [0].name = "popSize"; params [0].val = popSize;
    params [1].name = "Re";      params [1].val = Re;
    params [2].name = "Power";   params [2].val = Power;
    params [3].name = "PP1";     params [3].val = PP1;
  }

  void SetParams ()
  {
    popSize = (int)params [0].val;
    Re      = (int)params [1].val;
    Power   = params      [2].val;
    PP1     = params      [3].val;

    // Проверка корректности параметров
    if (Re < 0) Re = 0;
    if (PP1 < 0.0) PP1 = 0.0;
    if (PP1 > 1.0) PP1 = 1.0;
    if (Power < 0.1) Power = 0.1;
  }

  bool Init (const double &rangeMinP  [],  // минимальные значения
             const double &rangeMaxP  [],  // максимальные значения
             const double &rangeStepP [],  // шаг изменения
             const int     epochsP = 0);   // количество эпох

  void Moving   ();
  void Revision ();

  //------------------------------------------------------------------
  int    Re;           // эффективный радиус поиска
  double Power;        // степень кривой сходимости
  double PP1;          // фактор сходимости первой итерации

  private: //---------------------------------------------------------
  double PP;           // предопределенная вероятность для текущей итерации
  int    currentEpoch; // текущая эпоха
  int    totalEpochs;  // общее количество эпох
  double coeffA;       // динамический коэффициент для выбора позиций

  S_Coordinate coordData []; // данные по координатам (альтернативы и AF)

  void CalculatePP ();
  void CalculateAccumulativeFitness ();
  void ResetAFforBestLocation ();
  void SelectNextLocations    ();
  int  FindNearestAlternative (int coord, double value);
  void CalculateCoefficientA  ();
};
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Инициализация
bool C_AO_DEA::Init (const double &rangeMinP  [],
                     const double &rangeMaxP  [],
                     const double &rangeStepP [],
                     const int epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //------------------------------------------------------------------
  currentEpoch = 0;
  totalEpochs  = epochsP;

  // Инициализация структуры данных для координат
  ArrayResize (coordData, coords);

  // Создаем альтернативы для каждой координаты
  for (int c = 0; c < coords; c++)
  {
    if (rangeStep [c] != 0)
    {
      coordData [c].count = (int)((rangeMax [c] - rangeMin [c]) / rangeStep [c]) + 1;
    }
    else
    {
      coordData [c].count = 500;
    }

    // Проверяем, что Re не слишком большой для количества альтернатив
    if (Re > coordData [c].count / 4) Re = coordData [c].count / 4;

    ArrayResize (coordData [c].alt, coordData [c].count);

    // Заполняем альтернативы
    for (int i = 0; i < coordData [c].count; i++)
    {
      if (rangeStep [c] != 0)
      {
        coordData [c].alt [i].value = rangeMin [c] + i * rangeStep [c];
      }
      else
      {
        coordData [c].alt [i].value = rangeMin [c] + (rangeMax [c] - rangeMin [c]) * i / (coordData [c].count - 1);
      }
      coordData [c].alt [i].AF = 0.0;
    }
  }

  return true;
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Основной шаг алгоритма (согласно Algorithm 1)
void C_AO_DEA::Moving ()
{
  // Начальная инициализация
  if (!revision)
  {
    for (int p = 0; p < popSize; p++)
    {
      for (int c = 0; c < coords; c++)
      {
        a [p].c  [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        a [p].c  [c] = u.SeInDiSp (a [p].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
        a [p].cB [c] = a [p].c [c];
      }
    }

    revision = true;
    return;
  }

  // Увеличиваем счетчик эпох
  currentEpoch++;

  // Шаги алгоритма DEA согласно Algorithm 1:

  // 1. Вычисляем PP для текущей итерации
  CalculatePP ();

  // 2. Рассчитываем динамический коэффициент a
  CalculateCoefficientA ();

  // 3. Вычисляем накопленную пригодность
  CalculateAccumulativeFitness ();

  // 4. Находим лучшую локацию и сбрасываем для нее AF
  ResetAFforBestLocation ();

  // 5. Выбираем следующие позиции
  SelectNextLocations ();
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Вычисление предопределенной вероятности (согласно формуле 5)
void C_AO_DEA::CalculatePP ()
{
  if (totalEpochs <= 1)
  {
    PP = PP1;
    return;
  }

  // PP = PP1 + (1 - PP1) * ((Loop^Power - 1) / (LoopsNumber^Power - 1))
  double iterPower  = MathPow ((double)currentEpoch, Power) - 1.0;
  double totalPower = MathPow ((double)totalEpochs,  Power) - 1.0;

  if (totalPower != 0)
  {
    PP = PP1 + (1.0 - PP1) * iterPower / totalPower;
  }
  else
  {
    PP = PP1;
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Расчет динамического коэффициента a
void C_AO_DEA::CalculateCoefficientA ()
{
  double sumFitness = 0.0;

  for (int i = 0; i < popSize; i++)
  {
    sumFitness += fB - a [i].f;
  }

  coeffA = (fB - fW) / sumFitness;
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Поиск ближайшей альтернативы
int C_AO_DEA::FindNearestAlternative (int coord, double value)
{
  int nearest = 0;
  double minDist = DBL_MAX;

  for (int i = 0; i < coordData [coord].count; i++)
  {
    double dist = MathAbs (value - coordData [coord].alt [i].value);
    if (dist < minDist)
    {
      minDist = dist;
      nearest = i;
    }
  }

  return nearest;
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Вычисление накопленной пригодности (согласно Algorithm 2)
void C_AO_DEA::CalculateAccumulativeFitness ()
{
  // Очищаем накопленную пригодность для всех альтернатив
  for (int c = 0; c < coords; c++)
  {
    for (int i = 0; i < coordData [c].count; i++)
    {
      coordData [c].alt [i].AF = 0.0;
    }
  }

  double rangeFF = fB - fW;
  if (rangeFF == 0.0) rangeFF = DBL_EPSILON;

  // Для каждого агента (дельфина)
  for (int i = 0; i < popSize; i++)
  {
    // Нормализуем fitness для данного агента
    double normalizedFitness = (a [i].f - fW) / rangeFF;

    for (int c = 0; c < coords; c++)
    {
      // Находим ближайшую альтернативу для текущей координаты
      int nearestAlt = FindNearestAlternative (c, a [i].c [c]);

      // Обновляем накопленную пригодность в радиусе Re
      // Согласно Algorithm 2: AF(A+k)j = (1/Re) * (Re - |k|) * fitness(i) + AF(A+k)j
      for (int k = -Re; k <= Re; k++)
      {
        int altIndex = nearestAlt + k;

        // Проверка границ с учетом отражения (reflective characteristic)
        if (altIndex < 0)
        {
          altIndex = -altIndex; // отражение от нижней границы
        }
        else if (altIndex >= coordData [c].count)
        {
          altIndex = 2 * (coordData [c].count - 1) - altIndex; // отражение от верхней границы
        }

        if (altIndex >= 0 && altIndex < coordData [c].count)
        {
          double weight = (1.0 / (double)(Re + 1)) * (Re - MathAbs (k) + 1);
          coordData [c].alt [altIndex].AF += weight * normalizedFitness;
        }
      }
    }
  }

  // Добавляем малое значение epsilon ко всем AF для избежания нулевых вероятностей
  double epsilon = 0.0001;
  for (int c = 0; c < coords; c++)
  {
    for (int i = 0; i < coordData [c].count; i++)
    {
      coordData [c].alt [i].AF += epsilon;
    }
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Сброс AF для лучшей локации (согласно Algorithm 3)
void C_AO_DEA::ResetAFforBestLocation ()
{
  // Находим индекс лучшего решения
  int bestIndex = 0;
  double bestFitness = a [0].f;

  // Ищем решение с максимальным fitness (т.к. мы всегда максимизируем нормализованный fitness)
  for (int i = 1; i < popSize; i++)
  {
    if (a [i].f > bestFitness)
    {
      bestFitness = a [i].f;
      bestIndex = i;
    }
  }

  // Обнуляем AF для ВСЕХ альтернатив, соответствующих координатам лучшего решения
  for (int c = 0; c < coords; c++)
  {
    // Находим ближайшую альтернативу к координате лучшего решения
    int nearestAlt = FindNearestAlternative (c, a [bestIndex].c [c]);

    // Обнуляем AF только для этой альтернативы
    if (nearestAlt >= 0 && nearestAlt < coordData [c].count)
    {
      coordData [c].alt [nearestAlt].AF = 0.0;
    }
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Выбор следующих позиций на основе вероятностей
void C_AO_DEA::SelectNextLocations ()
{
  // Сначала находим индекс лучшего решения
  int bestIndex = 0;
  double bestFitness = a [0].f;

  for (int i = 1; i < popSize; i++)
  {
    if (a [i].f > bestFitness)
    {
      bestFitness = a [i].f;
      bestIndex = i;
    }
  }

  for (int i = 0; i < popSize; i++)
  {
    for (int c = 0; c < coords; c++)
    {
      // Для лучшей позиции применяем PP
      if (i == bestIndex && u.RNDprobab () < PP)
      {
        // Сохраняем текущее значение координаты лучшего решения с вероятностью PP
        continue;
      }

      // Выбираем на основе накопленной пригодности
      double totalAF = 0.0;
      for (int alt = 0; alt < coordData [c].count; alt++)
      {
        totalAF += coordData [c].alt [alt].AF;
      }

      if (totalAF > DBL_EPSILON) // Проверяем, что есть ненулевые AF
      {
        // Выбор альтернативы на основе рулетки
        double rnd = u.RNDprobab () * totalAF;
        double cumSum = 0.0;

        for (int alt = 0; alt < coordData [c].count; alt++)
        {
          cumSum += coordData [c].alt [alt].AF;
          if (cumSum >= rnd)
          {
            a [i].c [c] = coordData [c].alt [alt].value;
            break;
          }
        }
      }
      else
      {
        // Если все AF практически нулевые, используем случайный выбор
        // с динамическим коэффициентом coeffA для вероятности локального поиска
        if (u.RNDprobab () < coeffA) // Используем динамический коэффициент
        {
          // Локальный поиск - остаемся рядом с текущей позицией
          int currentAlt = FindNearestAlternative (c, a [i].c [c]);
          int shift = u.RNDminusOne (2 * Re + 1) - Re; // случайное смещение в пределах Re
          int newAlt = currentAlt + shift;

          // Проверка границ
          if (newAlt < 0) newAlt = 0;
          if (newAlt >= coordData [c].count) newAlt = coordData [c].count - 1;

          a [i].c [c] = coordData [c].alt [newAlt].value;
        }
        else
        {
          // Глобальный поиск - полностью случайный выбор
          int randAlt = u.RNDminusOne (coordData [c].count);
          a [i].c [c] = coordData [c].alt [randAlt].value;
        }
      }

      // Проверка границ и дискретизация
      a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Обновление лучшего решения
void C_AO_DEA::Revision ()
{
  int bestIND = -1;
  fW = fB;

  // Находим лучшее и худшее решения в текущей популяции
  for (int i = 0; i < popSize; i++)
  {
    // Обновляем лучшее решение
    if (a [i].f > fB)
    {
      fB = a [i].f;
      bestIND = i;
    }

    // Обновляем худшее решение
    if (a [i].f < fW)
    {
      fW = a [i].f;
    }
  }

  // Копируем координаты лучшего решения
  if (bestIND != -1)
  {
    ArrayCopy (cB, a [bestIND].c, 0, 0, WHOLE_ARRAY);
  }
}
//————————————————————————————————————————————————————————————————————
