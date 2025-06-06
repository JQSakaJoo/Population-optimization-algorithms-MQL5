//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_BBO |
//|                                            Copyright 2007-2025, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/18354

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
class C_AO_BBO : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_BBO () { }
  C_AO_BBO ()
  {
    ao_name = "BBO";
    ao_desc = "Biogeography-Based Optimization";
    ao_link = "https://www.mql5.com/ru/articles/18354";

    popSize        = 50;     // размер популяции (количество хабитатов)
    immigrationMax = 1.0;    // максимальная скорость иммиграции
    emigrationMax  = 1.0;    // максимальная скорость эмиграции
    mutationProb   = 0.5;    // вероятность мутации
    elitismCount   = 2;      // количество элитных решений
    speciesMax     = 50;     // максимальное количество видов

    ArrayResize (params, 6);

    params [0].name = "popSize";        params [0].val = popSize;
    params [1].name = "immigrationMax"; params [1].val = immigrationMax;
    params [2].name = "emigrationMax";  params [2].val = emigrationMax;
    params [3].name = "mutationProb";   params [3].val = mutationProb;
    params [4].name = "elitismCount";   params [4].val = elitismCount;
    params [5].name = "speciesMax";     params [5].val = speciesMax;
  }

  void SetParams ()
  {
    popSize        = (int)params [0].val;
    immigrationMax = params      [1].val;
    emigrationMax  = params      [2].val;
    mutationProb   = params      [3].val;
    elitismCount   = (int)params [4].val;
    speciesMax     = (int)params [5].val;
  }

  bool Init (const double &rangeMinP  [],  // минимальные значения
             const double &rangeMaxP  [],  // максимальные значения
             const double &rangeStepP [],  // шаг изменения
             const int     epochsP = 0);   // количество эпох

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  double immigrationMax;    // максимальная скорость иммиграции
  double emigrationMax;     // максимальная скорость эмиграции
  double mutationProb;      // вероятность мутации
  int    elitismCount;      // количество элитных решений
  int    speciesMax;        // максимальное количество видов

  private: //-------------------------------------------------------------------
  struct S_HabitatData
  {
      int    speciesCount;     // количество видов в хабитате
      double immigration;      // скорость иммиграции
      double emigration;       // скорость эмиграции
      double probability;      // вероятность существования
  };

  S_HabitatData habitatData   [];  // данные для каждого хабитата
  double        probabilities [];  // вероятности для подсчета мутаций

  // Вспомогательные методы
  void   InitializePopulation ();
  void   CalculateRates       ();
  void   Migration            ();
  void   Mutation             ();
  double CalculateProbability (int speciesCount);
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_BBO::Init (const double &rangeMinP  [],  // минимальные значения
                     const double &rangeMaxP  [],  // максимальные значения
                     const double &rangeStepP [],  // шаг изменения
                     const int     epochsP = 0)    // количество эпох
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  // Инициализация массивов для каждого хабитата
  ArrayResize (habitatData,   popSize);
  ArrayResize (probabilities, speciesMax + 1);

  // Расчет вероятностей для различного количества видов
  double sum = 0.0;
  for (int i = 0; i <= speciesMax; i++)
  {
    probabilities [i] = CalculateProbability (i);
    sum += probabilities [i];
  }

  // Нормализация вероятностей
  if (sum > 0)
  {
    for (int i = 0; i <= speciesMax; i++)
    {
      probabilities [i] /= sum;
    }
  }

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Основной метод оптимизации                                                 |
//+----------------------------------------------------------------------------+
void C_AO_BBO::Moving ()
{
  // Первая итерация - инициализация начальной популяции
  if (!revision)
  {
    InitializePopulation ();
    revision = true;
    return;
  }

  // Основной процесс оптимизации
  // 1. Сортировка популяции по HSI (fitness)
  static S_AO_Agent aTemp []; ArrayResize (aTemp, popSize);
  u.Sorting (a, aTemp, popSize);

  // 2. Расчет скоростей иммиграции и эмиграции
  CalculateRates ();

  // 3. Миграция (обмен SIV между хабитатами)
  Migration ();

  // 4. Мутация на основе вероятностей
  Mutation ();

  // 5. Сохранение состояния для следующей итерации
  for (int i = 0; i < popSize; i++)
  {
    a [i].fP = a [i].f;
  }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Обновление лучшего решения                                                 |
//+----------------------------------------------------------------------------+
void C_AO_BBO::Revision ()
{
  // Поиск лучшего решения в текущей популяции
  for (int i = 0; i < popSize; i++)
  {
    // Обновление лучшего решения
    if (a [i].f > fB)
    {
      fB = a [i].f;
      ArrayCopy (cB, a [i].c, 0, 0, WHOLE_ARRAY);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Инициализация начальной популяции                                          |
//+----------------------------------------------------------------------------+
void C_AO_BBO::InitializePopulation ()
{
  // Инициализация начальной популяции равномерно по всему пространству
  for (int i = 0; i < popSize; i++)
  {
    for (int c = 0; c < coords; c++)
    {
      // Генерация случайных координат в допустимых пределах
      a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
      // Округление до ближайшего допустимого шага
      a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }

    // Инициализация данных хабитата
    habitatData [i].speciesCount = 0;
    habitatData [i].immigration  = 0.0;
    habitatData [i].emigration   = 0.0;
    habitatData [i].probability  = 0.0;
  }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Расчет скоростей иммиграции и эмиграции                                    |
//+----------------------------------------------------------------------------+
void C_AO_BBO::CalculateRates ()
{
  // Для линейной модели миграции
  for (int i = 0; i < popSize; i++)
  {
    // Количество видов обратно пропорционально рангу (лучшие решения имеют больше видов)
    habitatData [i].speciesCount = speciesMax - (i * speciesMax / popSize);

    // Скорость иммиграции уменьшается с увеличением количества видов
    habitatData [i].immigration = immigrationMax * (1.0 - (double)habitatData [i].speciesCount / speciesMax);

    // Скорость эмиграции увеличивается с увеличением количества видов
    habitatData [i].emigration = emigrationMax * (double)habitatData [i].speciesCount / speciesMax;

    // Вероятность существования хабитата
    if (habitatData [i].speciesCount <= speciesMax)
    {
      habitatData [i].probability = probabilities [habitatData [i].speciesCount];
    }
    else
    {
      habitatData [i].probability = 0.0;
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Миграция (обмен SIV между хабитатами)                                      |
//+----------------------------------------------------------------------------+
void C_AO_BBO::Migration ()
{
  for (int i = 0; i < popSize; i++)
  {
    // Пропускаем элитные решения
    if (i < elitismCount) continue;

    // Определяем, будет ли хабитат модифицирован
    if (u.RNDprobab () < habitatData [i].immigration)
    {
      // Для каждой координаты (SIV)
      for (int c = 0; c < coords; c++)
      {
        // Определяем, будет ли эта координата модифицирована
        if (u.RNDprobab () < habitatData [i].immigration)
        {
          // Выбор источника миграции на основе скоростей эмиграции
          double sumEmigration = 0.0;
          for (int j = 0; j < popSize; j++)
          {
            if (j != i) sumEmigration += habitatData [j].emigration;
          }

          if (sumEmigration > 0)
          {
            // Рулеточная селекция источника
            double roulette = u.RNDprobab () * sumEmigration;
            double cumSum = 0.0;

            for (int j = 0; j < popSize; j++)
            {
              if (j != i)
              {
                cumSum += habitatData [j].emigration;
                if (roulette <= cumSum)
                {
                  // Копирование SIV из хабитата j в хабитат i
                  a [i].c [c] = a [j].c [c];
                  break;
                }
              }
            }
          }
        }
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Мутация на основе вероятностей                                             |
//+----------------------------------------------------------------------------+
void C_AO_BBO::Mutation ()
{
  for (int i = 0; i < popSize; i++)
  {
    // Пропускаем элитные решения
    if (i < elitismCount) continue;

    // Скорость мутации обратно пропорциональна вероятности существования
    double mutationRate = mutationProb * (1.0 - habitatData [i].probability);

    if (u.RNDprobab () < mutationRate)
    {
      // Выбираем случайную координату для мутации
      int mutateCoord = MathRand () % coords;

      // Генерируем новое значение для выбранной координаты
      a [i].c [mutateCoord] = u.RNDfromCI (rangeMin [mutateCoord], rangeMax [mutateCoord]);
      a [i].c [mutateCoord] = u.SeInDiSp (a [i].c [mutateCoord],
                                          rangeMin [mutateCoord],
                                          rangeMax [mutateCoord],
                                          rangeStep [mutateCoord]);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Расчет вероятности для заданного количества видов                          |
//+----------------------------------------------------------------------------+
double C_AO_BBO::CalculateProbability (int speciesCount)
{
  // Упрощенная модель вероятности
  // Максимальная вероятность в середине диапазона (равновесие)
  int equilibrium = speciesMax / 2;
  double distance = MathAbs (speciesCount - equilibrium);
  double probability = MathExp (-distance * distance / (2.0 * equilibrium * equilibrium));

  return probability;
}
//——————————————————————————————————————————————————————————————————————————————