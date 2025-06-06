//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_FBA |
//|                                            Copyright 2007-2025, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/17497

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
class C_AO_FBA : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_FBA () { }
  C_AO_FBA ()
  {
    ao_name = "FBA";
    ao_desc = "Fractal-Based Algorithm";
    ao_link = "https://www.mql5.com/ru/articles/17458";

    popSize = 50;  // размер популяции
    P1      = 60;  // процент перспективных точек
    P2      = 30;  // процент перспективных подпространств
    P3      = 0.8; // процент точек для случайной модификации
    m_value = 10;  // число интервалов для разделения каждого измерения

    ArrayResize (params, 5);

    params [0].name = "popSize"; params [0].val = popSize;
    params [1].name = "P1";      params [1].val = P1;
    params [2].name = "P2";      params [2].val = P2;
    params [3].name = "P3";      params [3].val = P3;
    params [4].name = "m_value"; params [4].val = m_value;
  }

  void SetParams ()
  {
    popSize = (int)params [0].val;
    P1      = (int)params [1].val;
    P2      = (int)params [2].val;
    P3      =      params [3].val;
    m_value = (int)params [4].val;
  }

  bool Init (const double &rangeMinP  [],  // минимальные значения
             const double &rangeMaxP  [],  // максимальные значения
             const double &rangeStepP [],  // шаг изменения
             const int     epochsP = 0);   // количество эпох

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  int    P1;           // процент перспективных точек
  int    P2;           // процент перспективных подпространств
  double P3;           // доля точек для случайной модификации
  int    m_value;      // число интервалов для разделения каждого измерения

  private: //-------------------------------------------------------------------

  // Структура для представления подпространства
  struct S_Subspace
  {
      double min [];        // минимальные границы подпространства
      double max [];        // максимальные границы подпространства
      double promisingRank; // ранг перспективности (нормализованное значение)
      bool   isPromising;   // флаг перспективности
      int    parentIndex;   // индекс родительского подпространства (-1 для корневых)
      int    level;         // уровень в иерархии (0 для исходного пространства)

      void Init (int coords)
      {
        ArrayResize (min, coords);
        ArrayResize (max, coords);
        promisingRank = 0.0;
        isPromising   = false;
        parentIndex   = -1;
        level         = 0;
      }
  };

  S_Subspace subspaces     []; // массив подпространств

  // Вспомогательные методы
  bool   IsPointInSubspace              (const double &point [], const S_Subspace &subspace);
  void   CreateInitialSpacePartitioning ();
  void   DivideSubspace                 (int subspaceIndex);
  void   IdentifyPromisingPoints        (int &promisingIndices []);
  void   CalculateSubspaceRanks         (const int &promisingIndices []);
  void   SelectPromisingSubspaces       ();
  void   DividePromisingSubspaces       ();
  void   GenerateNewPopulation          ();
  void   MutatePoints                   ();
  void   SortByFitness                  (double &values [], int &indices [], int size, bool ascending = false);
  void   QuickSort                      (double &values [], int &indices [], int low, int high, bool ascending);
  int    Partition                      (double &values [], int &indices [], int low, int high, bool ascending);
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_FBA::Init (const double &rangeMinP  [],  // минимальные значения
                     const double &rangeMaxP  [],  // максимальные значения
                     const double &rangeStepP [],  // шаг изменения
                     const int     epochsP = 0)    // количество эпох
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  // Создаем начальное разделение пространства поиска
  CreateInitialSpacePartitioning ();

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Основной метод оптимизации                                                 |
//+----------------------------------------------------------------------------+
void C_AO_FBA::Moving ()
{
  // Первая итерация - инициализация начальной популяции
  if (!revision)
  {
    // Инициализация начальной популяции равномерно по всему пространству
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

  // Основной процесс оптимизации

  // 1. Выявление перспективных точек (P1% точек с лучшими значениями функции)
  int promisingIndices [];
  IdentifyPromisingPoints (promisingIndices);

  // 2. Расчет рангов перспективности для каждого подпространства
  CalculateSubspaceRanks (promisingIndices);

  // 3. Выбор P2% самых перспективных подпространств
  SelectPromisingSubspaces ();

  // 4. Разделение перспективных подпространств на более мелкие
  DividePromisingSubspaces ();

  // 5. Генерация новых точек с учетом рангов перспективности
  GenerateNewPopulation ();

  // 6. Случайная модификация (мутация)
  MutatePoints ();
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Обновление лучшего решения                                                 |
//+----------------------------------------------------------------------------+
void C_AO_FBA::Revision ()
{
  // Поиск лучшего решения
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
//| Создание начального разделения пространства                                |
//+----------------------------------------------------------------------------+
void C_AO_FBA::CreateInitialSpacePartitioning ()
{
  // Создаем начальное разделение пространства
  int totalSubspaces = (int)MathPow (m_value, coords);

  // При очень большой размерности ограничиваем количество подпространств
  if (totalSubspaces > 10000) totalSubspaces = 10000;

  ArrayResize (subspaces, totalSubspaces);

  // Инициализируем все подпространства
  for (int i = 0; i < totalSubspaces; i++)
  {
    subspaces [i].Init (coords);
    subspaces [i].level = 0; // Начальный уровень
  }

  // Разделяем начальное пространство на равные подпространства
  int index = 0;

  // В зависимости от размерности пространства выбираем метод разделения
  if (coords == 1)
  {
    // Одномерный случай
    double intervalSize = (rangeMax [0] - rangeMin [0]) / m_value;

    for (int i = 0; i < m_value && index < totalSubspaces; i++)
    {
      subspaces [index].min [0] = rangeMin [0] + i * intervalSize;
      subspaces [index].max [0] = rangeMin [0] + (i + 1) * intervalSize;
      index++;
    }
  }
  else
    if (coords == 2)
    {
      // Двумерный случай
      double intervalSize0 = (rangeMax [0] - rangeMin [0]) / m_value;
      double intervalSize1 = (rangeMax [1] - rangeMin [1]) / m_value;

      for (int i = 0; i < m_value && index < totalSubspaces; i++)
      {
        for (int j = 0; j < m_value && index < totalSubspaces; j++)
        {
          subspaces [index].min [0] = rangeMin [0] + i * intervalSize0;
          subspaces [index].max [0] = rangeMin [0] + (i + 1) * intervalSize0;
          subspaces [index].min [1] = rangeMin [1] + j * intervalSize1;
          subspaces [index].max [1] = rangeMin [1] + (j + 1) * intervalSize1;
          index++;
        }
      }
    }
    else
    {
      // Многомерный случай - используем итеративный подход
      int indices [];
      ArrayResize (indices, coords);
      for (int i = 0; i < coords; i++) indices [i] = 0;

      while (index < totalSubspaces)
      {
        // Вычисляем границы текущего подпространства
        for (int c = 0; c < coords; c++)
        {
          double intervalSize = (rangeMax [c] - rangeMin [c]) / m_value;
          subspaces [index].min [c] = rangeMin [c] + indices [c] * intervalSize;
          subspaces [index].max [c] = rangeMin [c] + (indices [c] + 1) * intervalSize;
        }

        // Переходим к следующему подпространству
        int c = coords - 1;
        while (c >= 0)
        {
          indices [c]++;
          if (indices [c] < m_value) break;
          indices [c] = 0;
          c--;
        }

        // Если завершили полный цикл, выходим
        if (c < 0) break;

        index++;
      }
    }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Определение принадлежности точки подпространству                           |
//+----------------------------------------------------------------------------+
bool C_AO_FBA::IsPointInSubspace (const double &point [], const S_Subspace &subspace)
{
  // Проверяем, находится ли точка в указанном подпространстве
  for (int c = 0; c < coords; c++)
  {
    if (point [c] < subspace.min [c] || point [c] >= subspace.max [c])
    {
      return false;
    }
  }

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Выявление перспективных точек                                              |
//+----------------------------------------------------------------------------+
void C_AO_FBA::IdentifyPromisingPoints (int &promisingIndices [])
{
  // Выбираем P1% точек с лучшими значениями функции

  // Создаем массивы для сортировки
  double values  [];
  int    indices [];

  ArrayResize (values,  popSize);
  ArrayResize (indices, popSize);

  // Заполняем массивы
  for (int i = 0; i < popSize; i++)
  {
    values  [i] = a [i].f;
    indices [i] = i;
  }

  // Сортируем по убыванию (для задачи максимизации)
  SortByFitness (values, indices, popSize);

  // Выбираем P1% лучших точек
  int numPromisingPoints = (int)MathRound (popSize * P1 / 100.0);
  numPromisingPoints = MathMax (1, MathMin (numPromisingPoints, popSize));

  ArrayResize (promisingIndices, numPromisingPoints);

  for (int i = 0; i < numPromisingPoints; i++)
  {
    promisingIndices [i] = indices [i];
  }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Расчет рангов перспективности подпространств                               |
//+----------------------------------------------------------------------------+
void C_AO_FBA::CalculateSubspaceRanks (const int &promisingIndices [])
{
  // Сбрасываем ранги подпространств
  for (int i = 0; i < ArraySize (subspaces); i++)
  {
    subspaces [i].promisingRank = 0.0;
  }

  // Подсчитываем перспективные точки в каждом подпространстве
  for (int i = 0; i < ArraySize (promisingIndices); i++)
  {
    int pointIndex = promisingIndices [i];

    for (int j = 0; j < ArraySize (subspaces); j++)
    {
      if (IsPointInSubspace (a [pointIndex].c, subspaces [j]))
      {
        subspaces [j].promisingRank++;
        break; // Точка может находиться только в одном подпространстве
      }
    }
  }

  // Нормализуем ранги перспективности согласно статье
  // PromisingRank = Number of promising points in s / Total promising points
  int totalPromisingPoints = ArraySize (promisingIndices);
  if (totalPromisingPoints > 0)
  {
    for (int i = 0; i < ArraySize (subspaces); i++)
    {
      subspaces [i].promisingRank = subspaces [i].promisingRank / totalPromisingPoints;
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Выбор перспективных подпространств                                         |
//+----------------------------------------------------------------------------+
void C_AO_FBA::SelectPromisingSubspaces ()
{
  // Выбираем P2% подпространств с наивысшими рангами как перспективные

  // Создаем массивы для сортировки
  double ranks [];
  int indices [];

  int numSubspaces = ArraySize (subspaces);
  ArrayResize (ranks, numSubspaces);
  ArrayResize (indices, numSubspaces);

  // Заполняем массивы
  for (int i = 0; i < numSubspaces; i++)
  {
    ranks [i] = subspaces [i].promisingRank;
    indices [i] = i;
    // Сбрасываем флаг перспективности
    subspaces [i].isPromising = false;
  }

  // Сортируем по убыванию рангов
  SortByFitness (ranks, indices, numSubspaces);

  // Выбираем P2% самых перспективных подпространств
  int numPromisingSubspaces = (int)MathRound (numSubspaces * P2 / 100.0);
  numPromisingSubspaces = MathMax (1, MathMin (numPromisingSubspaces, numSubspaces));

  // Отмечаем перспективные подпространства
  for (int i = 0; i < numPromisingSubspaces && i < ArraySize (indices); i++)
  {
    // Защита от выхода за пределы массива
    if (indices [i] >= 0 && indices [i] < ArraySize (subspaces))
    {
      subspaces [indices [i]].isPromising = true;
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Разделение перспективных подпространств                                    |
//+----------------------------------------------------------------------------+
void C_AO_FBA::DividePromisingSubspaces ()
{
  // Собираем индексы перспективных подпространств
  int promisingIndices [];
  int numPromising = 0;

  for (int i = 0; i < ArraySize (subspaces); i++)
  {
    if (subspaces [i].isPromising)
    {
      numPromising++;
      ArrayResize (promisingIndices, numPromising);
      promisingIndices [numPromising - 1] = i;
    }
  }

  // Делим каждое перспективное подпространство
  for (int i = 0; i < numPromising; i++)
  {
    DivideSubspace (promisingIndices [i]);
  }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Разделение конкретного подпространства                                     |
//+----------------------------------------------------------------------------+
void C_AO_FBA::DivideSubspace (int subspaceIndex)
{
  // Делим указанное подпространство на m_value^coords подпространств

  S_Subspace parent = subspaces [subspaceIndex];

  // Ограничение на максимальное количество подпространств
  if (ArraySize (subspaces) > 10000) return;

  // Для каждого измерения делим на m_value частей
  int totalNewSubspaces = (int)MathPow (m_value, coords);
  int currentSize = ArraySize (subspaces);
  ArrayResize (subspaces, currentSize + totalNewSubspaces);

  // Создаем новые подпространства
  int newIndex = currentSize;
  int indices [];
  ArrayResize (indices, coords);
  for (int i = 0; i < coords; i++) indices [i] = 0;

  for (int idx = 0; idx < totalNewSubspaces && newIndex < ArraySize (subspaces); idx++)
  {
    subspaces [newIndex].Init (coords);
    subspaces [newIndex].level = parent.level + 1;
    subspaces [newIndex].parentIndex = subspaceIndex;

    // Вычисляем границы текущего подпространства
    for (int c = 0; c < coords; c++)
    {
      double intervalSize = (parent.max [c] - parent.min [c]) / m_value;
      subspaces [newIndex].min [c] = parent.min [c] + indices [c] * intervalSize;
      subspaces [newIndex].max [c] = parent.min [c] + (indices [c] + 1) * intervalSize;
    }

    // Переходим к следующему подпространству
    int c = coords - 1;
    while (c >= 0)
    {
      indices [c]++;
      if (indices [c] < m_value) break;
      indices [c] = 0;
      c--;
    }

    // Если завершили полный цикл, выходим
    if (c < 0) break;

    newIndex++;
  }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Генерация новой популяции                                                  |
//+----------------------------------------------------------------------------+
void C_AO_FBA::GenerateNewPopulation ()
{
  // Вычисляем сумму рангов всех подпространств
  double totalRank = 0.0;
  for (int i = 0; i < ArraySize (subspaces); i++)
  {
    totalRank += subspaces [i].promisingRank;
  }

  // Если все ранги равны 0, установим равномерное распределение
  if (totalRank <= 0.0001) // Проверка на приближенное равенство к нулю
  {
    for (int i = 0; i < ArraySize (subspaces); i++)
    {
      subspaces [i].promisingRank = 1.0;
    }
    totalRank = ArraySize (subspaces);
  }

  int points = 0;

  for (int i = 0; i < ArraySize (subspaces) && points < popSize; i++)
  {
    // Вычисляем количество точек для этого подпространства согласно формуле
    int pointsToGenerate = (int)MathRound ((subspaces [i].promisingRank / totalRank) * popSize);

    // Ограничение, чтобы не выйти за пределы
    pointsToGenerate = MathMin (pointsToGenerate, popSize - points);

    // Генерируем точки в этом подпространстве
    for (int j = 0; j < pointsToGenerate; j++)
    {
      // Создаем новую точку равномерно в пределах подпространства
      for (int c = 0; c < coords; c++)
      {
        a [points].c [c] = u.RNDfromCI (subspaces [i].min [c], subspaces [i].max [c]);
        a [points].c [c] = u.SeInDiSp (a [points].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }

      points++;
      if (points >= popSize) break;
    }
  }

  // Если не все точки были сгенерированы, заполняем оставшиеся равномерно по всему пространству
  while (points < popSize)
  {
    for (int c = 0; c < coords; c++)
    {
      a [points].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
      a [points].c [c] = u.SeInDiSp (a [points].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
    points++;
  }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Мутация точек в популяции                                                  |
//+----------------------------------------------------------------------------+
void C_AO_FBA::MutatePoints ()
{
  // Добавляем гауссовский шум к P3% случайно выбранных точек из новой популяции
  /*
  int numToMutate = (int)MathRound (popSize * P3 / 100.0);
  numToMutate = MathMax (1, MathMin (numToMutate, popSize));

  for (int i = 0; i < numToMutate; i++)
  {
    int index = u.RNDminusOne (popSize);

    // Добавляем шум к каждой координате
    for (int c = 0; c < coords; c++)
    {
      // Стандартное отклонение 10% от диапазона
      //double stdDev = (rangeMax [c] - rangeMin [c]) * 0.1;

      // Гауссовский шум с использованием метода Бокса-Мюллера
      //double noise = NormalRandom (0.0, stdDev);

      // Добавляем шум и ограничиваем значение
      a [index].c [c] += noise;
      a [index].c [c] = u.SeInDiSp (a [index].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }
  */

  for (int p = 0; p < popSize; p++)
  {
    for (int c = 0; c < coords; c++)
    {
      if (u.RNDprobab () < P3)
      {
        a [p].c [c] = u.PowerDistribution (cB [c], rangeMin [c], rangeMax [c], 20);
        a [p].c [c] = u.SeInDiSp (a [p].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

/*
//+----------------------------------------------------------------------------+
//| Генерация нормально распределенного случайного числа                       |
//+----------------------------------------------------------------------------+
double C_AO_FBA::NormalRandom (double mean, double stdDev)
{
  // Генерация случайного числа из нормального распределения (метод Бокса-Мюллера)

  // Генерируем два равномерно распределенных числа
  double u1 = u.RNDfromCI (0.0001, 0.9999); // Избегаем 0
  double u2 = u.RNDfromCI (0.0001, 0.9999); // Избегаем 0

  // Применяем преобразование Бокса-Мюллера для получения стандартного нормального распределения
  double z = MathSqrt (-2.0 * MathLog (u1)) * MathCos (2.0 * M_PI * u2);

  // Масштабируем и смещаем согласно требуемому среднему и стандартному отклонению
  return mean + z * stdDev;
}
//——————————————————————————————————————————————————————————————————————————————
*/
//+----------------------------------------------------------------------------+
//| Сортировка по значению фитнес-функции                                      |
//+----------------------------------------------------------------------------+
void C_AO_FBA::SortByFitness (double &values [], int &indices [], int size, bool ascending = false)
{
  if (size > 1) QuickSort (values, indices, 0, size - 1, ascending);
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Алгоритм быстрой сортировки                                                |
//+----------------------------------------------------------------------------+
void C_AO_FBA::QuickSort (double &values [], int &indices [], int low, int high, bool ascending)
{
  if (low < high)
  {
    int pi = Partition (values, indices, low, high, ascending);

    QuickSort (values, indices, low, pi - 1, ascending);
    QuickSort (values, indices, pi + 1, high, ascending);
  }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Функция разделения для QuickSort                                           |
//+----------------------------------------------------------------------------+
int C_AO_FBA::Partition (double &values [], int &indices [], int low, int high, bool ascending)
{
  double pivot = values [indices [high]];
  int i = low - 1;

  for (int j = low; j < high; j++)
  {
    bool condition = ascending ? (values [indices [j]] < pivot) : (values [indices [j]] > pivot);
    if (condition)
    {
      i++;
      // Обмен значениями
      int temp = indices [i];
      indices [i] = indices [j];
      indices [j] = temp;
    }
  }

  // Обмен значениями
  int temp = indices [i + 1];
  indices [i + 1] = indices [high];
  indices [high] = temp;

  return i + 1;
}
//——————————————————————————————————————————————————————————————————————————————