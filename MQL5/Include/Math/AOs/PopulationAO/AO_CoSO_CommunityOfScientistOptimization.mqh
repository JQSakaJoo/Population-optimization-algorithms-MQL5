//+——————————————————————————————————————————————————————————————————+
//|                                                        C_AO_CoSO |
//|                                  Copyright 2007-2025, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/18886

#include "#C_AO.mqh"

//————————————————————————————————————————————————————————————————————
// Структура для журнала
struct S_Journal_Entry
{
    double fitness;
    double decision [];

    void Init (int coords)
    {
      ArrayResize (decision, coords);
      fitness = -DBL_MAX;
    }
};

// Структура журнала
struct S_Journal
{
    S_Journal_Entry entries [];
    int length;
    int maxLength;

    void Init (int maxLen, int coords)
    {
      maxLength = maxLen;
      length = 0;

      ArrayResize (entries, maxLen);

      for (int i = 0; i < maxLen; i++)
      {
        entries [i].Init (coords);
        entries [i].fitness = -DBL_MAX; // Инициализируем минимальным значением
      }
    }

    void Add (double fit, const double &coord [])
    {
      // Быстрая проверка - если хуже всех и журнал полон, не добавляем
      if (length >= maxLength && fit <= entries [length - 1].fitness) return;

      int insertPos = length;

      // Находим позицию для вставки (бинарный поиск)
      if (length > 0)
      {
        int left = 0;
        int right = length - 1;

        while (left <= right)
        {
          int mid = (left + right) / 2;
          if (entries [mid].fitness < fit) right = mid - 1;
          else left = mid + 1;
        }
        insertPos = left;
      }

      // Если вставляем в конец и журнал полон
      if (insertPos >= maxLength) return;

      // Сдвигаем элементы
      if (length < maxLength) length++;

      for (int i = length - 1; i > insertPos; i--)
      {
        entries [i].fitness = entries [i - 1].fitness;
        ArrayCopy (entries [i].decision, entries [i - 1].decision, 0, 0, WHOLE_ARRAY);
      }

      // Вставляем новый элемент
      entries [insertPos].fitness = fit;
      ArrayCopy (entries [insertPos].decision, coord, 0, 0, WHOLE_ARRAY);
    }
};

// Расширенная структура исследователя
struct S_Researcher
{
    double x   [];   // текущая позиция
    double v   [];   // направление движения
    double b   [];   // личный лучший результат
    double rho [];   // вероятности публикации в журналах
    double s;        // стратегия управления средствами
    int    m;        // количество средств
    double f;        // fitness текущей позиции
    double fb;       // fitness лучшей позиции
    bool   alive;    // флаг активности

    void Init (int coords, int journalsNum)
    {
      if (ArraySize (x) != coords)
      {
        ArrayResize (x, coords);
        ArrayResize (v, coords);
        ArrayResize (b, coords);
      }
      if (ArraySize (rho) != journalsNum)
      {
        ArrayResize (rho, 0);
        ArrayResize (rho, journalsNum);
      }

      f  = -DBL_MAX;
      fb = -DBL_MAX;
      m  = 0;
      s  = 0.5;
      alive = true;
    }
};

//————————————————————————————————————————————————————————————————————
class C_AO_CoSO : public C_AO
{
  public: //----------------------------------------------------------
  ~C_AO_CoSO () { }
  C_AO_CoSO ()
  {
    ao_name = "CoSO";
    ao_desc = "Community of Scientist Optimization";
    ao_link = "https://www.mql5.com/ru/articles/18886";

    popSize      = 10;     // начальный размер популяции
    totalFunds   = 150;    // общее количество средств
    journalsNum  = 3;      // количество журналов
    journalLen   = 10;     // длина журнала
    omega        = 0.7;    // параметр инерции

    ArrayResize (params, 5);

    params [0].name = "popSize";     params [0].val = popSize;
    params [1].name = "totalFunds";  params [1].val = totalFunds;
    params [2].name = "journalsNum"; params [2].val = journalsNum;
    params [3].name = "journalLen";  params [3].val = journalLen;
    params [4].name = "omega";       params [4].val = omega;

    //----------------------------------------------------------------
    phi1         = 1.5;    // когнитивный параметр
    phi2         = 1.5;    // социальный параметр
    omegaMin     = 0.2;    // минимальный процент аутсайдеров
    omegaMax     = 0.5;    // максимальный процент аутсайдеров
    epsilonPlus  = 0.2;    // шаг увеличения разнообразия
    epsilonMinus = 0.1;    // шаг уменьшения разнообразия
  }

  void SetParams ()
  {
    popSize     = (int)params [0].val;
    totalFunds  = (int)params [1].val;
    journalsNum = (int)params [2].val;
    journalLen  = (int)params [3].val;
    omega       = params      [4].val;
  }

  bool Init (const double &rangeMinP  [],
             const double &rangeMaxP  [],
             const double &rangeStepP [],
             const int     epochsP = 0);

  void Moving   ();
  void Revision ();

  //------------------------------------------------------------------
  int    totalFunds;   // общее количество средств
  int    journalsNum;  // количество журналов
  int    journalLen;   // длина журнала
  double omega;        // параметр инерции
  double phi1;         // когнитивный параметр
  double phi2;         // социальный параметр
  double omegaMin;     // минимальный процент аутсайдеров
  double omegaMax;     // максимальный процент аутсайдеров
  double epsilonPlus;  // шаг увеличения разнообразия
  double epsilonMinus; // шаг уменьшения разнообразия

  private: //---------------------------------------------------------
  S_Researcher researchers [];  // массив исследователей
  S_Journal    journals    [];  // массив журналов
  double       omegaCurrent;    // текущий процент аутсайдеров
  double       sigma0;          // начальное стандартное отклонение
  int          actualPopSize;   // текущий размер популяции
  int          maxPopSize;      // максимальный размер популяции
  double       socialComponent []; // кэш для социального компонента

  struct S_GlobalReport
  {
      double fitness;
      int    index;
  };

  S_GlobalReport globalReport [];

  // Методы алгоритма
  void   UpdateDirection   (int idx);
  void   SubmitToJournal   (int idx);
  void   AssignFunds       (int availableFunds);
  void   HireResearchers   (int idx);
  void   CreateOutsiders   (int outsiderFunds);
  double ComputeStdDev     ();
  void   UpdateOmega       ();
  int    SelectJournal     (const double &probs []);
  void   NormalizeProbabilities (double &probs []);
  void   CompactPopulation ();
};
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
bool C_AO_CoSO::Init (const double &rangeMinP  [],
                      const double &rangeMaxP  [],
                      const double &rangeStepP [],
                      const int     epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //------------------------------------------------------------------
  // Инициализация переменных состояния
  actualPopSize = 0;
  maxPopSize    = 0;
  sigma0        = 0;
  omegaCurrent  = omegaMin + (omegaMax - omegaMin) / 2.0;

  // Проверка корректности параметров
  if (totalFunds < popSize) totalFunds = popSize;
  if (journalsNum < 1) journalsNum = 1;
  if (journalLen < 1) journalLen = 1;
  if (omega < 0.0) omega = 0.0;
  if (omega > 1.0) omega = 1.0;

  // Полная очистка от предыдущих запусков
  ArrayResize (researchers,     0);
  ArrayResize (journals,        0);
  ArrayResize (globalReport,    0);
  ArrayResize (socialComponent, 0);

  // Сброс параметров
  actualPopSize = 0;
  maxPopSize    = 0;
  sigma0        = 0;

  // Инициализация журналов
  ArrayResize (journals, journalsNum);
  for (int i = 0; i < journalsNum; i++)
  {
    journals [i].Init (journalLen, coords);
  }

  // Инициализация массива исследователей с запасом
  maxPopSize = MathMin (popSize * 3, 300); // Ограничиваем максимальный размер
  ArrayResize (researchers, maxPopSize);
  ArrayResize (socialComponent, coords);

  actualPopSize = popSize;
  int fundsPerResearcher = totalFunds / popSize;

  for (int i = 0; i < maxPopSize; i++)
  {
    researchers [i].Init (coords, journalsNum);
    researchers [i].alive = (i < popSize);

    if (i < popSize)
    {
      researchers [i].m = fundsPerResearcher;
      researchers [i].s = u.RNDprobab ();

      // Инициализация позиции
      for (int c = 0; c < coords; c++)
      {
        researchers [i].x [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        researchers [i].x [c] = u.SeInDiSp  (researchers [i].x [c], rangeMin [c], rangeMax [c], rangeStep [c]);
        researchers [i].b [c] = researchers [i].x [c];
        researchers [i].v [c] = u.GaussDistribution (0.0, -0.01, 0.01, 1);
      }

      // Инициализация вероятностей журналов
      for (int j = 0; j < journalsNum; j++)
      {
        researchers [i].rho [j] = u.RNDprobab ();
      }

      NormalizeProbabilities (researchers [i].rho);
    }
  }

  // Вычисляем начальное стандартное отклонение
  sigma0 = ComputeStdDev ();

  if (sigma0 == 0) sigma0 = 1.0; // Защита от деления на ноль
  omegaCurrent = omegaMin + (omegaMax - omegaMin) / 2.0;

  // Копируем исследователей в стандартный массив агентов
  for (int i = 0; i < popSize; i++)
  {
    ArrayCopy (a [i].c, researchers [i].x, 0, 0, WHOLE_ARRAY);
  }

  return true;
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
void C_AO_CoSO::Moving ()
{
  if (!revision)
  {
    revision = true;
    return;
  }

  //--- Основные шаги CoSO:

  // 1. Обновление fitness исследователей из массива агентов
  int aSize = ArraySize (a);
  for (int i = 0, j = 0; i < actualPopSize && j < aSize; i++)
  {
    if (!researchers [i].alive) continue;

    researchers [i].f = a [j].f;

    // Обновление личного лучшего
    if (researchers [i].f > researchers [i].fb)
    {
      researchers [i].fb = researchers [i].f;
      ArrayCopy (researchers [i].b, researchers [i].x, 0, 0, WHOLE_ARRAY);
    }
    j++;
  }

  // 2. Подача результатов в журналы
  for (int i = 0; i < actualPopSize; i++)
  {
    if (researchers [i].alive) SubmitToJournal (i);
  }

  // 3. Сбор глобального отчета и подсчет доступных средств
  int availableFunds = 0;
  int reportSize = 0;

  // Предварительный подсчет размера отчета
  for (int i = 0; i < actualPopSize; i++)
  {
    if (!researchers [i].alive) continue;

    researchers [i].m--;  // Тратим 1 единицу средств за итерацию
    availableFunds++;

    if (researchers [i].m > 0) reportSize++;
    else researchers [i].alive = false;
  }

  // Заполнение глобального отчета
  ArrayResize (globalReport, reportSize);
  int idx = 0;

  for (int i = 0; i < actualPopSize && idx < reportSize; i++)
  {
    if (researchers [i].alive && researchers [i].m > 0)
    {
      globalReport [idx].fitness = researchers [i].f;
      globalReport [idx].index = i;
      idx++;
    }
  }

  // 4. Быстрая сортировка глобального отчета
  for (int i = 0; i < reportSize - 1; i++)
  {
    for (int j = i + 1; j < reportSize; j++)
    {
      if (globalReport [i].fitness < globalReport [j].fitness)
      {
        S_GlobalReport temp = globalReport [i];
        globalReport [i] = globalReport [j];
        globalReport [j] = temp;
      }
    }
  }

  // 5. Распределение средств
  AssignFunds (availableFunds);

  // 6. Найм новых исследователей существующими
  for (int i = 0; i < actualPopSize; i++)
  {
    if (researchers [i].alive && researchers [i].m > 1) HireResearchers (i);
  }

  // 7. Обновление направления и позиции для каждого исследователя
  for (int i = 0; i < actualPopSize; i++)
  {
    if (!researchers [i].alive) continue;

    UpdateDirection (i);

    // Обновление позиции
    for (int c = 0; c < coords; c++)
    {
      researchers [i].x [c] += researchers [i].v [c];

      // Контроль границ
      if (researchers [i].x [c] < rangeMin [c]) researchers [i].x [c] = rangeMin [c];
      if (researchers [i].x [c] > rangeMax [c]) researchers [i].x [c] = rangeMax [c];

      researchers [i].x [c] = u.SeInDiSp (researchers [i].x [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }

  // 8. Обновление параметра разнообразия
  UpdateOmega ();

  // 9. Компактификация популяции
  CompactPopulation ();

  // 10. Копируем позиции в массив агентов для вычисления fitness
  ArrayResize (a, actualPopSize);
  idx = 0;
  for (int i = 0; i < maxPopSize && idx < actualPopSize; i++)
  {
    if (researchers [i].alive)
    {
      a [idx].Init (coords);
      ArrayCopy (a [idx].c, researchers [i].x, 0, 0, WHOLE_ARRAY);
      idx++;
    }
  }

  popSize = actualPopSize;  // Обновляем размер популяции
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
void C_AO_CoSO::UpdateDirection (int idx)
{
  double beta1 = u.RNDprobab ();
  double beta2 = u.RNDprobab ();

  // Социальный компонент
  ArrayInitialize (socialComponent, 0);

  for (int j = 0; j < journalsNum; j++)
  {
    if (journals [j].length > 0)
    {
      int entryIdx = u.RNDminusOne (journals [j].length);

      for (int c = 0; c < coords; c++)
      {
        socialComponent [c] += researchers [idx].rho [j] *
                               (journals [j].entries [entryIdx].decision [c] - researchers [idx].x [c]);
      }
    }
  }

  // Обновление направления
  for (int c = 0; c < coords; c++)
  {
    researchers [idx].v [c] = omega * researchers [idx].v [c] +
                              phi1 * beta1 * (researchers [idx].b [c] - researchers [idx].x [c]) +
                              phi2 * beta2 * socialComponent [c];
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
void C_AO_CoSO::SubmitToJournal (int idx)
{
  int journalIdx = SelectJournal (researchers [idx].rho);
  journals [journalIdx].Add (researchers [idx].f, researchers [idx].x);
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
int C_AO_CoSO::SelectJournal (const double &probs [])
{
  double rnd = u.RNDprobab ();
  double cumSum = 0;

  int probsSize = ArraySize (probs);
  for (int i = 0; i < probsSize; i++)
  {
    cumSum += probs [i];
    if (rnd <= cumSum) return i;
  }

  return probsSize - 1;
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
void C_AO_CoSO::AssignFunds (int availableFunds)
{
  // Средства для аутсайдеров
  int outsiderFunds = (int)(availableFunds * omegaCurrent);
  int existingFunds = availableFunds - outsiderFunds;

  int reportSize = ArraySize (globalReport);

  // Распределение средств существующим исследователям
  if (reportSize > 0)
  {
    int totalRank = reportSize * (reportSize + 1) / 2;

    for (int f = 0; f < existingFunds; f++)
    {
      double rnd = u.RNDprobab () * totalRank;
      double cumSum = 0;

      for (int i = 0; i < reportSize; i++)
      {
        cumSum += reportSize - i;
        if (rnd <= cumSum)
        {
          researchers [globalReport [i].index].m++;
          break;
        }
      }
    }
  }

  // Создание аутсайдеров
  CreateOutsiders (outsiderFunds);
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
void C_AO_CoSO::CreateOutsiders (int outsiderFunds)
{
  if (outsiderFunds <= 0) return;

  // Ограничиваем количество новых аутсайдеров, особенно если популяция уже большая
  int maxNew = (actualPopSize > 100) ? 2 : 5;
  int newResearchers = (int)(u.RNDfromCI (1, MathMin (outsiderFunds, maxNew)));
  int fundsPerNew = outsiderFunds / newResearchers;

  for (int i = 0; i < newResearchers; i++)
  {
    // Находим свободное место
    int idx = -1;
    for (int j = 0; j < actualPopSize; j++)
    {
      if (!researchers [j].alive)
      {
        idx = j;
        break;
      }
    }

    if (idx == -1 && actualPopSize < maxPopSize)
    {
      idx = actualPopSize;
    }

    if (idx == -1) // Нет места, расширяем массив
    {
      if (actualPopSize >= maxPopSize)
      {
        // Ограничиваем рост популяции
        int newMaxSize = MathMin (maxPopSize + 50, 500);
        if (newMaxSize == maxPopSize) continue; // Достигнут лимит, пропускаем создание

        maxPopSize = newMaxSize;
        ArrayResize (researchers, maxPopSize);
        for (int j = actualPopSize; j < maxPopSize; j++)
        {
          researchers [j].Init (coords, journalsNum);
          researchers [j].alive = false;
        }
        idx = actualPopSize;
      }
    }

    if (idx == -1) continue; // Не удалось создать

    researchers [idx].alive = true;
    researchers [idx].m = fundsPerNew;
    researchers [idx].s = u.RNDprobab ();

    // Случайная инициализация
    for (int c = 0; c < coords; c++)
    {
      researchers [idx].x [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
      researchers [idx].x [c] = u.SeInDiSp (researchers [idx].x [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      researchers [idx].b [c] = researchers [idx].x [c];
      researchers [idx].v [c] = u.GaussDistribution (0.0, -0.01, 0.01, 1);
    }

    // Инициализация вероятностей журналов
    for (int j = 0; j < journalsNum; j++)
    {
      researchers [idx].rho [j] = u.RNDprobab ();
    }
    NormalizeProbabilities (researchers [idx].rho);

    if (idx >= actualPopSize) actualPopSize = idx + 1;
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
void C_AO_CoSO::HireResearchers (int idx)
{
  if (researchers [idx].m <= 1) return;

  int keepFunds = (int)(researchers [idx].m * researchers [idx].s);
  int hireFunds = researchers [idx].m - keepFunds;
  researchers [idx].m = keepFunds;

  if (hireFunds <= 0) return;

  // Ограничиваем количество нанимаемых, особенно при большой популяции
  int maxNew = (actualPopSize > 100) ? 1 : 3;
  int newCount = (int)(u.RNDfromCI (1, MathMin (hireFunds, maxNew)));
  int fundsPerNew = hireFunds / newCount;

  for (int i = 0; i < newCount; i++)
  {
    // Находим свободное место
    int newIdx = -1;
    for (int j = 0; j < actualPopSize; j++)
    {
      if (!researchers [j].alive)
      {
        newIdx = j;
        break;
      }
    }

    if (newIdx == -1 && actualPopSize < maxPopSize)
    {
      newIdx = actualPopSize;
    }

    if (newIdx == -1) // Нет места
    {
      if (actualPopSize >= maxPopSize || actualPopSize >= 500) continue; // Пропускаем создание при достижении лимита
    }

    if (newIdx == -1) continue; // Не удалось найти место

    researchers [newIdx].alive = true;
    researchers [newIdx].m = fundsPerNew;
    researchers [newIdx].s = u.GaussDistribution (researchers [idx].s, 0, 1, 1);
    if (researchers [newIdx].s < 0) researchers [newIdx].s = 0;
    if (researchers [newIdx].s > 1) researchers [newIdx].s = 1;

    // Наследование от супервизора
    ArrayCopy (researchers [newIdx].b, researchers [idx].b, 0, 0, WHOLE_ARRAY);

    // Позиция около супервизора
    for (int c = 0; c < coords; c++)
    {
      researchers [newIdx].x [c] = researchers [idx].x [c] + u.GaussDistribution (0.0, -0.01, 0.01, 1);

      // Контроль границ
      if (researchers [newIdx].x [c] < rangeMin [c]) researchers [newIdx].x [c] = rangeMin [c];
      if (researchers [newIdx].x [c] > rangeMax [c]) researchers [newIdx].x [c] = rangeMax [c];

      researchers [newIdx].x [c] = u.SeInDiSp (researchers [newIdx].x [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      researchers [newIdx].v [c] = u.GaussDistribution (0.0, -0.01, 0.01, 1);
    }

    // Возмущенные вероятности журналов
    for (int j = 0; j < journalsNum; j++)
    {
      researchers [newIdx].rho [j] = u.GaussDistribution (researchers [idx].rho [j], 0, 1, 1);
      if (researchers [newIdx].rho [j] < 0) researchers [newIdx].rho [j] = 0;
    }
    NormalizeProbabilities (researchers [newIdx].rho);

    if (newIdx >= actualPopSize) actualPopSize = newIdx + 1;
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
double C_AO_CoSO::ComputeStdDev ()
{
  if (actualPopSize == 0) return 0;

  double mean = 0;
  int count = 0;

  for (int i = 0; i < actualPopSize; i++)
  {
    if (researchers [i].alive)
    {
      mean += researchers [i].f;
      count++;
    }
  }

  if (count == 0) return 0;
  mean /= count;

  double variance = 0;
  for (int i = 0; i < actualPopSize; i++)
  {
    if (researchers [i].alive)
    {
      variance += MathPow (researchers [i].f - mean, 2);
    }
  }
  variance /= count;

  return MathSqrt (variance);
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
void C_AO_CoSO::UpdateOmega ()
{
  double currentSigma = ComputeStdDev ();

  if (currentSigma < sigma0)
  {
    // Увеличиваем долю аутсайдеров при сходимости
    omegaCurrent += (omegaMax - omegaMin) / 2.0 * epsilonPlus;
  }
  else
  {
    // Уменьшаем долю аутсайдеров
    omegaCurrent -= (omegaMax - omegaMin) / 2.0 * epsilonMinus;
  }

  // Ограничиваем диапазон
  if (omegaCurrent < omegaMin) omegaCurrent = omegaMin;
  if (omegaCurrent > omegaMax) omegaCurrent = omegaMax;
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
void C_AO_CoSO::NormalizeProbabilities (double &probs [])
{
  double sum = 0;
  int size = ArraySize (probs);

  for (int i = 0; i < size; i++)
  {
    sum += probs [i];
  }

  if (sum > 0)
  {
    for (int i = 0; i < size; i++)
    {
      probs [i] /= sum;
    }
  }
  else
  {
    // Равномерное распределение
    double val = 1.0 / size;
    for (int i = 0; i < size; i++)
    {
      probs [i] = val;
    }
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
void C_AO_CoSO::CompactPopulation ()
{
  // Подсчитываем живых исследователей
  int aliveCount = 0;
  for (int i = 0; i < actualPopSize; i++)
  {
    if (researchers [i].alive) aliveCount++;
  }

  // Если слишком много мертвых, компактифицируем
  if (aliveCount < actualPopSize * 0.75 || actualPopSize > 200)
  {
    int newIdx = 0;
    for (int i = 0; i < actualPopSize; i++)
    {
      if (researchers [i].alive)
      {
        if (i != newIdx)
        {
          // Копируем живого исследователя на новое место
          researchers [newIdx] = researchers [i];
          researchers [i].alive = false;
        }
        newIdx++;
      }
    }
    actualPopSize = aliveCount;

    // Если популяция все еще слишком большая, ограничиваем
    if (actualPopSize > 150)
    {
      // Сортируем по fitness и оставляем лучших
      for (int i = 0; i < actualPopSize - 1; i++)
      {
        for (int j = i + 1; j < actualPopSize; j++)
        {
          if (researchers [i].f < researchers [j].f)
          {
            S_Researcher temp = researchers [i];
            researchers [i] = researchers [j];
            researchers [j] = temp;
          }
        }
      }

      // Убиваем худших
      for (int i = 150; i < actualPopSize; i++)
      {
        researchers [i].alive = false;
        researchers [i].m = 0;
      }
      actualPopSize = 150;
    }
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
void C_AO_CoSO::Revision ()
{
  int bestIND = -1;
  int aSize = ArraySize (a);

  for (int i = 0; i < aSize; i++)
  {
    if (a [i].f > fB)
    {
      fB = a [i].f;
      bestIND = i;
    }
  }

  if (bestIND != -1)
  {
    ArrayCopy (cB, a [bestIND].c, 0, 0, WHOLE_ARRAY);
  }
}
//————————————————————————————————————————————————————————————————————