//+——————————————————————————————————————————————————————————————————+
//|                                            C_AO_BSA_Backtracking |
//|                                  Copyright 2007-2025, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/18568

#include "#C_AO.mqh"

//————————————————————————————————————————————————————————————————————
class C_AO_BSA_Backtracking : public C_AO
{
  public: //----------------------------------------------------------
  ~C_AO_BSA_Backtracking () { }
  C_AO_BSA_Backtracking ()
  {
    ao_name = "BSA";
    ao_desc = "Backtracking Search Algorithm";
    ao_link = "https://www.mql5.com/ru/articles/18568";

    popSize = 10;     // размер популяции
    mixrate = 1.0;    // параметр кроссовера

    ArrayResize (params, 2);

    params [0].name = "popSize"; params [0].val = popSize;
    params [1].name = "mixrate"; params [1].val = mixrate;
  }

  void SetParams ()
  {
    popSize = (int)params [0].val;
    mixrate = params      [1].val;

    // Проверка корректности параметров
    //if (popSize < 2) popSize = 2;
    if (mixrate < 0.0) mixrate = 0.0;
    if (mixrate > 1.0) mixrate = 1.0;
  }

  bool Init (const double &rangeMinP  [],  // минимальные значения
             const double &rangeMaxP  [],  // максимальные значения
             const double &rangeStepP [],  // шаг изменения
             const int     epochsP = 0);   // количество эпох

  void Moving   ();
  void Revision ();

  //------------------------------------------------------------------
  double mixrate;        // параметр кроссовера

  private: //---------------------------------------------------------
  S_AO_Agent oldP [];    // историческая популяция
  S_AO_Agent M    [];    // мутантная популяция (Mutant)
  S_AO_Agent T    [];    // пробная популяция (Trial)

  double F;              // фактор амплитуды для мутации
  bool   needSelection;  // флаг необходимости выполнения Selection-II
  double prevFitness []; // массив для хранения предыдущих fitness

  // Вспомогательные структуры для кроссовера
  struct S_Map
  {
      int val [];      // бинарная карта для кроссовера

      void Init (int size)
      {
        ArrayResize (val, size);
        ArrayInitialize (val, 0);
      }
  };

  S_Map map [];        // массив бинарных карт для каждого агента

  // Методы алгоритма
  void SelectionI        ();
  void Mutation          ();
  void Crossover         ();
  void BoundaryControl   (S_AO_Agent &agent);
  void ShufflePopulation (S_AO_Agent &pop []);
};
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Инициализация
bool C_AO_BSA_Backtracking::Init (const double &rangeMinP  [],
                                  const double &rangeMaxP  [],
                                  const double &rangeStepP [],
                                  const int epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //------------------------------------------------------------------
  // Инициализация дополнительных структур BSA
  ArrayResize (oldP, popSize);
  ArrayResize (M,    popSize);
  ArrayResize (T,    popSize);
  ArrayResize (map,  popSize);
  ArrayResize (prevFitness, popSize);

  needSelection = false;

  for (int i = 0; i < popSize; i++)
  {
    oldP [i].Init (coords);
    M    [i].Init (coords);
    T    [i].Init (coords);
    map  [i].Init (coords);
  }

  // Инициализация исторической популяции oldP
  for (int p = 0; p < popSize; p++)
  {
    for (int c = 0; c < coords; c++)
    {
      oldP [p].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
      oldP [p].c [c] = u.SeInDiSp (oldP [p].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }

  return true;
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Основной шаг алгоритма
void C_AO_BSA_Backtracking::Moving ()
{
  // Начальная инициализация популяции
  if (!revision)
  {
    for (int p = 0; p < popSize; p++)
    {
      for (int c = 0; c < coords; c++)
      {
        a [p].c  [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        a [p].c  [c] = u.SeInDiSp (a [p].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }

    revision      = true;
    needSelection = false;
    return;
  }

  // Если нужно выполнить жадный отбор после вычисления fitness
  if (needSelection)
  {
    // Selection-II: Жадный отбор
    for (int i = 0; i < popSize; i++)
    {
      // Если текущее решение (из T) хуже предыдущего, возвращаем предыдущее
      if (a [i].f < prevFitness [i])
      {
        ArrayCopy (a [i].c, a [i].cP, 0, 0, WHOLE_ARRAY);
        a [i].f = prevFitness [i];
      }
    }

    needSelection = false;
  }

  //--- Основные шаги BSA:

  // Сохраняем текущие fitness перед генерацией новой популяции
  for (int i = 0; i < popSize; i++)
  {
    prevFitness [i] = a [i].f;
    ArrayCopy (a [i].cP, a [i].c, 0, 0, WHOLE_ARRAY);
  }

  // 1. Selection-I
  SelectionI ();

  // 2. Mutation
  Mutation ();

  // 3. Crossover
  Crossover ();

  // 4. Копируем пробную популяцию T в основную популяцию a для вычисления fitness
  for (int i = 0; i < popSize; i++)
  {
    ArrayCopy (a [i].c, T [i].c, 0, 0, WHOLE_ARRAY);
  }

  // Устанавливаем флаг для выполнения Selection-II после вычисления fitness
  needSelection = true;
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Selection-I: выбор исторической популяции
void C_AO_BSA_Backtracking::SelectionI ()
{
  // С вероятностью 50% обновляем историческую популяцию
  if (u.RNDprobab () < 0.5) // эквивалент if (a < b) где a,b ~ U(0,1)
  {
    // Копируем текущую популяцию в историческую
    for (int i = 0; i < popSize; i++)
    {
      ArrayCopy (oldP [i].c, a [i].c, 0, 0, WHOLE_ARRAY);
      oldP [i].f = a [i].f;
    }
  }

  // Перемешиваем историческую популяцию
  ShufflePopulation (oldP);
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Перемешивание популяции
void C_AO_BSA_Backtracking::ShufflePopulation (S_AO_Agent &pop [])
{
  for (int i = popSize - 1; i > 0; i--)
  {
    int j = u.RNDminusOne (i + 1);

    // Обмен местами элементов i и j
    S_AO_Agent temp;
    temp.Init (coords);

    ArrayCopy (temp.c, pop [i].c, 0, 0, WHOLE_ARRAY);
    temp.f = pop [i].f;

    ArrayCopy (pop [i].c, pop [j].c, 0, 0, WHOLE_ARRAY);
    pop [i].f = pop [j].f;

    ArrayCopy (pop [j].c, temp.c, 0, 0, WHOLE_ARRAY);
    pop [j].f = temp.f;
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Mutation: генерация мутантной популяции
void C_AO_BSA_Backtracking::Mutation ()
{
  // Генерируем фактор амплитуды
  F = u.GaussDistribution (0.0, -3.0, 3.0, 2);

  // Применяем мутацию: M = P + F * (oldP - P)
  for (int i = 0; i < popSize; i++)
  {
    for (int j = 0; j < coords; j++)
    {
      M [i].c [j] = a [i].c [j] + F * (oldP [i].c [j] - a [i].c [j]);
    }
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Crossover: генерация пробной популяции
void C_AO_BSA_Backtracking::Crossover ()
{
  // Инициализируем пробную популяцию как копию мутантной
  for (int i = 0; i < popSize; i++)
  {
    ArrayCopy (T [i].c, M [i].c, 0, 0, WHOLE_ARRAY);
  }

  // Выбираем стратегию кроссовера
  if (u.RNDprobab () < 0.4)
  {
    //--- СТРАТЕГИЯ 1: Использование mixrate
    for (int i = 0; i < popSize; i++)
    {
      // Сброс карты
      ArrayInitialize (map [i].val, 0);

      // Определяем количество элементов для кроссовера
      int numElements = (int)MathCeil (mixrate * u.RNDprobab () * coords);

      // Генерируем уникальные индексы для кроссовера
      for (int n = 0; n < numElements; n++)
      {
        int idx;
        do
        {
          idx = u.RNDminusOne (coords);
        }
        while (map [i].val [idx] == 1); // пока не найдем неиспользованный индекс

        map [i].val [idx] = 1;
      }

      // Применяем кроссовер
      for (int j = 0; j < coords; j++)
      {
        if (map [i].val [j] == 1)
        {
          T [i].c [j] = a [i].c [j];
        }
      }
    }
  }
  else
  {
    //--- СТРАТЕГИЯ 2: Мутация только одного элемента
    for (int i = 0; i < popSize; i++)
    {
      // Выбираем один случайный элемент
      int randomIndex = u.RNDminusOne (coords);
      T [i].c [randomIndex] = a [i].c [randomIndex];
    }
  }

  // Контроль границ для всех агентов пробной популяции
  for (int i = 0; i < popSize; i++)
  {
    BoundaryControl (T [i]);
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Контроль границ
void C_AO_BSA_Backtracking::BoundaryControl (S_AO_Agent &agent)
{
  for (int j = 0; j < coords; j++)
  {
    if (agent.c [j] < rangeMin [j] || agent.c [j] > rangeMax [j])
    {
      // Выбор стратегии обработки границ
      if (u.RNDprobab () < 0.5)
      {
        // Случайная регенерация
        agent.c [j] = u.RNDfromCI (rangeMin [j], rangeMax [j]);
      }
      else
      {
        // Установка на границу
        if (agent.c [j] < rangeMin [j]) agent.c [j] = rangeMin [j];
        else agent.c [j] = rangeMax [j];
      }
    }

    // Дискретизация
    agent.c [j] = u.SeInDiSp (agent.c [j], rangeMin [j], rangeMax [j], rangeStep [j]);
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Selection-II и обновление лучшего решения
void C_AO_BSA_Backtracking::Revision ()
{
  int bestIND = -1;

  for (int i = 0; i < popSize; i++)
  {
    // Обновляем глобальное лучшее решение
    if (a [i].f > fB)
    {
      fB = a [i].f;
      bestIND = i;
    }
  }

  // Копируем координаты лучшего решения
  if (bestIND != -1)
  {
    ArrayCopy (cB, a [bestIND].c, 0, 0, WHOLE_ARRAY);
  }
}
//————————————————————————————————————————————————————————————————————
