//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_SRA |
//|                                            Copyright 2007-2025, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/17380

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
class C_AO_SRA : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_SRA () { }
  C_AO_SRA ()
  {
    ao_name = "SRA";
    ao_desc = "Successful Restaurateur Algorithm (joo)";
    ao_link = "https://www.mql5.com/ru/articles/17380";

    popSize            = 50;   // число агентов (размер "меню")
    temperature        = 1.0;  // начальная "температура" для управления исследованием
    coolingRate        = 0.98; // скорость охлаждения
    menuInnovationRate = 0.3;  // интенсивность кулинарных экспериментов

    ArrayResize (params, 4);

    params [0].name = "popSize";            params [0].val = popSize;
    params [1].name = "temperature";        params [1].val = temperature;
    params [2].name = "coolingRate";        params [2].val = coolingRate;
    params [3].name = "menuInnovationRate"; params [3].val = menuInnovationRate;
  }

  void SetParams ()
  {
    popSize            = (int)params [0].val;
    temperature        = params      [1].val;
    coolingRate        = params      [2].val;
    menuInnovationRate = params      [3].val;
  }

  bool Init (const double &rangeMinP  [],  // минимальные значения
             const double &rangeMaxP  [],  // максимальные значения
             const double &rangeStepP [],  // шаг изменения
             const int     epochsP = 0);   // количество эпох

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  double temperature;        // текущая "температура"
  double coolingRate;        // скорость охлаждения
  double menuInnovationRate; // интенсивность кулинарных экспериментов

  private: //-------------------------------------------------------------------
  S_AO_Agent menu  [];
  S_AO_Agent menuT [];
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
//--- Инициализация
bool C_AO_SRA::Init (const double &rangeMinP  [],
                     const double &rangeMaxP  [],
                     const double &rangeStepP [],
                     const int epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  ArrayResize (menu,  popSize * 2);
  ArrayResize (menuT, popSize * 2);

  for (int p = 0; p < popSize * 2; p++) menu [p].Init (coords);

  temperature = 1.0; // сброс температуры при инициализации

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
//--- Основной шаг алгоритма
void C_AO_SRA::Moving ()
{
  //----------------------------------------------------------------------------
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

  //----------------------------------------------------------------------------
  // Снижаем температуру
  temperature *= coolingRate;

  // Основной цикл по агентам популяции
  for (int p = 0; p < popSize; p++)
  {
    // Берем худший элемент из первой половины отсортированной популяции (с индексом popSize-1)
    // Помним, что в menu элементы отсортированы от лучшего к худшему
    ArrayCopy (a [p].c, menu [popSize - 1].c, 0, 0, WHOLE_ARRAY);

    // Решаем, создавать ли гибрид или экспериментировать с новым "блюдом"
    // Вероятность эксперимента зависит от температуры - в начале больше экспериментов
    if (u.RNDprobab () < (1.0 - menuInnovationRate * temperature))
    {
      // Выбираем "рецепт-донор" с вероятностью пропорциональной успешности блюда
      double r = u.RNDprobab ();
      r = pow (r, 2);                                         // Усиление предпочтения к лучшим блюдам
      int menuIND = (int)u.Scale (r, 0, 1.0, 0, popSize - 1); // Лучшие в начале массива

      // Для каждой координаты
      for (int c = 0; c < coords; c++)
      {
        // С вероятностью, зависящей от температуры, берем параметр из успешного блюда
        if (u.RNDprobab () < 0.8)
        {
          a [p].c [c] = menu [menuIND].c [c];
        }

        // Мутация с адаптивной вероятностью - чем дальше от лучшего решения и выше температура, тем больше мутаций
        double mutationRate = 0.1 + 0.4 * temperature * (double)(p) / popSize;
        if (u.RNDprobab () < mutationRate)
        {
          // Комбинация различных типов мутаций
          if (u.RNDprobab () < 0.5) a [p].c [c] = u.GaussDistribution (a [p].c [c], rangeMin [c], rangeMax [c], 2);
          else                      a [p].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]); // Иногда полностью новое значение

          // Убедимся, что значение в допустимых пределах
          a [p].c [c] = u.SeInDiSp (a [p].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
        }
      }
    }
    else // Создаем абсолютно новое "блюдо"
    {
      for (int c = 0; c < coords; c++)
      {
        // Вариация 1: полностью случайное значение
        if (u.RNDprobab () < 0.7)
        {
          a [p].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        }
        // Вариация 2: основано на лучшем найденном решении с большим отклонением
        else
        {
          a [p].c [c] = u.GaussDistribution (cB [c], rangeMin [c], rangeMax [c], 1);
        }

        a [p].c [c] = u.SeInDiSp (a [p].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }

    // Иногда добавляем элементы из лучшего решения напрямую (элитизм)
    if (u.RNDprobab () < 0.1)
    {
      int numEliteCoords = u.RNDintInRange (1, coords / 3); // Берем от 1 до 30% координат
      for (int i = 0; i < numEliteCoords; i++)
      {
        int c = u.RNDminusOne (coords);
        a [p].c [c] = cB [c]; // Берем значение из лучшего решения
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
//--- Обновление лучшего решения с учетом жадного выбора и вероятности принятия худших решений
void C_AO_SRA::Revision ()
{
  int bestIND = -1;

  // Находим лучшего агента в текущей популяции
  for (int p = 0; p < popSize; p++)
  {
    if (a [p].f > fB)
    {
      fB = a [p].f;
      bestIND = p;
    }
  }

  // Если нашли лучшее решение, обновляем cB
  if (bestIND != -1) ArrayCopy (cB, a [bestIND].c, 0, 0, WHOLE_ARRAY);

  // Добавляем текущий набор блюд в общее "меню"
  for (int p = 0; p < popSize; p++)
  {
    menu [popSize + p] = a [p];
  }

  // Сортируем всё "меню" от лучших к худшим решениям
  // После сортировки, первая половина menu будет содержать лучшие решения,
  // которые будут использоваться на следующей итерации
  u.Sorting (menu, menuT, popSize * 2);

  // Не позволяем температуре упасть ниже определенного порога
  if (temperature < 0.1) temperature = 0.1;
}
//——————————————————————————————————————————————————————————————————————————————
