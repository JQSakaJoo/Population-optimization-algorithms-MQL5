//+————————————————————————————————————————————————————————————————————————————+
//|                                                                  C_AO_CROm |
//|                                            Copyright 2007-2025, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/17760

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
class C_AO_CROm : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_CROm () { }
  C_AO_CROm ()
  {
    ao_name = "CROm";
    ao_desc = "Coral Reef Optimization M";
    ao_link = "https://www.mql5.com/ru/articles/17760";

    popSize     = 50;      // размер популяции
    reefRows    = 20;      // высота рифа
    reefCols    = 20;      // ширина рифа
    rho0        = 0.2;     // начальная занятость рифа
    Fb          = 0.99;    // доля кораллов для внешнего размножения
    Fa          = 0.01;    // доля кораллов для бесполого размножения
    Fd          = 0.8;     // доля кораллов для удаления
    Pd          = 0.9;     // вероятность депредации
    attemptsNum = 20;      // количество попыток для оседания личинок


    ArrayResize (params, 9);

    params [0].name = "popSize";     params [0].val = popSize;
    params [1].name = "reefRows";    params [1].val = reefRows;
    params [2].name = "reefCols";    params [2].val = reefCols;
    params [3].name = "rho0";        params [3].val = rho0;
    params [4].name = "Fb";          params [4].val = Fb;
    params [5].name = "Fa";          params [5].val = Fa;
    params [6].name = "Fd";          params [6].val = Fd;
    params [7].name = "Pd";          params [7].val = Pd;
    params [8].name = "attemptsNum"; params [8].val = attemptsNum;
  }

  void SetParams ()
  {
    popSize     = (int)params [0].val;
    reefRows    = (int)params [1].val;
    reefCols    = (int)params [2].val;
    rho0        = params      [3].val;
    Fb          = params      [4].val;
    Fa          = params      [5].val;
    Fd          = params      [6].val;
    Pd          = params      [7].val;
    attemptsNum = (int)params [8].val;
  }

  bool Init (const double &rangeMinP  [],  // минимальный диапазон поиска
             const double &rangeMaxP  [],  // максимальный диапазон поиска
             const double &rangeStepP [],  // шаг поиска
             const int     epochsP = 0);   // количество эпох

  void Moving        ();
  void Revision      ();

  //----------------------------------------------------------------------------
  int                reefRows;      // высота рифа
  int                reefCols;      // ширина рифа
  double             rho0;          // начальная занятость рифа
  double             Fb;            // доля кораллов для внешнего размножения
  double             Fa;            // доля кораллов для бесполого размножения
  double             Fd;            // доля кораллов для удаления
  double             Pd;            // вероятность депредации
  int                attemptsNum;   // количество попыток для оседания личинок

  private: //-------------------------------------------------------------------
  int                totalReefSize; // общий размер рифа
  bool   occupied    [];   // флаги занятости клеток рифа
  int                reefIndices []; // индексы агентов в a[], соответствующие занятым клеткам

  // Вспомогательные методы
  void   InitReef            ();
  int    LarvaSettling       (S_AO_Agent &larva);
  void   BroadcastSpawning   (S_AO_Agent &larvae [], int &larvaCount);
  void   Brooding            (S_AO_Agent &larvae [], int &larvaCount);
  void   AsexualReproduction ();
  void   Depredation         ();
  int    GetReefCoordIndex   (int row, int col);
  void   SortAgentsByFitness (int &indices [], int &count);
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_CROm::Init (const double &rangeMinP  [],  // минимальный диапазон поиска
                      const double &rangeMaxP  [],  // максимальный диапазон поиска
                      const double &rangeStepP [],  // шаг поиска
                      const int     epochsP = 0)    // количество эпох
{
  // Стандартная инициализация родительского класса
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  // Рассчитываем общий размер рифа
  totalReefSize = reefRows * reefCols;

  // Количество начальных кораллов не должно превышать popSize
  int initialPopSize = (int)MathRound (rho0 * totalReefSize);
  if (initialPopSize > popSize) initialPopSize = popSize;

  // Инициализация массива занятости и индексов
  ArrayResize (occupied, totalReefSize);
  ArrayResize (reefIndices, totalReefSize);

  // Заполняем массивы начальными значениями
  for (int i = 0; i < totalReefSize; i++)
  {
    occupied    [i] = false;
    reefIndices [i] = -1;
  }

  // Инициализация рифа
  InitReef ();

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_CROm::InitReef ()
{
  // Количество начальных кораллов в рифе (исходя из rho0)
  int initialCorals = (int)MathRound (rho0 * totalReefSize);

  // Количество начальных кораллов не должно превышать размер популяции
  if (initialCorals > popSize) initialCorals = popSize;

  // Инициализируем initialCorals случайных позиций в рифе
  for (int i = 0; i < initialCorals; i++)
  {
    int pos;
    // Ищем свободную позицию
    do
    {
      pos = (int)MathFloor (u.RNDfromCI (0, totalReefSize));
      // Защита от выхода за пределы массива
      if (pos < 0) pos = 0;
      if (pos >= totalReefSize) pos = totalReefSize - 1;
    }
    while (occupied [pos]);

    // Создаем новый коралл на найденной позиции
    occupied [pos] = true;
    reefIndices [pos] = i;

    // Генерируем случайные координаты для нового коралла
    for (int c = 0; c < coords; c++)
    {
      double coordinate = u.RNDfromCI (rangeMin [c], rangeMax [c]);
      a [i].c [c] = u.SeInDiSp (coordinate, rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_CROm::Moving ()
{
  if (!revision)
  {
    // Первичная оценка всех кораллов в рифе
    for (int i = 0; i < totalReefSize; i++)
    {
      if (occupied [i])
      {
        int idx = reefIndices [i];
        if (idx >= 0 && idx < popSize)
        {
          // Расчет приспособленности не требует использования GetFitness()
          // так как он будет вычислен во внешнем коде (в FuncTests)
        }
      }
    }

    revision = true;
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_CROm::Revision ()
{
  // Обновление глобального наилучшего решения
  for (int i = 0; i < totalReefSize; i++)
  {
    if (occupied [i] && a [reefIndices [i]].f > fB)
    {
      fB = a [reefIndices [i]].f;
      ArrayCopy (cB, a [reefIndices [i]].c, 0, 0, WHOLE_ARRAY);
    }
  }

  // Формируем массив для хранения личинок
  S_AO_Agent larvae [];
  ArrayResize (larvae, totalReefSize * 2); // Выделяем с запасом
  int larvaCount = 0;

  // Этап 1: Broadcast Spawning (внешнее половое размножение)
  BroadcastSpawning (larvae, larvaCount);

  // Этап 2: Brooding (внутреннее половое размножение)
  Brooding (larvae, larvaCount);

  // Вычисляем функцию приспособленности для каждой личинки
  // (будет выполнено в внешнем коде в FuncTests)

  // Этап 3: Оседание личинок
  for (int i = 0; i < larvaCount; i++)
  {
    LarvaSettling (larvae [i]);
  }

  // Этап 4: Бесполое размножение
  AsexualReproduction ();

  // Этап 5: Депредация
  Depredation ();
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
int C_AO_CROm::LarvaSettling (S_AO_Agent &larva)
{
  // Пытаемся заселить личинку attemptsNum раз
  for (int attempt = 0; attempt < attemptsNum; attempt++)
  {
    // Выбираем случайную позицию в рифе
    int pos = (int)MathFloor (u.RNDfromCI (0, totalReefSize));

    // Проверяем, что позиция находится в пределах массива
    if (pos < 0 || pos >= totalReefSize) continue;

    // Если позиция свободна, заселяем личинку
    if (!occupied [pos])
    {
      // Ищем свободный индекс в массиве агентов
      int newIndex = -1;
      for (int i = 0; i < popSize; i++)
      {
        bool used = false;
        for (int j = 0; j < totalReefSize; j++)
        {
          if (reefIndices [j] == i)
          {
            used = true;
            break;
          }
        }

        if (!used)
        {
          newIndex = i;
          break;
        }
      }

      if (newIndex != -1)
      {
        // Копируем решение личинки
        ArrayCopy (a [newIndex].c, larva.c, 0, 0, WHOLE_ARRAY);
        a [newIndex].f = larva.f;

        // Обновляем информацию о рифе
        occupied [pos] = true;
        reefIndices [pos] = newIndex;

        return pos;
      }
    }
    // Если позиция занята, проверяем, лучше ли личинка текущего коралла
    else
      if (occupied [pos] && reefIndices [pos] >= 0 && reefIndices [pos] < popSize && larva.f > a [reefIndices [pos]].f)
      {
        // Личинка вытесняет существующий коралл
        ArrayCopy (a [reefIndices [pos]].c, larva.c, 0, 0, WHOLE_ARRAY);
        a [reefIndices [pos]].f = larva.f;

        return pos;
      }
  }

  // Если личинке не удалось заселиться, возвращаем -1
  return -1;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_CROm::BroadcastSpawning (S_AO_Agent &larvae [], int &larvaCount)
{
  // Находим все занятые позиции
  int occupiedIndices [];
  int occupiedCount = 0;

  for (int i = 0; i < totalReefSize; i++)
  {
    if (occupied [i])
    {
      ArrayResize (occupiedIndices, occupiedCount + 1);
      occupiedIndices [occupiedCount] = i;
      occupiedCount++;
    }
  }

  // Проверка на случай, если нет занятых позиций
  if (occupiedCount == 0) return;

  // Выбираем долю Fb для внешнего размножения
  int broadcastCount = (int)MathRound (Fb * occupiedCount);
  if (broadcastCount <= 0) broadcastCount = 1; // Минимум один коралл
  if (broadcastCount > occupiedCount) broadcastCount = occupiedCount;

  // Перемешиваем индексы
  for (int i = 0; i < occupiedCount; i++)
  {
    // Фиксируем проблему выхода за пределы массива
    int j = (int)MathFloor (u.RNDfromCI (0, occupiedCount));

    // Гарантируем, что j в пределах массива
    if (j >= 0 && j < occupiedCount && i < occupiedCount)
    {
      int temp = occupiedIndices [i];
      occupiedIndices [i] = occupiedIndices [j];
      occupiedIndices [j] = temp;
    }
  }

  // Образуем пары и создаем потомство
  for (int i = 0; i < broadcastCount - 1; i += 2)
  {
    if (i + 1 < broadcastCount) // Проверяем, что есть второй родитель
    {
      int idx1 = reefIndices [occupiedIndices [i]];
      int idx2 = reefIndices [occupiedIndices [i + 1]];

      if (idx1 >= 0 && idx1 < popSize && idx2 >= 0 && idx2 < popSize)
      {
        // Инициализируем личинку
        larvae [larvaCount].Init (coords);

        // Создаем новую личинку как результат скрещивания
        for (int c = 0; c < coords; c++)
        {
          // Простой метод скрещивания: среднее значение координат родителей с небольшой мутацией
          double value = (a [idx1].c [c] + a [idx2].c [c]) / 2.0 + u.RNDfromCI (-0.1, 0.1) * (rangeMax [c] - rangeMin [c]);
          larvae [larvaCount].c [c] = u.SeInDiSp (value, rangeMin [c], rangeMax [c], rangeStep [c]);
        }

        // Увеличиваем счетчик личинок
        larvaCount++;
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_CROm::Brooding (S_AO_Agent &larvae [], int &larvaCount)
{
  // Находим все занятые позиции
  int occupiedIndices [];
  int occupiedCount = 0;

  for (int i = 0; i < totalReefSize; i++)
  {
    if (occupied [i])
    {
      ArrayResize (occupiedIndices, occupiedCount + 1);
      occupiedIndices [occupiedCount] = i;
      occupiedCount++;
    }
  }

  // Проверка на случай, если нет занятых позиций
  if (occupiedCount == 0) return;

  // Количество кораллов для внутреннего размножения
  int broodingCount = (int)MathRound ((1.0 - Fb) * occupiedCount);
  if (broodingCount <= 0) broodingCount = 1; // Минимум один коралл
  if (broodingCount > occupiedCount) broodingCount = occupiedCount;

  // Перемешиваем индексы
  for (int i = 0; i < occupiedCount; i++)
  {
    // Фиксируем проблему выхода за пределы массива
    int j = (int)MathFloor (u.RNDfromCI (0, occupiedCount));

    // Гарантируем, что j в пределах массива
    if (j >= 0 && j < occupiedCount && i < occupiedCount)
    {
      int temp = occupiedIndices [i];
      occupiedIndices [i] = occupiedIndices [j];
      occupiedIndices [j] = temp;
    }
  }

  // Для каждого выбранного коралла создаем мутированную копию
  for (int i = 0; i < broodingCount; i++)
  {
    if (i < occupiedCount) // Проверка на выход за границы
    {
      int idx = reefIndices [occupiedIndices [i]];

      if (idx >= 0 && idx < popSize)
      {
        // Инициализируем личинку
        larvae [larvaCount].Init (coords);

        // Создаем новую личинку как мутацию исходного коралла
        for (int c = 0; c < coords; c++)
        {
          // Мутация: добавляем небольшое случайное отклонение
          double value = a [idx].c [c] + u.RNDfromCI (-0.2, 0.2) * (rangeMax [c] - rangeMin [c]);
          larvae [larvaCount].c [c] = u.SeInDiSp (value, rangeMin [c], rangeMax [c], rangeStep [c]);
        }

        // Увеличиваем счетчик личинок
        larvaCount++;
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_CROm::AsexualReproduction ()
{
  // Находим все занятые позиции и их индексы
  int occupiedIndices [];
  int occupiedCount = 0;

  for (int i = 0; i < totalReefSize; i++)
  {
    if (occupied [i])
    {
      ArrayResize (occupiedIndices, occupiedCount + 1);
      occupiedIndices [occupiedCount] = i;
      occupiedCount++;
    }
  }

  // Если нет занятых позиций, выходим
  if (occupiedCount == 0) return;

  // Сортируем индексы по приспособленности
  SortAgentsByFitness (occupiedIndices, occupiedCount);

  // Выбираем лучшие Fa% кораллов для бесполого размножения
  int budCount = (int)MathRound (Fa * occupiedCount);
  if (budCount <= 0) budCount = 1; // Минимум один коралл
  if (budCount > occupiedCount) budCount = occupiedCount;

  // Для каждого выбранного коралла создаем клон и пытаемся заселить его
  for (int i = 0; i < budCount; i++)
  {
    if (i < occupiedCount) // Проверка на выход за границы
    {
      int idx = reefIndices [occupiedIndices [i]];

      if (idx >= 0 && idx < popSize)
      {
        // Создаем новую личинку как точную копию исходного коралла
        S_AO_Agent clone;
        clone.Init (coords);
        ArrayCopy (clone.c, a [idx].c, 0, 0, WHOLE_ARRAY);
        clone.f = a [idx].f;

        // Пытаемся заселить клон
        LarvaSettling (clone);
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
/*oid C_AO_CRO::Depredation()
{
  // Применяем депредацию с вероятностью Pd
  if (u.RNDfromCI(0, 1) < Pd)
  {
    // Находим все занятые позиции и их индексы
    int occupiedIndices[];
    int occupiedCount = 0;

    for (int i = 0; i < totalReefSize; i++)
    {
      if (occupied[i])
      {
        ArrayResize(occupiedIndices, occupiedCount + 1);
        occupiedIndices[occupiedCount] = i;
        occupiedCount++;
      }
    }

    // Если нет занятых позиций, выходим
    if (occupiedCount == 0) return;

    // Сортируем индексы по приспособленности
    SortAgentsByFitness(occupiedIndices, occupiedCount);

    // Переворачиваем массив, чтобы худшие были первыми
    for (int i = 0; i < occupiedCount / 2; i++)
    {
      if(i < occupiedCount && (occupiedCount - 1 - i) < occupiedCount) // Проверка на выход за границы
      {
        int temp = occupiedIndices[i];
        occupiedIndices[i] = occupiedIndices[occupiedCount - 1 - i];
        occupiedIndices[occupiedCount - 1 - i] = temp;
      }
    }

    // Удаляем худшие Fd% кораллов
    int removeCount = (int)MathRound(Fd * occupiedCount);
    if (removeCount <= 0) removeCount = 1; // Минимум один коралл
    if (removeCount > occupiedCount) removeCount = occupiedCount; // Защита от переполнения

    for (int i = 0; i < removeCount; i++)
    {
      if(i < occupiedCount) // Проверка на выход за границы
      {
        int pos = occupiedIndices[i];
        if(pos >= 0 && pos < totalReefSize) // Дополнительная проверка
        {
          occupied[pos] = false;
          reefIndices[pos] = -1;
        }
      }
    }
  }
}*/

//——————————————————————————————————————————————————————————————————————————————
void C_AO_CROm::Depredation ()
{
  // Применяем депредацию с вероятностью Pd
  if (u.RNDfromCI (0, 1) < Pd)
  {
    // Находим все занятые позиции и их индексы
    int occupiedIndices [];
    int occupiedCount = 0;

    for (int i = 0; i < totalReefSize; i++)
    {
      if (occupied [i])
      {
        ArrayResize (occupiedIndices, occupiedCount + 1);
        occupiedIndices [occupiedCount] = i;
        occupiedCount++;
      }
    }

    // Если нет занятых позиций, выходим
    if (occupiedCount == 0) return;

    // Сортируем индексы по приспособленности (от лучших к худшим)
    SortAgentsByFitness (occupiedIndices, occupiedCount);

    // Определяем количество лучших решений, используемых для генерации новых
    int eliteCount = (int)MathRound (0.1 * occupiedCount); // Используем 10% лучших
    if (eliteCount <= 0) eliteCount = 1; // Минимум одно элитное решение
    if (eliteCount > occupiedCount) eliteCount = occupiedCount;

    // Определяем количество худших решений для замены
    int removeCount = (int)MathRound (Fd * occupiedCount);
    if (removeCount <= 0) removeCount = 1; // Минимум одно решение заменяем
    if (removeCount > occupiedCount - eliteCount) removeCount = occupiedCount - eliteCount; // Не удаляем элитные решения

    // Удаляем худшие решения и заменяем их новыми в окрестности лучших
    for (int i = 0; i < removeCount; i++)
    {
      // Индекс удаляемого решения (с конца отсортированного массива)
      int removeIndex = occupiedCount - 1 - i;
      if (removeIndex < 0 || removeIndex >= occupiedCount) continue;

      int posToRemove = occupiedIndices [removeIndex];
      if (posToRemove < 0 || posToRemove >= totalReefSize) continue;

      // Выбираем одно из элитных решений
      double power = 0.1; // Параметр степенного распределения
      double r = u.RNDfromCI (0, 1);
      int eliteIdx = (int)(eliteCount);
      if (eliteIdx >= eliteCount) eliteIdx = eliteCount - 1;

      int posBest = occupiedIndices [eliteIdx];
      if (posBest < 0 || posBest >= totalReefSize) continue;

      int bestAgentIdx = reefIndices [posBest];
      if (bestAgentIdx < 0 || bestAgentIdx >= popSize) continue;

      // Освобождаем позицию для нового решения
      occupied [posToRemove] = false;

      // Генерируем новое решение в окрестности выбранного элитного решения
      int newAgentIdx = reefIndices [posToRemove]; // Используем тот же индекс агента

      if (newAgentIdx >= 0 && newAgentIdx < popSize)
      {
        // Генерация через степенное распределение вокруг лучшего решения
        for (int c = 0; c < coords; c++)
        {
          // Определяем радиус поиска (можно адаптировать в зависимости от прогресса)
          double radius = 0.7 * (rangeMax [c] - rangeMin [c]); // 10% от диапазона

          // Генерируем значение по степенному закону
          double sign = u.RNDfromCI (0, 1) < 0.5 ? -1.0 : 1.0; // Случайный знак
          double deviation = sign * radius * MathPow (u.RNDfromCI (0, 1), 1.0 / power);

          // Новое значение в окрестности лучшего
          double newValue = a [bestAgentIdx].c [c] + deviation;

          // Ограничиваем значение в допустимом диапазоне
          a [newAgentIdx].c [c] = u.SeInDiSp (newValue, rangeMin [c], rangeMax [c], rangeStep [c]);
        }

        // Заселяем новое решение в риф
        occupied [posToRemove] = true;
        reefIndices [posToRemove] = newAgentIdx;
      }
    }
  }
}

//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
int C_AO_CROm::GetReefCoordIndex (int row, int col)
{
  // Проверка на выход за границы
  if (row < 0 || row >= reefRows || col < 0 || col >= reefCols) return -1;

  // Переводим двухмерную позицию в одномерный индекс
  return row * reefCols + col;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_CROm::SortAgentsByFitness (int &indices [], int &count)
{
  // Проверка на пустой массив
  if (count <= 0) return;

  // Сортировка пузырьком по убыванию приспособленности
  for (int i = 0; i < count - 1; i++)
  {
    for (int j = 0; j < count - i - 1; j++)
    {
      if (j + 1 < count) // Проверка, чтобы j+1 не выходил за пределы
      {
        int idx1 = reefIndices [indices [j]];
        int idx2 = reefIndices [indices [j + 1]];

        if (idx1 >= 0 && idx1 < popSize && idx2 >= 0 && idx2 < popSize) // Проверка индексов
        {
          if (a [idx1].f < a [idx2].f)
          {
            int temp = indices [j];
            indices [j] = indices [j + 1];
            indices [j + 1] = temp;
          }
        }
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————
