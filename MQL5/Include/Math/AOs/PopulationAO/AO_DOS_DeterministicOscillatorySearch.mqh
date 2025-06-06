//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_DOS |
//|                                            Copyright 2007-2025, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/18154

#include "#C_AO.mqh"

// Структура для хранения скорости частицы
struct S_DOS_Velocity
{
    int    slope; // Наклон частицы (-1: отрицательный, 0: неизвестный, 1: положительный)
    double v [];  // Компоненты скорости по каждому измерению


    void Init (int dims)
    {
      slope = 0;
      ArrayResize     (v, dims);
      ArrayInitialize (v, 0.0); // Быстрая инициализация нулями всего массива
    }

    // Проверка на нулевую скорость
    bool IsZero (double epsilon = 1e-10)
    {
      for (int i = 0; i < ArraySize (v); i++) if (MathAbs (v [i]) > epsilon) return false;
      return true;
    }
};

//——————————————————————————————————————————————————————————————————————————————
class C_AO_DOS : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_DOS () { }
  C_AO_DOS ()
  {
    ao_name = "DOS";
    ao_desc = "Deterministic Oscillatory Search";
    ao_link = "https://www.mql5.com/ru/articles/18154";

    // Установка параметров по умолчанию
    popSize        = 30;   // размер популяции
    movementFactor = 0.95; // фактор движения к лучшему решению

    // Создание и инициализация массива параметров
    ArrayResize (params, 2);
    params [0].name = "Population Size"; params [0].val = popSize;
    params [1].name = "Movement Factor"; params [1].val = movementFactor;
  }

  void SetParams ()
  {
    // Установка значений параметров с проверкой
    popSize        = (int)MathMax (5, params [0].val);             // Минимально 5 частиц для эффективности
    movementFactor = MathMax (0.1, MathMin (1.0, params [1].val)); // Ограничение от 0.1 до 1.0
  }

  bool Init (const double &rangeMinP  [],  // минимальные значения
             const double &rangeMaxP  [],  // максимальные значения
             const double &rangeStepP [],  // шаг изменения
             const int     epochsP = 0);

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  double movementFactor;          // фактор движения к лучшему решению

  S_DOS_Velocity velocities [];  // Массив структур скоростей частиц

  private: //-------------------------------------------------------------------
  void InitializeParticles     ();
  void ProcessParticleMovement (int particleIndex);
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_DOS::Init (const double &rangeMinP  [],  // минимальные значения
                     const double &rangeMaxP  [],  // максимальные значения
                     const double &rangeStepP [],  // шаг изменения
                     const int     epochsP = 0)    // количество эпох
{
  // Стандартная инициализация C_AO
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  // Выделение памяти под массивы
  ArrayResize (velocities, popSize);

  // Инициализация скоростей для каждого измерения
  for (int i = 0; i < popSize; i++) velocities [i].Init (coords);

  // Инициализация позиций частиц детерминистическим способом
  InitializeParticles ();

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Инициализация частиц с использованием комбинаторного детерминистического   |
//| подхода                                                                    |
//+----------------------------------------------------------------------------+
void C_AO_DOS::InitializeParticles ()
{
  // Равномерное распределение начальных позиций
  for (int i = 0; i < popSize; i++)
  {
    for (int d = 0; d < coords; d++)
    {
      a [i].c [d] = u.RNDfromCI (rangeMin [d], rangeMax [d]);
      a [i].c [d] = u.SeInDiSp (a [i].c [d], rangeMin [d], rangeMax [d], rangeStep [d]);
    }

    // Инициализация состояния частицы
    velocities [i].slope = 0; // Наклон неизвестен
  }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Основной метод движения частиц                                             |
//+----------------------------------------------------------------------------+
void C_AO_DOS::Moving ()
{
  // Обработка всех частиц
  for (int i = 0; i < popSize; i++)
  {
    // Сохранение значения фитнеса
    a [i].fP = a [i].f;

    // Вычисление новых координат на основе скорости
    for (int d = 0; d < coords; d++)
    {
      // Обновление позиции
      a [i].c [d] += velocities [i].v [d];

      // Округление до ближайшего допустимого шага
      a [i].c [d] = u.SeInDiSp (a [i].c [d], rangeMin [d], rangeMax [d], rangeStep [d]);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Метод обновления фитнес-функций                                            |
//+----------------------------------------------------------------------------+
void C_AO_DOS::Revision ()
{
  // Обработка каждой частицы
  for (int i = 0; i < popSize; i++)
  {
    // Обновление лучшего решения, если текущее решение лучше
    if (a [i].f > fB)
    {
      fB = a [i].f;
      ArrayCopy (cB, a [i].c, 0, 0, WHOLE_ARRAY);
    }

    // Обработка движения частицы на основе изменения фитнеса
    ProcessParticleMovement (i);
  }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Обработка движения частицы после обновления фитнеса                        |
//+----------------------------------------------------------------------------+
void C_AO_DOS::ProcessParticleMovement (int particleIndex)
{
  // Локальные переменные для оптимизации доступа
  double currentFitness  = a [particleIndex].f;
  double previousFitness = a [particleIndex].fP;
  int    currentSlope    = velocities [particleIndex].slope;

  // Сравнение фитнесов для определения направления движения
  double fitnessDiff = currentFitness - previousFitness;

  // Обработка наклона в соответствии с текущим состоянием
  if (currentSlope == 0) // Неизвестный наклон
  {
    // Определяем наклон на основе изменения фитнеса
    velocities [particleIndex].slope = (fitnessDiff > 0) ? 1 : (fitnessDiff < 0) ? -1 : 0;
  }
  else
    if (currentSlope == 1 && fitnessDiff < 0) // Положительный наклон и ухудшение фитнеса
    {
      // Меняем направление и уменьшаем скорость
      for (int d = 0; d < coords; d++) velocities [particleIndex].v [d] *= -0.5; // Оптимизированная форма деления на 2

      velocities [particleIndex].slope  = -1; // Меняем наклон на отрицательный
    }
    else
      if (currentSlope == -1 && fitnessDiff < 0) // Отрицательный наклон и ухудшение фитнеса
      {
        // Применяем механизм роения - движение к глобальному оптимуму
        for (int d = 0; d < coords; d++) velocities [particleIndex].v [d] += (cB [d] - a [particleIndex].c [d]) * movementFactor;

        velocities [particleIndex].slope = 0; // Сбрасываем наклон как неизвестный
      }

  // Проверка на нулевую скорость с использованием метода структуры
  if (velocities [particleIndex].IsZero ())
  {
    // Инициализируем скорость движением к глобальному оптимуму
    for (int d = 0; d < coords; d++) velocities [particleIndex].v [d] = (cB [d] - a [particleIndex].c [d]) * movementFactor;

    // Сбрасываем наклон
    velocities [particleIndex].slope  = 0;
  }
}
//——————————————————————————————————————————————————————————————————————————————