//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_SOA |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/16364

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
class C_AO_SOA : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_SOA () { }
  C_AO_SOA ()
  {
    ao_name = "SOA";
    ao_desc = "Simple Optimization Algorithm";
    ao_link = "https://www.mql5.com/ru/articles/16364";

    popSize = 50;   // Размер популяции
    minT    = 0.1;  // Минимальное значение T
    maxT    = 0.9;  // Максимальное значение T
    θ       = 10;   // Параметр θ

    ArrayResize (params, 4); // Изменение размера массива параметров

    // Инициализация параметров
    params [0].name = "popSize"; params [0].val = popSize;
    params [1].name = "minT";    params [1].val = minT;
    params [2].name = "maxT";    params [2].val = maxT;
    params [3].name = "θ";       params [3].val = θ;
  }

  void SetParams () // Метод для установки параметров
  {
    popSize = (int)params [0].val; // Установка размера популяции
    minT    = params      [1].val; // Установка минимального T
    maxT    = params      [2].val; // Установка максимального T
    θ       = params      [3].val; // Установка параметра θ
  }

  bool Init (const double &rangeMinP  [], // Минимальный диапазон поиска
             const double &rangeMaxP  [], // Максимальный диапазон поиска
             const double &rangeStepP [], // Шаг поиска
             const int     epochsP = 0);  // Количество эпох

  void Moving   (); // Метод перемещения частиц
  void Revision (); // Метод ревизии

  //----------------------------------------------------------------------------
  double minT; // Минимальное значение T
  double maxT; // Максимальное значение T
  double θ;    // Параметр θ

  private: //-------------------------------------------------------------------
  int epochs;    // Общее количество эпох
  int epochNow;  // Текущая эпоха
  double ϵ;      // Параметр для предотвращения деления на ноль
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_SOA::Init (const double &rangeMinP  [], // Минимальный диапазон поиска
                     const double &rangeMaxP  [], // Максимальный диапазон поиска
                     const double &rangeStepP [], // Шаг поиска
                     const int     epochsP = 0)   // Количество эпох
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false; // Инициализация стандартных параметров

  //----------------------------------------------------------------------------
  epochs   = epochsP;     // Установка общего количества эпох
  epochNow = 0;           // Инициализация текущей эпохи
  ϵ        = DBL_EPSILON; // Установка значения ϵ

  return true;            // Возвращаем true при успешной инициализации
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
// Метод перемещения частиц
void C_AO_SOA::Moving ()
{
  epochNow++; // Увеличиваем номер текущей эпохи

  // Начальное случайное позиционирование
  if (!revision) // Если еще не было ревизии
  {
    for (int i = 0; i < popSize; i++) // Для каждой частицы
    {
      for (int c = 0; c < coords; c++) // Для каждой координаты
      {
        a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);                             // Генерация случайной позиции
        a [i].c [c] = u.SeInDiSp  (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]); // Приведение к дискретным значениям
      }
    }

    revision = true; // Установка флага ревизии
    return;          // Выход из метода
  }

  //----------------------------------------------------------------------------
  double MoAc = minT + epochNow * ((maxT - minT) / epochs); // Вычисление значения MoAc
  double MoPr = 1.0 - pow (epochNow / epochs, (1.0 / θ));   // Вычисление значения MoPr
  double best = 0.0;                                        // Переменная для хранения лучшего значения

  // Фаза исследования с использованием операторов Деления (D) и Умножения (M)
  for (int i = 0; i < popSize; i++) // Для каждой частицы
  {
    for (int c = 0; c < coords; c++) // Для каждой координаты
    {
      // Вероятностное обновление позиции частицы
      if (u.RNDbool () < MoAc) a [i].c [c] = cB [c];                                        // Установка на лучшее значение
      else
        if (u.RNDbool () < MoPr) a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);    // Генерация новой случайной позиции

      a [i].c [c] = u.SeInDiSp  (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);   // Приведение к дискретным значениям
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_SOA::Revision ()
{
  int ind = -1;                     // Индекс для хранения лучшей частицы

  for (int i = 0; i < popSize; i++) // Для каждой частицы
  {
    if (a [i].f > fB)               // Если значение функции лучше, чем текущее лучшее
    {
      fB = a [i].f;                 // Обновление лучшего значения функции
      ind = i;                      // Сохранение индекса лучшей частицы
    }
  }

  if (ind != -1) ArrayCopy (cB, a [ind].c, 0, 0, WHOLE_ARRAY); // Копирование координат лучшей частицы
}
//——————————————————————————————————————————————————————————————————————————————