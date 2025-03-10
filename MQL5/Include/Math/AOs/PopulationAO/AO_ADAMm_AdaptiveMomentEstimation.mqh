//+————————————————————————————————————————————————————————————————————————————+
//|                                                                 C_AO_ADAMm |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/16443

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
// Структура для хранения градиентов и моментов
struct S_Gradients
{
    double g [];  // Градиенты
    double m [];  // Векторы первого момента
    double v [];  // Векторы второго момента

    // Метод инициализации градиентов
    void Init (int coords)
    {
      ArrayResize (g, coords);
      ArrayResize (m, coords);
      ArrayResize (v, coords);

      ArrayInitialize (g, 0.0); // Инициализация градиентов нулями
      ArrayInitialize (m, 0.0); // Инициализация первого момента нулями
      ArrayInitialize (v, 0.0); // Инициализация второго момента нулями
    }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
class C_AO_ADAMm : public C_AO
{
  public: //--------------------------------------------------------------------

  // Деструктор класса
  ~C_AO_ADAMm () { }

  // Конструктор класса
  C_AO_ADAMm ()
  {
    ao_name = "ADAMm";                                  // Название алгоритма
    ao_desc = "Adaptive Moment Estimation M";           // Описание алгоритма
    ao_link = "https://www.mql5.com/ru/articles/16443"; // Ссылка на статью

    popSize           = 100;    // Размер популяции
    hybridsPercentage = 0.5;    // Процент гибридов в популяции
    hybridsResistance = 10;     // Устойчивость гибридов к изменениям
    alpha             = 0.001;  // Коэффициент обучения
    beta1             = 0.9;    // Коэффициент экспоненциального затухания для первого момента
    beta2             = 0.999;  // Коэффициент экспоненциального затухания для второго момента
    epsilon           = 0.1;    // Маленькая константа для предотвращения деления на ноль

    // Инициализация массива параметров
    ArrayResize (params, 7);
    params [0].name = "popSize";           params [0].val = popSize;
    params [1].name = "hybridsPercentage"; params [1].val = hybridsPercentage;
    params [2].name = "hybridsResistance"; params [2].val = hybridsResistance;
    params [3].name = "alpha";             params [3].val = alpha;
    params [4].name = "beta1";             params [4].val = beta1;
    params [5].name = "beta2";             params [5].val = beta2;
    params [6].name = "epsilon";           params [6].val = epsilon;
  }

  // Метод для установки параметров
  void SetParams ()
  {
    popSize           = (int)params [0].val;   // Установка размера популяции
    hybridsPercentage = params      [1].val;   // Установка процента гибридов в популяции
    hybridsResistance = params      [2].val;   // Установка устойчивости гибридов к изменениям
    alpha             = params      [3].val;   // Установка коэффициента обучения
    beta1             = params      [4].val;   // Установка beta1
    beta2             = params      [5].val;   // Установка beta2
    epsilon           = params      [6].val;   // Установка epsilon
  }

  // Метод инициализации
  bool Init (const double &rangeMinP  [],  // минимальный диапазон поиска
             const double &rangeMaxP  [],  // максимальный диапазон поиска
             const double &rangeStepP [],  // шаг поиска
             const int     epochsP = 0);   // количество эпох

  void Moving   (); // Метод для перемещения
  void Revision (); // Метод для ревизии

  //----------------------------------------------------------------------------
  double hybridsPercentage;  // Процент гибридов в популяции
  double hybridsResistance;  // Устойчивость гибридов к изменениям
  double alpha;              // Коэффициент обучения
  double beta1;              // Коэффициент экспоненциального затухания для первого момента
  double beta2;              // Коэффициент экспоненциального затухания для второго момента
  double epsilon;            // Маленькая константа

  S_Gradients grad []; // Массив градиентов

  private: //-------------------------------------------------------------------
  int step;          // Шаг итерации
  int t;             // Счетчик итераций
  int hybridsNumber; // Число гибридов в популяции
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_ADAMm::Init (const double &rangeMinP  [],
                       const double &rangeMaxP  [],
                       const double &rangeStepP [],
                       const int     epochsP = 0)
{
  // Стандартная инициализация
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  step          = 0;                                        // Сброс счетчика шагов
  t             = 1;                                        // Сброс счетчика итераций
  hybridsNumber = int(popSize * hybridsPercentage);         // Расчет числа гибридов в популяции
  if (hybridsNumber > popSize) hybridsNumber = popSize;     // Корректировка

  ArrayResize (grad, popSize);                              // Изменение размера массива градиентов
  for (int i = 0; i < popSize; i++) grad [i].Init (coords); // Инициализация градиентов для каждого индивида

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_ADAMm::Moving ()
{
  //----------------------------------------------------------------------------
  if (step < 2) // Если шаг меньше 2
  {
    for (int i = 0; i < popSize; i++)
    {
      a [i].fP = a [i].f; // Сохранение предыдущего значения функции

      for (int c = 0; c < coords; c++)
      {
        a [i].cP [c] = a [i].c [c]; // Сохранение предыдущего значения координат

        // Генерация новых координат случайным образом
        a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        // Приведение новых координат к допустимым значениям
        a [i].c [c] = u.SeInDiSp  (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }

    step++; // Увеличение счетчика шагов
    return; // Выход из метода
  }

  //----------------------------------------------------------------------------
  double ΔF, ΔX; // Изменения функции и координат
  double cNew;

  for (int i = 0; i < popSize; i++)
  {
    ΔF = a [i].f - a [i].fP;           // Вычисление изменения функции

    for (int c = 0; c < coords; c++)
    {
      ΔX = a [i].c [c] - a [i].cP [c]; // Вычисление изменения координат

      if (ΔX == 0.0) ΔX = epsilon;     // Если изменение равно нулю, установить его на epsilon

      grad [i].g [c] = ΔF / ΔX;        // Вычисление градиента
    }
  }

  // Обновление параметров с использованием алгоритма ADAM
  for (int i = 0; i < popSize; i++)
  {
    // Сохранение предыдущего значения функции
    a [i].fP = a [i].f;

    for (int c = 0; c < coords; c++)
    {
      // Сохранение предыдущего значения координат
      a [i].cP [c] = a [i].c [c];

      if (i >= popSize - hybridsNumber)
      {
        double pr = u.RNDprobab ();
        pr *= pr;

        int ind = (int)u.Scale (pr, 0, 1, 0, popSize - 1);

        cNew = u.PowerDistribution (a [ind].c [c], rangeMin [c], rangeMax [c], hybridsResistance);
      }
      else
      {
        // Обновление смещенной оценки первого момента
        grad [i].m [c] = beta1 * grad [i].m [c] + (1.0 - beta1) * grad [i].g [c];

        // Обновление смещенной оценки второго момента
        grad [i].v [c] = beta2 * grad [i].v [c] + (1.0 - beta2) * grad [i].g [c] * grad [i].g [c];

        // Вычисление скорректированной оценки первого момента
        double m_hat = grad [i].m [c] / (1.0 - MathPow (beta1, t));

        // Вычисление скорректированной оценки второго момента
        double v_hat = grad [i].v [c] / (1.0 - MathPow (beta2, t));

        // Обновление координат
        //a [i].c [c] = a [i].c [c] + (alpha * m_hat / (MathSqrt (v_hat) + epsilon));
        cNew = a [i].c [c] + (alpha * m_hat / (MathSqrt (v_hat) + epsilon));
      }


      // Убедитесь, что координаты остаются в допустимых границах
      a [i].c [c] = u.SeInDiSp (cNew, rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }

  t++; // Увеличение счетчика итераций
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_ADAMm::Revision ()
{
  int ind = -1;       // Индекс лучшего индивида
  for (int i = 0; i < popSize; i++)
  {
    if (a [i].f > fB) // Если текущее значение функции больше лучшего
    {
      fB = a [i].f;   // Обновление лучшего значения функции
      ind = i;        // Сохранение индекса лучшего индивида
    }
  }

  if (ind != -1) ArrayCopy (cB, a [ind].c, 0, 0, WHOLE_ARRAY); // Копирование координат лучшего индивида

  //----------------------------------------------------------------------------
  S_AO_Agent aT [];
  ArrayResize (aT, popSize);
  u.Sorting (a, aT, popSize);
}
//——————————————————————————————————————————————————————————————————————————————
