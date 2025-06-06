//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_CFO |
//|                                            Copyright 2007-2025, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/17167

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
//--- Структура для пробы CFO
struct S_CFO_Agent : public S_AO_Agent
{
    double a []; // вектор ускорения

    void Init (int coords)
    {
      ArrayResize (c, coords);   // координаты
      ArrayResize (a, coords);   // ускорение
      ArrayInitialize (a, 0.0);  // обнуляем ускорения
      f = -DBL_MAX;              // значение фитнес-функции
    }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
//--- Основной класс алгоритма CFO
class C_AO_CFO : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_CFO () { }
  C_AO_CFO ()
  {
    ao_name = "CFO";
    ao_desc = "Central Force Optimization";
    ao_link = "https://www.mql5.com/ru/articles/17167";

    popSize     = 30;          // число проб
    g           = 1.0;         // гравитационная постоянная
    alpha       = 0.1;         // степень для массы
    beta        = 0.1;         // степень для расстояния
    initialFrep = 0.9;         // начальный фактор репозиционирования
    finalFrep   = 0.1;         // конечный фактор репозиционирования
    noiseFactor = 1.0;         // фактор случайного шума

    frep        = initialFrep; // текущий фактор репозиционирования

    ArrayResize (params, 7);
    params [0].name = "popSize";     params [0].val = popSize;
    params [1].name = "g";           params [1].val = g;
    params [2].name = "alpha";       params [2].val = alpha;
    params [3].name = "beta";        params [3].val = beta;
    params [4].name = "initialFrep"; params [4].val = initialFrep;
    params [5].name = "finalFrep";   params [5].val = finalFrep;
    params [6].name = "noiseFactor"; params [6].val = noiseFactor;
  }

  void SetParams ()
  {
    popSize     = (int)MathMax (1, params [0].val);
    g           = params [1].val;
    alpha       = params [2].val;
    beta        = params [3].val;
    initialFrep = params [4].val;
    finalFrep   = params [5].val;
    noiseFactor = params [6].val;

    frep        = initialFrep;
  }

  bool Init (const double &rangeMinP  [],  // минимальные значения
             const double &rangeMaxP  [],  // максимальные значения
             const double &rangeStepP [],  // шаг изменения
             const int     epochsP = 0);   // количество эпох

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  double g;              // гравитационная постоянная
  double alpha;          // степень для массы
  double beta;           // степень для расстояния
  double initialFrep;    // начальный фактор репозиционирования
  double finalFrep;      // конечный фактор репозиционирования
  double noiseFactor;    // фактор случайного шума

  S_CFO_Agent probe [];  // массив проб

  private: //-------------------------------------------------------------------
  int    epochs;         // общее число эпох
  int    epochNow;       // текущая эпоха
  double frep;           // фактор репозиционирования

  void   InitialDistribution      ();
  void   UpdateRepFactor          ();
  void   CalculateAccelerations   ();
  void   UpdatePositions          ();
  double CalculateDistanceSquared (const double &x1 [], const double &x2 []);
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
//--- Инициализация
bool C_AO_CFO::Init (const double &rangeMinP  [], // минимальные значения
                     const double &rangeMaxP  [], // максимальные значения
                     const double &rangeStepP [], // шаг изменения
                     const int     epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  epochs   = epochsP;
  epochNow = 0;

  // Инициализация проб
  ArrayResize (probe, popSize);
  for (int i = 0; i < popSize; i++) probe [i].Init (coords);

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
//--- Основной шаг алгоритма
void C_AO_CFO::Moving ()
{
  epochNow++;

  // Начальная инициализация
  if (!revision)
  {
    InitialDistribution ();
    revision = true;
    return;
  }

  //----------------------------------------------------------------------------
  // Копируем значения фитнес-функции из базового класса
  for (int p = 0; p < popSize; p++)
  {
    probe [p].f = a [p].f;
  }

  // Обновляем параметр репозиционирования
  UpdateRepFactor ();

  // Основной цикл CFO
  CalculateAccelerations ();
  UpdatePositions ();

  // Синхронизируем позиции с базовым классом
  for (int p = 0; p < popSize; p++)
  {
    ArrayCopy (a [p].c, probe [p].c);
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
//--- Начальное распределение проб
void C_AO_CFO::InitialDistribution ()
{
  for (int p = 0; p < popSize; p++)
  {
    ArrayInitialize (probe [p].a, 0.0);
    probe [p].f = -DBL_MAX;

    // Случайное распределение
    for (int c = 0; c < coords; c++)
    {
      probe [p].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
      probe [p].c [c] = u.SeInDiSp (probe [p].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }

    ArrayCopy (a [p].c, probe [p].c);
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
//--- Обновление фактора репозиционирования
void C_AO_CFO::UpdateRepFactor ()
{
  // Линейное уменьшение frep от начального к конечному значению
  if (epochs > 0) frep = initialFrep - (initialFrep - finalFrep) * epochNow / epochs;
  else frep = initialFrep;

  // Ограничение значения
  frep = MathMax (0.0, MathMin (1.0, frep));
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
//--- Вычисление ускорений
void C_AO_CFO::CalculateAccelerations ()
{
  for (int p = 0; p < popSize; p++)
  {
    // Обнуляем ускорение для текущей пробы
    ArrayInitialize (probe [p].a, 0.0);

    // Суммируем влияние всех других проб
    for (int k = 0; k < popSize; k++)
    {
      if (k == p) continue;

      // Разница масс (фитнес-значений)
      double massDiff = probe [k].f - probe [p].f;

      // Проверяем условие единичной функции U(...)
      if (massDiff > 0) // Строгое условие для единичной функции
      {
        // Вычисляем расстояние между пробами
        double distSquared = CalculateDistanceSquared (probe [k].c, probe [p].c);
        if (distSquared < DBL_EPSILON) continue;

        double distance = MathSqrt (distSquared);

        for (int c = 0; c < coords; c++)
        {
          // Направление силы
          double direction = (probe [k].c [c] - probe [p].c [c]) / distance;

          // Формула ускорения из статьи
          probe [p].a [c] += g * MathPow (massDiff, alpha) * direction / MathPow (distance, beta);
        }
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
//--- Обновление позиций
void C_AO_CFO::UpdatePositions ()
{
  // Коэффициент случайного шума, уменьшается с ростом номера эпохи
  double currentNoiseFactor = noiseFactor;
  if (epochs > 0) currentNoiseFactor *= (1.0 - (double)epochNow / epochs);

  for (int p = 0; p < popSize; p++)
  {
    for (int c = 0; c < coords; c++)
    {
      // Обновление позиции по формуле 13 из статьи (упрощенный вариант)
      probe [p].c [c] += 0.5 * probe [p].a [c];

      // Добавление небольшого случайного смещения непосредственно к позиции
      probe [p].c [c] += currentNoiseFactor * g * u.RNDfromCI (-1.0, 1.0);

      // Репозиционирование при выходе за границы
      if (probe [p].c [c] < rangeMin [c]) probe [p].c [c] = rangeMin [c];
      if (probe [p].c [c] > rangeMax [c]) probe [p].c [c] = rangeMax [c];

      // Дискретизация если задан шаг
      probe [p].c [c] = u.SeInDiSp (probe [p].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
//--- Вычисление расстояния (возвращает квадрат расстояния для оптимизации)
double C_AO_CFO::CalculateDistanceSquared (const double &x1 [], const double &x2 [])
{
  double sum = 0.0;
  for (int i = 0; i < coords; i++)
  {
    double diff = x1 [i] - x2 [i];
    sum += diff * diff;
  }
  return sum;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
//--- Обновление лучшего решения
void C_AO_CFO::Revision ()
{
  for (int p = 0; p < popSize; p++)
  {
    if (a [p].f > fB)
    {
      fB = a [p].f;
      ArrayCopy (cB, a [p].c, 0, 0, WHOLE_ARRAY);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————
