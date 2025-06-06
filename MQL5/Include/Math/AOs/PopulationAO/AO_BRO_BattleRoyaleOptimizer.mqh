//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_BRO |
//|                                            Copyright 2007-2025, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/17497

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
class C_AO_BRO : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_BRO () { }
  C_AO_BRO ()
  {
    ao_name = "BRO";
    ao_desc = "Battle Royale Optimizer";
    ao_link = "https://www.mql5.com/ru/articles/17688";

    popSize   = 50;    // размер популяции
    maxDamage = 3;     // максимальный порог повреждений

    ArrayResize (params, 2);

    params [0].name = "popSize";   params [0].val = popSize;
    params [1].name = "maxDamage"; params [1].val = maxDamage;
  }

  void SetParams ()
  {
    popSize   = (int)params [0].val;
    maxDamage = (int)params [1].val;
  }

  bool Init (const double &rangeMinP  [], // минимальный диапазон поиска
             const double &rangeMaxP  [], // максимальный диапазон поиска
             const double &rangeStepP [], // шаг поиска
             const int     epochsP = 0);  // количество эпох

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  int maxDamage;    // максимальный порог повреждений

  private: //-------------------------------------------------------------------
  int    delta;      // интервал для сужения пространства поиска
  int    damages []; // количество повреждений для каждого решения
  int    epoch;      // текущая эпоха
  int    epochs;     // максимальное количество эпох

  // Вспомогательные методы
  int    FindNearestNeighbor (int index);
  double CalculateDistance (int idx1, int idx2);
  void   CalculateStandardDeviations (double &sdValues []);
  void   ShrinkSearchSpace ();
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_BRO::Init (const double &rangeMinP  [],  // минимальный диапазон поиска
                     const double &rangeMaxP  [],  // максимальный диапазон поиска
                     const double &rangeStepP [],  // шаг поиска
                     const int     epochsP = 0)    // количество эпох
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  // Инициализация счетчиков повреждений для каждого решения
  ArrayResize (damages, popSize);
  ArrayInitialize (damages, 0);

  // Установка эпох
  epochs = epochsP;
  epoch  = 0;

  // Вычисление начального delta для сужения пространства поиска
  delta = (int)MathFloor (epochs / MathLog10 (epochs));
  if (delta <= 0) delta = 1;

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BRO::Moving ()
{
  if (!revision)
  {
    // Инициализация популяции случайными решениями
    for (int i = 0; i < popSize; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        double coordinate = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        a [i].c [c] = u.SeInDiSp (coordinate, rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }

    revision = true;
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BRO::Revision ()
{
  epoch++;

  // Обновление глобального наилучшего решения
  for (int i = 0; i < popSize; i++)
  {
    if (a [i].f > fB)
    {
      fB = a [i].f;
      ArrayCopy (cB, a [i].c, 0, 0, WHOLE_ARRAY);
    }
  }

  // Сравнение каждого решения с его ближайшим соседом и обновление счетчиков повреждений
  for (int i = 0; i < popSize; i++)
  {
    int neighbor = FindNearestNeighbor (i);

    if (neighbor != -1)
    {
      if (a [i].f >= a [neighbor].f)
      {
        // Решение i побеждает
        damages [i] = 0;
        damages [neighbor]++;

        // Проигравший (сосед) движется к наилучшему решению
        for (int c = 0; c < coords; c++)
        {
          double r = u.RNDfromCI (0, 1);
          a [neighbor].c [c] = a [neighbor].c [c] + r * (cB [c] - a [neighbor].c [c]);
          a [neighbor].c [c] = u.SeInDiSp (a [neighbor].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
        }
      }
      else
      {
        // Решение i проигрывает
        damages [i]++;
        damages [neighbor] = 0;

        // Проигравший (i) движется к наилучшему решению
        for (int c = 0; c < coords; c++)
        {
          double r = u.RNDfromCI (0, 1);
          a [i].c [c] = a [i].c [c] + r * (cB [c] - a [i].c [c]);
          a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
        }
      }
    }
  }

  // Проверка, достигло ли какое-либо решение максимального повреждения, и его замена
  for (int i = 0; i < popSize; i++)
  {
    if (damages [i] >= maxDamage)
    {
      // Сброс счетчика повреждений
      damages [i] = 0;

      // Генерация нового случайного решения
      for (int c = 0; c < coords; c++)
      {
        double coordinate = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        a [i].c [c] = u.SeInDiSp (coordinate, rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }
  }

  // Периодическое сужение пространства поиска
  if (epochs > 0 && epoch % delta == 0)
  {
    ShrinkSearchSpace ();
    // Обновление delta
    delta = delta + (int)MathRound (delta / 2);
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
int C_AO_BRO::FindNearestNeighbor (int index)
{
  double minDistance = DBL_MAX;
  int nearestIndex = -1;

  for (int i = 0; i < popSize; i++)
  {
    if (i == index) continue;

    double distance = CalculateDistance (index, i);
    if (distance < minDistance)
    {
      minDistance = distance;
      nearestIndex = i;
    }
  }

  return nearestIndex;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_BRO::CalculateDistance (int idx1, int idx2)
{
  double distanceSum = 0.0;

  for (int c = 0; c < coords; c++)
  {
    double diff = a [idx1].c [c] - a [idx2].c [c];
    distanceSum += diff * diff;
  }

  return MathSqrt (distanceSum);
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BRO::CalculateStandardDeviations (double &sdValues [])
{
  ArrayResize (sdValues, coords);

  for (int c = 0; c < coords; c++)
  {
    double sum = 0.0;
    double mean = 0.0;

    // Вычисление среднего
    for (int i = 0; i < popSize; i++) mean += a [i].c [c];

    mean /= popSize;

    // Вычисление суммы квадратов отклонений
    for (int i = 0; i < popSize; i++)
    {
      double diff = a [i].c [c] - mean;
      sum += diff * diff;
    }

    sdValues [c] = MathSqrt (sum / popSize);
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BRO::ShrinkSearchSpace ()
{
  double sdValues [];
  CalculateStandardDeviations (sdValues);

  for (int c = 0; c < coords; c++)
  {
    // Новые границы центрированы вокруг наилучшего решения с шириной стандартного отклонения
    double newMin = cB [c] - sdValues [c];
    double newMax = cB [c] + sdValues [c];

    // Убедитесь, что новые границы находятся в пределах исходных ограничений
    if (newMin < rangeMin [c]) newMin = rangeMin [c];
    if (newMax > rangeMax [c]) newMax = rangeMax [c];

    // Обновление границ
    rangeMin [c] = newMin;
    rangeMax [c] = newMax;
  }
}
//——————————————————————————————————————————————————————————————————————————————
