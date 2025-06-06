//+————————————————————————————————————————————————————————————————————————————+
//|                                                                  C_AO_BBBC |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/16701

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
class C_AO_BBBC : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_BBBC () { }
  C_AO_BBBC ()
  {
    ao_name = "BBBC";
    ao_desc = "Big Bang - Big Crunch Algorithm";
    ao_link = "https://www.mql5.com/ru/articles/16701";

    popSize       = 50;
    bigBangPeriod = 3;

    ArrayResize (params, 2);
    params [0].name = "popSize";       params [0].val = popSize;
    params [1].name = "bigBangPeriod"; params [1].val = bigBangPeriod;
  }

  void SetParams ()
  {
    popSize       = (int)params [0].val;
    bigBangPeriod = (int)params [1].val;
  }

  bool Init (const double &rangeMinP  [],  // минимальный диапазон поиска
             const double &rangeMaxP  [],  // максимальный диапазон поиска
             const double &rangeStepP [],  // шаг поиска
             const int epochsP = 0);       // количество эпох

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  int bigBangPeriod;       // периодичность большого взрыва


  private: //-------------------------------------------------------------------
  int epochs;              // общее число эпох
  int epochNow;            // текущая эпоха
  double centerMass [];    // центр масс
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_BBBC::Init (const double &rangeMinP  [],
                      const double &rangeMaxP  [],
                      const double &rangeStepP [],
                      const int epochsP = 0)
{
  // Инициализация базового класса
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  epochs   = epochsP;
  epochNow = 0;

  // Выделение памяти для массивов
  ArrayResize (centerMass, coords);

  return true;
}
//——————————————————————————————————————————————————————————————————————————————
/*
//——————————————————————————————————————————————————————————————————————————————
void C_AO_BBBC::Moving ()
{
  epochNow++;

  // Начальная инициализация (Big Bang)
  if (!revision)
  {
    for (int i = 0; i < popSize; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        // Генерация случайных начальных позиций
        a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        // Приведение к дискретной сетке поиска
        a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }
    revision = true;
    return;
  }

  //----------------------------------------------------------------------------
  // Фаза Big Crunch - большой коллапс
  if (epochNow % bigBangPeriod != 0)
  {
    for (int c = 0; c < coords; c++)
    {
      double numerator = 0;
      double denominator = 0;

      for (int i = 0; i < popSize; i++)
      {
        // Расчет веса как обратной величины от значения фитнес-функции
        double fitness = MathMax (MathAbs (a [i].f), 1e-10);
        double weight = 1.0 / fitness;

        // Суммирование для вычисления центра масс по формуле
        // xc = (Σ(1/fi)xi) / (Σ(1/fi))
        numerator += weight * a [i].c [c];
        denominator += weight;
      }

      // Определение координаты центра масс
      centerMass [c] = denominator > 1e-10 ? numerator / denominator : cB [c];
    }

    for (int i = 0; i < popSize; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        double r = u.GaussDistribution (0, -1.0, 1.0, 1);

        // Генерация новой точки по формуле
        // xnew = xc + r*xmax/k
        double newPoint = centerMass [c] + r * rangeMax [c] / epochNow;

        // Ограничение в пределах допустимой области и приведение к сетке
        a [i].c [c] = u.SeInDiSp (newPoint, rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }
  }
  //----------------------------------------------------------------------------
  // Фаза Big Bang - большой взрыв
  else
  {
    for (int i = 0; i < popSize; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        a [i].c [c] = u.GaussDistribution (cB [c], rangeMin [c], rangeMax [c], 8);
        a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————
*/

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BBBC::Moving ()
{
  epochNow++;

  // Начальная инициализация (Big Bang)
  if (!revision)
  {
    for (int i = 0; i < popSize; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        // Генерация случайных начальных позиций
        a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        // Приведение к дискретной сетке поиска
        a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }
    revision = true;
    return;
  }

  //--------------------------------------------------------------------------
  for (int i = 0; i < popSize; i++)
  {
    //Фаза Big Crunch - большой коллапс
    if (epochNow % bigBangPeriod != 0)
    {
      for (int c = 0; c < coords; c++)
      {
        // Вычисление размера пространства поиска для текущей координаты
        double range = rangeMax [c] - rangeMin [c];

        // Генерация случайного числа в диапазоне [-1, 1]
        double r = u.GaussDistribution (0, -1.0, 1.0, 1);

        // Генерация новой точки по формуле
        // xnew = xc + r*(xmax - xmin)/(k)
        double newPoint = cB [c] + r * range / epochNow;

        // Ограничение в пределах допустимой области и приведение к сетке
        a [i].c [c] = u.SeInDiSp (newPoint, rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }
    // Фаза Big Bang - большой взрыв
    else
    {
      for (int c = 0; c < coords; c++)
      {
        a [i].c [c] = u.GaussDistribution (cB [c], rangeMin [c], rangeMax [c], 8);
        a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BBBC::Revision ()
{
  int bestInd = -1;

  // Поиск лучшего решения в текущей популяции
  for (int i = 0; i < popSize; i++)
  {
    if (a [i].f > fB)
    {
      fB = a [i].f;
      bestInd = i;
    }
  }

  // Обновление лучшего известного решения
  if (bestInd != -1) ArrayCopy (cB, a [bestInd].c, 0, 0, WHOLE_ARRAY);
}
//——————————————————————————————————————————————————————————————————————————————