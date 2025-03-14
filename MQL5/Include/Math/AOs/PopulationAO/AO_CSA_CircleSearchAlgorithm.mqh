//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_CSA |
//|                                            Copyright 2007-2025, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/17143

#include "#C_AO.mqh"



//——————————————————————————————————————————————————————————————————————————————
class C_AO_CSA : public C_AO
{
  public: //--------------------------------------------------------------------
  C_AO_CSA ()
  {
    ao_name = "CSA";
    ao_desc = "Circle Search Algorithm";
    ao_link = "https://www.mql5.com/ru/articles/17143";

    popSize = 50;     // размер популяции
    constC  = 0.8;    // оптимальное значение для фазы исследования
    w       = M_PI;   // начальное значение w
    aParam  = M_PI;   // начальное значение a
    p       = 1.0;    // начальное значение p
    theta   = 0;      // начальное значение угла

    ArrayResize (params, 2);
    params [0].name = "popSize";     params [0].val = popSize;
    params [1].name = "constC";      params [1].val = constC;
  }

  void SetParams ()
  {
    popSize = (int)params [0].val;
    constC  = params      [1].val;
  }

  bool Init (const double &rangeMinP  [],  // минимальные значения
             const double &rangeMaxP  [],  // максимальные значения
             const double &rangeStepP [],  // шаг изменения
             const int     epochsP = 0);   // количество эпох

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  double constC;      // константа для определения фазы поиска [0,1]

  private: //-------------------------------------------------------------------
  int epochs;         // максимальное число итераций
  int epochNow;       // текущая итерация

  // Параметры для CSA
  double w;           // параметр для вычисления угла
  double aParam;      // параметр a из формулы (8)
  double p;           // параметр p из формулы (9)
  double theta;       // угол поиска

  double CalculateW ();
  double CalculateA ();
  double CalculateP ();
  double CalculateTheta (double currentW, double currentP);
  bool IsExplorationPhase ();
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_CSA::Init (const double &rangeMinP  [],
                     const double &rangeMaxP  [],
                     const double &rangeStepP [],
                     const int     epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  epochs   = epochsP;
  epochNow = 0;
  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_CSA::Moving ()
{
  epochNow++;

  //----------------------------------------------------------------------------
  if (!revision)
  {
    for (int i = 0; i < popSize; i++)
    {
      for (int j = 0; j < coords; j++)
      {
        a [i].c [j] = u.RNDfromCI (rangeMin [j], rangeMax [j]);
        a [i].c [j] = u.SeInDiSp (a [i].c [j], rangeMin [j], rangeMax [j], rangeStep [j]);
      }
    }
    revision = true;
    return;
  }

  //----------------------------------------------------------------------------
  w      = CalculateW ();    // Обновляем w по линейному закону
  aParam = CalculateA ();    // Обновляем a
  p      = CalculateP ();    // Обновляем p

  for (int i = 0; i < popSize; i++)
  {
    theta = CalculateTheta (w, p);

    for (int j = 0; j < coords; j++)
    {
      a [i].c [j] = cB [j] + u.RNDprobab () * (cB [j] - a [i].c  [j]) * tan (theta);
      a [i].c [j] = u.SeInDiSp (a [i].c [j], rangeMin [j], rangeMax [j], rangeStep [j]);
    }

  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_CSA::Revision ()
{
  for (int i = 0; i < popSize; i++)
  {
    // Обновляем лучшее глобальное решение
    if (a [i].f > fB)
    {
      fB = a [i].f;
      ArrayCopy (cB, a [i].c, 0, 0, WHOLE_ARRAY);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_CSA::CalculateW ()
{
  // Линейное уменьшение w от начального значения (M_PI) до 0
  return M_PI * (1.0 - (double)epochNow / epochs);
  //return w * u.RNDprobab () - w;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_CSA::CalculateA ()
{
  return M_PI - M_PI * MathPow ((double)epochNow / epochs, 2);
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_CSA::CalculateP ()
{
  return 1.0 - 0.9 * MathPow ((double)epochNow / epochs, 0.5);
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_CSA::CalculateTheta (double currentW, double currentP)
{
  // Используем параметр aParam для регулировки угла
  if (IsExplorationPhase ()) return currentW * u.RNDprobab ();
  else return currentW * currentP;

}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_CSA::IsExplorationPhase ()
{
  // Исследование в первой части итераций (constC обычно 0.8)
  return (epochNow <= constC * epochs);
}
//——————————————————————————————————————————————————————————————————————————————