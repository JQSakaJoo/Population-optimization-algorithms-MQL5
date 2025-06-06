//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_CGO |
//|                                            Copyright 2007-2025, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/17047

#include "#C_AO.mqh"

class C_AO_CGO : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_CGO () { }
  C_AO_CGO ()
  {
    ao_name = "CGO";
    ao_desc = "Chaos Game Optimization";
    ao_link = "https://www.mql5.com/ru/articles/17047";

    popSize = 50;

    ArrayResize (params, 1);
    params [0].name = "popSize"; params [0].val = popSize;
  }

  void SetParams ()
  {
    popSize = (int)params [0].val;
  }

  bool Init (const double &rangeMinP  [],  // минимальные значения
             const double &rangeMaxP  [],  // максимальные значения
             const double &rangeStepP [],  // шаг изменения
             const int     epochsP = 0);   // количество эпох

  void Moving ();
  void Revision ();

  private: //-------------------------------------------------------------------
  double GetAlpha ();
  void   GenerateNewSolution (int seedIndex, double &meanGroup []);
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_CGO::Init (const double &rangeMinP  [], // минимальные значения
                     const double &rangeMaxP  [], // максимальные значения
                     const double &rangeStepP [], // шаг изменения
                     const int     epochsP = 0)   // количество эпох
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_CGO::Moving ()
{
  //----------------------------------------------------------------------------
  if (!revision)
  {
    for (int i = 0; i < popSize; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }

    revision = true;
    return;
  }

  //----------------------------------------------------------------------------
  for (int i = 0; i < popSize; i++)
  {
    int randGroupSize = u.RNDminusOne (popSize) + 1;
    double meanGroup [];
    ArrayResize (meanGroup, coords);
    ArrayInitialize (meanGroup, 0);

    int randIndices [];
    ArrayResize (randIndices, randGroupSize);

    for (int j = 0; j < randGroupSize; j++) randIndices [j] = u.RNDminusOne (popSize);

    for (int j = 0; j < randGroupSize; j++)
    {
      for (int c = 0; c < coords; c++)
      {
        meanGroup [c] += a [randIndices [j]].c [c];
      }
    }
    for (int c = 0; c < coords; c++) meanGroup [c] /= randGroupSize;

    GenerateNewSolution (i, meanGroup);
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_CGO::Revision ()
{
  for (int i = 0; i < popSize; i++)
  {
    if (a [i].f > fB)
    {
      fB = a [i].f;
      ArrayCopy (cB, a [i].c, 0, 0, WHOLE_ARRAY);
    }
  }

  //static S_AO_Agent aT []; ArrayResize (aT, popSize);
  //u.Sorting (a, aT, popSize);
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_CGO::GetAlpha ()
{
  int Ir = u.RNDminusOne (2);

  switch (u.RNDminusOne (4))
  {
    case 0: return u.RNDfromCI (0, 1);
    case 1: return 2 * u.RNDfromCI (0, 1) - 1;
    case 2: return Ir * u.RNDfromCI (0, 1) + 1;
    case 3: return Ir * u.RNDfromCI (0, 1) + (1 - Ir);
  }
  return 0;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_CGO::GenerateNewSolution (int seedIndex, double &meanGroup [])
{
  double alpha = GetAlpha ();
  int    beta  = u.RNDminusOne (2) + 1;
  int    gamma = u.RNDminusOne (2) + 1;

  int formula = u.RNDminusOne (4);

  for (int c = 0; c < coords; c++)
  {
    double newPos = 0;

    switch (formula)
    {
      case 0:
        newPos = a [seedIndex].c [c] + alpha * (beta * cB [c] - gamma * meanGroup [c]);
        break;

      case 1:
        newPos = cB [c] + alpha * (beta * meanGroup [c] - gamma * a [seedIndex].c [c]);
        break;

      case 2:
        newPos = meanGroup [c] + alpha * (beta * cB [c] - gamma * a [seedIndex].c [c]);
        break;

      case 3:
        newPos = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        break;
    }

    a [seedIndex].c [c] = u.SeInDiSp (newPos, rangeMin [c], rangeMax [c], rangeStep [c]);
  }
}
//——————————————————————————————————————————————————————————————————————————————