//+————————————————————————————————————————————————————————————————————————————+
//|                                                                  C_AO_AMOm |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/15543

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
class C_AO_AMOm : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_AMOm () { }
  C_AO_AMOm ()
  {
    ao_name = "AMOm";
    ao_desc = "Animal Migration Optimization M";
    ao_link = "https://www.mql5.com/ru/articles/15543";

    popSize               = 50;   // Размер популяции
    deviation             = 8;
    neighborsNumberOnSide = 10;

    ArrayResize (params, 3);

    params [0].name = "popSize";               params [0].val = popSize;

    params [1].name = "deviation";             params [1].val = deviation;
    params [2].name = "neighborsNumberOnSide"; params [2].val = neighborsNumberOnSide;
  }

  void SetParams ()
  {
    popSize               = (int)params [0].val;

    deviation             = params      [1].val;
    neighborsNumberOnSide = (int)params [2].val;
  }

  bool Init (const double &rangeMinP  [],
             const double &rangeMaxP  [],
             const double &rangeStepP [],
             const int     epochsP = 0);

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  double deviation;
  int    neighborsNumberOnSide;

  S_AO_Agent population []; // Популяция животных
  S_AO_Agent pTemp      []; // Временная популяция животных

  private: //-------------------------------------------------------------------
  int   GetNeighborsIndex (int i);
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_AMOm::Init (const double &rangeMinP  [],
                      const double &rangeMaxP  [],
                      const double &rangeStepP [],
                      const int     epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  ArrayResize (population, popSize * 2);
  ArrayResize (pTemp,      popSize * 2);

  for (int i = 0; i < popSize * 2; i++) population [i].Init (coords);

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_AMOm::Moving ()
{
  //----------------------------------------------------------------------------
  if (!revision)
  {
    for (int i = 0; i < popSize; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        a [i].c [c] = u.SeInDiSp  (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }

    revision = true;
    return;
  }

  //----------------------------------------------------------------------------
  int    ind1    = 0;
  int    ind2    = 0;
  double dist    = 0.0;
  double x       = 0.0;
  double min     = 0.0;
  double max     = 0.0;
  double prob    = 0.0;

  for (int i = 0; i < popSize; i++)
  {
    prob = 1.0 - (1.0 / (i + 1));
    
    for (int c = 0; c < coords; c++)
    {
      //------------------------------------------------------------------------
      ind1 = GetNeighborsIndex (i);

      dist = fabs (population [ind1].c [c] - a [i].c [c]);

      x    = population [ind1].c [c];
      min  = x - dist;
      max  = x + dist;

      if (min < rangeMin [c]) min = rangeMin [c];
      if (max > rangeMax [c]) max = rangeMax [c];

      a [i].c [c] = u.GaussDistribution (x, min, max, deviation);

      //------------------------------------------------------------------------
      if (u.RNDprobab() < prob)
      {
        if (u.RNDprobab() <= 0.01)
        {
          ind1 = u.RNDminusOne (popSize);
          ind2 = u.RNDminusOne (popSize);

          a [i].c [c] = (population [ind1].c [c] + population [ind2].c [c]) * 0.5;
        }
      }
      
      a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_AMOm::Revision ()
{
  //----------------------------------------------------------------------------
  int ind = -1;

  for (int i = 0; i < popSize; i++)
  {
    if (a [i].f > fB)
    {
      fB  = a [i].f;
      ind = i;
    }
  }

  if (ind != -1) ArrayCopy (cB, a [ind].c, 0, 0, WHOLE_ARRAY);

  //----------------------------------------------------------------------------
  for (int i = 0; i < popSize; i++) population [popSize + i] = a [i];

  u.Sorting (population, pTemp, popSize * 2);
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
int C_AO_AMOm::GetNeighborsIndex (int i)
{
  int Ncount = neighborsNumberOnSide;
  int N = Ncount * 2 + 1;
  int neighborIndex;

  // Выбор случайного соседа с учетом границ массива
  if (i < Ncount)
  {
    // Для первых Ncount элементов
    neighborIndex = MathRand () % N;
  }
  else
  {
    if (i >= popSize - Ncount)
    {
      // Для последних Ncount элементов
      neighborIndex = (popSize - N) + MathRand () % N;
    }
    else
    {
      // Для всех остальных элементов
      neighborIndex = i - Ncount + MathRand () % N;
    }
  }

  return neighborIndex;
}
//——————————————————————————————————————————————————————————————————————————————