//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_ABO |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/16024

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
struct S_Buffalo
{
    double w [];
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
class C_AO_ABO : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_ABO () { }
  C_AO_ABO ()
  {
    ao_name = "ABO";
    ao_desc = "African Buffalo Optimization";
    ao_link = "https://www.mql5.com/ru/articles/16024";

    popSize = 50;    // размер популяции

    lp1     = 1.0;   // фактор обучения 1
    lp2     = 0.1;   // фактор обучения 2
    lambda  = 0.9;   // лямбда для уравнения движения

    ArrayResize (params, 4);

    params [0].name = "popSize"; params [0].val = popSize;
    params [1].name = "lp1";     params [1].val = lp1;
    params [2].name = "lp2";     params [2].val = lp2;
    params [3].name = "lambda";  params [3].val = lambda;
  }

  void SetParams ()
  {
    popSize = (int)params [0].val;
    lp1     = params      [1].val;
    lp2     = params      [2].val;
    lambda  = params      [3].val;
  }

  bool Init (const double &rangeMinP  [], //minimum search range
             const double &rangeMaxP  [], //maximum search range
             const double &rangeStepP [], //step search
             const int     epochsP = 0);  //number of epochs

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  double lp1;    // фактор обучения 1
  double lp2;    // фактор обучения 2
  double lambda; // лямбда для уравнения движения

  private: //-------------------------------------------------------------------
  S_Buffalo b [];
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_ABO::Init (const double &rangeMinP [],
                     const double &rangeMaxP [],
                     const double &rangeStepP [],
                     const int epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  ArrayResize (b, popSize);
  for (int i = 0; i < popSize; i++)
  {
    ArrayResize(b [i].w, coords);
    ArrayInitialize (b [i].w, 0.0);
  }

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_ABO::Moving ()
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
  double w  = 0.0;
  double m  = 0.0;
  double r1 = 0.0;
  double r2 = 0.0;
  double bg = 0.0;
  double bp = 0.0;

  for (int i = 0; i < popSize; i++)
  {
    for (int c = 0; c < coords; c++)
    {
      /*
      r1 = u.RNDfromCI (0, lp1);
      r2 = u.RNDfromCI (0, lp2);

      bg = cB [c];
      bp = a [i].cB [c];

      m = a [i].c [c];
      w = b [i].w [c];

      b [i].w [c] = w + r1 * (bg - m) + r2 * (bp - m);

      m = lambda * (m + b [i].w [c]);

      a [i].c [c] = u.SeInDiSp (m, rangeMin [c], rangeMax [c], rangeStep [c]);
      */
      
      r1 = u.RNDfromCI (-lp1, lp1);
      r2 = u.RNDfromCI (-lp2, lp2);

      bg = cB [c];
      bp = a [i].cB [c];
      
      m = a [i].c [c];

      m = m + r1 * (bg - m) + r2 * (bp - m);

      a [i].c [c] = u.SeInDiSp (m, rangeMin [c], rangeMax [c], rangeStep [c]);
      
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_ABO::Revision ()
{
  int ind = -1;

  for (int i = 0; i < popSize; i++)
  {
    if (a [i].f > fB)
    {
      fB = a [i].f;
      ind = i;
    }
  }

  if (ind != -1) ArrayCopy (cB, a [ind].c, 0, 0, WHOLE_ARRAY);

  //----------------------------------------------------------------------------
  for (int i = 0; i < popSize; i++)
  {
    if (a [i].f > a [i].fB)
    {
      a [i].fB = a [i].f;
      ArrayCopy (a [i].cB, a [i].c, 0, 0, WHOLE_ARRAY);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————