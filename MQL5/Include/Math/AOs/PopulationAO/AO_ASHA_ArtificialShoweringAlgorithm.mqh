//+————————————————————————————————————————————————————————————————————————————+
//|                                                                  C_AO_ASHA |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/15980

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
class C_AO_ASHA : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_ASHA () { }
  C_AO_ASHA ()
  {
    ao_name = "ASHA";
    ao_desc = "Artificial Showering Algorithm";
    ao_link = "https://www.mql5.com/ru/articles/15980";

    popSize       = 100;  //population size

    F             = 0.3;  //water flow velocity
    δ             = 2;    //resistance level(infiltration threshold)
    β             = 0.8;  //parameter that controls the rate of change in probability
    ρ0            = 0.1;  //initial probability

    ArrayResize (params, 5);

    params [0].name = "popSize"; params [0].val = popSize;
    params [1].name = "F";       params [1].val = F;
    params [2].name = "δ";       params [2].val = δ;
    params [3].name = "β";       params [3].val = β;
    params [4].name = "ρ0";      params [4].val = ρ0;

  }

  void SetParams ()
  {
    popSize = (int)params [0].val;
    F       = params      [1].val;
    δ       = (int)params [2].val;
    β       = params      [3].val;
    ρ0      = params      [4].val;
  }

  bool Init (const double &rangeMinP  [], //minimum search range
             const double &rangeMaxP  [], //maximum search range
             const double &rangeStepP [], //step search
             const int     epochsP = 0);  //number of epochs

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  double F;  //water flow velocity
  int    δ;  //resistance level(infiltration threshold)
  double β;  //parameter that controls the rate of change in probability
  double ρ0; //initial probability

  private: //-------------------------------------------------------------------
  S_AO_Agent aT [];
  int  epochs;
  int  epochNow;
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_ASHA::Init (const double &rangeMinP  [],
                      const double &rangeMaxP  [],
                      const double &rangeStepP [],
                      const int     epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  epochs   = epochsP;
  epochNow = 0;

  ArrayResize (aT, popSize);

  return true;
}
//——————————————————————————————————————————————————————————————————————————————


//——————————————————————————————————————————————————————————————————————————————
void C_AO_ASHA::Moving ()
{
  epochNow++;

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
  double xOld    = 0.0;
  double xNew    = 0.0;
  double xLower  = 0.0;
  double xLowest = 0.0;
  double ρ       = MathMax (β * (epochs - epochNow) / epochs, ρ0);
  double inf     = 0.0;
  int    ind     = 0;
  double rnd     = 0.0;

  for (int i = 0; i < popSize; i++)
  {
    inf = u.Scale (a [i].cnt, 0, δ, 0, 1);
    inf = inf * inf * inf * inf;

    rnd = u.RNDprobab ();

    for (int c = 0; c < coords; c++)
    {
      ind = (int)u.RNDintInRange (0, i - 1);
      
      if (i < 1)
      {
        if (rnd < inf)
        {
          a [i].c [c] = u.GaussDistribution (cB [c], rangeMin [c], rangeMax [c], 8);
        }
      }
      else
      {
        if (rnd < inf)
        {
          a [i].c [c] = u.GaussDistribution (a [ind].cB [c], rangeMin [c], rangeMax [c], 8);
        }
        else
        {
          xOld = a [i].c [c];

          if (u.RNDprobab () < ρ)
          {
            xLower = a [ind].cB [c];

            xNew = xOld + F * (u.RNDprobab () * (xLower - xOld));
          }
          else
          {
            xLowest = cB [c];

            xNew = xOld + F * (u.RNDprobab () * (xLowest - xOld));
          }

          a [i].c [c] = xNew;
        }
      }

      a [i].c [c] = u.SeInDiSp  (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_ASHA::Revision ()
{
  //----------------------------------------------------------------------------
  int ind = -1;

  for (int i = 0; i < popSize; i++)
  {
    if (a [i].f > fB)
    {
      fB = a [i].f;
      ind = i;
    }

    if (a [i].f > a [i].fB)
    {
      a [i].fB = a [i].f;
      ArrayCopy (a [i].cB, a [i].c, 0, 0, WHOLE_ARRAY);
      a [i].cnt = 0;
    }
    else
    {
      a [i].cnt++;
    }
  }

  if (ind != -1) ArrayCopy (cB, a [ind].c, 0, 0, WHOLE_ARRAY);

  //----------------------------------------------------------------------------
  u.Sorting_fB (a, aT, popSize);
}
//——————————————————————————————————————————————————————————————————————————————