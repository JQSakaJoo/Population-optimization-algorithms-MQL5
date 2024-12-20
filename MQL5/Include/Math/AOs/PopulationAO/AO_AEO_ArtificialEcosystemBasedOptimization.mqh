//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_AEO |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/16058

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
class C_AO_AEO : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_AEO () { }
  C_AO_AEO ()
  {
    ao_name = "AEO";
    ao_desc = "Artificial Ecosystem-based Optimization Algorithm";
    ao_link = "https://www.mql5.com/ru/articles/16058";


    popSize    = 50;       // population size
    levisPower = 10.0;

    ArrayResize (params, 2);

    params [0].name = "popSize";    params [0].val = popSize;
    params [1].name = "levisPower"; params [1].val = levisPower;
  }

  void SetParams ()
  {
    popSize    = (int)params [0].val;
    levisPower = params      [1].val;
  }

  bool Init (const double &rangeMinP  [],  // minimum search range
             const double &rangeMaxP  [],  // maximum search range
             const double &rangeStepP [],  // step search
             const int     epochsP = 0);   // number of epochs

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  double levisPower;


  private: //-------------------------------------------------------------------
  int  epochs;
  int  epochNow;
  int  consModel; // consumption model;
  S_AO_Agent aT [];
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_AEO::Init (const double &rangeMinP [],
                     const double &rangeMaxP [],
                     const double &rangeStepP [],
                     const int epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  epochs    = epochsP;
  epochNow  = 0;
  consModel = 0;
  ArrayResize (aT, popSize);

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_AEO::Moving ()
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
        a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }
    revision = true;
    return;
  }

  //----------------------------------------------------------------------------
  double α = (1.0 - (double)epochNow / epochs);

  double Xi   = 0.0;
  double Xb   = 0.0;
  double Xr   = 0.0;
  double Xj   = 0.0;
  double C    = 0.0;
  int    j    = 0;
  double r    = 0.0;

  // Production ---------------------------------------------------------------- Производство
  // X(t + 1) = (1 - α) * Xb(t) + α * Xrnd (t)
  // α = (1 - t / T) * rnd [0.0; 1.0]
  // Xrnd = rnd [Xmin; Xmax]
  if (consModel == 0)
  {
    for (int i = 0; i < popSize; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        Xb = cB [c];
        Xr = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        //Xi = Xb + α * (Xr - Xb);
        Xi = Xb + α * (Xb - Xr);

        a [i].c [c] = u.SeInDiSp (Xi, rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }

    consModel++;
    return;
  }

  // Consumption --------------------------------------------------------------- Потребление
  if (consModel == 1)
  {
    for (int i = 0; i < popSize; i++)
    {
      if (i > 1)
      {
        for (int c = 0; c < coords; c++)
        {
          r = u.RNDprobab ();

          // Herbivore behavior ------------------------------------------------ Травоядный
          //Xi (t + 1) = Xi (t) + C * (Xi (t) - Xb (t));
          if (r < 0.333)
          {
            C  = u.LevyFlightDistribution (levisPower);
            Xb = cB [c];
            //Xi = a [i].c [c];
            Xi = a [i].cB [c];

            //Xi = Xi + C * (Xi - Xb);
            Xi = Xi + C * (Xb - Xi);
          }
          else
          {
            // Carnivore behavior ---------------------------------------------- Плотоядный
            //Xi (t + 1) = Xi (t) + C * (Xi (t) - Xj (t));
            if (r < 0.667)
            {
              C  = u.LevyFlightDistribution (levisPower);
              j  = u.RNDminusOne (i);
              //Xj = a [j].c [c];
              Xj = a [j].cB [c];
              //Xi = a [i].c [c];
              Xi = a [i].cB [c];

              //Xi = Xi + C * (Xi - Xj);
              Xi = Xi + C * (Xj - Xi);
            }
            // Omnivore behavior ----------------------------------------------- Всеядный
            //Xi (t + 1) = Xi (t) + C * r2 * (Xi (t) - Xb (t)) +
            //                    (1 - r2) * (Xi (t) - Xj (t));
            else
            {
              C  = u.LevyFlightDistribution (levisPower);
              Xb = cB [c];
              j  = u.RNDminusOne (i);
              //Xj = a [j].c [c];
              Xj = a [j].cB [c];
              //Xi = a [i].c [c];
              Xi = a [i].cB [c];
              r = u.RNDprobab ();

              //Xi = Xi + C * r * (Xi - Xb) +
              //     (1 - r) * (Xi - Xj);
              Xi = Xi + C * r * (Xb - Xi) +
                   (1 - r) * (Xj - Xi);
            }
          }

          a [i].c [c] = u.SeInDiSp (Xi, rangeMin [c], rangeMax [c], rangeStep [c]);
        }
      }
    }

    consModel++;
    return;
  }

  // Decomposition -------------------------------------------------------------
  double D = 0.0;
  double h = 0.0;

  for (int i = 0; i < popSize; i++)
  {
    D = 3 * u.RNDprobab ();
    h = u.RNDprobab () * (u.RNDprobab () < 0.5 ? 1 : -1);
    C = u.LevyFlightDistribution (levisPower);
    j = u.RNDminusOne (popSize);

    for (int c = 0; c < coords; c++)
    {
      double x = a [i].cB [c] + D * (C * a [i].cB [c] - h * a [j].c [c]);
      a [i].c [c] = u.SeInDiSp (x, rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }

  consModel = 0;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_AEO::Revision ()
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

  u.Sorting_fB (a, aT, popSize);
}
//——————————————————————————————————————————————————————————————————————————————
