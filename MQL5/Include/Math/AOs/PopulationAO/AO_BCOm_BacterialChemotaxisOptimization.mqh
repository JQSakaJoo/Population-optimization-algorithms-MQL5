//+————————————————————————————————————————————————————————————————————————————+
//|                                                                  C_AO_BCOm |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/15711

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
struct S_BCO_Bacterium
{
    double fPrev;       // предыдущее значение целевой функции
    double fHistory []; // история значений целевой функции

    void Init (int coords, int historySize)
    {
      ArrayResize     (fHistory, historySize);
      ArrayInitialize (fHistory, 0.0);
      fPrev = -DBL_MAX;
    }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
class C_AO_BCOm : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_BCOm () { }
  C_AO_BCOm ()
  {
    ao_name = "BCOm";
    ao_desc = "Bacterial Chemotaxis Optimization M";
    ao_link = "https://www.mql5.com/ru/articles/15711";

    popSize = 50;     // размер популяции (количество бактерий)
    hs      = 10;     // количество шагов для вычисления среднего изменения

    ArrayResize (params, 2);

    params [0].name = "popSize";       params [0].val = popSize;
    params [1].name = "hs";            params [1].val = hs;
  }

  void SetParams ()
  {
    popSize       = (int)params [0].val;
    hs            = (int)params [1].val;
  }

  bool Init (const double &rangeMinP  [], //minimum search range
             const double &rangeMaxP  [], //maximum search range
             const double &rangeStepP [], //step search
             const int     epochsP = 0);  //number of epochs

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  int hs;

  S_BCO_Bacterium bacterium [];

  private: //-------------------------------------------------------------------
  double allowedDispl [];  //допустимые смещения бактерии

  double CalculateAverageDeltaFpr (int bacteriumIndex);
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_BCOm::Init (const double &rangeMinP  [], //minimum search range
                     const double &rangeMaxP  [], //maximum search range
                     const double &rangeStepP [], //step search
                     const int     epochsP = 0)   //number of epochs
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  ArrayResize (bacterium, popSize);
  for (int i = 0; i < popSize; i++) bacterium [i].Init (coords, hs);

  ArrayResize (allowedDispl, coords);
  for (int c = 0; c < coords; c++) allowedDispl [c] = rangeMax [c] - rangeMin [c];

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BCOm::Moving ()
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
  double x    = 0.0;
  double xMin = 0.0;
  double xMax = 0.0;
  double Δ    = 0.0;
  double ΔAve = 0.0;
  double dist = 0.0;

  for (int i = 0; i < popSize; i++)
  {
    ΔAve = CalculateAverageDeltaFpr (i) + DBL_EPSILON;

    Δ = fabs (a [i].f - bacterium [i].fPrev) / ΔAve;

    Δ = 1.0 - Δ;

    if (Δ < 0.0001) Δ = 0.0001;

    for (int c = 0; c < coords; c++)
    {
      if (u.RNDprobab () < 0.5)
      {
        dist = allowedDispl [c] * Δ;

        x    = a [i].c [c];
        xMin = x - dist;
        xMax = x + dist;

        x = u.GaussDistribution (x, xMin, xMax, 8);

        if (x > rangeMax [c]) x = u.RNDfromCI (xMin, rangeMax [c]);
        if (x < rangeMin [c]) x = u.RNDfromCI (rangeMin [c], xMax);

        a [i].c [c] = u.SeInDiSp (x, rangeMin [c], rangeMax [c], rangeStep [c]);
      }
      else a [i].c [c] = cB [c];
    }

    bacterium [i].fPrev = a [i].f;
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BCOm::Revision ()
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

  if (ind != -1)
  {
    ArrayCopy (cB, a [ind].c, 0, 0, WHOLE_ARRAY);
  }

  //----------------------------------------------------------------------------
  // Обновление истории значений целевой функции
  for (int i = 0; i < popSize; i++)
  {
    for (int j = 1; j < hs; j++)
    {
      bacterium [i].fHistory [j - 1] = bacterium [i].fHistory [j];
    }

    bacterium [i].fHistory [hs - 1] = a [i].f;
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_BCOm::CalculateAverageDeltaFpr (int bacteriumIndex)
{
  double sum = 0;

  for (int i = 1; i < hs; i++)
  {
    sum += bacterium [bacteriumIndex].fHistory [i] - bacterium [bacteriumIndex].fHistory [i - 1];
  }

  return sum / (hs - 1);
}
//——————————————————————————————————————————————————————————————————————————————

/*
#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
struct S_BCO_Bacterium
{
    double fPrev;        // previous value of the objective function
    double fHistory [];  // history of objective function values

    void Init (int coords, int historySize)
    {
      ArrayResize     (fHistory, historySize);
      ArrayInitialize (fHistory, 0.0);
      fPrev = -DBL_MAX;
    }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
class C_AO_BCO : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_BCO () { }
  C_AO_BCO ()
  {
    ao_name = "BCO";
    ao_desc = "Bacterial Chemotaxis Optimization";
    ao_link = "https://www.mql5.com/ru/articles/15711";

    popSize = 50;     // population size (number of bacteria)
    hs      = 10;     // number of steps to calculate average change
    T0      = 1.0;    // initial temperature
    b       = 0.5;    // parameter b
    tau_c   = 1.0;    // parameter tau_c
    v       = 1.0;    // velocity

    ArrayResize (params, 6);

    params [0].name = "popSize"; params [0].val = popSize;
    params [1].name = "hs";      params [1].val = hs;
    params [2].name = "T0";      params [2].val = T0;
    params [3].name = "b";       params [3].val = b;
    params [4].name = "tau_c";   params [4].val = tau_c;
    params [5].name = "v";       params [5].val = v;
  }

  void SetParams ()
  {
    popSize = (int)params [0].val;
    hs      = (int)params [1].val;
    T0      = params      [2].val;
    b       = params      [3].val;
    tau_c   = params      [4].val;
    v       = params      [5].val;
  }

  bool Init (const double &rangeMinP  [], //minimum search range
             const double &rangeMaxP  [], //maximum search range
             const double &rangeStepP [], //step search
             const int     epochsP = 0);  //number of epochs

  void Moving ();
  void Revision ();

  //----------------------------------------------------------------------------
  int    hs;
  double T0;
  double b;
  double tau_c;
  double v;

  S_BCO_Bacterium bacterium [];

  private: //-------------------------------------------------------------------
  double CalculateAverageDeltaFpr            (int bacteriumIndex);
  void   CalculateNewAngles                  (double &angles []);
  void   CalculateNewPosition                (double l, const double &angles [], double &new_position [], const double &current_position []);
  double GenerateFromExponentialDistribution (double T);
  double CalculateCorrectedB                 (double f_tr, double f_pr, int bacteriumIndex);
  double CalculateRotationAngle              ();
  int    RNDdir                              ();
};

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_BCO::Init (const double &rangeMinP  [], //minimum search range
                     const double &rangeMaxP  [], //maximum search range
                     const double &rangeStepP [], //step search
                     const int     epochsP = 0)   //number of epochs
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  ArrayResize (bacterium, popSize);
  for (int i = 0; i < popSize; i++) bacterium [i].Init (coords, hs);

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BCO::Moving ()
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
    double f_tr = a [i].f;
    double f_pr = bacterium [i].fPrev;

    double T;

    if (f_tr <= f_pr)
    {
      T = T0;
    }
    else
    {
      double b_corr = CalculateCorrectedB (f_tr, f_pr, i);

      T = T0 * exp (b_corr * (f_tr - f_pr));
    }

    double tau = GenerateFromExponentialDistribution (T);

    double new_angles [];
    ArrayResize (new_angles, coords - 1);
    CalculateNewAngles (new_angles);

    double l = v * tau;
    double new_position [];
    ArrayResize (new_position, coords);
    CalculateNewPosition (l, new_angles, new_position, a [i].c);

    for (int c = 0; c < coords; c++)
    {
      a [i].c [c] = u.SeInDiSp (new_position [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }

    bacterium [i].fPrev = a [i].f;
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BCO::Revision ()
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

  if (ind != -1)
  {
    ArrayCopy (cB, a [ind].c, 0, 0, WHOLE_ARRAY);
  }

  for (int i = 0; i < popSize; i++)
  {
    for (int j = 1; j < hs; j++)
    {
      bacterium [i].fHistory [j - 1] = bacterium [i].fHistory [j];
    }

    bacterium [i].fHistory [hs - 1] = a [i].f;
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_BCO::CalculateAverageDeltaFpr (int bacteriumIndex)
{
  double sum = 0;

  for (int i = 1; i < hs; i++)
  {
    sum += bacterium [bacteriumIndex].fHistory [i] - bacterium [bacteriumIndex].fHistory [i - 1];
  }

  return sum / (hs - 1);
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BCO::CalculateNewAngles (double &angles [])
{
  for (int i = 0; i < coords - 1; i++)
  {
    double theta = CalculateRotationAngle ();
    double mu    = 62 * (1 - MathCos (theta)) * M_PI / 180.0;
    double sigma = 26 * (1 - MathCos (theta)) * M_PI / 180.0;

    angles [i] = u.GaussDistribution (mu, mu - M_PI, mu + M_PI, 8);
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BCO::CalculateNewPosition (double l, const double &angles [], double &new_position [], const double &current_position [])
{
  new_position [0] = current_position [0] + l * (rangeMax [0] - rangeMin [0]);

  for (int k = 0; k < coords - 1; k++)
  {
    new_position [0] *= MathCos (angles [k]);
  }

  new_position [0] *= RNDdir ();

  for (int i = 1; i < coords - 1; i++)
  {
    new_position [i] = current_position [i] + l * MathSin (angles [i - 1]) * (rangeMax [0] - rangeMin [0]);

    for (int k = i; k < coords - 1; k++)
    {
      new_position [i] *= MathCos (angles [k]);
    }

    new_position [i] *= RNDdir ();
  }

  if (coords > 1)
  {
    new_position [coords - 1] = current_position [coords - 1] + l * MathSin (angles [coords - 2]);
  }

}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
int C_AO_BCO::RNDdir ()
{
  if (u.RNDbool () < 0.5) return -1;

  return 1;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_BCO::GenerateFromExponentialDistribution (double T)
{
  return -T * MathLog (u.RNDprobab ());
}

double C_AO_BCO::CalculateCorrectedB (double f_tr, double f_pr, int bacteriumIndex)
{
  double delta_f_tr = f_tr - f_pr;
  double delta_f_pr = CalculateAverageDeltaFpr (bacteriumIndex);

  if (MathAbs (delta_f_pr) < DBL_EPSILON)
  {
    return b;
  }
  else
  {
    return b * (1 / (MathAbs (delta_f_tr / delta_f_pr) + 1) + 1 / (MathAbs (f_pr) + 1));
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_BCO::CalculateRotationAngle ()
{
  double cos_theta = MathExp (-tau_c / T0);
  return MathArccos (cos_theta);
}
//——————————————————————————————————————————————————————————————————————————————
*/