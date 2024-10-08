//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_AAm |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/15782

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
struct S_AA_Agent
{
    double cPrev []; // previous coordinates
    double fPrev;    // previous fitness

    void Init (int coords)
    {
      ArrayResize (cPrev, coords);
      fPrev = -DBL_MAX;
    }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
class C_AO_AAm : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_AAm () { }
  C_AO_AAm ()
  {
    ao_name = "AAm";
    ao_desc = "Archery Algorithm M";
    ao_link = "https://www.mql5.com/ru/articles/15782";

    popSize   = 50;    // population size
    inhProbab = 0.3;

    ArrayResize (params, 2);

    params [0].name = "popSize";   params [0].val = popSize;
    params [1].name = "inhProbab"; params [1].val = inhProbab;
  }

  void SetParams ()
  {
    popSize   = (int)params [0].val;
    inhProbab = params      [1].val;
  }

  bool Init (const double &rangeMinP  [], // minimum search range
             const double &rangeMaxP  [], // maximum search range
             const double &rangeStepP [], // step search
             const int     epochsP = 0);  // number of epochs

  void Moving ();
  void Revision ();

  //----------------------------------------------------------------------------
  double  inhProbab; //probability of inheritance

  S_AA_Agent agent [];

  private: //-------------------------------------------------------------------
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_AAm::Init (const double &rangeMinP  [],
                     const double &rangeMaxP  [],
                     const double &rangeStepP [],
                     const int     epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  ArrayResize (agent, popSize);
  for (int i = 0; i < popSize; i++) agent [i].Init (coords);

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_AAm::Moving ()
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

  //-------------------------------------------------------------------------
  // Calculate probability vector P and cumulative probability C
  double P [], C [];
  ArrayResize (P, popSize);
  ArrayResize (C, popSize);
  double F_worst = DBL_MAX; // a[ArrayMaximum(a, WHOLE_ARRAY, 0, popSize)].f;
  double sum = 0;

  for (int i = 0; i < popSize; i++)
  {
    if (a [i].f < F_worst) F_worst = a [i].f;
  }

  for (int i = 0; i < popSize; i++)
  {
    P [i] =  a [i].f - F_worst; ////F_worst - a[i].f;
    sum += P [i];
  }

  for (int i = 0; i < popSize; i++)
  {
    P [i] /= sum;
    C [i] = (i == 0) ? P [i] : C [i - 1] + P [i];
  }

  double x;
  
  double maxFF = fB;
  double minFF = DBL_MAX;
  double prob1;
  double prob2;
  
  for (int i = 0; i < popSize; i++)
  {
    if (a [i].f < minFF) minFF = a [i].f;
  } 

  for (int i = 0; i < popSize; i++)
  {
    for (int c = 0; c < coords; c++)
    {
      // Select archer (k) using cumulative probability
      int k = 0;
      double r = u.RNDprobab ();
      while (k < popSize - 1 && C [k] < r) k++;

      if (u.RNDbool () < inhProbab)
      {
        x = a [k].c [c];
      }
      else
      {
        
        // Update position using Eq. (5) and (6)
        //double I   = MathRound (1 + u.RNDprobab ());
        double rnd = u.GaussDistribution (0, -1, 1, 8);
        /*
        if (a [k].f > a [i].f)
        {
          x = agent [i].cPrev [c] + rnd * (a [k].c [c] - I * agent [i].cPrev [c]);
        }
        else
        {
          x = agent [i].cPrev [c] + rnd * (agent [i].cPrev [c] - I * a [k].c [c]);
        }
        */
        
        prob1 = u.Scale (a [i].f, minFF, maxFF, 0, 1);
        prob2 = u.Scale (a [k].f, minFF, maxFF, 0, 1);
        
        x = agent [i].cPrev [c] + rnd * (a [k].c [c] - agent [i].cPrev [c]) * (1 - prob1 - prob2); 
        //x = agent [i].cPrev [c] + rnd * (a [k].c [c] - agent [i].cPrev [c]) * (2 - prob1 - prob2);  
      }

      a [i].c [c] = u.SeInDiSp (x, rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_AAm::Revision ()
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
    if (a [i].f > agent [i].fPrev)
    {
      agent [i].fPrev = a [i].f;
      ArrayCopy (agent [i].cPrev, a [i].c, 0, 0, WHOLE_ARRAY);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————
