//+————————————————————————————————————————————————————————————————————————————+
//|                                                                  C_AO_ABHA |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/15511

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
struct S_ASO_Member
{
    double pPrev [];     // Previous position
    double pBest  [];    // Personal best position
    double pBestFitness; // Personal best fitness


    void Init (int coords)
    {
      ArrayResize (pBest, coords);
      ArrayResize (pPrev, coords);
      pBestFitness = -DBL_MAX;
    }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
class C_AO_ASO : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_ASO () { }
  C_AO_ASO ()
  {
    ao_name = "ASO";
    ao_desc = "Anarchy Society Optimization";
    ao_link = "https://www.mql5.com/ru/articles/15511";

    popSize     = 50;     // Population size
    anarchyProb = 0.01;    // Probability of anarchic behavior

    omega       = 0.7;    // Inertia weight
    lambda1     = 1.5;    // Acceleration coefficient for P-best
    lambda2     = 1.5;    // Acceleration coefficient for G-best

    alpha       = 0.5;    // Parameter for FI calculation
    theta       = 0.1;    // Parameter for EI calculation
    delta       = 0.1;    // Parameter for II calculation

    ArrayResize (params, 8);

    params [0].name = "popSize";     params [0].val = popSize;
    params [1].name = "anarchyProb"; params [1].val = anarchyProb;

    params [2].name = "omega";       params [2].val = omega;
    params [3].name = "lambda1";     params [3].val = lambda1;
    params [4].name = "lambda2";     params [4].val = lambda2;

    params [5].name = "alpha";       params [5].val = alpha;
    params [6].name = "theta";       params [6].val = theta;
    params [7].name = "delta";       params [7].val = delta;
  }

  void SetParams ()
  {
    popSize     = (int)params [0].val;
    anarchyProb = params      [1].val;

    omega       = params      [2].val;
    lambda1     = params      [3].val;
    lambda2     = params      [4].val;

    alpha       = params      [5].val;
    theta       = params      [6].val;
    delta       = params      [7].val;
  }

  bool Init (const double &rangeMinP  [],
             const double &rangeMaxP  [],
             const double &rangeStepP [],
             const int     epochsP = 0);

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  double anarchyProb; // Probability of anarchic behavior
  double omega;       // Inertia weight
  double lambda1;     // Acceleration coefficient for P-best
  double lambda2;     // Acceleration coefficient for G-best
  double alpha;       // Parameter for FI calculation
  double theta;       // Parameter for EI calculation
  double delta;       // Parameter for II calculation

  S_ASO_Member member []; // Vector of society members

  private: //-------------------------------------------------------------------

  double CalculateFI (int memberIndex);
  double CalculateEI (int memberIndex);
  double CalculateII (int memberIndex);
  void   CurrentMP   (S_AO_Agent &agent, S_ASO_Member &memb, int coordInd);
  void   SocietyMP   (S_AO_Agent &agent, int coordInd);
  void   PastMP      (S_AO_Agent &agent, S_ASO_Member &memb, int coordInd);
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_ASO::Init (const double &rangeMinP  [],
                     const double &rangeMaxP  [],
                     const double &rangeStepP [],
                     const int     epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  ArrayResize (member, popSize);
  for (int i = 0; i < popSize; i++) member [i].Init (coords);

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_ASO::Moving ()
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

        member [i].pPrev [c] = a [i].c [c];
      }
    }

    revision = true;
    return;
  }

  //----------------------------------------------------------------------------
  double fi  = 0.0; //индекс недовольства
  double ei  = 0.0; //индекс внешней нерегулярности
  double ii  = 0.0; //индекс внутренней нерегулярности
  double rnd = 0.0;

  for (int i = 0; i < popSize; i++)
  {
    fi = CalculateFI (i);
    ei = CalculateEI (i);
    ii = CalculateII (i);

    for (int c = 0; c < coords; c++)
    {
      member [i].pPrev [c] = a [i].c [c];
      rnd = u.RNDprobab ();

      if (u.RNDprobab () < anarchyProb) a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
      else
      {
        if (rnd > fi) CurrentMP (a [i], member [i], c);
        else
        {
          if (rnd < ei) SocietyMP (a [i], c);
          else
          {
            if (rnd < ii) PastMP (a [i], member [i], c);
          }
        }
      }
    }

    for (int c = 0; c < coords; c++)
    {
      a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_ASO::Revision ()
{
  int ind = -1;

  for (int i = 0; i < popSize; i++)
  {
    if (a [i].f > fB)
    {
      fB = a [i].f;
      ind = i;
    }

    if (a [i].f > member [i].pBestFitness)
    {
      member [i].pBestFitness = a [i].f;
      ArrayCopy (member [i].pBest, a [i].c, 0, 0, WHOLE_ARRAY);
    }
  }

  if (ind != -1) ArrayCopy (cB, a [ind].c, 0, 0, WHOLE_ARRAY);
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_ASO::CalculateFI (int memberIndex)
{
  double currentFitness      = a      [memberIndex].f;
  double personalBestFitness = member [memberIndex].pBestFitness;
  double globalBestFitness   = fB;

  //1 - 0.9 * (800-x)/(1000-x)
  return 1 - alpha * (personalBestFitness - currentFitness) / (globalBestFitness - currentFitness);
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_ASO::CalculateEI (int memberIndex)
{
  double currentFitness    = a [memberIndex].f;
  double globalBestFitness = fB;

  //1-exp(-(10000-x)/(10000*0.9))
  return 1 - MathExp (-(globalBestFitness - currentFitness) / (globalBestFitness * theta));
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_ASO::CalculateII (int memberIndex)
{
  double currentFitness      = a      [memberIndex].f;
  double personalBestFitness = member [memberIndex].pBestFitness;

  //1-exp(-(10000-x)/(10000*0.9))
  return 1 - MathExp (-(personalBestFitness - currentFitness) / (personalBestFitness * delta));
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_ASO::CurrentMP (S_AO_Agent &agent, S_ASO_Member &memb, int coordInd)
{
  double r1 = u.RNDprobab ();
  double r2 = u.RNDprobab ();

  double velocity = omega   *      (agent.c    [coordInd] - memb.pBest [coordInd]) +
                    lambda1 * r1 * (memb.pBest [coordInd] - agent.c    [coordInd]) +
                    lambda2 * r2 * (cB         [coordInd] - agent.c    [coordInd]);

  agent.c [coordInd] += velocity;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_ASO::SocietyMP (S_AO_Agent &agent, int coordInd)
{
  int otherMember = u.RNDminusOne (popSize);

  agent.c [coordInd] = u.RNDprobab () < 0.5 ? cB [coordInd] : member [otherMember].pBest [coordInd];
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_ASO::PastMP (S_AO_Agent &agent, S_ASO_Member &memb, int coordInd)
{
  agent.c [coordInd] = u.RNDprobab () < 0.5 ? memb.pBest [coordInd] :
                                              //memb.pPrev [coordInd];
                                              u.GaussDistribution (agent.c [coordInd], rangeMin [coordInd], rangeMax [coordInd], 8);
}
//——————————————————————————————————————————————————————————————————————————————