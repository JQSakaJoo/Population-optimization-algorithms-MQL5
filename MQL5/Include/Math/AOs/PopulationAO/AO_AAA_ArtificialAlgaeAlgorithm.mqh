//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_AAA |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/15565

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
struct S_AAA_Agent
{
    double energy;
    int    hunger;
    double size;
    double friction;

    void Init ()
    {
      energy   = 1.0;
      hunger   = 0;
      size     = 1.0;
      friction = 0.0;
    }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
class C_AO_AAA : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_AAA () { }
  C_AO_AAA ()
  {
    ao_name = "AAA";
    ao_desc = "Algae Adaptive Algorithm";
    ao_link = "https://www.mql5.com/ru/articles/15565";

    popSize                = 200;

    adaptationProbability  = 0.2;
    energyLoss             = 0.05;
    maxGrowthRate          = 0.1;
    halfSaturationConstant = 1.0;

    ArrayResize (params, 5);

    params [0].name = "popSize";                params [0].val = popSize;

    params [1].name = "adaptationProbability";  params [1].val = adaptationProbability;
    params [2].name = "energyLoss";             params [2].val = energyLoss;
    params [3].name = "maxGrowthRate";          params [3].val = maxGrowthRate;
    params [4].name = "halfSaturationConstant"; params [4].val = halfSaturationConstant;
  }

  void SetParams ()
  {
    popSize                = (int)params [0].val;
    adaptationProbability  = params      [1].val;
    energyLoss             = params      [2].val;
    maxGrowthRate          = params      [3].val;
    halfSaturationConstant = params      [4].val;
  }

  bool Init (const double &rangeMinP  [],
             const double &rangeMaxP  [],
             const double &rangeStepP [],
             const int     epochsP = 0);

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  double adaptationProbability;
  double energyLoss;
  double maxGrowthRate;
  double halfSaturationConstant;

  S_AAA_Agent agent [];

  private: //-------------------------------------------------------------------
  void   EvolutionProcess    ();
  void   AdaptationProcess   ();
  double CalculateEnergy     (int index);
  int    TournamentSelection ();
  double fMin, fMax;
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_AAA::Init (const double &rangeMinP  [],
                     const double &rangeMaxP  [],
                     const double &rangeStepP [],
                     const int     epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  ArrayResize (agent, popSize);

  fMin = -DBL_MAX;
  fMax =  DBL_MAX;

  for (int i = 0; i < popSize; i++)
  {
    agent [i].Init ();

    for (int c = 0; c < coords; c++)
    {
      a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
      a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_AAA::Moving ()
{
  //----------------------------------------------------------------------------
  if (!revision)
  {
    revision = true;
    return;
  }

  //----------------------------------------------------------------------------
  for (int i = 0; i < popSize; i++)
  {
    int variant = 0;

    int j = TournamentSelection ();

    for (int c = 0; c < coords; c++)
    {
      double α = u.RNDfromCI (0.0, 2 * M_PI);
      double β = u.RNDfromCI (0.0, 2 * M_PI);
      double ρ = u.RNDfromCI (-1.0, 1.0);

      if (variant == 0)
      {
        a [i].c [c] += (a [j].c [c] - a [i].c [c]) * agent [i].friction * MathCos (α);
      }
      if (variant == 1)
      {
        a [i].c [c] += (a [j].c [c] - a [i].c [c]) * agent [i].friction * MathSin (β);
      }
      if (variant == 2)
      {
        a [i].c [c] += (a [j].c [c] - a [i].c [c]) * agent [i].friction * ρ;
      }

      variant++;

      if (variant > 2) variant = 0;

      a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_AAA::Revision ()
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

    agent [i].size = a [i].f;

    if (a [i].f < fMin) fMin = a [i].f;
    if (a [i].f > fMax) fMax = a [i].f;
  }

  if (ind != -1) ArrayCopy (cB, a [ind].c, 0, 0, WHOLE_ARRAY);

  //----------------------------------------------------------------------------
  for (int i = 0; i < popSize; i++)
  {
    agent [i].energy   = CalculateEnergy (i);
    agent [i].friction = u.Scale (a [i].f, fMin, fMax, 0.1, 1.0, false);

    agent [i].energy -= energyLoss;

    double newEnergy = CalculateEnergy (i);

    if (newEnergy > agent [i].energy)
    {
      agent [i].energy += energyLoss / 2;
      agent [i].hunger = 0;
    }
    else
    {
      agent [i].hunger++;
    }

    double growthRate = maxGrowthRate * agent [i].size / (agent [i].size + halfSaturationConstant);
    agent [i].size *= (1 + growthRate);
  }

  //----------------------------------------------------------------------------
  EvolutionProcess  ();
  AdaptationProcess ();
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_AAA::EvolutionProcess ()
{
  int smallestIndex = 0;

  for (int i = 1; i < popSize; i++)
  {
    if (agent [i].size < agent [smallestIndex].size) smallestIndex = i;
  }

  int m = 0;

  for (int c = 0; c < coords; c++)
  {
    m = u.RNDminusOne (popSize);
    a [smallestIndex].c [c] = a [m].c [c];
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_AAA::AdaptationProcess ()
{
  int starvingIndex = 0;

  for (int i = 1; i < popSize; i++)
  {
    if (agent [i].hunger > agent [starvingIndex].hunger) starvingIndex = i;
  }

  if (u.RNDprobab () < adaptationProbability)
  {
    int biggestIndex = 0;

    for (int i = 1; i < popSize; i++)
    {
      if (agent [i].size > agent [biggestIndex].size) biggestIndex = i;
    }

    for (int j = 0; j < coords; j++)
    {
      a [starvingIndex].c [j] += (a [biggestIndex].c [j] - a [starvingIndex].c [j]) * u.RNDprobab ();
      a [starvingIndex].c [j] = u.SeInDiSp (a [starvingIndex].c [j], rangeMin [j], rangeMax [j], rangeStep [j]);
    }

    agent [starvingIndex].size = a [starvingIndex].f;
    agent [starvingIndex].hunger = 0;
    agent [starvingIndex].energy = 1.0;
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_AAA::CalculateEnergy (int index)
{
  double colony_size              = agent [index].size;
  double max_growth_rate          = maxGrowthRate;
  double half_saturation_constant = halfSaturationConstant;

  // Используем нормализованное значение фитнес-функции
  double nutrient_concentration = (a [index].f - fMin) / (fMax - fMin + 1e-10);

  double current_growth_rate = agent [index].energy;

  double growth_rate = max_growth_rate * nutrient_concentration / (half_saturation_constant + current_growth_rate) * colony_size;

  double energy = growth_rate - energyLoss;

  if (energy < 0) energy = 0;

  return energy;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
int C_AO_AAA::TournamentSelection ()
{
  int candidate1 = u.RNDminusOne (popSize);
  int candidate2 = u.RNDminusOne (popSize);

  return (a [candidate1].f > a [candidate2].f) ? candidate1 : candidate2;
}
//——————————————————————————————————————————————————————————————————————————————