//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_BSA |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/14491

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
struct S_BSA_Agent
{
    double cBest []; //best coordinates
    double fBest;    //best fitness

    void Init (int coords)
    {
      ArrayResize (cBest, coords);
      fBest = -DBL_MAX;
    }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
class C_AO_BSA : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_BSA () { }
  C_AO_BSA ()
  {
    ao_name = "BSA";
    ao_desc = "Bird Swarm Algorithm";
    ao_link = "https://www.mql5.com/ru/articles/14491";

    popSize        = 20;  //population size

    flyingProb     = 0.8;  //Flight probability
    producerProb   = 0.25; //Producer probability
    foragingProb   = 0.55; //Foraging probability
    a1             = 0.6;  //a1 constant [0...2]
    a2             = 0.05; //a2 constant [0...2]
    C              = 0.05; //Cognitive coefficient
    S              = 1.1;  //Social coefficient
    FL             = 1.75; //FL constant [0...2]
    producerPower  = 7.05; //Producer power
    scroungerPower = 2.60; //Scrounger power

    ArrayResize (params, 11);

    params [0].name = "popSize";         params [0].val  = popSize;

    params [1].name  = "flyingProb";     params [1].val  = flyingProb;
    params [2].name  = "producerProb";   params [2].val  = producerProb;
    params [3].name  = "foragingProb";   params [3].val  = foragingProb;
    params [4].name  = "a1";             params [4].val  = a1;
    params [5].name  = "a2";             params [5].val  = a2;
    params [6].name  = "C";              params [6].val  = C;
    params [7].name  = "S";              params [7].val  = S;
    params [8].name  = "FL";             params [8].val  = FL;
    params [9].name  = "producerPower";  params [9].val  = producerPower;
    params [10].name = "scroungerPower"; params [10].val = scroungerPower;
  }

  void SetParams ()
  {
    popSize        = (int)params [0].val;

    flyingProb     = params [1].val;
    producerProb   = params [2].val;
    foragingProb   = params [3].val;
    a1             = params [4].val;
    a2             = params [5].val;
    C              = params [6].val;
    S              = params [7].val;
    FL             = params [8].val;
    producerPower  = params [9].val;
    scroungerPower = params [10].val;
  }

  bool Init (const double &rangeMinP  [], //minimum search range
             const double &rangeMaxP  [], //maximum search range
             const double &rangeStepP [], //step search
             const int     epochsP = 0);  //number of epochs

  void Moving    ();
  void Revision  ();
  void Injection (const int popPos, const int coordPos, const double value);

  //----------------------------------------------------------------------------
  double flyingProb;      //Flight probability
  double producerProb;    //Producer probability
  double foragingProb;    //Foraging probability
  double a1;              //a1 constant [0...2]
  double a2;              //a2 constant [0...2]
  double C;               //Cognitive coefficient
  double S;               //Social coefficient
  double FL;              //FL constant [0...2]
  double producerPower;   //Producer power
  double scroungerPower;  //Scrounger power

  S_BSA_Agent agent [];

  private: //-------------------------------------------------------------------
  double mean [];  //represents the element of the average position of the whole bird’s swarm
  double N;
  double e;        //epsilon

  void BirdProducer  (int pos);
  void BirdScrounger (int pos);
  void BirdForaging  (int pos);
  void BirdVigilance (int pos);
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_BSA::Init (const double &rangeMinP  [], //minimum search range
                     const double &rangeMaxP  [], //maximum search range
                     const double &rangeStepP [], //step search
                     const int     epochsP = 0)   //number of epochs
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  ArrayResize (agent, popSize);
  for (int i = 0; i < popSize; i++) agent [i].Init (coords);

  ArrayResize (mean, coords);

  N = popSize;
  e = DBL_MIN;

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BSA::Moving ()
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
  for (int i = 0; i < popSize; i++)
  {
    //bird is flying------------------------------------------------------------
    if (u.RNDprobab () < flyingProb)
    {
      //bird producer
      if (u.RNDprobab () < producerProb) BirdProducer  (i); //bird is looking for a new place to eat
      //bird is not a producer
      else                               BirdScrounger (i); //scrounger follows the  producer
    }
    //bird is not flying--------------------------------------------------------
    else
    {
      //bird foraging
      if (u.RNDprobab () < foragingProb) BirdForaging  (i); //bird feeds
      //bird is not foraging
      else                               BirdVigilance (i); //bird vigilance
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BSA::Revision ()
{
  //----------------------------------------------------------------------------
  int ind = -1;

  for (int i = 0; i < popSize; i++)
  {
    if (a [i].f > fB) ind = i;
  }

  if (ind != -1)
  {
    fB = a [ind].f;
    ArrayCopy (cB, a [ind].c, 0, 0, WHOLE_ARRAY);
  }

  //----------------------------------------------------------------------------
  for (int i = 0; i < popSize; i++)
  {
    if (a [i].f > agent [i].fBest)
    {
      agent [i].fBest = a [i].f;
      ArrayCopy (agent [i].cBest, a [i].c, 0, 0, WHOLE_ARRAY);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void  C_AO_BSA::BirdProducer  (int pos)
{
  double x = 0.0; //bird position

  for (int c = 0; c < coords; c++)
  {
    x = a [pos].c [c];
    x = u.GaussDistribution (x, rangeMin [c], rangeMax [c], producerPower);

    a [pos].c [c] = u.SeInDiSp (x, rangeMin [c], rangeMax [c], rangeStep [c]);
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void  C_AO_BSA::BirdScrounger (int pos)
{
  int    K  = 0;   //position of a randomly selected bird in a swarm
  double x  = 0.0; //best bird position
  double xK = 0.0; //current best position of a randomly selected bird in a swarm

  for (int c = 0; c < coords; c++)
  {
    do K = u.RNDminusOne (popSize);
    while (K == pos);

    x  = agent [pos].cBest [c];
    xK = agent   [K].cBest [c];

    x = x + (xK - x) * FL * u.GaussDistribution (0, -1.0, 1.0, scroungerPower);

    a [pos].c [c] = u.SeInDiSp (x, rangeMin [c], rangeMax [c], rangeStep [c]);
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void  C_AO_BSA::BirdForaging  (int pos)
{
  double x  = 0.0; //current bird position
  double p  = 0.0; //best bird position
  double g  = 0.0; //best global position
  double r1 = 0.0; //uniform random number [0.0 ... 1.0]
  double r2 = 0.0; //uniform random number [0.0 ... 1.0]

  for (int c = 0; c < coords; c++)
  {
    x = a     [pos].c     [c];
    p = agent [pos].cBest [c];
    g = cB                [c];

    r1 = u.RNDprobab ();
    r2 = u.RNDprobab ();

    x = x + (p - x) * C * r1 + (g - x) * S * r2;

    a [pos].c [c] = u.SeInDiSp (x, rangeMin [c], rangeMax [c], rangeStep [c]);
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void  C_AO_BSA::BirdVigilance (int pos)
{
  int    K      = 0;   //position of a randomly selected bird in a swarm
  double sumFit = 0.0; //best birds fitness sum
  double pFitK  = 0.0; //best fitness of a randomly selected bird
  double pFit   = 0.0; //best bird fitness
  double A1     = 0.0;
  double A2     = 0.0;
  double r1     = 0.0; //uniform random number [ 0.0 ... 1.0]
  double r2     = 0.0; //uniform random number [-1.0 ... 1.0]
  double x      = 0.0; //best bird position
  double xK     = 0.0; //best position of a randomly selected bird in a swarm

  ArrayInitialize (mean, 0.0);

  for (int i = 0; i < popSize; i++) sumFit += agent [i].fBest;

  for (int c = 0; c < coords; c++)
  {
    for (int i = 0; i < popSize; i++) mean [c] += a [i].c [c];

    mean [c] /= popSize;
  }

  do K = u.RNDminusOne (popSize);
  while (K == pos);

  pFit  = agent [pos].fBest;
  pFitK = agent   [K].fBest;

  A1 = a1 * exp (-pFit * N / (sumFit + e));
  A2 = a2 * exp (((pFit - pFitK) / (fabs (pFitK - pFit) + e)) * (N * pFitK / (sumFit + e)));

  for (int c = 0; c < coords; c++)
  {
    r1 = u.RNDprobab ();
    r2 = u.RNDfromCI (-1, 1);

    x  = agent [pos].cBest [c];
    xK = agent   [K].cBest [c];

    x = x + A1 * (mean [c] - x) * r1 + A2 * (xK - x) * r2;

    a [pos].c [c] = u.SeInDiSp (x, rangeMin [c], rangeMax [c], rangeStep [c]);
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BSA::Injection (const int popPos, const int coordPos, const double value)
{
  if (popPos   < 0 || popPos   >= popSize) return;
  if (coordPos < 0 || coordPos >= coords) return;

  if (value < rangeMin [coordPos])
  {
    a [popPos].c [coordPos] = rangeMin [coordPos];
  }

  if (value > rangeMax [coordPos])
  {
    a [popPos].c [coordPos] = rangeMax [coordPos];
  }

  a [popPos].c [coordPos] = u.SeInDiSp (value, rangeMin [coordPos], rangeMax [coordPos], rangeStep [coordPos]);
}
//——————————————————————————————————————————————————————————————————————————————