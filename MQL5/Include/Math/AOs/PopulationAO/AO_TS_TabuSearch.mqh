//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_TSm |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/15654

#include "#C_AO.mqh"

/*
Каждая координата разбивается на указанное количество секторов равной длины.
Популяция: 50
Количество секторов на каждой координате: 100

1. Инициализировать популяцию случайными числами в диапазоне каждой координаты.
2. Вычислить приспособленность особи.
3. Проверить каждую особь:
   a) если приспособленность особи улучшилась, то прибавить счетчики в соответсвующем секторе в белом листе.
   b) если приспособленность особи ухудшилась, то прибавить счетчики в соответсвующем секторе в черном листе.
4. Для каждой особи сгенерировать новые координаты в соответсвующих секторах пропорционально из вероятности (на основе счетчиков по секторам в белом листе).
5. Для каждой особи проверить сгенерированные координаты по черному листу, вычислить вероятность,
   если она выпадает на сектор черного списка, то выбрать случайно сектор и сгенерировать новую координату.
6. Повторить с п.2
*/

//——————————————————————————————————————————————————————————————————————————————
struct S_TSmSector
{
    int sector [];
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
struct S_TSmAgent
{
    S_TSmSector blacklist []; //черный лист по секторам каждой координаты
    S_TSmSector whitelist []; //белый лист по секторам каждой координаты

    double fPrev;             //предыдущая приспособленность

    void Init (int coords, int sectorsPerCord)
    {
      ArrayResize (blacklist, coords);
      ArrayResize (whitelist, coords);

      for (int i = 0; i < coords; i++)
      {
        ArrayResize (blacklist [i].sector, sectorsPerCord);
        ArrayResize (whitelist [i].sector, sectorsPerCord);

        ArrayInitialize (blacklist [i].sector, 0);
        ArrayInitialize (whitelist [i].sector, 0);
      }

      fPrev = -DBL_MAX;
    }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
class C_AO_TSm : public C_AO
{
  public: //--------------------------------------------------------------------
  C_AO_TSm ()
  {
    ao_name = "TSm";
    ao_desc = "Tabu Search M";
    ao_link = "https://www.mql5.com/ru/articles/15654";

    popSize         = 50;
    sectorsPerCoord = 100;
    bestProbab      = 0.8;

    ArrayResize (params, 3);

    params [0].name = "popSize";         params [0].val = popSize;
    params [1].name = "sectorsPerCoord"; params [1].val = sectorsPerCoord;
    params [2].name = "bestProbab";      params [2].val = bestProbab;
  }

  void SetParams ()
  {
    popSize         = (int)params [0].val;
    sectorsPerCoord = (int)params [1].val;
    bestProbab      = params      [2].val;
  }

  bool Init (const double &rangeMinP  [], //minimum search range
             const double &rangeMaxP  [], //maximum search range
             const double &rangeStepP [], //step search
             const int     epochsP = 0);

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  int    sectorsPerCoord;
  double bestProbab;

  S_TSmAgent agents [];

  private: //-------------------------------------------------------------------
  void   InitializePopulation      ();
  void   UpdateLists               ();
  void   GenerateNewCoordinates    ();
  int    GetSectorIndex            (double coord, int dimension);
  int    ChooseSectorFromWhiteList (int agentIndex, int dimension);
  double GenerateCoordInSector     (int sectorIndex, int dimension);
  bool   IsInBlackList             (int agentIndex, int dimension, int sectorIndex);

};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_TSm::Init (const double &rangeMinP  [], //minimum search range
                     const double &rangeMaxP  [], //maximum search range
                     const double &rangeStepP [], //step search
                     const int     epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  ArrayResize (agents, popSize);

  for (int i = 0; i < popSize; i++) agents [i].Init (coords, sectorsPerCoord);

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_TSm::Moving ()
{
  //----------------------------------------------------------------------------
  if (!revision)
  {
    InitializePopulation ();
    revision = true;
    return;
  }

  //----------------------------------------------------------------------------
  GenerateNewCoordinates ();
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_TSm::Revision ()
{
  //----------------------------------------------------------------------------
  for (int i = 0; i < popSize; i++)
  {
    if (a [i].f > fB)
    {
      fB = a [i].f;
      ArrayCopy (cB, a [i].c);
    }
  }

  //----------------------------------------------------------------------------
  UpdateLists ();
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_TSm::InitializePopulation ()
{
  for (int i = 0; i < popSize; i++)
  {
    for (int c = 0; c < coords; c++)
    {
      a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
      a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
    agents [i].fPrev = -DBL_MAX;
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_TSm::UpdateLists ()
{
  for (int i = 0; i < popSize; i++)
  {
    for (int c = 0; c < coords; c++)
    {
      int sectorIndex = GetSectorIndex (a [i].c [c], c);

      if (a [i].f > agents [i].fPrev)
      {
        agents [i].whitelist [c].sector [sectorIndex]++;
      }
      else
        if (a [i].f < agents [i].fPrev)
        {
          agents [i].blacklist [c].sector [sectorIndex]++;
        }
    }
    agents [i].fPrev = a [i].f;
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_TSm::GenerateNewCoordinates ()
{
  for (int i = 0; i < popSize; i++)
  {
    for (int c = 0; c < coords; c++)
    {
      if (u.RNDprobab () < bestProbab)
      {
        a [i].c [c] = cB [c];
      }
      else
      {
        int sectorIndex = ChooseSectorFromWhiteList (i, c);
        double newCoord = GenerateCoordInSector (sectorIndex, c);

        if (IsInBlackList (i, c, sectorIndex))
        {
          sectorIndex = u.RNDminusOne (sectorsPerCoord);
          newCoord = GenerateCoordInSector (sectorIndex, c);
        }

        newCoord = u.SeInDiSp (newCoord, rangeMin [c], rangeMax [c], rangeStep [c]);

        a [i].c [c] = newCoord;
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
int C_AO_TSm::GetSectorIndex (double coord, int dimension)
{
  if (rangeMax [dimension] == rangeMin [dimension]) return 0;

  double sL =  (rangeMax [dimension] - rangeMin [dimension]) / sectorsPerCoord;

  int ind = (int)MathFloor ((coord - rangeMin [dimension]) / sL);

  // Особая обработка для максимального значения
  if (coord == rangeMax [dimension]) return sectorsPerCoord - 1;

  if (ind >= sectorsPerCoord) return sectorsPerCoord - 1;
  if (ind < 0) return 0;

  return ind;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
int C_AO_TSm::ChooseSectorFromWhiteList (int agentIndex, int dimension)
{
  int totalCount = 0;

  for (int s = 0; s < sectorsPerCoord; s++)
  {
    totalCount += agents [agentIndex].whitelist [dimension].sector [s];
  }

  if (totalCount == 0)
  {
    int randomSector = u.RNDminusOne (sectorsPerCoord);
    return randomSector;
  }

  int randomValue = u.RNDminusOne (totalCount);
  int cumulativeCount = 0;

  for (int s = 0; s < sectorsPerCoord; s++)
  {
    cumulativeCount += agents [agentIndex].whitelist [dimension].sector [s];

    if (randomValue <= cumulativeCount)
    {
      return s;
    }
  }

  return sectorsPerCoord - 1;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_TSm::GenerateCoordInSector (int sectorIndex, int dimension)
{
  double sectorSize  = (rangeMax [dimension] - rangeMin [dimension]) / sectorsPerCoord;
  double sectorStart = rangeMin [dimension] + sectorIndex * sectorSize;
  double sectorEnd   = sectorStart + sectorSize;

  double newCoord = u.RNDfromCI (sectorStart, sectorEnd);

  return newCoord;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_TSm::IsInBlackList (int agentIndex, int dimension, int sectorIndex)
{
  int blackCount = agents [agentIndex].blacklist [dimension].sector [sectorIndex];
  int whiteCount = agents [agentIndex].whitelist [dimension].sector [sectorIndex];
  int totalCount = blackCount + whiteCount;

  if (totalCount == 0) return false;

  double blackProbability = (double)blackCount / totalCount;
  return u.RNDprobab () < blackProbability;

  /*
  int totalCount = 0;

  for (int s = 0; s < sectorsPerCoord; s++)
  {
    totalCount += agents [agentIndex].blacklist [dimension].sector [s];
  }

  if (totalCount == 0)
  {
    return false;
  }

  int randomValue = u.RNDminusOne (totalCount);
  int cumulativeCount = 0;

  for (int s = 0; s < sectorsPerCoord; s++)
  {
    cumulativeCount += agents [agentIndex].blacklist [dimension].sector [s];

    if (randomValue <= cumulativeCount)
    {
      return true;
    }
  }

  return false;
  */
}
//——————————————————————————————————————————————————————————————————————————————
