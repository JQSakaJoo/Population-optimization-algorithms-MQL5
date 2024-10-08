//+————————————————————————————————————————————————————————————————————————————+
//|                                                                  C_AO_ACMO |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/15921

#include "#C_AO.mqh"

/*
1. создать случайные капли по пространству
2. внести информацию о влажности и давлении по регионам согласно выпавшим в регионах каплях
3.

*/

//——————————————————————————————————————————————————————————————————————————————
// Region structure
struct S_ACMO_Region
{
    double humidity; //humidity in the region
    double pressure; //pressure in the region
    double centre;   //the center of the region
    double x;        //point of highest pressure in the region

    void Init ()
    {
      humidity = -DBL_MAX;
      pressure = 0;
    }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
struct S_ACMO_Area
{
    S_ACMO_Region regions [];
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
// Cloud structure
struct S_ACMO_Cloud
{
    double center       [];    // cloud center
    double entropy      [];    // entropy
    double entropyStart [];    // initial entropy
    double hyperEntropy;       // hyperEntropy
    int    regionIndex  [];    // index of regions
    double averageHumidity;    // average humidity by regions
    double droplets;           // droplets

    void Init (int coords)
    {
      ArrayResize (center,       coords);
      ArrayResize (entropy,      coords);
      ArrayResize (entropyStart, coords);
      ArrayResize (regionIndex,  coords);
      droplets = 0.0;
    }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
class C_AO_ACMO : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_ACMO () { }
  C_AO_ACMO ()
  {
    ao_name = "ACMO";
    ao_desc = "Atmospheric Cloud Model Optimization";
    ao_link = "https://www.mql5.com/ru/articles/15921";

    popSize       = 50;    //population size

    cloudsNumber  = 4;     // Number of clouds
    regionsNumber = 10;    // Number of regions per dimension  (M)
    dMin          = 0.2;   // Minimum number of drops relative to the average number of drops in the clouds (dN)
    EnM0          = 0.2;   // Initial value of entropy
    HeM0          = 2.0;   // Initial value of hyperentropy
    λ             = 0.9;   // Threshold factor (threshold of the rainiest regions)
    γ             = 0.9;   // Weaken rate

    ArrayResize (params, 8);

    params [0].name = "popSize";       params [0].val = popSize;

    params [1].name = "cloudsNumber";  params [1].val = cloudsNumber;
    params [2].name = "regionsNumber"; params [2].val = regionsNumber;
    params [3].name = "dMin";          params [3].val = dMin;
    params [4].name = "EnM0";          params [4].val = EnM0;
    params [5].name = "HeM0";          params [5].val = HeM0;
    params [6].name = "λ";             params [6].val = λ;
    params [7].name = "γ";             params [7].val = γ;
  }

  void SetParams ()
  {
    popSize       = (int)params [0].val;

    cloudsNumber  = (int)params [1].val;
    regionsNumber = (int)params [2].val;
    dMin          = params      [3].val;
    EnM0          = params      [4].val;
    HeM0          = params      [5].val;
    λ             = params      [6].val;
    γ             = params      [7].val;
  }

  bool Init (const double &rangeMinP  [], //minimum search range
             const double &rangeMaxP  [], //maximum search range
             const double &rangeStepP [], //step search
             const int     epochsP = 0);  //number of epochs

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  int    cloudsNumber;  // Number of clouds
  int    regionsNumber; // Number of regions per dimension (M)
  double dMin;          // Minimum number of drops relative to the average number of drops in the clouds (dN)
  double EnM0;          // Initial value of entropy
  double HeM0;          // Initial value of hyperentropy
  double λ;             // Threshold factor
  double γ;             // Weaken rate


  S_ACMO_Area   areas  [];
  S_ACMO_Cloud  clouds [];

  private: //-------------------------------------------------------------------
  int    epochs;
  int    epochNow;
  int    dTotal;         // Maximum total number of droplets (N)
  double entropy [];     // Entropy
  double minGp;          // Minimum global pressure


  void   MoveClouds                 (bool &rev);
  int    GetRegionIndex             (double point, int ind);
  void   RainProcess                (bool &rev);
  void   DropletsDistribution       (double &clouds [], int &droplets []);

  void   UpdateRegionProperties     ();

  void   GenerateClouds             ();
  double CalculateHumidityThreshold ();
  void   CalculateNewEntropy        (S_ACMO_Cloud &cl, int t);
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_ACMO::Init (const double &rangeMinP  [],
                      const double &rangeMaxP  [],
                      const double &rangeStepP [],
                      const int     epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  epochs   = epochsP;
  epochNow = 0;

  //----------------------------------------------------------------------------
  dTotal       = popSize;
  dMin         = dMin * (popSize / (double)cloudsNumber);

  ArrayResize (entropy, coords);
  ArrayResize (areas,   coords);

  for (int c = 0; c < coords; c++)
  {
    entropy [c] = (rangeMax [c] - rangeMin [c]) / regionsNumber;

    ArrayResize (areas [c].regions, regionsNumber);

    for (int r = 0; r < regionsNumber; r++)
    {
      areas [c].regions [r].Init ();
      areas [c].regions [r].centre = rangeMin [c] + entropy [c] * (r + 0.5);
      areas [c].regions [r].centre = u.SeInDiSp (areas [c].regions [r].centre, rangeMin [c], rangeMax [c], rangeStep [c]);
      areas [c].regions [r].x      = areas [c].regions [r].centre;
    }
  }

  ArrayResize (clouds, cloudsNumber);
  for (int i = 0; i < cloudsNumber; i++) clouds [i].Init (coords);

  minGp = DBL_MAX;

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

/*
//------------------------------------------------------------------------------
1. На первой эпохе случайное размещение облаков:
   EnCk = EnM0;
   HeCk = HeM0;
//------------------------------------------------------------------------------
1.1 Движение облаков в сторону регионов с меньшим давлением:
   β = deltaP / normP
   d = Tck.x - Cck.c
   VC = β * d;
   Ck = Ck + VC

   изменение поличества капель после движения:
   nk = nk × (1 - γ)

   изменение энтропии и гиперэнтропии:
   α = ΔP / ΔPmax;
   EnCk = EnCk * (1 + α)
   HeCk = HeCk * (1 - α)
//------------------------------------------------------------------------------
2. Процесс дождя, выпадение капель:
   распределение капель между облаками пропорционально влажности региона
   увеличение количества капель к существующим в облаках
//------------------------------------------------------------------------------
3. Расчет фитнес-функции для капель
//------------------------------------------------------------------------------
4. Обновление глобального решения и минимального давления в регионах, где выпадали капли
//------------------------------------------------------------------------------
5. Проверка на распад облаков и создание новых в замен распавшихся в регионах больше порога:
   правило распада вследствии расширения больше допустимого значения (разрыв облака):
   En > 5 * EnM0_t
   правило распада при содержании влаги ниже критического значения (высушение облака):
   dCk < dMin

   пороговое значение, выше которого регионы могут образовать облака:
   HT = H_min + λ * (H_max - H_min);
//------------------------------------------------------------------------------
6. Энтропия и гиперэнтропия для новых облаков расчитать:
   En = EnM0 / (1+2.72^(-(8-16*(t/maxT))))
   He = HeM0 / (1+2.72^((8-16*(t/maxT))))
*/
//——————————————————————————————————————————————————————————————————————————————
void C_AO_ACMO::Moving ()
{
  MoveClouds  (revision);
  RainProcess (revision);
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_ACMO::MoveClouds (bool &rev)
{
  //----------------------------------------------------------------------------
  if (!rev)
  {
    //creating clouds with random centers---------------------------------------
    for (int i = 0; i < cloudsNumber; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        clouds [i].center [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        clouds [i].center [c] = u.SeInDiSp  (clouds [i].center [c], rangeMin [c], rangeMax [c], rangeStep [c]);

        clouds [i].regionIndex  [c] = GetRegionIndex (clouds [i].center [c], c);
        clouds [i].entropy      [c] = entropy [c] * EnM0;
        clouds [i].entropyStart [c] = clouds [i].entropy [c];
      }

      if (i != 0) clouds [i].hyperEntropy = HeM0;
      else        clouds [i].hyperEntropy = 8;
    }

    return;
  }

  //----------------------------------------------------------------------------
  //search for the region with the lowest pressure------------------------------
  int targetRegion = 0;

  int lHind []; //lowest humidity index
  ArrayResize     (lHind, coords);
  ArrayInitialize (lHind, 0);

  double minP;
  double maxP;
  double normP [];
  ArrayResize (normP, coords);

  for (int c = 0; c < coords; c++)
  {
    minP =  DBL_MAX;
    maxP = -DBL_MAX;

    for (int r = 0; r < regionsNumber; r++)
    {
      if (areas [c].regions [r].pressure < areas [c].regions [lHind [c]].pressure)
      {
        lHind [c] = r;
      }

      if (areas [c].regions [r].pressure < minP) minP = areas [c].regions [r].pressure;
      if (areas [c].regions [r].pressure > maxP) maxP = areas [c].regions [r].pressure;
    }

    normP [c] = maxP - minP + DBL_EPSILON;
  }

  //moving the cloud to a region with less pressure-----------------------------
  int    clRegIND = 0;
  double deltaP   = 0.0;
  double α        = 0.0; // Entropy factor
  double β        = 0.0; // Atmospheric pressure factor
  double VC       = 0.0; // Cloud velocity
  double d        = 0.0; // Cloud direction

  for (int i = 0; i < cloudsNumber; i++)
  {
    //create an artificial cloud in the region with the highest humidity,
    //even if it is the region with the highest pressure.
    if (i == 0)
    {
      for (int c = 0; c < coords; c++)
      {
        clouds [i].regionIndex [c] = GetRegionIndex (cB [c], c);
        clouds [i].center      [c] = cB [c];
      }
    }
    else
    {
      for (int c = 0; c < coords; c++)
      {
        //find a region with lower pressure-------------------------------------
        if (clouds [i].regionIndex [c] == lHind [c]) continue;

        clRegIND = clouds [i].regionIndex [c];

        do targetRegion = u.RNDminusOne (regionsNumber);
        while (areas [c].regions [clRegIND].pressure < areas [c].regions [targetRegion].pressure);

        //------------------------------------------------------------------------
        deltaP = areas [c].regions [clRegIND].pressure - areas [c].regions [targetRegion].pressure;

        β = deltaP / normP [c];
        d = areas [c].regions [targetRegion].x - areas [c].regions [clRegIND].centre;
        //d = clouds [i].entropy [c];

        VC = β * d;

        clouds [i].center      [c] += VC;
        clouds [i].center      [c] = u.SeInDiSp (clouds [i].center [c], rangeMin [c], rangeMax [c], rangeStep [c]);

        clouds [i].regionIndex [c] = GetRegionIndex (clouds [i].center [c], c);

        α = β;
        clouds [i].entropy [c] *=(1 + α);
      }

      clouds [i].droplets     *=(1 - γ);
      clouds [i].hyperEntropy *=(1 + α);
      if (clouds [i].hyperEntropy > 8) clouds [i].hyperEntropy = 8;
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
int C_AO_ACMO::GetRegionIndex (double point, int ind)
{
  if (point <= rangeMin [ind]) return 0;
  if (point >= rangeMax [ind]) return regionsNumber - 1;

  int regPos = (int)((point - rangeMin [ind]) / entropy [ind]);

  return regPos;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_ACMO::RainProcess (bool &rev)
{
  //to shed drops from every cloud----------------------------------------------
  double cloud [];
  int    drops [];
  ArrayResize (cloud, cloudsNumber);
  ArrayResize (drops, cloudsNumber);

  if (!rev)
  {
    ArrayInitialize (cloud, 1.0);
  }
  else
  {
    ArrayInitialize (cloud, 0.0);

    double humidity;

    for (int i = 0; i < cloudsNumber; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        humidity = areas [c].regions [clouds [i].regionIndex [c]].humidity;
        if (humidity != -DBL_MAX) cloud [i] += humidity;
        else                      cloud [i] += minGp;
      }
    }
  }

  DropletsDistribution (cloud, drops);
  //ArrayPrint (drops);

  double dist   = 0.0;
  double centre = 0.0;
  double xMin   = 0.0;
  double xMax   = 0.0;
  double x      = 0.0;
  int    dCNT   = 0;

  for (int i = 0; i < cloudsNumber; i++)
  {
    for (int dr = 0; dr < drops [i]; dr++)
    {
      for (int c = 0; c < coords; c++)
      {
        dist   = clouds [i].entropy [c];
        centre = clouds [i].center  [c];
        xMin   = centre - dist;
        xMax   = centre + dist;

        x = u.GaussDistribution (centre, xMin, xMax, clouds [i].hyperEntropy);

        if (x < rangeMin [c]) x = u.RNDfromCI (rangeMin [c], centre);
        if (x > rangeMax [c]) x = u.RNDfromCI (centre, rangeMax [c]);

        x = u.SeInDiSp (x, rangeMin [c], rangeMax [c], rangeStep [c]);

        int p = u.RNDminusOne (popSize);

        if (a [p].f > a [dCNT].f)
        {
          if (u.RNDprobab () < 0.95) a [dCNT].c [c] = a [p].c [c];
        }
        else
        {
          a [dCNT].c [c] = x;
        }
      }

      dCNT++;
    }

    clouds [i].droplets += drops [i];
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_ACMO::DropletsDistribution (double &cloud [], int &droplets [])
{
  double minHumidity    = DBL_MAX;
  int    indMinHumidity = -1;
  double totalHumidity  = DBL_EPSILON; //total amount of humidity in all clouds

  for (int i = 0; i < ArraySize (cloud); i++)
  {
    totalHumidity += cloud [i];

    if (cloud [i] < minHumidity)
    {
      minHumidity = cloud [i];
      indMinHumidity = i;
    }
  }

  if (totalHumidity == 0.0) totalHumidity = DBL_EPSILON;

  // Filling the droplets array in proportion to the value in clouds
  for (int i = 0; i < ArraySize (clouds); i++)
  {
    droplets [i] = int((cloud [i] / totalHumidity) * popSize); //proportional distribution of droplets
  }

  // Distribute the remaining drops, if any
  int totalDrops = 0;

  for (int i = 0; i < ArraySize (droplets); i++)
  {
    totalDrops += droplets [i];
  }

  // If not all drops are distributed, add the remaining drops to the element with the lowest humidity
  int remainingDrops = popSize - totalDrops;

  if (remainingDrops > 0)
  {
    droplets [indMinHumidity] += remainingDrops; //add the remaining drops to the lightest cloud
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_ACMO::Revision ()
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

    if (a [i].f < minGp) minGp = a [i].f;
  }

  if (ind != -1) ArrayCopy (cB, a [ind].c, 0, 0, WHOLE_ARRAY);

  //----------------------------------------------------------------------------
  UpdateRegionProperties (); //updating humidity and pressure in the regions
  GenerateClouds         (); //disappearance of clouds and the creation of new ones

  revision = true;
  epochNow++;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_ACMO::UpdateRegionProperties ()
{
  int regionIndex = 0;

  for (int i = 0; i < dTotal; i++)
  {
    for (int c = 0; c < coords; c++)
    {
      regionIndex = GetRegionIndex (a [i].c [c], c);

      if (a [i].f > areas [c].regions [regionIndex].humidity)
      {
        areas [c].regions [regionIndex].humidity  = a [i].f;
        areas [c].regions [regionIndex].pressure += 1.0;
        areas [c].regions [regionIndex].x         = a [i].c [c];
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_ACMO::GenerateClouds ()
{
  //Collecting statistics of regions capable of creating clouds-----------------
  double Ht = CalculateHumidityThreshold ();

  struct S_Areas
  {
      int regsIND []; //index of the potential region
  };

  S_Areas ar [];
  ArrayResize (ar, coords);

  int sizePr = 0;

  for (int i = 0; i < coords; i++)
  {
    for (int r = 0; r < regionsNumber; r++)
    {
      if (areas [i].regions [r].humidity > Ht)
      {
        sizePr = ArraySize (ar [i].regsIND);
        sizePr++;
        ArrayResize (ar [i].regsIND, sizePr, coords);
        ar [i].regsIND [sizePr - 1] = r;
      }
    }
  }

  //Check the conditions for cloud decay----------------------------------------
  bool cloudDecay = false;

  for (int i = 0; i < cloudsNumber; i++)
  {
    if (i != 0)
    {
      cloudDecay = false;

      //checking the cloud for too much entropy---------------------------------
      for (int c = 0; c < coords; c++)
      {
        if (clouds [i].entropy [c] > 5 * clouds [i].entropyStart [c])
        {
          cloudDecay = true;
          break;
        }
      }

      //checking the cloud for decay--------------------------------------------
      if (!cloudDecay)
      {
        if (clouds [i].droplets < dMin)
        {
          cloudDecay = true;
        }
      }

      //if the cloud has decayed------------------------------------------------
      int regIND = 0;

      if (cloudDecay)
      {
        //creating a cloud in a very humid region-------------------------------
        for (int c = 0; c < coords; c++)
        {
          regIND = u.RNDminusOne (ArraySize (ar [c].regsIND));
          regIND = ar [c].regsIND [regIND];

          clouds [i].center [c] = areas [c].regions [regIND].x;

          clouds [i].regionIndex [c] = regIND;
        }

        CalculateNewEntropy (clouds [i], epochNow);
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_ACMO::CalculateHumidityThreshold ()
{
  double H_max = fB;
  double H_min = DBL_MAX;

  for (int c = 0; c < coords; c++)
  {
    for (int r = 0; r < regionsNumber; r++)
    {
      if (areas [c].regions [r].humidity != -DBL_MAX)
      {
        if (areas [c].regions [r].humidity < H_min)
        {
          H_min = areas [c].regions [r].humidity;
        }
      }
    }
  }

  return H_min + λ * (H_max - H_min);
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_ACMO::CalculateNewEntropy (S_ACMO_Cloud &cl, int t)
{
  //----------------------------------------------------------------------------
  //En: 1/(1+2.72^(-(8-16*(t/maxT))))
  for (int c = 0; c < coords; c++)
  {
    cl.entropy      [c] = entropy [c] * EnM0 / (1.0 + pow (M_E, (-(8.0 - 16.0 * (t / epochs)))));
    //cl.entropy      [c] = entropy [c] * EnM0 * (1-0.9999*(t / epochs));

    cl.entropyStart [c] = cl.entropy [c] = entropy [c];
  }

  //----------------------------------------------------------------------------
  //He: 1/(1+2.72^((8-16*(t/maxT))))
  cl.hyperEntropy = 1.0 / (1.0 + pow (M_E, ((8.0 - 16.0 * (t / epochs)))));

  cl.hyperEntropy = u.Scale (cl.hyperEntropy, 0.0, 8.0, HeM0, 8.0);
}
//——————————————————————————————————————————————————————————————————————————————
