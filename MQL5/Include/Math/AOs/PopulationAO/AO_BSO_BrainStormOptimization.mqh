//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_BSO |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/14622

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
struct S_BSO_Agent
{
    double c     []; //coordinates
    double f;        //fitness
    int    label;    //cluster membership label

    void Init (int coords)
    {
      ArrayResize (c,     coords);
      f     = -DBL_MAX;
      label = -1;
    }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
struct S_Clusters
{
    double centroid [];  //cluster centroid
    double f;            //centroid fitness
    int    count;        //number of points in the cluster
    int    ideasList []; //list of ideas

    void Init (int coords)
    {
      ArrayResize (centroid, coords);
      f = -DBL_MAX;
      ArrayResize (ideasList, 0, 100);
    }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
class S_BSO_KMeans
{
  public: //--------------------------------------------------------------------

  void KMeansInit (S_BSO_Agent &data [], int dataSizeClust, S_Clusters &clust [])
  {
    for (int i = 0; i < ArraySize (clust); i++)
    {
      int ind = MathRand () % dataSizeClust;
      ArrayCopy (clust [i].centroid, data [ind].c, 0, 0, WHOLE_ARRAY);
    }
  }

  double VectorDistance (double &v1 [], double &v2 [])
  {
    double distance = 0.0;
    for (int i = 0; i < ArraySize (v1); i++)
    {
      distance += (v1 [i] - v2 [i]) * (v1 [i] - v2 [i]);
    }
    return MathSqrt (distance);
  }

  void KMeans (S_BSO_Agent &data [], int dataSizeClust, S_Clusters &clust [])
  {
    bool changed   = true;
    int  nClusters = ArraySize (clust);
    int  cnt       = 0;

    while (changed && cnt < 100)
    {
      cnt++;
      changed = false;

      // Назначение точек данных к ближайшему центроиду
      for (int d = 0; d < dataSizeClust; d++)
      {
        int    closest_centroid = -1;
        double closest_distance = DBL_MAX;

        if (data [d].f != -DBL_MAX)
        {
          for (int cl = 0; cl < nClusters; cl++)
          {
            double distance = VectorDistance (data [d].c, clust [cl].centroid);

            if (distance < closest_distance)
            {
              closest_distance = distance;
              closest_centroid = cl;
            }
          }

          if (data [d].label != closest_centroid)
          {
            data [d].label = closest_centroid;
            changed = true;
          }
        }
        else
        {
          data [d].label = -1;
        }
      }


      // Обновление центроидов
      double sum_c [];
      ArrayResize (sum_c, ArraySize (data [0].c));

      for (int cl = 0; cl < nClusters; cl++)
      {
        ArrayInitialize (sum_c, 0.0);

        clust [cl].count = 0;
        ArrayResize (clust [cl].ideasList, 0);

        for (int d = 0; d < dataSizeClust; d++)
        {
          if (data [d].label == cl)
          {
            for (int k = 0; k < ArraySize (data [d].c); k++)
            {
              sum_c [k] += data [d].c [k];
            }

            clust [cl].count++;
            ArrayResize (clust [cl].ideasList, clust [cl].count);
            clust [cl].ideasList [clust [cl].count - 1] = d;
          }
        }

        if (clust [cl].count > 0)
        {
          for (int k = 0; k < ArraySize (sum_c); k++)
          {
            clust [cl].centroid [k] = sum_c [k] / clust [cl].count;
          }
        }
      }
    }
  }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
class C_AO_BSO : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_BSO () { }
  C_AO_BSO ()
  {
    ao_name = "BSO";
    ao_desc = "Brain Storm Optimization";
    ao_link = "https://www.mql5.com/ru/articles/14622";

    popSize        = 25;   //population size

    parentPopSize  = 50;   //parent population size;
    clustersNumb   = 5;    //number of clusters
    p_Replace      = 0.1;  //replace probability
    p_One          = 0.5;  //probability of choosing one
    p_One_center   = 0.3;  //probability of choosing one center
    p_Two_center   = 0.2;  //probability of choosing two centers
    k_Mutation     = 20.0; //mutation coefficient
    distribCoeff   = 1.0;  //distribution coefficient

    ArrayResize (params, 9);

    params [0].name = "popSize";       params [0].val  = popSize;

    params [1].name = "parentPopSize"; params [1].val  = parentPopSize;
    params [2].name = "clustersNumb";  params [2].val  = clustersNumb;
    params [3].name = "p_Replace";     params [3].val  = p_Replace;
    params [4].name = "p_One";         params [4].val  = p_One;
    params [5].name = "p_One_center";  params [5].val  = p_One_center;
    params [6].name = "p_Two_center";  params [6].val  = p_Two_center;
    params [7].name = "k_Mutation";    params [7].val  = k_Mutation;
    params [8].name = "distribCoeff";  params [8].val  = distribCoeff;
  }

  void SetParams ()
  {
    popSize       = (int)params [0].val;

    parentPopSize = (int)params [1].val;
    clustersNumb  = (int)params [2].val;
    p_Replace     = params      [3].val;
    p_One         = params      [4].val;
    p_One_center  = params      [5].val;
    p_Two_center  = params      [6].val;
    k_Mutation    = params      [7].val;
    distribCoeff  = params      [8].val;
  }

  bool Init (const double &rangeMinP  [], //minimum search range
             const double &rangeMaxP  [], //maximum search range
             const double &rangeStepP [], //step search
             const int     epochsP = 0);  //number of epochs

  void Moving    ();
  void Revision  ();
  void Injection (const int popPos, const int coordPos, const double value);

  //----------------------------------------------------------------------------
  int    parentPopSize; //parent population size;
  int    clustersNumb;  //number of clusters
  double p_Replace;     //replace probability
  double p_One;         //probability of choosing one
  double p_One_center;  //probability of choosing one center
  double p_Two_center;  //probability of choosing two centers
  double k_Mutation;    //mutation coefficient
  double distribCoeff;  //distribution coefficient

  S_BSO_Agent  agent   [];
  S_BSO_Agent  parents [];

  S_Clusters   clusters [];
  S_BSO_KMeans km;

  private: //-------------------------------------------------------------------
  S_BSO_Agent  parentsTemp [];
  int          epochs;
  int          epochsNow;
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_BSO::Init (const double &rangeMinP  [], //minimum search range
                     const double &rangeMaxP  [], //maximum search range
                     const double &rangeStepP [], //step search
                     const int     epochsP = 0)   //number of epochs
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  ArrayResize (agent, popSize);
  for (int i = 0; i < popSize; i++) agent [i].Init (coords);

  ArrayResize (clusters, clustersNumb);
  for (int i = 0; i < clustersNumb; i++) clusters [i].Init (coords);

  ArrayResize (parents,     parentPopSize + popSize);
  ArrayResize (parentsTemp, parentPopSize + popSize);

  for (int i = 0; i < parentPopSize + popSize; i++)
  {
    parents     [i].Init (coords);
    parentsTemp [i].Init (coords);
  }

  epochs    = epochsP;
  epochsNow = 0;

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

/*
1. Инициализация популяции из n индивидов, m кластеров и максимального числа итераций gmax.
2. Оценка приспособленности

Цикл итераций до достижения максимального числа итераций gmax.
  Кластеризация: Индивиды группируются в m кластеров в зависимости от их приспособленности.
  Установить лучшее решение в кластере как центр кластера.

  Если Preplace
      генерируется новый индивид, который заменяет выбранный центр кластера (Из центра кластера)

  Если Pone
      выбирается индивид из одного кластера.

      Если Pone_center
          выбирается центр кластера для мутации,
      иначе
          случайный индивид из этого кластера.
  Иначе
      выбираются индивиды из двух кластеров.

      Если Ptwo_center,
          то два центра кластера объединяются и мутируют
      иначе
          случайно выбираются два индивида из каждого выбранного кластера, которые затем объединяются и мутируют.

  Мутация:
      Полученный индивид подвергается мутации с помощью гауссовой мутации

  Вычисляется его приспособленность.

  Отбор:
      отбор, в результате которого в популяции остаются только наилучшие индивиды.
*/

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BSO::Moving ()
{
  epochsNow++;

  //----------------------------------------------------------------------------
  if (!revision)
  {
    for (int i = 0; i < popSize; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        a [i].c [c] = u.SeInDiSp  (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);

        agent [i].c [c] = a [i].c [c];
      }
    }

    return;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  int    cIndx_1    = 0;  //индекс в списке непустых кластеров
  int    iIndx_1    = 0;  //индекс в списке идей в кластере
  int    cIndx_2    = 0;  //индекс в списке непустых кластеров
  int    iIndx_2    = 0;  //индекс в списке идей в кластере
  double min        = 0.0;
  double max        = 0.0;
  double dist       = 0.0;
  double val        = 0.0;
  double X1         = 0.0;
  double X2         = 0.0;
  int    clListSize = 0;
  int    clustList [];
  ArrayResize (clustList, 0, clustersNumb);

  //----------------------------------------------------------------------------
  //составим список непустых кластеров
  for (int cl = 0; cl < clustersNumb; cl++)
  {
    if (clusters [cl].count > 0)
    {
      clListSize++;
      ArrayResize (clustList, clListSize);
      clustList [clListSize - 1] = cl;
    }
  }

  for (int i = 0; i < popSize; i++)
  {
    //==========================================================================
    //генерация новой идеи, которая заменяет выбранный центр кластера (смещение центра кластера)
    if (u.RNDprobab () < p_Replace)
    {
      cIndx_1 = u.RNDminusOne (clListSize);

      for (int c = 0; c < coords; c++)
      {
        val = clusters [clustList [cIndx_1]].centroid [c];

        dist = (rangeMax [c] - rangeMin [c]) * 0.8;

        min = val - dist; if (min < rangeMin [c]) min = rangeMin [c];
        max = val + dist; if (max > rangeMax [c]) max = rangeMax [c];

        val = u.GaussDistribution (val, min, max, 3);
        val = u.SeInDiSp  (val, rangeMin [c], rangeMax [c], rangeStep [c]);

        clusters [clustList [cIndx_1]].centroid [c] = val;
      }
    }

    //==========================================================================
    //выбирается идея из одного кластера
    if (u.RNDprobab () < p_One)
    {
      cIndx_1 = u.RNDminusOne (clListSize);

      //------------------------------------------------------------------------
      if (u.RNDprobab () < p_One_center) //выбирается центр кластера
      {
        for (int c = 0; c < coords; c++)
        {
          a [i].c [c] = clusters [clustList [cIndx_1]].centroid [c];
        }
      }
      //------------------------------------------------------------------------
      else                               //случайная идея из этого кластера
      {
        iIndx_1 = u.RNDminusOne (clusters [clustList [cIndx_1]].count);

        for (int c = 0; c < coords; c++)
        {
          a [i].c [c] = parents [clusters [clustList [cIndx_1]].ideasList [iIndx_1]].c [c];
        }
      }
    }
    //==========================================================================
    //выбираются идеи из двух кластеров
    else
    {
      if (clListSize == 1)
      {
        cIndx_1 = 0;
        cIndx_2 = 0;
      }
      else
      {
        if (clListSize == 2)
        {
          cIndx_1 = 0;
          cIndx_2 = 1;
        }
        else
        {
          cIndx_1 = u.RNDminusOne (clListSize);

          do
          {
            cIndx_2 = u.RNDminusOne (clListSize);
          }
          while (cIndx_1 == cIndx_2);
        }
      }

      //------------------------------------------------------------------------
      if (u.RNDprobab () < p_Two_center) //выбрали два центра кластеров
      {
        for (int c = 0; c < coords; c++)
        {
          X1 = clusters [clustList [cIndx_1]].centroid [c];
          X2 = clusters [clustList [cIndx_2]].centroid [c];

          a [i].c [c] = u.RNDfromCI (X1, X2);
        }
      }
      //------------------------------------------------------------------------
      else //две идеи из двух выбранных кластеров
      {
        iIndx_1 = u.RNDminusOne (clusters [clustList [cIndx_1]].count);
        iIndx_2 = u.RNDminusOne (clusters [clustList [cIndx_2]].count);

        for (int c = 0; c < coords; c++)
        {
          X1 = parents [clusters [clustList [cIndx_1]].ideasList [iIndx_1]].c [c];
          X2 = parents [clusters [clustList [cIndx_2]].ideasList [iIndx_2]].c [c];

          a [i].c [c] = u.RNDfromCI (X1, X2);
        }
      }
    }

    //==========================================================================
    //Мутация
    for (int c = 0; c < coords; c++)
    {
      int x = (int)u.Scale (epochsNow, 1, epochs, 1, 200);

      double ξ = (1.0 / (1.0 + exp (-((100 - x) / k_Mutation))));// * u.RNDprobab ();

      double dist = (rangeMax [c] - rangeMin [c]) * distribCoeff * ξ;
      double min = a [i].c [c] - dist; if (min < rangeMin [c]) min = rangeMin [c];
      double max = a [i].c [c] + dist; if (max > rangeMax [c]) max = rangeMax [c];

      val = a [i].c [c];

      a [i].c [c] = u.GaussDistribution (val, min, max, 8);
    }

    //Сохраним агента-----------------------------------------------------------
    for (int c = 0; c < coords; c++)
    {
      val = u.SeInDiSp  (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);

      a     [i].c [c] = val;
      agent [i].c [c] = val;
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BSO::Revision ()
{
  //получить приспособленность--------------------------------------------------
  for (int i = 0; i < popSize; i++)
  {
    agent [i].f = a [i].f;
  }

  //перенести новые идеи в популяцию--------------------------------------------
  for (int i = parentPopSize; i < parentPopSize + popSize; i++)
  {
    parents [i] = agent [i - parentPopSize];
  }

  //отсортировать родительскую популяцию----------------------------------------
  u.Sorting (parents, parentsTemp, parentPopSize + popSize);

  if (parents [0].f > fB)
  {
    fB = parents [0].f;
    ArrayCopy (cB, parents [0].c, 0, 0, WHOLE_ARRAY);
  }


  //выполнить кластеризацию-----------------------------------------------------
  if (!revision)
  {
    km.KMeansInit (parents, parentPopSize, clusters);
    revision = true;
  }

  km.KMeansInit (parents, parentPopSize, clusters);
  km.KMeans     (parents, parentPopSize, clusters);

  //Назначить лучшее решение кластера центром кластера--------------------------
  for (int cl = 0; cl < clustersNumb; cl++)
  {
    clusters [cl].f = -DBL_MAX;

    if (clusters [cl].count > 0)
    {
      for (int p = 0; p < parentPopSize; p++)
      {
        if (parents [p].label == cl)
        {
          if (parents [p].f > clusters [cl].f)
          {
            clusters [cl].f = parents [p].f;
            ArrayCopy (clusters [cl].centroid, parents [p].c, 0, 0, WHOLE_ARRAY);
          }
        }
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BSO::Injection (const int popPos, const int coordPos, const double value)
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
