//+————————————————————————————————————————————————————————————————————————————+
//|                                                                  C_AO_TSEA |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/14789

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
struct S_TSEA_Agent
{
    double c     [];     //coordinates
    double f;            //fitness
    int    label;        //cluster membership label
    int    labelClustV;  //clusters vertically
    //int    labelClustH;  //clusters horizontally
    double minDist;      //minimum distance to the nearest centroid

    void Init (int coords)
    {
      ArrayResize (c,     coords);
      f           = -DBL_MAX;
      label       = -1;
      labelClustV = -1;
      minDist     = DBL_MAX;
    }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
struct S_TSEA_horizontal
{
    //double cB [];
    int indBest;
    S_TSEA_Agent agent [];
};

struct S_TSEA_vertical
{
    S_TSEA_horizontal cell [];
};

//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
struct S_T_Cluster
{
    double centroid [];  //cluster centroid
    double f;            //centroid fitness
    int    count;        //number of points in the cluster
    int    ideasList []; //list of ideas

    void Init (int coords)
    {
      ArrayResize (centroid, coords);
      f     = -DBL_MAX;
      count = 0;
      ArrayResize (ideasList, 0, 100);
    }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
class C_TSEA_clusters
{
  public: //--------------------------------------------------------------------

  void KMeansInit (S_TSEA_Agent &data [], int dataSizeClust, S_T_Cluster &clust [])
  {
    for (int i = 0; i < ArraySize (clust); i++)
    {
      int ind = MathRand () % dataSizeClust;
      ArrayCopy (clust [i].centroid, data [ind].c, 0, 0, WHOLE_ARRAY);
    }
  }

  void KMeansPlusPlusInit (S_TSEA_Agent &data [], int dataSizeClust, S_T_Cluster &clust [])
  {
    // Choose the first centroid randomly
    int ind = MathRand () % dataSizeClust;
    ArrayCopy (clust [0].centroid, data [ind].c, 0, 0, WHOLE_ARRAY);

    for (int i = 1; i < ArraySize (clust); i++)
    {
      double sum = 0;

      // Compute the distance from each data point to the nearest centroid
      for (int j = 0; j < dataSizeClust; j++)
      {
        double minDist = DBL_MAX;

        for (int k = 0; k < i; k++)
        {
          double dist = VectorDistance (data [j].c, clust [k].centroid);

          if (dist < minDist)
          {
            minDist = dist;
          }
        }

        data [j].minDist = minDist;
        sum += minDist;
      }

      // Choose the next centroid with a probability proportional to the distance
      double randomValue    = ((double)rand () / 32767) * sum; // Generate a random value in the range [0, sum)
      double partialSum     = 0;
      bool   centroidChosen = false;

      for (int j = 0; j < dataSizeClust; j++)
      {
        partialSum += data [j].minDist;
        if (randomValue <= partialSum)
        {
          ArrayCopy (clust [i].centroid, data [j].c, 0, 0, WHOLE_ARRAY);
          centroidChosen = true;
          break;
        }
      }


      // If no point was assigned to the centroid, reassign it to the farthest point
      if (!centroidChosen)
      {
        double maxDist = -DBL_MAX;
        int farthestPointIndex = 0;

        for (int j = 0; j < dataSizeClust; j++)
        {
          if (data [j].minDist > maxDist)
          {
            maxDist = data [j].minDist;
            farthestPointIndex = j;
          }
        }
        Print ("Центроид ", i, " пустой");
        ArrayCopy (clust [i].centroid, data [farthestPointIndex].c, 0, 0, WHOLE_ARRAY);
      }
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

  void KMeans (S_TSEA_Agent &data [], int dataSizeClust, S_T_Cluster &clust [])
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


  //----------------------------------------------------------------------------
  struct DistanceIndex
  {
      double distance;
      int index;
  };

  void BubbleSort (DistanceIndex &arr [], int start, int end)
  {
    for (int i = start; i < end; i++)
    {
      for (int j = start; j < end - i; j++)
      {
        if (arr [j].distance > arr [j + 1].distance)
        {
          DistanceIndex temp = arr [j];
          arr [j] = arr [j + 1];
          arr [j + 1] = temp;
        }
      }
    }
  }

  int KNN (S_TSEA_Agent &data [], S_TSEA_Agent &point, int k_neighbors, int n_clusters)
  {
    int n = ArraySize (data);
    DistanceIndex distances_indices [];

    // Вычисление расстояний от точки до всех других точек
    for (int i = 0; i < n; i++)
    {
      DistanceIndex dist;
      dist.distance = VectorDistance (point.c, data [i].c);
      dist.index = i;
      ArrayResize (distances_indices, n);
      distances_indices [i] = dist;
    }

    // Сортировка расстояний
    BubbleSort (distances_indices, 0, n - 1);

    // Определение кластера для точки
    int votes [];
    ArrayResize (votes, n_clusters);
    ArrayInitialize (votes, 0);

    for (int j = 0; j < k_neighbors; j++)
    {
      int label = data [distances_indices [j].index].label;

      if (label != -1 && label < n_clusters)
      {
        votes [label]++;
      }
    }

    int max_votes = 0;
    int max_votes_cluster = -1;

    for (int j = 0; j < n_clusters; j++)
    {
      if (votes [j] > max_votes)
      {
        max_votes = votes [j];
        max_votes_cluster = j;
      }
    }

    return max_votes_cluster;
  }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
class C_AO_TSEA : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_TSEA () { }
  C_AO_TSEA ()
  {
    ao_name = "TSEA";
    ao_desc = "Turtle Shell Evolution Algorithm";
    ao_link = "https://www.mql5.com/ru/articles/14789";

    popSize         = 100; //population size

    vClusters       = 3;   //number of vertical clusters
    hClusters       = 20;  //number of horizontal clusters
    neighbNumb      = 5;   //number of nearest neighbors
    maxAgentsInCell = 3;   //max agents in cell

    ArrayResize (params, 5);

    params [0].name = "popSize";          params [0].val  = popSize;

    params [1].name = "vClusters";        params [1].val  = vClusters;
    params [2].name = "hClusters";        params [2].val  = hClusters;
    params [3].name = "neighbNumb";       params [3].val  = neighbNumb;
    params [4].name = "maxAgentsInCell";  params [4].val  = maxAgentsInCell;
  }

  void SetParams ()
  {
    popSize         = (int)params [0].val;

    vClusters       = (int)params [1].val;
    hClusters       = (int)params [2].val;
    neighbNumb      = (int)params [3].val;
    maxAgentsInCell = (int)params [4].val;
  }

  bool Init (const double &rangeMinP  [], //minimum search range
             const double &rangeMaxP  [], //maximum search range
             const double &rangeStepP [], //step search
             const int     epochsP = 0);  //number of epochs

  void Moving    ();
  void Revision  ();
  void Injection (const int popPos, const int coordPos, const double value);

  //----------------------------------------------------------------------------
  int    vClusters;      //number of vertical clusters
  int    hClusters;      //number of horizontal clusters
  int    neighbNumb;     //number of nearest neighbors
  int    maxAgentsInCell;

  S_TSEA_Agent    agent   [];
  S_TSEA_vertical cell    [];

  S_T_Cluster      clusters [];
  C_TSEA_clusters km;

  private: //-------------------------------------------------------------------
  double minFval;
  double stepF;

  int    epochs;
  int    epochsNow;
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_TSEA::Init (const double &rangeMinP  [], //minimum search range
                      const double &rangeMaxP  [], //maximum search range
                      const double &rangeStepP [], //step search
                      const int     epochsP = 0)   //number of epochs
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  ArrayResize (agent, popSize);
  for (int i = 0; i < popSize; i++) agent [i].Init (coords);

  ArrayResize (clusters, hClusters);
  for (int i = 0; i < hClusters; i++) clusters [i].Init (coords);

  ArrayResize (cell, vClusters);

  for (int i = 0; i < vClusters; i++)
  {
    ArrayResize (cell [i].cell, hClusters);

    for (int c = 0; c < hClusters; c++) ArrayResize (cell [i].cell [c].agent, 0, maxAgentsInCell);
  }

  minFval   = DBL_MAX;
  stepF     = 0.0;
  epochs    = epochsP;
  epochsNow = 0;

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_TSEA::Moving ()
{
  epochsNow++;

  //----------------------------------------------------------------------------
  //1. Сгенерировать случайные особи в популяцию
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
  int    vPos = 0;
  int    hPos = 0;
  int    pos  = 0;
  int    size = 0;
  double val  = 0.0;
  double rnd  = 0.0;
  double min  = 0.0;
  double max  = 0.0;

  for (int v = 0; v < vClusters; v++)
  {
    for (int h = 0; h < hClusters; h++)
    {
      size = ArraySize (cell [v].cell [h].agent);

      if (size > 0)
      {
        max = -DBL_MAX;
        pos = -1;

        for (int c = 0; c < size; c++)
        {
          if (cell [v].cell [h].agent [c].f > max)
          {
            max = cell [v].cell [h].agent [c].f;
            pos = c;
            cell [v].cell [h].indBest = c;
          }
        }
      }
    }
  }

  for (int i = 0; i < popSize; i++)
  {
    if (u.RNDprobab () < 0.8)
    {
      while (true)
      {
        rnd = u.RNDprobab ();
        rnd = (-rnd * rnd + 1.0) * vClusters;

        vPos = (int)rnd;
        if (vPos > vClusters - 1) vPos = vClusters - 1;

        hPos = u.RNDminusOne (hClusters);

        size = ArraySize (cell [vPos].cell [hPos].agent);

        if (size > 0) break;
      }

      pos = u.RNDminusOne (size);

      if (u.RNDprobab () < 0.5) pos = cell [vPos].cell [hPos].indBest;

      for (int c = 0; c < coords; c++)
      {
        if (u.RNDprobab () < 0.6) val = cell [vPos].cell [hPos].agent [pos].c [c];
        else                      val = cB [c];

        double dist = (rangeMax [c] - rangeMin [c]) * 0.1;
        min = val - dist; if (min < rangeMin [c]) min = rangeMin [c];
        max = val + dist; if (max > rangeMax [c]) max = rangeMax [c];

        val = u.PowerDistribution (val, min, max, 30);

        a [i].c [c] = u.SeInDiSp  (val, rangeMin [c], rangeMax [c], rangeStep [c]);

        agent [i].c [c] = a [i].c [c];
      }
    }
    else
    {
      int size2 = 0;
      int hPos2 = 0;
      int pos2  = 0;

      while (true)
      {
        rnd = u.RNDprobab ();
        rnd = (-rnd * rnd + 1.0) * vClusters;

        vPos = (int)rnd;
        if (vPos > vClusters - 1) vPos = vClusters - 1;

        hPos = u.RNDminusOne (hClusters);
        size = ArraySize (cell [vPos].cell [hPos].agent);

        hPos2 = u.RNDminusOne (hClusters);
        size2 = ArraySize (cell [vPos].cell [hPos2].agent);

        if (size > 0 && size2 > 0) break;
      }

      pos  = u.RNDminusOne (size);
      pos2 = u.RNDminusOne (size2);

      for (int c = 0; c < coords; c++)
      {
        val = (cell [vPos].cell [hPos ].agent [pos ].c [c] +
               cell [vPos].cell [hPos2].agent [pos2].c [c]) * 0.5;

        a [i].c [c] = u.SeInDiSp  (val, rangeMin [c], rangeMax [c], rangeStep [c]);

        agent [i].c [c] = a [i].c [c];
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_TSEA::Revision ()
{
  //получить приспособленность--------------------------------------------------
  int pos = -1;

  for (int i = 0; i < popSize; i++)
  {
    agent [i].f = a [i].f;

    if (a [i].f > fB)
    {
      fB = a [i].f;
      pos = i;
    }

    if (a [i].f < minFval) minFval = a [i].f;
  }

  if (pos != -1) ArrayCopy (cB, a [pos].c, 0, 0, WHOLE_ARRAY);

  stepF = (fB - minFval) / vClusters;

  //Разметка по вертикали дочерней популяции------------------------------------
  for (int i = 0; i < popSize; i++)
  {
    if (agent [i].f == fB) agent [i].labelClustV = vClusters - 1;
    else
    {
      agent [i].labelClustV = int((agent [i].f - minFval) / stepF);
      if (agent [i].labelClustV > vClusters - 1) agent [i].labelClustV = vClusters - 1;
    }
  }

  //----------------------------------------------------------------------------
  if (!revision)
  {
    km.KMeansPlusPlusInit (agent, popSize, clusters);
    km.KMeans             (agent, popSize, clusters);

    revision = true;
  }
  //----------------------------------------------------------------------------
  else
  {
    static S_TSEA_Agent data [];
    ArrayResize (data, 0, 1000);
    int size = 0;

    for (int v = 0; v < vClusters; v++)
    {
      for (int h = 0; h < hClusters; h++)
      {
        for (int c = 0; c < ArraySize (cell [v].cell [h].agent); c++)
        {
          size++;
          ArrayResize (data, size);

          data [size - 1] = cell [v].cell [h].agent [c];
        }
      }
    }

    for (int i = 0; i < popSize; i++)
    {
      agent [i].label = km.KNN (data, agent [i], neighbNumb, hClusters);
    }

    if (epochsNow % 50 == 0)
    {
      //km.KMeansPlusPlusInit (data, ArraySize (data), clusters);
      //km.KMeans             (data, ArraySize (data), clusters);

      for (int v = 0; v < vClusters; v++)
      {
        for (int h = 0; h < hClusters; h++)
        {
          ArrayResize (cell [v].cell [h].agent, 0);
        }
      }

      for (int i = 0; i < ArraySize (data); i++)
      {
        if (data [i].f == fB) data [i].labelClustV = vClusters - 1;
        else
        {
          data [i].labelClustV = int((data [i].f - minFval) / stepF);
          if (data [i].labelClustV > vClusters - 1) data [i].labelClustV = vClusters - 1;
        }

        int v = data [i].labelClustV;
        int h = data [i].label;

        int size = ArraySize (cell [v].cell [h].agent) + 1;
        ArrayResize (cell [v].cell [h].agent, size);

        cell [v].cell [h].agent [size - 1] = data [i];
      }
    }
  }

  //5, 10. Поместить популяцию в панцирь----------------------------------------
  for (int i = 0; i < popSize; i++)
  {
    int v = agent [i].labelClustV;
    int h = agent [i].label;

    int size = ArraySize (cell [v].cell [h].agent);
    int pos    = 0;
    int posMin = 0;
    int posMax = 0;

    if (size >= maxAgentsInCell)
    {
      double minF =  DBL_MAX;
      double maxF = -DBL_MAX;

      for (int c = 0; c < maxAgentsInCell; c++)
      {
        if (cell [v].cell [h].agent [c].f < minF)
        {
          minF = cell [v].cell [h].agent [c].f;
          posMin = c;
        }
        if (cell [v].cell [h].agent [c].f > maxF)
        {
          maxF = cell [v].cell [h].agent [c].f;
          posMax = c;
        }
      }

      if (v == 0)
      {
        if (agent [i].f < minF)
        {
          pos = posMin;
        }
        else
        {
          pos = posMax;
        }
      }
      else  pos = posMin;
    }
    else
    {
      ArrayResize (cell [v].cell [h].agent, size + 1);
      pos = size;
    }

    cell [v].cell [h].agent [pos] = agent [i];
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_TSEA::Injection (const int popPos, const int coordPos, const double value)
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
