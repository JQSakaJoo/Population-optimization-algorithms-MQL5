//+——————————————————————————————————————————————————————————————————+
//|                                                  C_AO_BCO_Bezier |
//|                                  Copyright 2007-2026, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

#include "#C_AO.mqh"

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class C_AO_BCO_Bezier : public C_AO
  {
public:
                    ~C_AO_BCO_Bezier() {}
                     C_AO_BCO_Bezier()
     {
      ao_name = "BCO(Bezier)";
      ao_desc = "Bezier Curve-based Optimization";
      ao_link = "https://www.mql5.com/ru/articles/24988";

      popSize = 20;     // размер популяции

      ArrayResize(params, 1);
      params [0].name = "popSize";
      params [0].val  = popSize;
     }

   void               SetParams()
     {
      popSize = (int)params [0].val;

      if(popSize < 3)
         popSize = 3;   // Eq.17 требует трёх разных агентов

      params [0].val = popSize;
     }

   bool               Init(const double &rangeMinP  [],
                           const double &rangeMaxP  [],
                           const double &rangeStepP [],
                           const int     epochsP = 0);

   void               Moving();
   void               Revision();

private:
   int                epochs;      // бюджет стенда в эпохах
   int                epochNow;    // выполненных эпох (номер итерации It)
   double             alpha0;      // Eq.10: знаменатель разброса alpha
   double             a0;          // Eq.19: огибающая фактора баланса A
   double             mean [];     // центр принятой популяции на начало итерации

   void               MakeCandidate(int i);
   void               Cubic(int i, double al);        // Eq.17-18
   void               Quadratic(int i, double al);    // Eq.15-16
   void               LinearBest(int i, double al);   // Eq.8-9
   void               LinearMean(int i, double al);   // Eq.11-12
   void               LinearCross(int i, double al);  // Eq.13-14

   void               PutCoord(int i, int c, double x);
   int                Other(int i, int j);
   double             RndOpen();
   double             Gauss();
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                              Init                                |
//+------------------------------------------------------------------+
bool C_AO_BCO_Bezier::Init(const double &rangeMinP  [],
                           const double &rangeMaxP  [],
                           const double &rangeStepP [],
                           const int     epochsP = 0)
  {
   if(!StandardInit(rangeMinP, rangeMaxP, rangeStepP))
      return false;

   epochs   = (epochsP > 0) ? epochsP : 1;
   epochNow = 0;
   alpha0   = 1.0;
   a0       = 1.0;

   ArrayResize(mean, coords);

   return true;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| SpaceBound оригинала: вышедшая координата — равномерно в диапазон|
//+------------------------------------------------------------------+
void C_AO_BCO_Bezier::PutCoord(int i, int c, double x)
  {
   if(x < rangeMin [c] || x > rangeMax [c])
      x = u.RNDfromCI(rangeMin [c], rangeMax [c]);

   a [i].c [c] = u.SeInDiSp(x, rangeMin [c], rangeMax [c], rangeStep [c]);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| случайный индекс агента, не равный i и j (j = -1 — не учитывать) |
//+------------------------------------------------------------------+
int C_AO_BCO_Bezier::Other(int i, int j)
  {
   int k;
   do
      k = u.RNDminusOne(popSize);
   while(k == i || k == j);
   return k;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| равномерное на (0; 1]                                            |
//+------------------------------------------------------------------+
double C_AO_BCO_Bezier::RndOpen()
  {
   double r;
   do
      r = u.RNDprobab();
   while(r <= 0.0);
   return r;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| N(0,1), Box-Muller                                               |
//+------------------------------------------------------------------+
double C_AO_BCO_Bezier::Gauss()
  {
   return sqrt(-2.0 * log(RndOpen())) * cos(2.0 * M_PI * u.RNDprobab());
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   MakeCandidate — выбор оператора (Eq.10, Eq.19).                |
//|   A > 1:  50% кубическая / 50% квадратичная;                     |
//|   A <= 1: 50% к лучшему / 25% к центру / 25% покоординатная.     |
//+------------------------------------------------------------------+
void C_AO_BCO_Bezier::MakeCandidate(int i)
  {
   double al = 1.0 + (2.0 * u.RNDprobab() - 1.0) / alpha0;    // Eq.10
   double A  = (0.4 + 2.0 * log(1.0 / RndOpen())) * a0;       // Eq.19

   if(A > 1.0)
     {
      if(u.RNDprobab() > 0.5)
         Cubic(i, al);
      else
         Quadratic(i, al);
     }
   else
     {
      if(u.RNDprobab() > 0.5)
         LinearBest(i, al);
      else
        {
         if(u.RNDprobab() > 0.5)
            LinearMean(i, al);
         else
            LinearCross(i, al);
        }
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   Cubic — кубическая кривая, исследование (Eq.17-18).            |
//|   P0 = x_i, P1 = x_j, P2 = x_k, P3 — по одной из трёх формул.    |
//+------------------------------------------------------------------+
void C_AO_BCO_Bezier::Cubic(int i, double al)
  {
   int    j    = Other(i, -1);
   int    k    = i;
   double ur   = 0.0;
   double r    = u.RNDprobab();
   double p3   = 0.0;
   int    mode = 0;

   if(r < 0.8)
     {
      mode = 0;
      k  = Other(i, j);
     }
   else
      if(r > 0.9)
        {
         mode = 1;
         k  = Other(i, j);
        }
      else
        {
         mode = 2;
        }

   double b  = 1.0 - al;
   double w0 = b * b * b;
   double w1 = 3.0 * b * b * al;
   double w2 = 3.0 * b * al * al;
   double w3 = al * al * al;

   for(int c = 0; c < coords; c++)
     {
      switch(mode)
        {
         case 0:
            p3 = mean [c] + 2.0 * Gauss() * (a [k].cB [c] - a [j].cB [c]);
            break;
         case 1:
            p3 = a [j].cB [c] + 2.0 * Gauss() * (a [k].cB [c] - a [j].cB [c]);
            break;
         default:
            p3 = a [j].cB [c] + 2.0 * Gauss() * (rangeMin [c] + u.RNDprobab() * (rangeMax [c] - rangeMin [c]));
            break;
        }

      PutCoord(i, c, w0 * a [i].cB [c] + w1 * a [j].cB [c] + w2 * a [k].cB [c] + w3 * p3);
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   Quadratic — квадратичная кривая, выход из ловушек (Eq.15-16).  |
//|   P0 = x_i, P1 = x_j, P2 = x_i + randn*(xbest|M - x_j).          |
//+------------------------------------------------------------------+
void C_AO_BCO_Bezier::Quadratic(int i, double al)
  {
   int    j      = Other(i, -1);
   double n      = Gauss();
   bool   toBest = u.RNDprobab() > 0.5;

   double b  = 1.0 - al;
   double w0 = b * b;
   double w1 = 2.0 * b * al;
   double w2 = al * al;
   double ref, p2;

   for(int c = 0; c < coords; c++)
     {
      ref = toBest ? cB [c] : mean [c];
      p2  = a [i].cB [c] + n * (ref - a [j].cB [c]);

      PutCoord(i, c, w0 * a [i].cB [c] + w1 * a [j].cB [c] + w2 * p2);
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   LinearBest — линейная кривая к лучшему (Eq.8-9).               |
//+------------------------------------------------------------------+
void C_AO_BCO_Bezier::LinearBest(int i, double al)
  {
   int  j    = Other(i, -1);
   bool diff = u.RNDprobab() < 0.5;
   double p1;

   for(int c = 0; c < coords; c++)
     {
      if(diff)
         p1 = cB [c] + (a [j].cB [c] - a [i].cB [c]) * 0.5;
      else
         p1 = cB [c] - (a [j].cB [c] + a [i].cB [c]) * 0.5;

      PutCoord(i, c, (1.0 - al) * a [i].cB [c] + al * p1);
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   LinearMean — линейная кривая к центру популяции (Eq.11-12).    |
//+------------------------------------------------------------------+
void C_AO_BCO_Bezier::LinearMean(int i, double al)
  {
   double n   = Gauss();
   bool   rel = u.RNDprobab() > 0.5;
   double p1;

   for(int c = 0; c < coords; c++)
     {
      if(rel)
         p1 = mean [c] + n * (a [i].cB [c] - mean [c]);
      else
         p1 = mean [c] + n * mean [c];

      PutCoord(i, c, (1.0 - al) * a [i].cB [c] + al * p1);
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   LinearCross — покоординатная интерполяция (Eq.13-14).          |
//|   Партнёр k выбирается заново для каждой координаты.             |
//+------------------------------------------------------------------+
void C_AO_BCO_Bezier::LinearCross(int i, double al)
  {
   double p1;

   for(int c = 0; c < coords; c++)
     {
      if(u.RNDprobab() > 0.5)
        {
         int k = Other(i, -1);
         p1 = (1.0 - al) * a [i].cB [c] + al * a [k].cB [c];
        }
      else
         p1 = a [i].cB [c];

      PutCoord(i, c, (1.0 - al) * a [i].cB [c] + al * p1);
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                            Moving                                |
//+------------------------------------------------------------------+
void C_AO_BCO_Bezier::Moving()
  {
//--- первый прогон: стартовая популяция
   if(!revision)
     {
      for(int i = 0; i < popSize; i++)
         for(int c = 0; c < coords; c++)
            a [i].c [c] = u.SeInDiSp(u.RNDfromCI(rangeMin [c], rangeMax [c]), rangeMin [c], rangeMax [c], rangeStep [c]);
      return;
     }

//--- расписания: It = epochNow (1..epochs-1), MaxIt = epochs-1
   int maxIt = epochs - 1;
   if(maxIt < 1)
      maxIt = 1;

   double tau = (double)epochNow / (double)maxIt;
   if(tau > 1.0)
      tau = 1.0;

   alpha0 = 1.0 + 80.0 * tau + 0.02 * pow(10.0 * tau, 3.0);   // Eq.10
   a0     = sin(M_PI * 0.5 * (1.0 - tau));                    // Eq.19

//--- центр принятой популяции (фиксируется на итерацию)
   for(int c = 0; c < coords; c++)
     {
      mean [c] = 0.0;
      for(int i = 0; i < popSize; i++)
         mean [c] += a [i].cB [c];
      mean [c] /= (double)popSize;
     }

   for(int i = 0; i < popSize; i++)
      MakeCandidate(i);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                           Revision                               |
//+------------------------------------------------------------------+
void C_AO_BCO_Bezier::Revision()
  {
//--- обновление лучших решений
   for(int i = 0; i < popSize; i++)
     {
      if(a [i].f > fB)
        {
         fB = a [i].f;
         ArrayCopy(cB, a [i].c, 0, 0, coords);
        }
      if(a [i].f > a [i].fB)
        {
         a [i].fB = a [i].f;
         ArrayCopy(a [i].cB, a [i].c, 0, 0, coords);
        }
     }

   epochNow++;

//--- первый проход: стартовая популяция становится принятой
   if(!revision)
     {
      for(int i = 0; i < popSize; i++)
        {
         ArrayCopy(a [i].cB, a [i].c, 0, 0, coords);
         a [i].fB = a [i].f;
        }

      revision = true;
      return;
     }
  }
//+------------------------------------------------------------------+
