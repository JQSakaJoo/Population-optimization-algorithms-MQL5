//+——————————————————————————————————————————————————————————————————+
//|                                                         C_AO_FOX |
//|                                  Copyright 2007-2026, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

#include "#C_AO.mqh"

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class C_AO_FOX : public C_AO
  {
public:
                    ~C_AO_FOX() {}
                     C_AO_FOX()
     {
      ao_name = "FOX";
      ao_desc = "Fox-inspired Optimization Algorithm";
      ao_link = "https://www.mql5.com/ru/articles/25488";

      popSize  = 30;     // размер популяции
      pExploit = 0.5;    // вероятность ветки прыжка (r >= 0.5)
      pDir     = 0.18;   // порог направления прыжка (p > 0.18 -> c1)
      c1       = 0.18;   // коэффициент прыжка на северо-восток
      c2       = 0.82;   // коэффициент прыжка в обратную сторону
      fixMinT  = false;  // false — MinT = +inf на старте, как у авторов

      ArrayResize(params, 6);
      params [0].name = "popSize";
      params [0].val = popSize;
      params [1].name = "pExploit";
      params [1].val = pExploit;
      params [2].name = "pDir";
      params [2].val = pDir;
      params [3].name = "c1";
      params [3].val = c1;
      params [4].name = "c2";
      params [4].val = c2;
      params [5].name = "fixMinT";
      params [5].val = fixMinT ? 1.0 : 0.0;
     }

   void               SetParams()
     {
      popSize  = (int)params [0].val;
      pExploit = params [1].val;
      pDir     = params [2].val;
      c1       = params [3].val;
      c2       = params [4].val;
      fixMinT  = params [5].val != 0.0;

      if(popSize < 1)
         popSize = 1;
      if(pExploit < 0.0)
         pExploit = 0.0;
      if(pExploit > 1.0)
         pExploit = 1.0;
      if(pDir < 0.0)
         pDir = 0.0;
      if(pDir > 1.0)
         pDir = 1.0;

      params [0].val = popSize;
      params [1].val = pExploit;
      params [2].val = pDir;
      params [3].val = c1;
      params [4].val = c2;
      params [5].val = fixMinT ? 1.0 : 0.0;
     }

   bool               Init(const double &rangeMinP  [],
                           const double &rangeMaxP  [],
                           const double &rangeStepP [],
                           const int     epochsP = 0);

   void               Moving();
   void               Revision();

   //----------------------------------------------------------------
   double             pExploit;   // доля шагов-прыжков
   double             pDir;       // порог выбора c1 / c2
   double             c1;
   double             c2;
   bool               fixMinT;    // true — MinT стартует с 0.5 вместо +inf

private:
   int                epochs;     // Max_iter
   int                epochNow;   // l
   double             minT;       // MinT, глобальный за весь прогон
   double             tm [];      // Time(i,:) — буфер текущего агента

   void               Jump(int i, double c);
   void               Walk(int i, double aCoef);

   double             Gauss();
   double             RndOpen();
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                              Init                                |
//+------------------------------------------------------------------+
bool C_AO_FOX::Init(const double &rangeMinP  [],
                    const double &rangeMaxP  [],
                    const double &rangeStepP [],
                    const int     epochsP = 0)
  {
   if(!StandardInit(rangeMinP, rangeMaxP, rangeStepP))
      return false;

   epochs   = epochsP > 0 ? epochsP : 1;
   epochNow = 0;
   minT     = fixMinT ? 0.5 : DBL_MAX;   // 0.5 = E[mean(U(0,1))]

   ArrayResize(tm, coords);
   return true;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| равномерное на (0; 1]                                            |
//+------------------------------------------------------------------+
double C_AO_FOX::RndOpen()
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
double C_AO_FOX::Gauss()
  {
   return sqrt(-2.0 * log(RndOpen())) * cos(2.0 * M_PI * u.RNDprobab());
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| прыжок: X = 0.5 * proj_T(cB) * Jump * c                          |
//| MATLAB: sps = Best_pos / Time — матричное деление, скаляр        |
//|         sps = <cB,T> / <T,T>                                     |
//+------------------------------------------------------------------+
void C_AO_FOX::Jump(int i, double c)
  {
   double dotBT = 0.0, dotTT = 0.0, sumT = 0.0;

   for(int k = 0; k < coords; k++)
     {
      tm [k] = u.RNDprobab();
      dotBT += cB [k] * tm [k];
      dotTT += tm [k] * tm [k];
      sumT  += tm [k];
     }

   double sps  = dotTT > 0.0 ? dotBT / dotTT : 0.0;
   double tt   = sumT / coords;
   double t    = tt / 2.0;
   double jump = 0.5 * 9.81 * t * t;
   double x;

   for(int k = 0; k < coords; k++)
     {
      x = 0.5 * sps * tm [k] * jump * c;   // Distance_Fox_Rat * Jump * c
      a [i].c [k] = u.SeInDiSp(x, rangeMin [k], rangeMax [k], rangeStep [k]);
     }

   if(tt < minT)
      minT = tt;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| случайное блуждание: X = cB + randn * MinT * a                   |
//| при MinT = +inf координата уходит на границу (как у авторов)     |
//+------------------------------------------------------------------+
void C_AO_FOX::Walk(int i, double aCoef)
  {
   double x;

   for(int k = 0; k < coords; k++)
     {
      //x = cB [k] + Gauss() * minT * aCoef;
      x = u.PowerDistribution (cB [k], rangeMin [k], rangeMax [k], 20);
      a [i].c [k] = u.SeInDiSp(x, rangeMin [k], rangeMax [k], rangeStep [k]);
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                            Moving                                |
//+------------------------------------------------------------------+
void C_AO_FOX::Moving()
  {
//--- первый прогон: стартовая популяция
   if(!revision)
     {
      for(int i = 0; i < popSize; i++)
         for(int c = 0; c < coords; c++)
            a [i].c [c] = u.SeInDiSp(u.RNDfromCI(rangeMin [c], rangeMax [c]), rangeMin [c], rangeMax [c], rangeStep [c]);
      return;
     }

//--- одна итерация FOX = одна эпоха стенда
   double aCoef = 2.0 * (1.0 - (double)epochNow / epochs);

   for(int i = 0; i < popSize; i++)
     {
      if(u.RNDprobab() < pExploit)                // r >= 0.5
         Jump(i, u.RNDprobab() > pDir ? c1 : c2);
      else
         Walk(i, aCoef);
     }

   epochNow++;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                           Revision                               |
//| у FOX нет личной памяти и отбора: X перезаписывается целиком,    |
//| сохраняется только глобальный лучший                             |
//+------------------------------------------------------------------+
void C_AO_FOX::Revision()
  {
   for(int i = 0; i < popSize; i++)
     {
      if(a [i].f > fB)
        {
         fB = a [i].f;
         ArrayCopy(cB, a [i].c, 0, 0, coords);
        }
     }

   if(!revision)
      revision = true;
  }
//+------------------------------------------------------------------+