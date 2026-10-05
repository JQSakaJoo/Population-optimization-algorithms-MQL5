//+——————————————————————————————————————————————————————————————————+
//|                                                         C_AO_FPA |
//|                                  Copyright 2007-2026, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

#include "#C_AO.mqh"

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class C_AO_FPA : public C_AO
  {
public:
                    ~C_AO_FPA() {}
                     C_AO_FPA()
     {
      ao_name = "FPA";
      ao_desc = "Flower Pollination Algorithm";
      ao_link = "https://www.mql5.com/ru/articles/25138";

      popSize   = 10;     // размер популяции
      pSwitch   = 0.2;    // вероятность локального шага
      levyScale = 0.5;    // масштаб шага Леви
      fromCode  = false;  // база локального шага

      ArrayResize(params, 4);
      params [0].name = "popSize";
      params [0].val = popSize;
      params [1].name = "pSwitch";
      params [1].val = pSwitch;
      params [2].name = "levyScale";
      params [2].val = levyScale;
      params [3].name = "fromCode";
      params [3].val = fromCode ? 1.0 : 0.0;
     }

   void               SetParams()
     {
      popSize   = (int)params [0].val;
      pSwitch   = params [1].val;
      levyScale = params [2].val;
      fromCode  = params [3].val != 0.0;

      if(popSize < 2)
         popSize = 2;   // локальному шагу нужны два разных агента
      if(pSwitch < 0.0)
         pSwitch = 0.0;
      if(pSwitch > 1.0)
         pSwitch = 1.0;
      if(levyScale < 0.0)
         levyScale = 0.0;

      params [0].val = popSize;
      params [1].val = pSwitch;
      params [2].val = levyScale;
      params [3].val = fromCode ? 1.0 : 0.0;
     }

   bool               Init(const double &rangeMinP  [],
                           const double &rangeMaxP  [],
                           const double &rangeStepP [],
                           const int     epochsP = 0);

   void               Moving();
   void               Revision();

   //----------------------------------------------------------------
   double             pSwitch;     // вероятность локального опыления
   double             levyScale;   // масштаб шага Леви
   bool               fromCode;    // база локального шага: true — a[i].c, false — a[i].cB

private:
   void               GlobalPollination(int i);
   void               LocalPollination(int i);

   double             Levy();
   double             Gauss();
   double             RndOpen();
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                              Init                                |
//+------------------------------------------------------------------+
bool C_AO_FPA::Init(const double &rangeMinP  [],
                    const double &rangeMaxP  [],
                    const double &rangeStepP [],
                    const int     epochsP = 0)
  {
   if(!StandardInit(rangeMinP, rangeMaxP, rangeStepP))
      return false;

   return true;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| равномерное на (0; 1]                                            |
//+------------------------------------------------------------------+
double C_AO_FPA::RndOpen()
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
double C_AO_FPA::Gauss()
  {
   return sqrt(-2.0 * log(RndOpen())) * cos(2.0 * M_PI * u.RNDprobab());
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| шаг Леви, алгоритм Мантенья: u / |v|^(1/beta)                    |
//| v = 0 возможно (RndOpen = 1 даёт нулевой радиус) — перевыборка   |
//+------------------------------------------------------------------+
double C_AO_FPA::Levy()
  {
   double uu = Gauss() * 0.6965745025576968;
   double v;
   do
      v = Gauss();
   while(v == 0.0);

   return uu / pow(fabs(v), 1.0 / 1.5);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| глобальное опыление: s_i = x_i + levyScale * L (.) (x_i - g)     |
//+------------------------------------------------------------------+
void C_AO_FPA::GlobalPollination(int i)
  {
   double x;

   for(int c = 0; c < coords; c++)
     {
      x = a [i].cB [c] + levyScale * Levy() * (a [i].cB [c] - cB [c]);
      a [i].c [c] = u.SeInDiSp(x, rangeMin [c], rangeMax [c], rangeStep [c]);
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| локальное опыление: s_i = base_i + eps * (x_j - x_k)             |
//+------------------------------------------------------------------+
void C_AO_FPA::LocalPollination(int i)
  {
   double eps = u.RNDprobab();

   int j = u.RNDminusOne(popSize);
   int k;
   do
      k = u.RNDminusOne(popSize);
   while(k == j);

   double base, x;

   for(int c = 0; c < coords; c++)
     {
      base = fromCode ? a [i].c [c] : a [i].cB [c];
      x    = base + eps * (a [j].cB [c] - a [k].cB [c]);
      a [i].c [c] = u.SeInDiSp(x, rangeMin [c], rangeMax [c], rangeStep [c]);
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                            Moving                                |
//+------------------------------------------------------------------+
void C_AO_FPA::Moving()
  {
//--- первый прогон: стартовая популяция
   if(!revision)
     {
      for(int i = 0; i < popSize; i++)
         for(int c = 0; c < coords; c++)
            a [i].c [c] = u.SeInDiSp(u.RNDfromCI(rangeMin [c], rangeMax [c]), rangeMin [c], rangeMax [c], rangeStep [c]);
      return;
     }

//--- одна итерация FPA = одна эпоха стенда
   for(int i = 0; i < popSize; i++)
     {
      if(u.RNDprobab() > pSwitch)
         GlobalPollination(i);
      else
         LocalPollination(i);
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                           Revision                               |
//+------------------------------------------------------------------+
void C_AO_FPA::Revision()
  {
//--- сохранить лучшее решение
   for(int i = 0; i < popSize; i++)
     {
      if(a [i].f > fB)
        {
         fB = a [i].f;
         ArrayCopy(cB, a [i].c, 0, 0, coords);
        }
     }

//--- первый проход: a[i].fB = -DBL_MAX, стартовая популяция принята целиком
   if(!revision)
      revision = true;
  }
//+------------------------------------------------------------------+
