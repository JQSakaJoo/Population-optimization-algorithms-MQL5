//+——————————————————————————————————————————————————————————————————+
//|                                                        C_AO_OOAm |
//|                                  Copyright 2007-2026, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

#include "#C_AO.mqh"


//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class C_AO_OOAm : public C_AO
  {
public:
                    ~C_AO_OOAm() {}
                     C_AO_OOAm() : P_MIN(1.0), P_MAX(10.0)
     {
      ao_name = "OOAm";
      ao_desc = "Osprey Optimization Algorithm M";
      ao_link = "https://www.mql5.com/ru/articles/24779";

      popSize  = 50;    // размер популяции
      pertProb = 0.1;   // фаза 1: вероятность возмущения координаты
      copyProb = 0.5;   // фаза 0: вероятность взять координату рыбы как есть

      ArrayResize(params, 3);
      params [0].name = "popSize";
      params [0].val  = popSize;
      params [1].name = "pertProb";
      params [1].val  = pertProb;
      params [2].name = "copyProb";
      params [2].val  = copyProb;
     }

   void               SetParams()
     {
      popSize  = (int)params [0].val;
      pertProb =      params [1].val;
      copyProb =      params [2].val;

      if(popSize < 1)
         popSize = 1;

      if(pertProb < 0.0)
         pertProb = 0.0;
      if(pertProb > 1.0)
         pertProb = 1.0;

      if(copyProb < 0.0)
         copyProb = 0.0;
      if(copyProb > 1.0)
         copyProb = 1.0;

      params [0].val = popSize;
      params [1].val = pertProb;
      params [2].val = copyProb;
     }

   bool               Init(const double &rangeMinP  [],
                           const double &rangeMaxP  [],
                           const double &rangeStepP [],
                           const int     epochsP = 0);

   void               Moving();
   void               Revision();

   //--- видимые параметры
   double             pertProb;
   double             copyProb;

private:
   //--- расписание степени PowerDistribution (найдено, зашито)
   const double       P_MIN;       // старт: равномерно до границы
   const double       P_MAX;       // финиш: медиана шага 0.5^P_MAX от расстояния до границы

   int                phase;       // 0 — разведка; 1 — эксплуатация
   int                epochs;      // бюджет стенда в эпохах
   int                epochNow;    // выполненных эпох (расписание p)
   int                idx [];      // буфер индексов агентов лучше текущего

   void               MakeExploration(int i);
   void               MakeExploitation(int i);
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                              Init                                |
//+------------------------------------------------------------------+
bool C_AO_OOAm::Init(const double &rangeMinP  [],
                     const double &rangeMaxP  [],
                     const double &rangeStepP [],
                     const int     epochsP = 0)
  {
   if(!StandardInit(rangeMinP, rangeMaxP, rangeStepP))
      return false;

   phase    = 0;
   epochs   = (epochsP > 0) ? epochsP : 1;
   epochNow = 0;

   ArrayResize(idx, popSize);

   return true;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   MakeExploration — фаза 0 (M1 + M4).                            |
//|   Рыба равновероятно из {агенты лучше i} U {xbest}.              |
//|   Координата: с вероятностью copyProb — копия координаты рыбы,   |
//|   иначе x_j + r_j*(fish_j - x_j).                                |
//+------------------------------------------------------------------+
void C_AO_OOAm::MakeExploration(int i)
  {
   int nBetter = 0;

   for(int k = 0; k < popSize; k++)
      if(a [k].fB > a [i].fB)
         idx [nBetter++] = k;

   int fish = u.RNDminusOne(nBetter + 1);   // == nBetter -> xbest
   double f, x;

   for(int c = 0; c < coords; c++)
     {
      f = (fish == nBetter) ? cB [c] : a [idx [fish]].cB [c];

      if(u.RNDprobab() < copyProb)
         x = f;
      else
         x = a [i].cB [c] + u.RNDprobab() * (f - a [i].cB [c]);

      a [i].c [c] = u.SeInDiSp(x, rangeMin [c], rangeMax [c], rangeStep [c]);
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   MakeExploitation — фаза 1 (M3).                                |
//|   Каждая координата с вероятностью pertProb получает шаг         |
//|   PowerDistribution со степенью p, растущей от P_MIN к P_MAX по  |
//|   доле бюджета; хотя бы одна координата возмущается всегда.      |
//+------------------------------------------------------------------+
void C_AO_OOAm::MakeExploitation(int i)
  {
   double frac = (double)epochNow / (double)epochs;
   if(frac > 1.0)
      frac = 1.0;

   double p     = P_MIN + (P_MAX - P_MIN) * frac;
   int    cLast = -1;

   for(int c = 0; c < coords; c++)
     {
      if(u.RNDprobab() < pertProb)
        {
         a [i].c [c] = u.SeInDiSp(u.PowerDistribution(a [i].cB [c], rangeMin [c], rangeMax [c], p),
                                  rangeMin [c], rangeMax [c], rangeStep [c]);
         cLast = c;
        }
      else
         a [i].c [c] = a [i].cB [c];
     }

   if(cLast < 0)
     {
      int c = u.RNDminusOne(coords);
      a [i].c [c] = u.SeInDiSp(u.PowerDistribution(a [i].cB [c], rangeMin [c], rangeMax [c], p),
                               rangeMin [c], rangeMax [c], rangeStep [c]);
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                            Moving                                |
//+------------------------------------------------------------------+
void C_AO_OOAm::Moving()
  {
//--- первый прогон: стартовая популяция
   if(!revision)
     {
      for(int i = 0; i < popSize; i++)
         for(int c = 0; c < coords; c++)
            a [i].c [c] = u.SeInDiSp(u.RNDfromCI(rangeMin [c], rangeMax [c]),
                                     rangeMin [c], rangeMax [c], rangeStep [c]);
      return;
     }

   if(phase == 0)
     {
      for(int i = 0; i < popSize; i++)
         MakeExploration(i);
     }
   else
     {
      for(int i = 0; i < popSize; i++)
         MakeExploitation(i);
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                           Revision                               |
//+------------------------------------------------------------------+
void C_AO_OOAm::Revision()
  {
//--- глобальный лучший
   for(int i = 0; i < popSize; i++)
     {
      if(a [i].f > fB)
        {
         fB = a [i].f;
         ArrayCopy(cB, a [i].c, 0, 0, coords);
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

      phase    = 0;
      revision = true;
      return;
     }

//--- жадная приёмка (строго лучше)
   for(int i = 0; i < popSize; i++)
     {
      if(a [i].f > a [i].fB)
        {
         a [i].fB = a [i].f;
         ArrayCopy(a [i].cB, a [i].c, 0, 0, coords);
        }
     }

   phase = 1 - phase;
  }
//+------------------------------------------------------------------+