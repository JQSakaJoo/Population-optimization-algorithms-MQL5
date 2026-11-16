//+------------------------------------------------------------------+
//|                                  Copyright 2007-2026, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//+------------------------------------------------------------------+

#include "#C_AO.mqh"

//+------------------------------------------------------------------+
//| дерево                                                           |
//+------------------------------------------------------------------+
struct S_FOA_Tree
  {
   double            c [];   // координаты
   double            f;      // приспособленность
   int               age;    // возраст
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
struct S_FOA_Pop
  {
   S_FOA_Tree        t [];
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class C_AO_FOA : public C_AO
  {
public:
                    ~C_AO_FOA() {}
                     C_AO_FOA()
     {
      ao_name = "FOA";
      ao_desc = "Forest Optimization Algorithm";
      ao_link = "https://www.mql5.com/ru/articles/25876"; // ссылка на статью

      popSize      = 10;     // порция вычислений ФФ за эпоху
      areaLimit    = 20;     // предел размера леса
      lifeTime     = 30;     // предельный возраст дерева
      transferRate = 0.15;   // доля кандидатов для глобального посева
      lscRatio     = 0.0;    // LSC как доля размерности
      gscRatio     = 0.01;   // GSC как доля размерности
      dxRatio      = 0.4;    // радиус локального посева как доля диапазона

      ArrayResize(params, 7);
      params [0].name = "popSize";
      params [0].val = popSize;
      params [1].name = "areaLimit";
      params [1].val = areaLimit;
      params [2].name = "lifeTime";
      params [2].val = lifeTime;
      params [3].name = "transferRate";
      params [3].val = transferRate;
      params [4].name = "lscRatio";
      params [4].val = lscRatio;
      params [5].name = "gscRatio";
      params [5].val = gscRatio;
      params [6].name = "dxRatio";
      params [6].val = dxRatio;
     }

   void               SetParams()
     {
      popSize      = (int)params [0].val;
      areaLimit    = (int)params [1].val;
      lifeTime     = (int)params [2].val;
      transferRate = params [3].val;
      lscRatio     = params [4].val;
      gscRatio     = params [5].val;
      dxRatio      = params [6].val;

      if(popSize < 1)
         popSize = 1;
      if(areaLimit < 1)
         areaLimit = 1;
      if(lifeTime < 1)
         lifeTime = 1;
      if(transferRate < 0.0)
         transferRate = 0.0;
      if(transferRate > 1.0)
         transferRate = 1.0;
      if(lscRatio < 0.0)
         lscRatio = 0.0;
      if(gscRatio < 0.0)
         gscRatio = 0.0;
      if(gscRatio > 1.0)
         gscRatio = 1.0;
      if(dxRatio < 0.0)
         dxRatio = 0.0;
      if(dxRatio > 1.0)
         dxRatio = 1.0;

      params [0].val = popSize;
      params [1].val = areaLimit;
      params [2].val = lifeTime;
      params [3].val = transferRate;
      params [4].val = lscRatio;
      params [5].val = gscRatio;
      params [6].val = dxRatio;
     }

   bool               Init(const double &rangeMinP  [],
                           const double &rangeMaxP  [],
                           const double &rangeStepP [],
                           const int     epochsP = 0);

   void               Moving();
   void               Revision();

   //----------------------------------------------------------------
   int                areaLimit;
   int                lifeTime;
   double             transferRate;
   double             lscRatio;
   double             gscRatio;
   double             dxRatio;

   int                iterations;   // число завершённых итераций FOA за прогон

private:
   S_FOA_Pop          buf [2];      // лес: текущий буфер и буфер для пересборки
   int                cur;
   int                fN;           // число оценённых деревьев леса в buf [cur]

   int                LSC;          // локальных потомков на дерево
   int                GSC;          // заменяемых координат при глобальном посеве
   double             dx [];        // радиус локального посева по координатам

   //--- текущая итерация: деревья лежат в buf [cur].t [fN + g]
   bool               batchOpen;
   int                bSeeds;       // готовые деревья в начале итерации (глобальный посев / начальный лес)
   int                bLocal;       // локальных потомков по канону: pN * LSC
   int                bTotal;       // всего вычислений в итерации, кратно popSize
   int                bDone;        // уже оценено
   int                par [];       // индексы родителей локального посева
   int                pN;

   double             srt [][2];    // сортировка: [f, индекс]
   int                cand [];      // кандидаты (индексы)

   void               Reserve(int b, int size);
   void               CopyTree(S_FOA_Tree &dst, const S_FOA_Tree &src);
   void               InitBatch();
   void               NextBatch();
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                              Init                                |
//+------------------------------------------------------------------+
bool C_AO_FOA::Init(const double &rangeMinP  [],
                    const double &rangeMaxP  [],
                    const double &rangeStepP [],
                    const int     epochsP = 0)
  {
   if(!StandardInit(rangeMinP, rangeMaxP, rangeStepP))
      return false;

   LSC = (int)MathFloor(lscRatio * coords);
   if(LSC < 1)
      LSC = 1;
   GSC = (int)MathFloor(gscRatio * coords);
   if(GSC < 1)
      GSC = 1;

   ArrayResize(dx, coords);
   for(int k = 0; k < coords; k++)
      dx [k] = dxRatio * (rangeMax [k] - rangeMin [k]);

   cur        = 0;
   iterations = 0;

   InitBatch();
   batchOpen = true;
   return true;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
void C_AO_FOA::Reserve(int b, int size)
  {
   if(ArraySize(buf [b].t) < size)
      ArrayResize(buf [b].t, size, 256);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
void C_AO_FOA::CopyTree(S_FOA_Tree &dst, const S_FOA_Tree &src)
  {
   if(ArraySize(dst.c) != coords)
      ArrayResize(dst.c, coords);
   ArrayCopy(dst.c, src.c, 0, 0, coords);
   dst.f   = src.f;
   dst.age = src.age;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| начальный лес: случайные деревья возраста 0,                     |
//| areaLimit с округлением вверх до целого числа порций             |
//+------------------------------------------------------------------+
void C_AO_FOA::InitBatch()
  {
   int n = ((areaLimit + popSize - 1) / popSize) * popSize;

   fN = 0;
   Reserve(cur, n);

   for(int i = 0; i < n; i++)
     {
      if(ArraySize(buf [cur].t [i].c) != coords)
         ArrayResize(buf [cur].t [i].c, coords);

      for(int k = 0; k < coords; k++)
         buf [cur].t [i].c [k] = u.SeInDiSp(u.RNDfromCI(rangeMin [k], rangeMax [k]),
                                            rangeMin [k], rangeMax [k], rangeStep [k]);
      buf [cur].t [i].f   = -DBL_MAX;
      buf [cur].t [i].age = 0;
     }

   bSeeds = n;
   bLocal = 0;
   bTotal = n;
   bDone  = 0;
   pN     = 0;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| одна итерация FOA над полностью оценённым лесом:                 |
//| ограничение популяции -> глобальный посев -> лучшее дерево ->    |
//| план локального посева следующей итерации                        |
//+------------------------------------------------------------------+
void C_AO_FOA::NextBatch()
  {
   int nxt = cur ^ 1;
   int all = fN;

//--- 1. ограничение популяции -------------------------------------
   int keptN = 0;
   int candN = 0;

   ArrayResize(srt, all);
   if(ArraySize(cand) < all)
      ArrayResize(cand, all);

   for(int i = 0; i < all; i++)
     {
      if(buf [cur].t [i].age < lifeTime)
        {
         srt [keptN][0] = buf [cur].t [i].f;
         srt [keptN][1] = i;
         keptN++;
        }
      else
         cand [candN++] = i;         // старые — в кандидаты, в порядке следования
     }

   if(keptN == 0)                    // недостижимо при lifeTime >= 1; страховка
     {
      for(int i = 0; i < all; i++)
        {
         srt [i][0] = buf [cur].t [i].f;
         srt [i][1] = i;
        }
      keptN = all;
      candN = 0;
     }

   if(keptN == 0)                    // лес пуст — начать заново
     {
      InitBatch();
      return;
     }

   ArrayResize(srt, keptN);
   ArraySort(srt);                   // по возрастанию f; лучший — последний

   int K = keptN < areaLimit ? keptN : areaLimit;

   for(int r = K; r < keptN; r++)    // лишние — в кандидаты после старых, лучшие раньше
      cand [candN++] = (int)srt [keptN - 1 - r][1];

//--- 2. глобальный посев ------------------------------------------
   int n = (int)MathFloor(candN * transferRate);

   Reserve(nxt, K + n);

   for(int r = 0; r < K; r++)        // выжившие, отсортированы: лучший в [0]
      CopyTree(buf [nxt].t [r], buf [cur].t [(int)srt [keptN - 1 - r][1]]);

   for(int s = 0; s < n; s++)        // первые n кандидатов
     {
      CopyTree(buf [nxt].t [K + s], buf [cur].t [cand [s]]);

      for(int g = 0; g < GSC; g++)   // свой набор координат у каждого дерева
        {
         int k = u.RNDminusOne(coords);
         buf [nxt].t [K + s].c [k] = u.SeInDiSp(u.RNDfromCI(rangeMin [k], rangeMax [k]),
                                                rangeMin [k], rangeMax [k], rangeStep [k]);
        }
      buf [nxt].t [K + s].f   = -DBL_MAX;
      buf [nxt].t [K + s].age = 1;   // возраст 0 при посеве, +1 на шаге 4
     }

//--- 3. лучшее дерево не стареет ----------------------------------
   buf [nxt].t [0].age = 0;

//--- 4. локальный посев: родители — все деревья возраста 0 --------
   if(ArraySize(par) < K + n)
      ArrayResize(par, K + n);
   pN = 0;

   for(int r = 0; r < K; r++)
      if(buf [nxt].t [r].age == 0)
         par [pN++] = r;
   for(int s = 0; s < n; s++)
      par [pN++] = K + s;

   for(int r = 0; r < K; r++)        // стареют все, кроме новых потомков
      buf [nxt].t [r].age++;

//--- план итерации ------------------------------------------------
   cur    = nxt;
   fN     = K;
   bSeeds = n;
   bLocal = pN * LSC;
   bTotal = ((bSeeds + bLocal + popSize - 1) / popSize) * popSize;
   bDone  = 0;

   iterations++;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                            Moving                                |
//| очередная порция итерации: сначала готовые деревья, затем        |
//| локальные потомки, создаваемые по мере выдачи                    |
//+------------------------------------------------------------------+
void C_AO_FOA::Moving()
  {
   if(!batchOpen)
     {
      NextBatch();
      batchOpen = true;
     }

   for(int i = 0; i < popSize; i++)
     {
      int g = bDone + i;

      if(g < bSeeds)
        {
         ArrayCopy(a [i].c, buf [cur].t [fN + g].c, 0, 0, coords);
         continue;
        }

      //--- локальный потомок: LSC на родителя, остаток порции — по кругу
      int h  = g - bSeeds;
      int pi = h < bLocal ? h / LSC : (h - bLocal) % pN;

      ArrayCopy(a [i].c, buf [cur].t [par [pi]].c, 0, 0, coords);

      int    k = u.RNDminusOne(coords);
      double x = a [i].c [k] + u.RNDfromCI(-dx [k], dx [k]);
      a [i].c [k] = u.SeInDiSp(x, rangeMin [k], rangeMax [k], rangeStep [k]);
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                           Revision                               |
//| оценённые деревья дописываются в лес; итерация закрывается,      |
//| когда оценены все её деревья                                     |
//+------------------------------------------------------------------+
void C_AO_FOA::Revision()
  {
   Reserve(cur, fN + bDone + popSize);

   for(int i = 0; i < popSize; i++)
     {
      int g   = bDone + i;
      int idx = fN + g;

      if(g >= bSeeds)                // локальный потомок: сохранить координаты
        {
         if(ArraySize(buf [cur].t [idx].c) != coords)
            ArrayResize(buf [cur].t [idx].c, coords);
         ArrayCopy(buf [cur].t [idx].c, a [i].c, 0, 0, coords);
         buf [cur].t [idx].age = 0;
        }
      buf [cur].t [idx].f = a [i].f;

      if(a [i].f > fB)
        {
         fB = a [i].f;
         ArrayCopy(cB, a [i].c, 0, 0, coords);
        }
     }

   bDone += popSize;

   if(bDone >= bTotal)
     {
      fN       += bTotal;
      batchOpen = false;
     }

   if(!revision)
      revision = true;
  }
//+------------------------------------------------------------------+