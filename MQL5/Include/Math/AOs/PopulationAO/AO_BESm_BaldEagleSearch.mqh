//+------------------------------------------------------------------+
//|                                                        C_AO_BESm |
//|                                  Copyright 2007-2026, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//+------------------------------------------------------------------+

#include "#C_AO.mqh"

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
struct S_BES_Coord
  {
   double            wid;     // ширина диапазона
   double            mean;    // центр популяции (нормированный), фиксируется на эпоху
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| BES                                                              |
//+------------------------------------------------------------------+
class C_AO_BESm : public C_AO
  {
public:
                    ~C_AO_BESm() {}
                     C_AO_BESm()
     {
      ao_name = "BESm";
      ao_desc = "Bald Eagle Search M";
      ao_link = "https://www.mql5.com/ru/articles/26062";

      popSize   = 10;     // размер популяции
      alpha     = 2.0;    // α ∈ (1.5, 2)
      angleK    = 10.0;   // a: число витков спирали, (5...10)
      radR      = 1.5;    // R: добавка к радиусу спирали, (0.5...2)
      c1        = 0.5;    // вес центра в пикировании
      c2        = 2.0;    // вес лучшего в пикировании
      rndProbab = 0.97;   // вероятность случайного положения

      ArrayResize(params, 7);
      params [0].name = "popSize";
      params [0].val = popSize;
      params [1].name = "alpha";
      params [1].val = alpha;
      params [2].name = "a";
      params [2].val = angleK;
      params [3].name = "R";
      params [3].val = radR;
      params [4].name = "c1";
      params [4].val = c1;
      params [5].name = "c2";
      params [5].val = c2;
      params [6].name = "rndProbab";
      params [6].val = rndProbab;
     }

   void               SetParams()
     {
      popSize   = (int)params [0].val;
      alpha     = params [1].val;
      angleK    = params [2].val;
      radR      = params [3].val;
      c1        = params [4].val;
      c2        = params [5].val;
      rndProbab = params [6].val;

      if(popSize < 2)
         popSize = 2;          // фазе поиска нужен партнёр
      if(alpha < 0.0)
         alpha = 0.0;
      if(angleK < 0.0)
         angleK = 0.0;
      if(radR < 0.0)
         radR = 0.0;

      //--- обратная запись, чтобы печать параметров совпадала с прогоном
      params [0].val = popSize;
      params [1].val = alpha;
      params [2].val = angleK;
      params [3].val = radR;
      params [4].val = c1;
      params [5].val = c2;
      params [6].val = rndProbab;
     }

   bool               Init(const double &rangeMinP  [],
                           const double &rangeMaxP  [],
                           const double &rangeStepP [],
                           const int     epochsP = 0);

   void               Moving();
   void               Revision();

   //----------------------------------------------------------------
   double             alpha;
   double             angleK;
   double             radR;
   double             c1;
   double             c2;
   double             rndProbab;

private:
   int                phase;      // 0 выбор, 1 поиск, 2 пикирование
   int                epochs;     // число эпох прогона (для BES_STAT)
   int                revCnt;     // счётчик вызовов Revision
   S_BES_Coord        cr   [];    // нормирование и центр, по координатам

   double             Norm(const double x, const int j)
     {
      if(u.RNDprobab() < rndProbab)
         return (cB [j] - x) / cr [j].wid;
      else
         return u.RNDfromCI(-1.0, 1.0);
     }
   double             Denorm(const double v, const int j) { return cB [j] - v * cr [j].wid; }

   void               CalcMean();
   int                Partner(const int i);
   void               Spiral(double &x, double &y);
   void               Swoop(double &x1, double &y1);
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                              Init                                |
//+------------------------------------------------------------------+
bool C_AO_BESm::Init(const double &rangeMinP  [],
                     const double &rangeMaxP  [],
                     const double &rangeStepP [],
                     const int     epochsP = 0)
  {
   if(!StandardInit(rangeMinP, rangeMaxP, rangeStepP))
      return false;

   ArrayResize(cr, coords);
   for(int j = 0; j < coords; j++)
     {
      cr [j].wid  = rangeMax [j] - rangeMin [j];
      cr [j].mean = 0.0;
     }

   phase  = 0;
   epochs = epochsP;
   revCnt = 0;

   return true;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| центр популяции по принятым позициям, в нормированных координатах|
//+------------------------------------------------------------------+
void C_AO_BESm::CalcMean()
  {
   double s;

   for(int j = 0; j < coords; j++)
     {
      s = 0.0;
      for(int i = 0; i < popSize; i++)
         s += a [i].cB [j];
      cr [j].mean = Norm(s / popSize, j);
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| случайный партнёр k != i                                         |
//+------------------------------------------------------------------+
int C_AO_BESm::Partner(const int i)
  {
   int k = u.RNDminusOne(popSize - 1);   // [0; popSize-2]
   if(k >= i)
      k++;
   return k;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| polr: θ = a·π·rand, r = θ + R·rand, dx = r·sinθ, dy = r·cosθ     |
//| нормировка по max|·| среди popSize свежих значений, берётся [0]  |
//+------------------------------------------------------------------+
void C_AO_BESm::Spiral(double &x, double &y)
  {
   double th, r, dx, dy;
   double dx0 = 0.0, dy0 = 0.0, mx = 0.0, my = 0.0;

   for(int k = 0; k < popSize; k++)
     {
      th = angleK * M_PI * u.RNDprobab();
      r  = th + radR * u.RNDprobab();
      dx = r * MathSin(th);
      dy = r * MathCos(th);

      if(k == 0)
        {
         dx0 = dx;
         dy0 = dy;
        }
      if(MathAbs(dx) > mx)
         mx = MathAbs(dx);
      if(MathAbs(dy) > my)
         my = MathAbs(dy);
     }

   x = mx > 0.0 ? dx0 / mx : 0.0;
   y = my > 0.0 ? dy0 / my : 0.0;
  }
//+------------------------------------------------------------------+

//+---------------------------------------------------------------------+
//| swoo_p: dx = θ·sinh θ, dy = θ·cosh θ, нормировка по максимуму       |
//| θ >= 0, обе функции монотонны => максимум при θmax. Отношение       |
//| считается без переполнения (тождественно прямому расчёту):          |
//| θ0 sinhθ0 / (θm sinhθm) = θ0/θm · e^(θ0-θm) · (1-e^-2θ0)/(1-e^-2θm) |
//| θ0 coshθ0 / (θm coshθm) = θ0/θm · e^(θ0-θm) · (1+e^-2θ0)/(1+e^-2θm) |
//+---------------------------------------------------------------------+
void C_AO_BESm::Swoop(double &x1, double &y1)
  {
   double th, th0 = 0.0, thM = 0.0;

   for(int k = 0; k < popSize; k++)
     {
      th = angleK * M_PI * u.RNDprobab();
      if(k == 0)
         th0 = th;
      if(th > thM)
         thM = th;
     }

   if(thM <= 0.0)
     {
      x1 = 0.0;
      y1 = 0.0;
      return;
     }

   double e  = (th0 / thM) * MathExp(th0 - thM);
   double p0 = MathExp(-2.0 * th0);
   double pM = MathExp(-2.0 * thM);

   x1 = e * (1.0 - p0) / (1.0 - pM);
   y1 = e * (1.0 + p0) / (1.0 + pM);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                            Moving                                |
//+------------------------------------------------------------------+
void C_AO_BESm::Moving()
  {
//--- первый прогон: стартовая популяция
   if(!revision)
     {
      for(int i = 0; i < popSize; i++)
         for(int j = 0; j < coords; j++)
            a [i].c [j] = u.SeInDiSp(u.RNDfromCI(rangeMin [j], rangeMax [j]), rangeMin [j], rangeMax [j], rangeStep [j]);
      return;
     }

   CalcMean();

   double xn, bn, kn, v;

//--- phase 0: выбор пространства
//    u_new = u_best + alpha·q·(u_mean − u_i),  q — по координатам
   if(phase == 0)
     {
      for(int i = 0; i < popSize; i++)
        {
         for(int j = 0; j < coords; j++)
           {
            bn = Norm(cB [j], j);
            xn = Norm(a [i].cB [j], j);
            v  = bn + alpha * u.RNDprobab() * (cr [j].mean - xn);
            a [i].c [j] = u.SeInDiSp(Denorm(v, j), rangeMin [j], rangeMax [j], rangeStep [j]);
           }
        }
      return;
     }

//--- phase 1: поиск по спирали
//    u_new = u_i + y·(u_i − u_k) + x·(u_i − u_mean),  x, y — скаляры на агента
   if(phase == 1)
     {
      double x, y;
      int    k;

      for(int i = 0; i < popSize; i++)
        {
         k = Partner(i);
         Spiral(x, y);

         for(int j = 0; j < coords; j++)
           {
            xn = Norm(a [i].cB [j], j);
            kn = Norm(a [k].cB [j], j);
            v  = xn + y * (xn - kn) + x * (xn - cr [j].mean);
            a [i].c [j] = u.SeInDiSp(Denorm(v, j), rangeMin [j], rangeMax [j], rangeStep [j]);
           }
        }
      return;
     }

//--- phase 2: пикирование
//    u_new = q·u_best + x1·(u_i − c1·u_mean) + y1·(u_i − c2·u_best),
//    x1, y1 — скаляры на агента, q — по координатам
   double x1, y1;

   for(int i = 0; i < popSize; i++)
     {
      Swoop(x1, y1);

      for(int j = 0; j < coords; j++)
        {
         bn = Norm(cB [j], j);
         xn = Norm(a [i].cB [j], j);
         v  = u.RNDprobab() * bn + x1 * (xn - c1 * cr [j].mean) + y1 * (xn - c2 * bn);
         a [i].c [j] = u.SeInDiSp(Denorm(v, j), rangeMin [j], rangeMax [j], rangeStep [j]);
        }
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                           Revision                               |
//| стенд максимизирует: у авторов newsol.cost < pop.cost -> f > fB  |
//+------------------------------------------------------------------+
void C_AO_BESm::Revision()
  {
   revCnt++;

//--- стартовая популяция: принимается целиком
   if(!revision)
     {
      for(int i = 0; i < popSize; i++)
        {
         a [i].fB = a [i].f;
         ArrayCopy(a [i].cB, a [i].c, 0, 0, coords);

         if(a [i].fB > fB)
           {
            fB = a [i].fB;
            ArrayCopy(cB, a [i].cB, 0, 0, coords);
           }
        }

      phase    = 0;
      revision = true;
      return;
     }

//--- жадный отбор против своей принятой позиции
   for(int i = 0; i < popSize; i++)
     {
      if(a [i].f > a [i].fB)
        {
         a [i].fB = a [i].f;
         ArrayCopy(a [i].cB, a [i].c, 0, 0, coords);

         if(a [i].fB > fB)
           {
            fB = a [i].fB;
            ArrayCopy(cB, a [i].cB, 0, 0, coords);
           }
        }
     }

   phase = (phase + 1) % 3;
  }
//+------------------------------------------------------------------+
