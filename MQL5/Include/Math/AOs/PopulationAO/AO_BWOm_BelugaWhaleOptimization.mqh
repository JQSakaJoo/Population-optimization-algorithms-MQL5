//+——————————————————————————————————————————————————————————————————+
//|                                                        C_AO_BWOm |
//|                                  Copyright 2007-2026, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

#include "#C_AO.mqh"

//+------------------------------------------------------------------+
//| Beluga Whale Optimization — МОДИФИЦИРОВАННАЯ версия.             |
//|                                                                  |
//| Общий принцип модификации: в оригинальной BWO каждая из трёх фаз |
//| использует случайные множители, разыгранные ОДИН раз на агента   |
//| (r1..r7 — скаляры), домножающие целые векторы. Любой такой       |
//| скаляр-на-вектор запирает кандидата в низкоразмерное подпрост-   |
//| ранство. На бенче из N копий одной функции оптимум осесимметри-  |
//| чен, и это подпространство совпадает с главной диагональю, на    |
//| которой оптимум и лежит — алгоритм «считывает диагональ» даром,  |
//| причём эффект УСИЛИВАЕТСЯ с размерностью: покоординатный шум     |
//| Леви усредняется при greedy-приёмке по среднему фитнесу.         |
//|                                                                  |
//| BWOm заменяет скаляры-на-агента на ПОКООРДИНАТНУЮ случайность во |
//| всех трёх фазах:                                                 |
//|   - разведка: r1, r2 на каждую пару; тяга к одноимённой паре     |
//|     соседа (без общего скаляра-опоры и кросс-перестановки);      |
//|   - эксплуатация: тяга к глобальному лучшему — покоординатная    |
//|     интерполяция (r3 на каждую координату); масштаб-скаляр       |
//|     r3*cB убран — именно он порождал луч origin->cB;             |
//|   - падение кита: r5, r6 покоординатные.                         |
//| Трёхфазная структура, балансовый фактор, вероятность падения     |
//| кита, полёт Леви и временное затухание C1/C2 сохранены.          |
//+------------------------------------------------------------------+
class C_AO_BWOm : public C_AO
 {
public:
                    ~C_AO_BWOm() {}
                     C_AO_BWOm()
   {
    ao_name = "BWOm";
    ao_desc = "Beluga Whale Optimization (modified)";
    ao_link = "https://www.mql5.com/ru/articles/22686";

    popSize = 20;
    KD      = 1.5;

    ArrayResize(params, 2);
    params [0].name = "popSize"; params [0].val  = popSize;
    params [1].name = "KD";      params [1].val  = KD;
   }

  void               SetParams()
   {
    popSize = (int)params [0].val;
    KD      = params      [1].val;

    //--- предохранители
    if(popSize < 1)
      popSize = 1;
    if(KD < 0.0)
      KD = 0.0;
   }

  bool               Init(const double &rangeMinP  [],
                          const double &rangeMaxP  [],
                          const double &rangeStepP [],
                          const int     epochsP = 0);

  void               Moving();
  void               Revision();

  //--- видимые параметры
  double             KD;           // интенсивность полёта Леви (оригинал: 0.05)

private:
  //--- snapshot для greedy acceptance
  double             snap_c  [];   // [popSize * coords] — плоский буфер
  double             snap_f  [];   // [popSize]

  //--- счётчик итераций (для отношения T/Max_it)
  int                epochs;       // всего эпох (из Init)
  int                epochsDen;    // знаменатель для T/Max_it
  int                epochNow;     // номер текущей «боевой» итерации

  //--- константа масштаба для полёта Леви (Mantegna, beta = 1.5)
  double             levySigma;

  //--- стандартное нормальное число (Box–Muller)
  double             Gauss();
  //--- один отсчёт полёта Леви по координате
  double             LevyStep();
 };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                              Init                                |
//+------------------------------------------------------------------+
bool C_AO_BWOm::Init(const double &rangeMinP  [],
                     const double &rangeMaxP  [],
                     const double &rangeStepP [],
                     const int     epochsP = 0)
 {
  if(!StandardInit(rangeMinP, rangeMaxP, rangeStepP))
    return false;

//--- начальная случайная популяция в cP[] (попадёт в c[] на первом Moving)
  for(int i = 0; i < popSize; i++)
    for(int c = 0; c < coords; c++)
      a [i].cP [c] = u.RNDfromCI(rangeMin [c], rangeMax [c]);

//--- внутренние буферы
  ArrayResize(snap_c, popSize * coords);
  ArrayResize(snap_f, popSize);

//--- счётчик итераций
  epochs    = epochsP;
  epochsDen = (epochs > 1) ? epochs - 1 : 1;
  epochNow  = 0;

//--- sigma для полёта Леви (Mantegna), beta = 1.5
//    значения гамма-функции посчитаны заранее, чтобы не зависеть
//    от наличия MathGamma:
//      Г(1 + beta)     = Г(2.5)  = 1.3293403881791
//      Г((1+beta)/2)   = Г(1.25) = 0.9064024770554
  double beta = 1.5;
  double g1   = 1.3293403881791;
  double g2   = 0.9064024770554;
  levySigma = MathPow((g1 * MathSin(M_PI * beta / 2.0)) /
                      (g2 * beta * MathPow(2.0, (beta - 1.0) / 2.0)),
                      1.0 / beta);              // ≈ 0.69655

  return true;
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|              Gauss — стандартное нормальное (Box–Muller)         |
//+------------------------------------------------------------------+
double C_AO_BWOm::Gauss()
 {
  double u1 = u.RNDfromCI(0.0, 1.0);
  double u2 = u.RNDfromCI(0.0, 1.0);
  if(u1 < 1e-15)
    u1 = 1e-15;
  return MathSqrt(-2.0 * MathLog(u1)) * MathCos(2.0 * M_PI * u2);
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|        LevyStep — один отсчёт полёта Леви (Mantegna)             |
//+------------------------------------------------------------------+
double C_AO_BWOm::LevyStep()
 {
  double uu = Gauss() * levySigma;
  double vv = Gauss();
  double av = MathAbs(vv);
  if(av < 1e-300)
    av = 1e-300;
  return KD * (uu / MathPow(av, 1.0 / 1.5));
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                            Moving                                |
//|                                                                  |
//|  Структура итерации:                                             |
//|    1. snapshot текущей популяции;                                |
//|    2. для каждого агента — балансовый фактор Bf и сосед RJ;      |
//|    3. выбор ветви:                                               |
//|         Bf <= WF   -> падение кита;                              |
//|         Bf >  0.5  -> разведка (честное парное плавание);        |
//|         иначе      -> эксплуатация (покоординатная тяга к cB);   |
//|    4. новый кандидат пишется в a[i].c.                           |
//|  Во всех ветвях случайность ПОКООРДИНАТНАЯ — нет скаляра на      |
//|  агента, домножающего целый вектор.                              |
//+------------------------------------------------------------------+
void C_AO_BWOm::Moving()
 {
//--- первый прогон: cP -> c, чтобы внешний цикл оценил FF
  if(!revision)
   {
    for(int i = 0; i < popSize; i++)
      for(int c = 0; c < coords; c++)
        a [i].c [c] = u.SeInDiSp(a [i].cP [c], rangeMin  [c], rangeMax  [c], rangeStep [c]);
    return;
   }

//--- номер текущей итерации и отношение T/Max_it
  epochNow++;
  double tRatio = (double)epochNow / (double)epochsDen;
  if(tRatio > 1.0)
    tRatio = 1.0;

//--- snapshot текущей популяции
  for(int i = 0; i < popSize; i++)
   {
    for(int c = 0; c < coords; c++)
      snap_c [i * coords + c] = a [i].c [c];
    snap_f [i] = a [i].f;
   }

//--- вероятность падения кита
  double WF = 0.1 - 0.05 * tRatio;

//--- движение
  for(int i = 0; i < popSize; i++)
   {
    //--- балансовый фактор Bf = B0 * (1 - T/2Tmax)
    double Bf = (1.0 - 0.5 * tRatio) * u.RNDprobab();

    //--- случайный сосед RJ != i
    int RJ = i;
    if(popSize > 1)
      while(RJ == i)
        RJ = u.RNDminusOne(popSize);

    //================================================================
    if(Bf <= WF)
     {
      //--- ФАЗА 3: ПАДЕНИЕ КИТА (whale fall) — модификация
      //    Оригинал: v = r5*xi - r6*xr + step, r5/r6 — скаляры на
      //    агента -> скелет заперт в плоскости (xi, xr), на стянутой
      //    диагонали — в линию. BWOm: r5, r6 покоординатные.
      //    r7 (магнитуда шага) оставлен на агента — равномерный
      //    масштаб шага не запирает кандидата.
      double r7    = u.RNDprobab();
      double C2    = 2.0 * (double)popSize * WF;
      double decay = MathExp(-C2 * tRatio);

      for(int c = 0; c < coords; c++)
       {
        double r5c  = u.RNDprobab();
        double r6c  = u.RNDprobab();
        double step = r7 * (rangeMax [c] - rangeMin [c]) * decay;
        double v    = r5c * snap_c [i  * coords + c] -
                      r6c * snap_c [RJ * coords + c] + step;
        a [i].c [c] = u.SeInDiSp(v, rangeMin [c], rangeMax [c], rangeStep [c]);
       }
     }
    //================================================================
    else
      if(Bf > 0.5)
       {
        //--- ФАЗА 1: РАЗВЕДКА — честное парное плавание (модификация)
        //    каждая пара координат (dA, dB) тянется к своей одноимённой
        //    паре соседа RJ; r1, r2 разыгрываются на КАЖДУЮ пару.
        //    Ни общего скаляра-опоры, ни кросс-перестановки измерений.

        //--- незатронутые координаты остаются прежними
        for(int c = 0; c < coords; c++)
          a [i].c [c] = snap_c [i * coords + c];

        int pairs = coords / 2;

        for(int j = 0; j < pairs; j++)
         {
          int dA = 2 * j;            // x-координата пары
          int dB = 2 * j + 1;        // y-координата пары

          double r1   = u.RNDprobab();
          double r2   = u.RNDprobab();
          double ang  = 2.0 * M_PI * r2;
          double sinA = MathSin(ang);
          double cosA = MathCos(ang);

          double xA   = snap_c [i  * coords + dA];
          double xB   = snap_c [i  * coords + dB];
          double refA = snap_c [RJ * coords + dA];
          double refB = snap_c [RJ * coords + dB];

          double vA = xA + (refA - xA) * (r1 + 1.0) * sinA;
          double vB = xB + (refB - xB) * (r1 + 1.0) * cosA;
          a [i].c [dA] = u.SeInDiSp(vA, rangeMin [dA], rangeMax [dA], rangeStep [dA]);
          a [i].c [dB] = u.SeInDiSp(vB, rangeMin [dB], rangeMax [dB], rangeStep [dB]);
         }
       }
      //================================================================
      else
       {
        //--- ФАЗА 2: ЭКСПЛУАТАЦИЯ (охота, полёт Леви) — модификация
        //    Оригинал: v = r3*cB - r4*xi + C1*LF*(xr-xi), r3 — скаляр
        //    на агента. Слагаемое r3*cB — масштабированный осесиммет-
        //    ричный глобальный лучший: кандидат скользит по лучу
        //    origin->cB (главная диагональ на бенче из копий функции).
        //    На больших размерностях покоординатный шум Леви усред-
        //    няется при greedy-приёмке по среднему фитнесу, и отбор
        //    оптимизирует фактически одномерное сканирование диагонали.
        //    BWOm: тяга к лучшему — покоординатная интерполяция (r3 на
        //    КАЖДУЮ координату), скаляр-масштаб r3*cB убран.
        //    C1 (затухание Леви) оставлен на агента: равномерный
        //    масштаб покоординатно-случайного вектора Леви не запирает
        //    кандидата в подпространство.
        double r4 = u.RNDprobab();
        double C1 = 2.0 * r4 * (1.0 - tRatio);

        for(int c = 0; c < coords; c++)
         {
          double r3c = u.RNDprobab();
          double xi  = snap_c [i  * coords + c];
          double xr  = snap_c [RJ * coords + c];
          double Lv  = LevyStep();
          double v   = xi + r3c * (cB [c] - xi) + C1 * Lv * (xr - xi);
          a [i].c [c] = u.SeInDiSp(v, rangeMin [c], rangeMax [c], rangeStep [c]);
         }
       }
   }
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                           Revision                               |
//|                                                                  |
//|  Стандартное обновление fB/cB + greedy acceptance: если новый    |
//|  фитнес хуже снапшота — откат координат и фитнеса агента.        |
//+------------------------------------------------------------------+
void C_AO_BWOm::Revision()
 {
//--- обновление личных и глобального лучших
  for(int i = 0; i < popSize; i++)
   {
    if(a [i].f > a [i].fB)
     {
      a [i].fB = a [i].f;
      ArrayCopy(a [i].cB, a [i].c, 0, 0, coords);
     }
    if(a [i].f > fB)
     {
      fB = a [i].f;
      ArrayCopy(cB, a [i].c, 0, 0, coords);
     }
   }

  if(!revision)
   {
    revision = true;
    return;
   }

//--- greedy acceptance: rollback если стало хуже
  for(int i = 0; i < popSize; i++)
   {
    if(a [i].f < snap_f [i])
     {
      for(int c = 0; c < coords; c++)
        a [i].c [c] = snap_c [i * coords + c];
      a [i].f = snap_f [i];
     }
   }
 }
//+------------------------------------------------------------------+