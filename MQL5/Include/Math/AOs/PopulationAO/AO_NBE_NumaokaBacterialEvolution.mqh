//+——————————————————————————————————————————————————————————————————+
//|                                                         C_AO_NBE |
//|                                  Copyright 2007-2026, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

#include "#C_AO.mqh"

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
struct S_Bacteroid
 {
  double c           [];     // генотип (валидные координаты)
  double             f;      // фитнес
  double             E;      // энергия

  void               Init(int coords)
   {
    ArrayResize(c, coords);
    f = -DBL_MAX;
    E = 0.0;
   }
 };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class C_AO_NBE : public C_AO
 {
public:
                    ~C_AO_NBE() {}
                     C_AO_NBE()
   {
    ao_name = "NBE";
    ao_desc = "Numaoka Bacterial Evolution";
    ao_link = "https://www.mql5.com/ru/articles/23267";

    popSize   = 30;     // число бактероидов
    stepSize  = 0.1;    // стартовый масштаб движения (доля диапазона)
    pPlasmid  = 0.1;    // вероятность плазмидного переноса на бактероид/эпоху
    transFrac = 0.2;    // длина плазмиды как доля coords (0..1]
    cost      = 0.5;    // метаболическая стоимость (приход в [0,1])
    pDeath    = 0.0005; // «средовая» вероятность гибели на бактероид/эпоху

    ArrayResize(params, 6);
    params [0].name = "popSize";   params [0].val = popSize;
    params [1].name = "stepSize";  params [1].val = stepSize;
    params [2].name = "pPlasmid";  params [2].val = pPlasmid;
    params [3].name = "transFrac"; params [3].val = transFrac;
    params [4].name = "cost";      params [4].val = cost;
    params [5].name = "pDeath";    params [5].val = pDeath;
   }

  void               SetParams()
   {
    popSize   = (int)params [0].val;
    stepSize  =      params [1].val;
    pPlasmid  =      params [2].val;
    transFrac =      params [3].val;
    cost      =      params [4].val;
    pDeath    =      params [5].val;

    //--- предохранители
    if(popSize   < 1)
      popSize   = 1;
    if(stepSize  <= 0.0)
      stepSize  = 0.001;
    if(pPlasmid  < 0.0)
      pPlasmid  = 0.0;
    if(pPlasmid  > 1.0)
      pPlasmid  = 1.0;
    if(transFrac <= 0.0)
      transFrac = 0.01;
    if(transFrac >  1.0)
      transFrac = 1.0;
    if(cost      < 0.0)
      cost      = 0.0;
    if(pDeath    < 0.0)
      pDeath    = 0.0;
    if(pDeath    > 1.0)
      pDeath    = 1.0;
   }

  bool               Init(const double &rangeMinP  [],
                          const double &rangeMaxP  [],
                          const double &rangeStepP [],
                          const int     epochsP = 0);

  void               Moving();
  void               Revision();

  //--- видимые параметры
  double             stepSize;    // стартовый масштаб движения
  double             pPlasmid;    // вероятность переноса
  double             transFrac;   // доля coords для плазмиды
  double             cost;        // метаболическая стоимость
  double             pDeath;      // средовая гибель

private:
  //--- данные (массивы структур)
  S_Bacteroid        bro    [];   // [popSize] — бактероиды
  int                actType [];  // [popSize] 0=движение, 1=перенос (на эпоху)
  int                allIdx [];   // [popSize] — 0..popSize-1
  int                living [];   // [popSize] — живые в текущем Revision

  int                transLen;    // длина плазмиды (из transFrac)

  //--- прогресс прогона (для отжига шага)
  int                epochsDen;
  int                epochNow;

  //--- константы
  double             E0;          // стартовая энергия
  double             Emax;        // потолок энергии
  double             STEP_MIN_RATIO;

  //--- вспомогательные
  double             Gauss();
  double             StepNow();
  int                RouletteByE(int &pool [], int n, int exclude);
  void               Respawn(int i);
 };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                              Init                                |
//+------------------------------------------------------------------+
bool C_AO_NBE::Init(const double &rangeMinP  [],
                    const double &rangeMaxP  [],
                    const double &rangeStepP [],
                    const int     epochsP = 0)
 {
  if(!StandardInit(rangeMinP, rangeMaxP, rangeStepP))
    return false;

//--- длина плазмиды
  transLen = (int)MathRound(transFrac * coords);
  if(transLen < 1)
    transLen = 1;
  if(transLen > coords)
    transLen = coords;

//--- прогресс / константы
  epochsDen      = (epochsP > 1) ? epochsP - 1 : 1;
  epochNow       = 0;
  E0             = 1.0;
  Emax           = 5.0;
  STEP_MIN_RATIO = 0.025;

//--- буферы
  ArrayResize(bro,     popSize);
  ArrayResize(actType, popSize);
  ArrayResize(allIdx,  popSize);
  ArrayResize(living,  popSize);

  for(int i = 0; i < popSize; i++)
   {
    bro [i].Init(coords);
    allIdx [i] = i;
   }

//--- стартовая популяция
  for(int i = 0; i < popSize; i++)
    Respawn(i);

  return true;
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   Respawn — случайный генотип, стартовая энергия                 |
//+------------------------------------------------------------------+
void C_AO_NBE::Respawn(int i)
 {
  for(int c = 0; c < coords; c++)
    bro [i].c [c] = u.SeInDiSp(u.RNDfromCI(rangeMin [c], rangeMax [c]),
                               rangeMin [c], rangeMax [c], rangeStep [c]);
  bro [i].f = -DBL_MAX;
  bro [i].E = E0;
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|         Gauss — стандартное нормальное (Box–Muller)              |
//+------------------------------------------------------------------+
double C_AO_NBE::Gauss()
 {
  double u1 = u.RNDfromCI(0.0, 1.0);
  double u2 = u.RNDfromCI(0.0, 1.0);
  if(u1 < 1e-15)
    u1 = 1e-15;
  return MathSqrt(-2.0 * MathLog(u1)) * MathCos(2.0 * M_PI * u2);
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   StepNow — масштаб движения с линейным отжигом                 |
//+------------------------------------------------------------------+
double C_AO_NBE::StepNow()
 {
  double t = (double)epochNow / (double)epochsDen;
  if(t > 1.0)
    t = 1.0;
  return stepSize * (STEP_MIN_RATIO + (1.0 - STEP_MIN_RATIO) * (1.0 - t));
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   RouletteByE — индекс из pool[0..n-1], взвешенный по энергии,   |
//|   с исключением exclude (или -1)                                 |
//+------------------------------------------------------------------+
int C_AO_NBE::RouletteByE(int &pool [], int n, int exclude)
 {
  double sum = 0.0;
  for(int k = 0; k < n; k++)
   {
    int idx = pool [k];
    if(idx == exclude)
      continue;
    double e = bro [idx].E;
    if(e < 1e-12)
      e = 1e-12;
    sum += e;
   }

  if(sum <= 0.0)
   {
    int idx = pool [u.RNDminusOne(n)];
    return idx;
   }

  double r   = u.RNDprobab() * sum;
  double acc = 0.0;
  for(int k = 0; k < n; k++)
   {
    int idx = pool [k];
    if(idx == exclude)
      continue;
    double e = bro [idx].E;
    if(e < 1e-12)
      e = 1e-12;
    acc += e;
    if(acc >= r)
      return idx;
   }
  return pool [n - 1];
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                            Moving                                |
//|   Каждый бактероид этой эпохой делает ЛИБО плазмидный перенос,   |
//|   ЛИБО хемотаксис-движение. Кандидат пишется в a[i].c.           |
//+------------------------------------------------------------------+
void C_AO_NBE::Moving()
 {
//--- первый прогон: генотипы -> слоты для первичной оценки
  if(!revision)
   {
    for(int i = 0; i < popSize; i++)
      for(int c = 0; c < coords; c++)
        a [i].c [c] = bro [i].c [c];
    return;
   }

  epochNow++;
  double step = StepNow();

  for(int i = 0; i < popSize; i++)
   {
    //--- стартуем кандидата с текущего генотипа
    for(int c = 0; c < coords; c++)
      a [i].c [c] = bro [i].c [c];

    if(popSize > 1 && u.RNDprobab() < pPlasmid)
     {
      //--- ПЛАЗМИДНЫЙ ПЕРЕНОС: блок от донора (рулетка по энергии)
      actType [i] = 1;

      int donor = RouletteByE(allIdx, popSize, i);
      int span  = coords - transLen + 1;
      if(span < 1)
        span = 1;
      int start = u.RNDminusOne(span);

      for(int t = 0; t < transLen; t++)
       {
        int c = start + t;
        if(c < coords)
          a [i].c [c] = bro [donor].c [c];
       }
     }
    else
     {
      //--- ХЕМОТАКСИС: покоординатное гауссово движение
      actType [i] = 0;

      for(int c = 0; c < coords; c++)
       {
        double rng = rangeMax [c] - rangeMin [c];
        double v   = bro [i].c [c] + Gauss() * step * rng;
        a [i].c [c] = u.SeInDiSp(v, rangeMin [c], rangeMax [c], rangeStep [c]);
       }
     }
   }
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                           Revision                               |
//|   Приёмка хода -> энергобаланс -> гибель -> деление.             |
//+------------------------------------------------------------------+
void C_AO_NBE::Revision()
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

//--- первичная фиксация
  if(!revision)
   {
    for(int i = 0; i < popSize; i++)
     {
      bro [i].f = a [i].f;
      bro [i].E = E0;
     }
    revision = true;
    return;
   }

//--- приёмка хода: перенос безусловно, движение — greedy
  for(int i = 0; i < popSize; i++)
   {
    if(actType [i] == 1)
     {
      for(int c = 0; c < coords; c++)
        bro [i].c [c] = a [i].c [c];
      bro [i].f = a [i].f;
     }
    else
     {
      if(a [i].f >= bro [i].f)
       {
        for(int c = 0; c < coords; c++)
          bro [i].c [c] = a [i].c [c];
        bro [i].f = a [i].f;
       }
      //--- иначе откат: bro[i] остаётся прежним
     }
   }

//--- энергобаланс: приход ~ популяционно-относительной адаптированности
  double fmin =  DBL_MAX;
  double fmax = -DBL_MAX;
  for(int i = 0; i < popSize; i++)
   {
    if(bro [i].f < fmin)
      fmin = bro [i].f;
    if(bro [i].f > fmax)
      fmax = bro [i].f;
   }
  double span = fmax - fmin;
  if(span < 1e-12)
    span = 1e-12;

  for(int i = 0; i < popSize; i++)
   {
    double gain = (bro [i].f - fmin) / span;   // [0,1]
    bro [i].E += gain - cost;
    if(bro [i].E > Emax)
      bro [i].E = Emax;
   }

//--- гибель: энергетическая (E<=0) или средовая (pDeath)
  int nLiving = 0;
  bool dead [];
  ArrayResize(dead, popSize);

  for(int i = 0; i < popSize; i++)
   {
    bool die = (bro [i].E <= 0.0) || (u.RNDprobab() < pDeath);
    dead [i] = die;
    if(!die)
     {
      living [nLiving] = i;
      nLiving++;
     }
   }

//--- заполнение освободившихся слотов делением (или массовый респаун)
  if(nLiving == 0)
   {
    for(int i = 0; i < popSize; i++)
      Respawn(i);
    return;
   }

  double step = StepNow();

  for(int i = 0; i < popSize; i++)
   {
    if(!dead [i])
      continue;

    //--- родитель — энергичный живой (рулетка по энергии)
    int p = RouletteByE(living, nLiving, -1);

    //--- деление: потомок = мутированная копия родителя рядом
    for(int c = 0; c < coords; c++)
     {
      double rng = rangeMax [c] - rangeMin [c];
      double v   = bro [p].c [c] + Gauss() * step * rng;
      bro [i].c [c] = u.SeInDiSp(v, rangeMin [c], rangeMax [c], rangeStep [c]);
     }

    //--- энергия делится пополам
    double half = bro [p].E * 0.5;
    bro [i].E = half;
    bro [p].E = bro [p].E - half;

    //--- фитнес потомка оценится на следующей эпохе; оценка-заглушка
    bro [i].f = bro [p].f;
   }
 }
//+------------------------------------------------------------------+