//+——————————————————————————————————————————————————————————————————+
//|                                                         C_AO_CVO |
//|                                  Copyright 2007-2026, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

#include "#C_AO.mqh"

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
struct S_CVO_Patient
 {
  double             c [];
  double             f;
 };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class C_AO_CVO : public C_AO
 {
public:
                    ~C_AO_CVO() {}
                     C_AO_CVO()
   {
    ao_name = "CVO";
    ao_desc = "Corona Virus Optimization";
    ao_link = "https://www.mql5.com/ru/articles/22887";

    popSize = 50;     // бюджет FE на эпоху (= кандидатов за проход)
    nPop    = 200;    // потолок инфекционной популяции
    basePop = 5;      // стартовая инфекционная популяция
    R0      = 5;      // потомков-восприимчивых на одного носителя
    rho     = 0.1;    // контактность: общий множитель масштаба λ

    ArrayResize(params, 5);
    params [0].name = "popSize"; params [0].val = popSize;
    params [1].name = "nPop";    params [1].val = nPop;
    params [2].name = "basePop"; params [2].val = basePop;
    params [3].name = "R0";      params [3].val = R0;
    params [4].name = "rho";     params [4].val = rho;
   }

  void               SetParams()
   {
    popSize = (int)params [0].val;
    nPop    = (int)params [1].val;
    basePop = (int)params [2].val;
    R0      = (int)params [3].val;
    rho     =      params [4].val;

    //--- предохранители
    if(popSize < 1)       popSize = 1;
    if(nPop    < 1)       nPop    = 1;
    if(R0      < 1)       R0      = 1;
    if(rho   < 0.0)       rho     = 0.0;
    if(basePop < 1)       basePop = 1;
    if(basePop > popSize) basePop = popSize; // пул засевается из 1-го батча
    if(basePop > nPop)    basePop = nPop;
   }

  bool               Init(const double &rangeMinP  [],
                          const double &rangeMaxP  [],
                          const double &rangeStepP [],
                          const int     epochsP = 0);

  void               Moving();
  void               Revision();

  //--- видимые параметры
  int                nPop;       // потолок инфекционной популяции
  int                basePop;    // стартовая инфекционная популяция
  int                R0;         // потомков на носителя
  double             rho;        // множитель масштаба λ

private:
  //--- инфекционный пул (массив структур: координаты + фитнес вместе)
  S_CVO_Patient      inf [];        // [nPop]
  int                infCount;      // текущий размер пула (<= nPop)

  //--- привязка слота-кандидата к носителю
  int                parentOf [];   // [popSize]

  //--- временный пул пересборки
  S_CVO_Patient      tmp [];        // [nPop + popSize]

  //--- стандартное нормальное (Box-Muller), Eq.3
  double             Gauss();
  //--- индексная сортировка по убыванию f
  void               ArgSortDesc(const double &f [], int &idx [], const int count);
 };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                              Init                                |
//+------------------------------------------------------------------+
bool C_AO_CVO::Init(const double &rangeMinP  [],
                    const double &rangeMaxP  [],
                    const double &rangeStepP [],
                    const int     epochsP = 0)
 {
  if(!StandardInit(rangeMinP, rangeMaxP, rangeStepP))
    return false;

//--- стартовый случайный батч в cP[] (оценится на первом проходе)
  for(int i = 0; i < popSize; i++)
    for(int c = 0; c < coords; c++)
      a [i].cP [c] = u.RNDfromCI(rangeMin [c], rangeMax [c]);

//--- инфекционный пул и временный пул (массивы структур)
  ArrayResize(inf, nPop);
  for(int i = 0; i < nPop; i++)
    ArrayResize(inf [i].c, coords);

  ArrayResize(tmp, nPop + popSize);
  for(int i = 0; i < nPop + popSize; i++)
    ArrayResize(tmp [i].c, coords);

  ArrayResize(parentOf, popSize);
  infCount = 0;

  return true;
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                            Moving                                |
//|                                                                  |
//|  Первый проход: случайный батч cP -> c (стенд его оценит, далее  |
//|  Revision засеет инфекционный пул лучшими basePop особями).      |
//|  Боевой проход: каждый слот-кандидат = восприимчивый, рождённый  |
//|  носителем p = (k/R0) % infCount возмущением λ (покоординатно).  |
//+------------------------------------------------------------------+
void C_AO_CVO::Moving()
 {
//--- первый прогон: cP -> c
  if(!revision)
   {
    for(int i = 0; i < popSize; i++)
      for(int c = 0; c < coords; c++)
        a [i].c [c] = u.SeInDiSp(a [i].cP [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    return;
   }

//--- генерация кандидатов-восприимчивых
  for(int k = 0; k < popSize; k++)
   {
    int denom = infCount > 0 ? infCount : 1; // защита от %0
    int p     = (k / R0) % denom;            // носитель данного слота
    parentOf [k] = p;

    for(int c = 0; c < coords; c++)
     {
      double phi  = Gauss();                                   // Eq.3, [D1] свой на координату
      double frac = (double)(infCount + k) / (double)nPop;     // I/N, Eq.1*Eq.2, [D3]
      if(frac > 1.0)
        frac = 1.0;
      double lam  = rho * phi * frac;        // абсолютный шаг, ТОЧНО по статье (Eq.1*2*3)
      double v    = inf [p].c [c] + lam * (rangeMax [c] - rangeMin [c]);
      a [k].c [c] = u.SeInDiSp(v, rangeMin [c], rangeMax [c], rangeStep [c]);
     }
   }
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                           Revision                               |
//|                                                                  |
//|  1. обновление личных и глобального лучших по кандидатам;        |
//|  2. первый проход: посев пула лучшими basePop особями;          |
//|  3. LocalBest эпохи; решение замена/рост пула;                  |
//|  4. приёмка восприимчивых (иммунитет >= носителя);              |
//|  5. отсечение пула до nPop по убыванию f (выздоровление).       |
//+------------------------------------------------------------------+
void C_AO_CVO::Revision()
 {
//--- сохранить прежний глобальный лучший ДО обновления (для решения замена/рост)
  double fB_prev = fB;

//--- 1) обновление лучших
  for(int i = 0; i < popSize; i++)
   {
    if(!MathIsValidNumber(a [i].f))    // нечисло не должно отравлять лучших/пул
      continue;
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

//--- 2) первый проход: посев инфекционного пула лучшими basePop
  if(!revision)
   {
    int    idx [];
    double af  [];
    ArrayResize(idx, popSize);
    ArrayResize(af,  popSize);
    for(int i = 0; i < popSize; i++) af [i] = a [i].f;

    ArgSortDesc(af, idx, popSize);

    infCount = basePop;                // basePop <= popSize (зажато в SetParams)
    for(int n = 0; n < infCount; n++)
     {
      int s = idx [n];
      ArrayCopy(inf [n].c, a [s].c, 0, 0, coords);
      inf [n].f = a [s].f;
     }

    revision = true;
    return;
   }

//--- 3) лучший кандидат эпохи и решение замена/рост
  double localBest = -DBL_MAX;
  for(int i = 0; i < popSize; i++)
    if(MathIsValidNumber(a [i].f) && a [i].f > localBest)
      localBest = a [i].f;

  bool improved = (localBest > fB_prev);

//--- 4) сборка обновлённого пула во временный пул (предвыделен в Init)
  int cnt = 0;

  if(!improved)
   {
    //--- рост: сначала весь прежний пул
    for(int n = 0; n < infCount; n++)
     {
      ArrayCopy(tmp [cnt].c, inf [n].c, 0, 0, coords);
      tmp [cnt].f = inf [n].f;
      cnt++;
     }
   }
//--- принятые восприимчивые: иммунитет не лучше носителя -> заразился
//    (читаем inf[p].f ДО перезаписи пула — наложения нет)
  for(int k = 0; k < popSize; k++)
   {
    int p = parentOf [k];
    if(MathIsValidNumber(a [k].f) && a [k].f >= inf [p].f)
     {
      ArrayCopy(tmp [cnt].c, a [k].c, 0, 0, coords);
      tmp [cnt].f = a [k].f;
      cnt++;
     }
   }

//--- защита: при замене с пустым набором сохраняем прежний пул
  if(cnt == 0)
    return;

//--- 5) отсечение до nPop (выздоровление/смерть слабейших)
  if(cnt > nPop)
   {
    int    idx [];
    double tf  [];
    ArrayResize(idx, cnt);
    ArrayResize(tf,  cnt);
    for(int n = 0; n < cnt; n++)
      tf [n] = tmp [n].f;

    ArgSortDesc(tf, idx, cnt);
    for(int n = 0; n < nPop; n++)
     {
      int s = idx [n];
      ArrayCopy(inf [n].c, tmp [s].c, 0, 0, coords);
      inf [n].f = tmp [s].f;
     }
    infCount = nPop;
   }
  else
   {
    for(int n = 0; n < cnt; n++)
     {
      ArrayCopy(inf [n].c, tmp [n].c, 0, 0, coords);
      inf [n].f = tmp [n].f;
     }
    infCount = cnt;
   }
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|              Gauss — стандартное нормальное (Box-Muller)         |
//+------------------------------------------------------------------+
double C_AO_CVO::Gauss()
 {
  double u1 = u.RNDfromCI(0.0, 1.0);
  double u2 = u.RNDfromCI(0.0, 1.0);
  if(u1 < 1e-15)
    u1 = 1e-15;                       // защита от ln(0)
  return MathSqrt(-2.0 * MathLog(u1)) * MathCos(2.0 * M_PI * u2);
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|     ArgSortDesc — порядок индексов по убыванию f                 |
//+------------------------------------------------------------------+
void C_AO_CVO::ArgSortDesc(const double &f [], int &idx [], const int count)
 {
  for(int i = 0; i < count; i++)
    idx [i] = i;

  for(int i = 0; i < count - 1; i++)
   {
    int b = i;
    for(int j = i + 1; j < count; j++)
      if(f [idx [j]] > f [idx [b]])
        b = j;
    if(b != i)
     {
      int t = idx [i];
      idx [i] = idx [b];
      idx [b] = t;
     }
   }
 }
//+------------------------------------------------------------------+
