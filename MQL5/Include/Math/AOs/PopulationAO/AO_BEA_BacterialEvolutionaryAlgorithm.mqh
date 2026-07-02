//+——————————————————————————————————————————————————————————————————+
//|                                                         C_AO_BEA |
//|                                  Copyright 2007-2026, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

#include "#C_AO.mqh"

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
struct S_BEA_Ind
 {
  double c           [];   // координаты (валидные, уже прогнаны через SeInDiSp)
  double             f;      // фитнес

  void               Init(int coords)
   {
    ArrayResize(c, coords);
    f = -DBL_MAX;
   }
 };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class C_AO_BEA : public C_AO
 {
public:
                    ~C_AO_BEA() {}
                     C_AO_BEA()
   {
    ao_name = "BEA";
    ao_desc = "Bacterial Evolutionary Algorithm";
    ao_link = "https://www.mql5.com/ru/articles/23232";

    popSize   = 10;     // число бактерий
    Nc        = 2;      // свежих мутантов сегмента на бактерию (>=1)
    Nseg      = 10;     // число сегментов хромосомы
    nTrans    = 2;      // эпох переноса генов на поколение
    transFrac = 0.2;    // длина переносимого блока как доля coords (0..1]

    ArrayResize(params, 5);
    params [0].name = "popSize";   params [0].val = popSize;
    params [1].name = "Nc";        params [1].val = Nc;
    params [2].name = "Nseg";      params [2].val = Nseg;
    params [3].name = "nTrans";    params [3].val = nTrans;
    params [4].name = "transFrac"; params [4].val = transFrac;
   }

  void               SetParams()
   {
    popSize   = (int)params [0].val;
    Nc        = (int)params [1].val;
    Nseg      = (int)params [2].val;
    nTrans    = (int)params [3].val;
    transFrac =      params [4].val;

    //--- предохранители
    if(popSize   < 1)
      popSize   = 1;
    if(Nc        < 1)
      Nc        = 1;     // число свежих мутантов на сегмент
    if(Nseg      < 1)
      Nseg      = 1;
    if(nTrans    < 0)
      nTrans    = 0;
    if(transFrac <= 0.0)
      transFrac = 0.01;
    if(transFrac >  1.0)
      transFrac = 1.0;
   }

  bool               Init(const double &rangeMinP  [],
                          const double &rangeMaxP  [],
                          const double &rangeStepP [],
                          const int     epochsP = 0);

  void               Moving();
  void               Revision();

  //--- видимые параметры
  int                Nc;          // клонов
  int                Nseg;        // сегментов
  int                nTrans;      // эпох переноса генов
  double             transFrac;   // доля coords для блока переноса

private:
  //--- фазы конечного автомата
  enum E_Phase { PH_MUT = 0, PH_TRANS = 1 };
  E_Phase            phase;
  int                segIdx;      // текущий сегмент в мутации
  int                cloneIdx;    // текущий клон (срез оценки) в мутации
  int                transIdx;    // текущая эпоха переноса

  //--- данные (массивы структур)
  S_BEA_Ind          bact  [];    // [popSize]       — закреплённые бактерии
  S_BEA_Ind          clone [];    // [popSize * Nc]  — рабочие клоны
  int                ord   [];    // [popSize]       — индексы для сортировки

  int                transLen;    // длина блока переноса (из transFrac)

  //--- вспомогательные
  void               SegBounds(int s, int &lo, int &hi);
  void               BuildClonesForSeg(int s);
  void               DoGeneTransfer();
 };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                              Init                                |
//+------------------------------------------------------------------+
bool C_AO_BEA::Init(const double &rangeMinP  [],
                    const double &rangeMaxP  [],
                    const double &rangeStepP [],
                    const int     epochsP = 0)
 {
  if(!StandardInit(rangeMinP, rangeMaxP, rangeStepP))
    return false;

//--- сегментов не больше, чем координат
  if(Nseg > coords)
    Nseg = coords;

//--- длина блока переноса
  transLen = (int)MathRound(transFrac * coords);
  if(transLen < 1)
    transLen = 1;
  if(transLen > coords)
    transLen = coords;

//--- буферы (массивы структур)
  ArrayResize(bact,  popSize);
  ArrayResize(clone, popSize * Nc);
  ArrayResize(ord,   popSize);

  for(int i = 0; i < popSize;      i++)
    bact  [i].Init(coords);
  for(int i = 0; i < popSize * Nc; i++)
    clone [i].Init(coords);

//--- стартовая случайная популяция бактерий
  for(int i = 0; i < popSize; i++)
   {
    for(int c = 0; c < coords; c++)
      bact [i].c [c] = u.SeInDiSp(u.RNDfromCI(rangeMin [c], rangeMax [c]),
                                  rangeMin [c], rangeMax [c], rangeStep [c]);
    bact [i].f = -DBL_MAX;
   }

//--- стартовое состояние автомата
  phase    = PH_MUT;
  segIdx   = 0;
  cloneIdx = 0;
  transIdx = 0;

  return true;
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   SegBounds — границы сегмента s в [lo, hi) с равномерным        |
//|   распределением остатка координат по сегментам                  |
//+------------------------------------------------------------------+
void C_AO_BEA::SegBounds(int s, int &lo, int &hi)
 {
  lo = (int)((long)s       * coords / Nseg);
  hi = (int)((long)(s + 1) * coords / Nseg);
  if(lo < 0)
    lo = 0;
  if(hi > coords)
    hi = coords;
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   BuildClonesForSeg — снять с каждой бактерии Nc клонов и во     |
//|   ВСЕХ клонах случайно регенерировать координаты сегмента s.     |
//|   Инкумбент не хранится отдельным клоном: его фитнес помнится в  |
//|   bact.f и сравнивается при закреплении (см. Revision).          |
//+------------------------------------------------------------------+
void C_AO_BEA::BuildClonesForSeg(int s)
 {
  int lo, hi;
  SegBounds(s, lo, hi);

  for(int b = 0; b < popSize; b++)
   {
    for(int j = 0; j < Nc; j++)
     {
      int idx = b * Nc + j;

      //--- копия бактерии вне сегмента
      for(int c = 0; c < coords; c++)
        clone [idx].c [c] = bact [b].c [c];

      //--- регенерация сегмента во ВСЕХ клонах. Инкумбент хранится
      //    отдельно (bact.f) и сравнивается на этапе закрепления —
      //    переоценивать его клоном не нужно, экономим вызовы FF.
      for(int c = lo; c < hi; c++)
        clone [idx].c [c] = u.SeInDiSp(u.RNDfromCI(rangeMin [c], rangeMax [c]),
                                       rangeMin [c], rangeMax [c], rangeStep [c]);

      clone [idx].f = -DBL_MAX;
     }
   }
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|   DoGeneTransfer — сортировка популяции и инфицирование «слабой» |
//|   половины блоками из «сильной». Модифицирует bact[] на месте.   |
//+------------------------------------------------------------------+
void C_AO_BEA::DoGeneTransfer()
 {
  if(popSize < 2)
    return;

//--- сортировка индексов по фитнесу бактерий (по убыванию)
  for(int i = 0; i < popSize; i++)
    ord [i] = i;

  for(int i = 0; i < popSize - 1; i++)
   {
    int bi = i;
    for(int j = i + 1; j < popSize; j++)
      if(bact [ord [j]].f > bact [ord [bi]].f)
        bi = j;
    int t = ord [i];
    ord [i] = ord [bi];
    ord [bi] = t;
   }

  int half = popSize / 2;
  if(half < 1)
    half = 1;

  int span = coords - transLen + 1;
  if(span < 1)
    span = 1;

//--- каждый приёмник «слабой» половины принимает блок от случайного
//    источника «сильной» половины
  for(int k = half; k < popSize; k++)
   {
    int dst   = ord [k];
    int src   = ord [u.RNDminusOne(half)];
    int start = u.RNDminusOne(span);

    for(int t = 0; t < transLen; t++)
     {
      int c = start + t;
      if(c < coords)
        bact [dst].c [c] = bact [src].c [c];
     }
   }
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                            Moving                                |
//|   Заполняет слоты a[i].c кандидатами для оценки этой эпохи.      |
//+------------------------------------------------------------------+
void C_AO_BEA::Moving()
 {
//--- первый прогон: бактерии -> слоты для первичной оценки FF
  if(!revision)
   {
    for(int i = 0; i < popSize; i++)
      for(int c = 0; c < coords; c++)
        a [i].c [c] = u.SeInDiSp(bact [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    return;
   }

  if(phase == PH_MUT)
   {
    //--- на старте сегмента строим клоны
    if(cloneIdx == 0)
      BuildClonesForSeg(segIdx);

    //--- грузим текущий клон каждой бактерии в слоты
    for(int b = 0; b < popSize; b++)
     {
      int idx = b * Nc + cloneIdx;
      for(int c = 0; c < coords; c++)
        a [b].c [c] = clone [idx].c [c];
     }
   }
  else // PH_TRANS
   {
    //--- инфицируем «слабую» половину и грузим бактерии в слоты
    DoGeneTransfer();

    for(int b = 0; b < popSize; b++)
      for(int c = 0; c < coords; c++)
        a [b].c [c] = bact [b].c [c];
   }
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                           Revision                               |
//|   Считывает a[i].f, обновляет cB/fB и продвигает автомат.        |
//+------------------------------------------------------------------+
void C_AO_BEA::Revision()
 {
//--- глобальный лучший (элитизм: лучшее не теряется при переносе)
  for(int i = 0; i < popSize; i++)
   {
    if(a [i].f > fB)
     {
      fB = a [i].f;
      ArrayCopy(cB, a [i].c, 0, 0, coords);
     }
   }

//--- первичная фиксация фитнеса бактерий
  if(!revision)
   {
    for(int i = 0; i < popSize; i++)
      bact [i].f = a [i].f;

    revision = true;
    phase    = PH_MUT;
    segIdx   = 0;
    cloneIdx = 0;
    transIdx = 0;
    return;
   }

  if(phase == PH_MUT)
   {
    //--- сохраняем фитнес оценённого среза клонов
    for(int b = 0; b < popSize; b++)
      clone [b * Nc + cloneIdx].f = a [b].f;

    cloneIdx++;

    //--- сегмент полностью оценён -> закрепляем лучший клон
    if(cloneIdx >= Nc)
     {
      for(int b = 0; b < popSize; b++)
       {
        int    best = 0;
        double bf   = clone [b * Nc + 0].f;
        for(int j = 1; j < Nc; j++)
         {
          double cf = clone [b * Nc + j].f;
          if(cf > bf)
           {
            bf = cf;
            best = j;
           }
         }

        //--- закрепляем лучший мутант ТОЛЬКО если он превзошёл
        //    инкумбента; иначе бактерия остаётся прежней (монотонность
        //    мутации сохраняется без переоценки инкумбента клоном)
        if(bf > bact [b].f)
         {
          int bi = b * Nc + best;
          for(int c = 0; c < coords; c++)
            bact [b].c [c] = clone [bi].c [c];
          bact [b].f = bf;
         }
       }

      cloneIdx = 0;
      segIdx++;

      //--- все сегменты пройдены -> перенос генов (или новое поколение)
      if(segIdx >= Nseg)
       {
        segIdx = 0;
        if(nTrans > 0)
         {
          phase = PH_TRANS;
          transIdx = 0;
         }
        else
          phase = PH_MUT;
       }
     }
   }
  else // PH_TRANS
   {
    //--- принимаем инфицированные генотипы (не-greedy)
    for(int b = 0; b < popSize; b++)
      bact [b].f = a [b].f;

    transIdx++;
    if(transIdx >= nTrans)
     {
      phase    = PH_MUT;
      segIdx   = 0;
      cloneIdx = 0;
     }
   }
 }
//+------------------------------------------------------------------+