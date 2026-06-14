//+——————————————————————————————————————————————————————————————————+
//|                                                   C_AO_COA_Coyote|
//|                                  Copyright 2007-2025, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//+——————————————————————————————————————————————————————————————————+
#include "#C_AO.mqh"
//+------------------------------------------------------------------+
//| Coyote Optimization Algorithm (COA).                             |
//| Pierezan J., Coelho L.S. "Coyote Optimization Algorithm: A new   |
//| metaheuristic for global optimization problems", IEEE CEC 2018.  |
//|                                                                  |
//| Структура: популяция разбита на nPacks стай по nCoy койотов.     |
//| За эпоху в каждой стае:                                          |
//|   1. определяется альфа (лучший по фитнесу, Eq. 5);              |
//|   2. считается «социальная тенденция» — покоординатная медиана   |
//|      стаи (Eq. 6);                                               |
//|   3. каждый койот делает социальный ход (Eq. 12):                |
//|        new = x + r1*(alpha - x_rc1) + r2*(tend - x_rc2),         |
//|      r1, r2 — СКАЛЯРЫ на агента (как в оригинале);               |
//|   4. рождается щенок (Eq. 7, Alg. 1): бинарный кроссовер двух    |
//|      случайных родителей + равномерный шум по части координат;   |
//|   5. принятие — только при строгом улучшении (Eq. 14).           |
//| Редкий обмен койотами между стаями (Eq. 4, p_leave) и счётчик    |
//| возраста (щенок обнуляет возраст слота).                         |
//|                                                                  |
//| АДАПТАЦИЯ К БАТЧЕВОЙ СХЕМЕ СТЕНДА (Moving/Revision):             |
//| 1. Оригинал последовательный: каждый новый койот оценивается     |
//|    сразу, плюс отдельная оценка щенка — nCoy+1 вызовов FF на     |
//|    стаю за «год». Здесь щенок занимает слот ХУДШЕГО койота стаи  |
//|    (при равенстве фитнеса — старшего по возрасту); социальный    |
//|    ход этого койота в текущую эпоху пропускается. Итог: ровно    |
//|    popSize вызовов FF за эпоху, бюджет расходуется честно.       |
//| 2. Приживание щенка в оригинале: замена старейшего среди тех,    |
//|    кто хуже щенка. Здесь слот выбран заранее (худший/старший),   |
//|    щенок принимается, если строго лучше него — поведенчески      |
//|    эквивалентное приближение в рамках greedy-acceptance.         |
//| 3. Все ходы эпохи считаются от снапшота популяции (в оригинале   |
//|    обновления внутри стаи частично последовательные).            |
//|                                                                  |
//| ВНИМАНИЕ (для off-diagonal проверки): в социальном ходе r1 и r2  |
//| — скаляры, домножающие целые векторы разностей. Это потенциаль-  |
//| ный кандидат на «диагональный чит» на бенче из копий функции —   |
//| проверить составным тестом с |x_norm - y_norm| > 0.3.            |
//+------------------------------------------------------------------+
class C_AO_COA_Coyote : public C_AO
 {
public:
                    ~C_AO_COA_Coyote() {}
                     C_AO_COA_Coyote()
   {
    ao_name = "COA(Coyote)";
    ao_desc = "Coyote Optimization Algorithm";
    ao_link = "https://www.mql5.com/ru/articles/23053";

    popSize = 20;
    nCoy    = 5;

    ArrayResize(params, 2);
    params [0].name = "popSize";
    params [0].val = popSize;
    params [1].name = "nCoy";
    params [1].val = nCoy;

    nPacks = popSize / nCoy;
   }

  void               SetParams()
   {
    popSize = (int)params [0].val;
    nCoy    = (int)params [1].val;

    nPacks = popSize / nCoy;
   }

  bool               Init(const double &rangeMinP  [],
                          const double &rangeMaxP  [],
                          const double &rangeStepP [],
                          const int     epochsP = 0);

  void               Moving();
  void               Revision();

  //--- видимые параметры
  int                nCoy;       // койотов в стае (оригинал: 5)

private:
  int                nPacks;     // число стай = popSize / nCoy
  double             pLeave;     // вероятность обмена между стаями: 0.005*nCoy^2
  double             probPup;    // (1 - Ps)/2, Ps = 1/coords

  int                packs  [];  // [nPacks * nCoy] — индексы агентов по стаям
  int                ages   [];  // [popSize] — возраст койота
  bool               isPup  [];  // [popSize] — слот щенка в текущую эпоху

  //--- snapshot для greedy acceptance
  double             snap_c [];  // [popSize * coords]
  double             snap_f [];  // [popSize]

  //--- рабочие буферы
  double             medBuf [];  // [nCoy]   — для медианы
  double             tend   [];  // [coords] — социальная тенденция стаи
  int                pdr    [];  // [coords] — перестановка измерений (щенок)
  int                mask   [];  // [coords] — 0: шум, 1: родитель 1, 2: родитель 2
 };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                              Init                                |
//+------------------------------------------------------------------+
bool C_AO_COA_Coyote::Init(const double &rangeMinP  [],
                    const double &rangeMaxP  [],
                    const double &rangeStepP [],
                    const int     epochsP = 0)
 {
  if(!StandardInit(rangeMinP, rangeMaxP, rangeStepP))
    return false;

//--- начальная случайная популяция в cP[] (Eq. 2)
  for(int i = 0; i < popSize; i++)
    for(int c = 0; c < coords; c++)
      a [i].cP [c] = u.RNDfromCI(rangeMin [c], rangeMax [c]);

//--- параметры, зависящие от nCoy и coords
  pLeave  = 0.005 * (double)nCoy * (double)nCoy;          // Eq. 4
  probPup = (coords > 0) ? (1.0 - 1.0 / (double)coords) * 0.5 : 0.0; // (1-Ps)/2

//--- буферы
  ArrayResize(snap_c, popSize * coords);
  ArrayResize(snap_f, popSize);
  ArrayResize(packs,  nPacks * nCoy);
  ArrayResize(ages,   popSize);
  ArrayResize(isPup,  popSize);
  ArrayResize(medBuf, nCoy);
  ArrayResize(tend,   coords);
  ArrayResize(pdr,    coords);
  ArrayResize(mask,   coords);

  ArrayInitialize(ages, 0);
  ArrayInitialize(isPup, false);

//--- случайное распределение койотов по стаям (Fisher–Yates)
  for(int i = 0; i < popSize; i++)
    packs [i] = i;
  for(int j = popSize - 1; j > 0; j--)
   {
    int r     = u.RNDminusOne(j + 1);
    int tmp   = packs [j];
    packs [j] = packs [r];
    packs [r] = tmp;
   }

  return true;
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                            Moving                                |
//|                                                                  |
//|  Структура эпохи:                                                |
//|    1. snapshot популяции;                                        |
//|    2. редкий обмен койотами между стаями (Eq. 4);                |
//|    3. для каждой стаи: альфа, тенденция-медиана, слот щенка;     |
//|    4. социальные ходы (Eq. 12) для всех, кроме слота щенка;      |
//|    5. в слот щенка — кроссовер двух родителей + шум (Alg. 1).    |
//+------------------------------------------------------------------+
void C_AO_COA_Coyote::Moving()
 {
//--- первый прогон: cP -> c, чтобы внешний цикл оценил FF
  if(!revision)
   {
    for(int i = 0; i < popSize; i++)
      for(int c = 0; c < coords; c++)
        a [i].c [c] = u.SeInDiSp(a [i].cP [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    return;
   }

//--- snapshot текущей популяции
  for(int i = 0; i < popSize; i++)
   {
    for(int c = 0; c < coords; c++)
      snap_c [i * coords + c] = a [i].c [c];
    snap_f [i] = a [i].f;
   }

//--- обмен койотами между стаями (Eq. 4)
  if(nPacks > 1 && u.RNDprobab() < pLeave)
   {
    int p1 = u.RNDminusOne(nPacks);
    int p2 = p1;
    while(p2 == p1)
      p2 = u.RNDminusOne(nPacks);

    int c1 = u.RNDminusOne(nCoy);
    int c2 = u.RNDminusOne(nCoy);

    int tmp = packs [p1 * nCoy + c1];
    packs [p1 * nCoy + c1] = packs [p2 * nCoy + c2];
    packs [p2 * nCoy + c2] = tmp;
   }

  ArrayInitialize(isPup, false);

//--- операции внутри каждой стаи
  for(int p = 0; p < nPacks; p++)
   {
    int base = p * nCoy;

    //--- альфа стаи (Eq. 5): максимум фитнеса (стенд максимизирует)
    int aIdx = packs [base];
    for(int k = 1; k < nCoy; k++)
     {
      int idx = packs [base + k];
      if(snap_f [idx] > snap_f [aIdx])
        aIdx = idx;
     }

    //--- социальная тенденция (Eq. 6): покоординатная медиана стаи
    for(int c = 0; c < coords; c++)
     {
      for(int k = 0; k < nCoy; k++)
        medBuf [k] = snap_c [packs [base + k] * coords + c];
      ArraySort(medBuf);
      tend [c] = (nCoy % 2 == 1) ? medBuf [nCoy / 2]
                 : 0.5 * (medBuf [nCoy / 2 - 1] + medBuf [nCoy / 2]);
     }

    //--- слот щенка: худший фитнес, при равенстве — старший возраст
    int wK = 0;
    for(int k = 1; k < nCoy; k++)
     {
      int idx  = packs [base + k];
      int wIdx = packs [base + wK];
      if(snap_f [idx] <  snap_f [wIdx] ||
         (snap_f [idx] == snap_f [wIdx] && ages [idx] > ages [wIdx]))
        wK = k;
     }
    int pupIdx = packs [base + wK];
    isPup [pupIdx] = true;

    //--- движение койотов стаи
    for(int k = 0; k < nCoy; k++)
     {
      int i = packs [base + k];

      if(i == pupIdx)
       {
        //--- РОЖДЕНИЕ ЩЕНКА (Eq. 7, Alg. 1) -------------------------
        //    два случайных различных родителя из стаи
        int k1 = u.RNDminusOne(nCoy);
        int k2 = k1;
        while(k2 == k1)
          k2 = u.RNDminusOne(nCoy);
        int par1 = packs [base + k1];
        int par2 = packs [base + k2];

        //--- перестановка измерений (Fisher–Yates)
        for(int c = 0; c < coords; c++)
          pdr [c] = c;

        for(int j = coords - 1; j > 0; j--)
         {
          int r   = u.RNDminusOne(j + 1);
          int tmp = pdr [j];
          pdr [j] = pdr [r];
          pdr [r] = tmp;
         }

        //--- маска: гарантированно по одной координате от каждого
        //    родителя, остальные — родитель 1 / родитель 2 / шум
        ArrayInitialize(mask, 0);
        mask [pdr [0]] = 1;
        if(coords > 1)
          mask [pdr [1]] = 2;
        for(int j = 2; j < coords; j++)
         {
          double r = u.RNDprobab();
          if(r < probPup)
            mask [pdr [j]] = 1;
          else
            if(r > 1.0 - probPup)
              mask [pdr [j]] = 2;
          // иначе 0 — равномерный шум по координате
         }

        for(int c = 0; c < coords; c++)
         {
          double v;
          if(mask [c] == 1)
            v = snap_c [par1 * coords + c];
          else
            if(mask [c] == 2)
              v = snap_c [par2 * coords + c];
            else
              v = u.RNDfromCI(rangeMin [c], rangeMax [c]);

          a [i].c [c] = u.SeInDiSp(v, rangeMin [c], rangeMax [c], rangeStep [c]);
         }
       }
      else
       {
        //--- СОЦИАЛЬНЫЙ ХОД (Eq. 12) --------------------------------
        //    rc1, rc2 — различные, не совпадающие с k
        int k1 = k;
        while(k1 == k)
          k1 = u.RNDminusOne(nCoy);
        int k2 = k;
        while(k2 == k || k2 == k1)
          k2 = u.RNDminusOne(nCoy);
        int rc1 = packs [base + k1];
        int rc2 = packs [base + k2];

        //--- r1, r2 — скаляры на агента, как в оригинале
        double r1 = u.RNDprobab();
        double r2 = u.RNDprobab();

        for(int c = 0; c < coords; c++)
         {
          double v = snap_c [i * coords + c] +
                     r1 * (snap_c [aIdx * coords + c] - snap_c [rc1 * coords + c]) +
                     r2 * (tend [c]                   - snap_c [rc2 * coords + c]);

          a [i].c [c] = u.SeInDiSp(v, rangeMin [c], rangeMax [c], rangeStep [c]);
         }
       }
     } // k
   } // p
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                           Revision                               |
//|                                                                  |
//|  Стандартное обновление fB/cB + адаптация (Eq. 14): кандидат     |
//|  принимается только при СТРОГОМ улучшении, иначе откат к         |
//|  снапшоту (касается и социальных ходов, и щенка). Возраст всех   |
//|  койотов растёт на 1 за эпоху; прижившийся щенок обнуляет        |
//|  возраст слота.                                                  |
//+------------------------------------------------------------------+
void C_AO_COA_Coyote::Revision()
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

//--- адаптация (Eq. 14): принятие только при строгом улучшении
  for(int i = 0; i < popSize; i++)
   {
    bool accepted = (a [i].f > snap_f [i]);

    if(!accepted)
     {
      for(int c = 0; c < coords; c++)
        a [i].c [c] = snap_c [i * coords + c];
      a [i].f = snap_f [i];
     }

    //--- возраст
    ages [i]++;
    if(isPup [i] && accepted)
      ages [i] = 0;
   }
 }
//+------------------------------------------------------------------+