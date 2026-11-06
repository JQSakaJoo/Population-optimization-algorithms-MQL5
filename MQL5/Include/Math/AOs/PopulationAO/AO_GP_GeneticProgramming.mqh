//+------------------------------------------------------------------+
//|                                  Copyright 2007-2026, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//+------------------------------------------------------------------+

#include "#C_AO.mqh"

// GP — Genetic Programming (Koza 1992; Koza 1994, Statistics and
// Computing 4(2):87-112).
//
// Перенос на стенд C_AO:
//   особь = программа из coords результирующих ветвей, по одной на
//   координату; входов у ветвей нет, терминал один — ERC из [-1; 1];
//   функции {+, -, *, %}; выход ветви через обёртку (wrapper)
//   отображается в диапазон координаты: clamp(v,-1,1) -> [min; max].
//
// Канон (статья 1994):
//   кроссовер 90% / репродукция 10%, мутации нет; точка кроссовера —
//   функция 90%, лист 10%; глубина начальных деревьев <= 6, после
//   кроссовера <= 17; отбор обоих родителей пропорционально
//   скорректированной приспособленности; без элитизма; результат —
//   лучший за прогон; ERC фиксируются при рождении; % при делении
//   на ноль возвращает 1; результат операции ограничен +-1e10.
//
// Из книги GP I (1992), в открытой статье не расписано:
//   ramped half-and-half с глубинами 2..initDepth, корень — функция;
//   a = 1/(1+s); потомок глубже maxDepth заменяется копией родителя.
//
// Не из статьи (решения реализации):
//   s_i = fMax_gen - f_i (стенд максимизирует, r_max неизвестен) —
//         отбор зависит от абсолютного масштаба f;
//   кроссовер с типизацией ветвей (branch typing, как для ADF в GP II):
//         точка в P1 выбирается по всей программе, точка в P2 — в той
//         же ветви, значит один кроссовер меняет одну координату;
//   правило Коза про |v| < 1e-10 не реализовано;
//   удаление дубликатов поколения 0 не реализовано: при coords >= 2
//         и вещественных ERC совпадение программ практически исключено.
//
// Одна эпоха стенда = одно поколение = popSize вычислений ФФ.
//+------------------------------------------------------------------+

#define GP_OP_ADD   0
#define GP_OP_SUB   1
#define GP_OP_MUL   2
#define GP_OP_DIV   3
#define GP_OP_ERC   4
#define GP_NFUNC    4        // все функции двухместные
#define GP_MIN_DEP  2        // нижняя глубина ramped half-and-half
#define GP_BIG      1.0e10   // ограничение результата операции

//+------------------------------------------------------------------+
//| узел дерева в префиксной записи                                  |
//+------------------------------------------------------------------+
struct S_GP_Node
  {
   int               op;     // GP_OP_*
   double            val;    // константа для GP_OP_ERC
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| ветвь = одно дерево = одна координата                            |
//| все функции двухместные: внутренних узлов (len-1)/2,             |
//| листьев (len+1)/2                                                |
//+------------------------------------------------------------------+
struct S_GP_Branch
  {
   S_GP_Node         n [];
   int               len;
   double            v;      // кэш значения ветви
   bool              ok;     // кэш действителен
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
struct S_GP_Prog
  {
   S_GP_Branch       br [];
   double            f;
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
struct S_GP_Pop
  {
   S_GP_Prog         p [];
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class C_AO_GP : public C_AO
  {
public:
                    ~C_AO_GP() {}
                     C_AO_GP()
     {
      ao_name = "GP";
      ao_desc = "Genetic Programming";
      ao_link = "https://www.mql5.com/ru/articles/25549";

      popSize   = 50;      // M
      pCross    = 0.9;     // доля кроссовера, остальное — репродукция
      pInternal = 0.9;     // вероятность точки кроссовера во внутреннем узле
      initDepth = 6;       // максимальная глубина начальных деревьев
      maxDepth  = 17;      // максимальная глубина после кроссовера
      fScale    = 10000.0; // множитель s в отборе; 1.0 — канон

      ArrayResize(params, 6);
      params [0].name = "popSize";
      params [0].val = popSize;
      params [1].name = "pCross";
      params [1].val = pCross;
      params [2].name = "pInternal";
      params [2].val = pInternal;
      params [3].name = "initDepth";
      params [3].val = initDepth;
      params [4].name = "maxDepth";
      params [4].val = maxDepth;
      params [5].name = "fScale";
      params [5].val = fScale;
     }

   void               SetParams()
     {
      popSize   = (int)params [0].val;
      pCross    = params [1].val;
      pInternal = params [2].val;
      initDepth = (int)params [3].val;
      maxDepth  = (int)params [4].val;
      fScale    = params [5].val;

      if(popSize < 2)
         popSize = 2;
      if(pCross < 0.0)
         pCross = 0.0;
      if(pCross > 1.0)
         pCross = 1.0;
      if(pInternal < 0.0)
         pInternal = 0.0;
      if(pInternal > 1.0)
         pInternal = 1.0;
      if(initDepth < GP_MIN_DEP)
         initDepth = GP_MIN_DEP;
      if(maxDepth < initDepth)
         maxDepth = initDepth;
      if(fScale < 0.0)
         fScale = 0.0;

      params [0].val = popSize;
      params [1].val = pCross;
      params [2].val = pInternal;
      params [3].val = initDepth;
      params [4].val = maxDepth;
      params [5].val = fScale;
     }

   bool               Init(const double &rangeMinP  [],
                           const double &rangeMaxP  [],
                           const double &rangeStepP [],
                           const int     epochsP = 0);

   void               Moving();
   void               Revision();

   //----------------------------------------------------------------
   double             pCross;
   double             pInternal;
   int                initDepth;
   int                maxDepth;
   double             fScale;     // a = 1/(1 + fScale*s); 1.0 — канон

   double             avgSize;    // средний размер ветви в узлах (рост деревьев)

private:
   S_GP_Pop           gen [2];    // текущее и следующее поколения
   int                cur;
   double             cum [];     // накопленная скорректированная приспособленность
   double             cumSum;

   void               Ramped();
   void               Gen(S_GP_Branch &b, int d, bool full, bool root);
   void               Push(S_GP_Branch &b, int op, double val);

   double             EvalRec(const S_GP_Branch &b, int &pos);
   double             BranchValue(S_GP_Branch &b);
   double             Wrap(double v, int k);

   int                SubtreeEnd(const S_GP_Branch &b, int i);
   int                DepthRec(const S_GP_Branch &b, int &pos);
   int                Depth(const S_GP_Branch &b);
   int                NthNode(const S_GP_Branch &b, int r, bool internal);
   void               PickInProgram(const S_GP_Prog &p, int &k, int &i);
   int                PickInBranch(const S_GP_Branch &b);

   void               CopyBranch(S_GP_Branch &dst, const S_GP_Branch &src);
   void               CopyProg(S_GP_Prog &dst, const S_GP_Prog &src);
   void               Splice(S_GP_Branch &dst,
                             const S_GP_Branch &a, int i, int ie,
                             const S_GP_Branch &b, int j, int je);

   void               PrepRoulette();
   int                Spin();
   void               Crossover(int p1, int p2, int c1, int c2);
   void               Breed();

   int                RndInt(int n);
  };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                              Init                                |
//+------------------------------------------------------------------+
bool C_AO_GP::Init(const double &rangeMinP  [],
                   const double &rangeMaxP  [],
                   const double &rangeStepP [],
                   const int     epochsP = 0)
  {
   if(!StandardInit(rangeMinP, rangeMaxP, rangeStepP))
      return false;

   cur     = 0;
   avgSize = 0.0;

   for(int g = 0; g < 2; g++)
     {
      ArrayResize(gen [g].p, popSize);
      for(int i = 0; i < popSize; i++)
        {
         ArrayResize(gen [g].p [i].br, coords);
         gen [g].p [i].f = -DBL_MAX;
         for(int k = 0; k < coords; k++)
           {
            gen [g].p [i].br [k].len = 0;
            gen [g].p [i].br [k].v   = 0.0;
            gen [g].p [i].br [k].ok  = false;
           }
        }
     }

   ArrayResize(cum, popSize);
   cumSum = 0.0;
   return true;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| случайное целое [0; n-1], 30 бит (MathRand даёт только 15)       |
//+------------------------------------------------------------------+
int C_AO_GP::RndInt(int n)
  {
   if(n <= 1)
      return 0;
   long r = ((long)MathRand() << 15) | (long)MathRand();
   return (int)(r % n);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| добавить узел в конец ветви                                      |
//+------------------------------------------------------------------+
void C_AO_GP::Push(S_GP_Branch &b, int op, double val)
  {
   if(ArraySize(b.n) <= b.len)
      ArrayResize(b.n, b.len + 1, 64);
   b.n [b.len].op  = op;
   b.n [b.len].val = val;
   b.len++;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| full / grow, корень всегда функция [GP I]                        |
//| grow: терминал с вероятностью |T|/(|T|+|F|) = 1/5                |
//+------------------------------------------------------------------+
void C_AO_GP::Gen(S_GP_Branch &b, int d, bool full, bool root)
  {
   if(d == 0 || (!full && !root && u.RNDprobab() < 1.0 / (GP_NFUNC + 1)))
     {
      Push(b, GP_OP_ERC, u.RNDfromCI(-1.0, 1.0));   // ERC фиксируется навсегда
      return;
     }

   Push(b, RndInt(GP_NFUNC), 0.0);
   Gen(b, d - 1, full, false);
   Gen(b, d - 1, full, false);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| ramped half-and-half: глубины GP_MIN_DEP..initDepth по кругу,    |
//| половина full, половина grow; все ветви особи — с её настройками |
//+------------------------------------------------------------------+
void C_AO_GP::Ramped()
  {
   int levels = initDepth - GP_MIN_DEP + 1;

   for(int i = 0; i < popSize; i++)
     {
      int  d    = GP_MIN_DEP + i % levels;
      bool full = ((i / levels) % 2) == 0;

      for(int k = 0; k < coords; k++)
        {
         gen [cur].p [i].br [k].len = 0;
         gen [cur].p [i].br [k].ok  = false;
         Gen(gen [cur].p [i].br [k], d, full, true);
        }
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| интерпретатор префиксной записи                                  |
//+------------------------------------------------------------------+
double C_AO_GP::EvalRec(const S_GP_Branch &b, int &pos)
  {
   int op = b.n [pos].op;

   if(op == GP_OP_ERC)
     {
      double c = b.n [pos].val;
      pos++;
      return c;
     }

   pos++;
   double x = EvalRec(b, pos);
   double y = EvalRec(b, pos);
   double v = 0.0;

   switch(op)
     {
      case GP_OP_ADD:
         v = x + y;
         break;
      case GP_OP_SUB:
         v = x - y;
         break;
      case GP_OP_MUL:
         v = x * y;
         break;
      case GP_OP_DIV:
         v = (y == 0.0) ? 1.0 : x / y;   // защищённое деление
         break;
     }

   if(v >  GP_BIG)
      v =  GP_BIG;
   if(v < -GP_BIG)
      v = -GP_BIG;
   return v;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| значение ветви с кэшем: при кроссовере меняется одна ветвь       |
//+------------------------------------------------------------------+
double C_AO_GP::BranchValue(S_GP_Branch &b)
  {
   if(!b.ok)
     {
      int pos = 0;
      b.v  = EvalRec(b, pos);
      b.ok = true;
     }
   return b.v;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| обёртка: насыщение на [-1; 1], линейно в диапазон координаты     |
//+------------------------------------------------------------------+
double C_AO_GP::Wrap(double v, int k)
  {
   if(v >  1.0)
      v =  1.0;
   if(v < -1.0)
      v = -1.0;
   double x = rangeMin [k] + (v + 1.0) * 0.5 * (rangeMax [k] - rangeMin [k]);
   return u.SeInDiSp(x, rangeMin [k], rangeMax [k], rangeStep [k]);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| индекс сразу за поддеревом с корнем i                            |
//+------------------------------------------------------------------+
int C_AO_GP::SubtreeEnd(const S_GP_Branch &b, int i)
  {
   int need = 1;
   while(need > 0)
     {
      need += (b.n [i].op == GP_OP_ERC ? 0 : 2) - 1;
      i++;
     }
   return i;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
int C_AO_GP::DepthRec(const S_GP_Branch &b, int &pos)
  {
   int op = b.n [pos].op;
   pos++;
   if(op == GP_OP_ERC)
      return 0;
   int d1 = DepthRec(b, pos);
   int d2 = DepthRec(b, pos);
   return 1 + (d1 > d2 ? d1 : d2);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
int C_AO_GP::Depth(const S_GP_Branch &b)
  {
   int pos = 0;
   return DepthRec(b, pos);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| r-й (с нуля) внутренний узел или лист ветви                      |
//+------------------------------------------------------------------+
int C_AO_GP::NthNode(const S_GP_Branch &b, int r, bool internal)
  {
   for(int i = 0; i < b.len; i++)
     {
      if((b.n [i].op != GP_OP_ERC) == internal)
        {
         if(r == 0)
            return i;
         r--;
        }
     }
   return 0;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| точка в P1 по всей программе: правило 90/10, затем равномерно    |
//| среди всех внутренних узлов (или листьев) всех ветвей            |
//+------------------------------------------------------------------+
void C_AO_GP::PickInProgram(const S_GP_Prog &p, int &k, int &i)
  {
   bool wantInt = u.RNDprobab() < pInternal;
   int  totInt  = 0;
   int  totLeaf = 0;

   for(int b = 0; b < coords; b++)
     {
      totInt  += (p.br [b].len - 1) / 2;
      totLeaf += (p.br [b].len + 1) / 2;
     }

   if(wantInt && totInt == 0)
      wantInt = false;

   int r = RndInt(wantInt ? totInt : totLeaf);

   for(int b = 0; b < coords; b++)
     {
      int nb = wantInt ? (p.br [b].len - 1) / 2 : (p.br [b].len + 1) / 2;
      if(r < nb)
        {
         k = b;
         i = NthNode(p.br [b], r, wantInt);
         return;
        }
      r -= nb;
     }

   k = coords - 1;
   i = 0;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| точка в P2 внутри заданной ветви, правило 90/10                  |
//+------------------------------------------------------------------+
int C_AO_GP::PickInBranch(const S_GP_Branch &b)
  {
   int  nInt    = (b.len - 1) / 2;
   bool wantInt = nInt > 0 && u.RNDprobab() < pInternal;
   int  n       = wantInt ? nInt : (b.len + 1) / 2;
   return NthNode(b, RndInt(n), wantInt);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
void C_AO_GP::CopyBranch(S_GP_Branch &dst, const S_GP_Branch &src)
  {
   if(ArraySize(dst.n) < src.len)
      ArrayResize(dst.n, src.len, 64);
   for(int i = 0; i < src.len; i++)
      dst.n [i] = src.n [i];
   dst.len = src.len;
   dst.v   = src.v;
   dst.ok  = src.ok;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
void C_AO_GP::CopyProg(S_GP_Prog &dst, const S_GP_Prog &src)
  {
   for(int k = 0; k < coords; k++)
      CopyBranch(dst.br [k], src.br [k]);
   dst.f = src.f;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| dst = a[0..i) + b[j..je) + a[ie..len)                            |
//| dst не должен совпадать с a или b                                |
//+------------------------------------------------------------------+
void C_AO_GP::Splice(S_GP_Branch &dst,
                     const S_GP_Branch &ab, int i, int ie,
                     const S_GP_Branch &bb, int j, int je)
  {
   int newLen = i + (je - j) + (ab.len - ie);
   if(ArraySize(dst.n) < newLen)
      ArrayResize(dst.n, newLen, 64);

   int w = 0;
   for(int s = 0;  s < i;      s++)
      dst.n [w++] = ab.n [s];
   for(int s = j;  s < je;     s++)
      dst.n [w++] = bb.n [s];
   for(int s = ie; s < ab.len; s++)
      dst.n [w++] = ab.n [s];

   dst.len = newLen;
   dst.ok  = false;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| a_i = 1/(1+fScale*s_i), s_i = fMax - f_i  (s_i — не из статьи)   |
//| fScale = 1 — канон; > 1 растягивает шкалу f стенда для отбора    |
//+------------------------------------------------------------------+
void C_AO_GP::PrepRoulette()
  {
   double fMax = -DBL_MAX;
   for(int i = 0; i < popSize; i++)
      if(gen [cur].p [i].f > fMax)
         fMax = gen [cur].p [i].f;

   cumSum = 0.0;
   for(int i = 0; i < popSize; i++)
     {
      cumSum += 1.0 / (1.0 + fScale * (fMax - gen [cur].p [i].f));
      cum [i] = cumSum;
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| пропорциональный отбор, бинарный поиск по накопленной сумме      |
//+------------------------------------------------------------------+
int C_AO_GP::Spin()
  {
   double r  = u.RNDprobab() * cumSum;
   int    lo = 0;
   int    hi = popSize - 1;

   while(lo < hi)
     {
      int mid = (lo + hi) / 2;
      if(cum [mid] < r)
         lo = mid + 1;
      else
         hi = mid;
     }
   return lo;
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| кроссовер поддеревьев с типизацией ветвей;                       |
//| c2 < 0 — второй потомок не нужен (последнее место в поколении)   |
//+------------------------------------------------------------------+
void C_AO_GP::Crossover(int p1, int p2, int c1, int c2)
  {
   int nxt = cur ^ 1;
   int k, i;

   PickInProgram(gen [cur].p [p1], k, i);
   int ie = SubtreeEnd(gen [cur].p [p1].br [k], i);
   int j  = PickInBranch(gen [cur].p [p2].br [k]);
   int je = SubtreeEnd(gen [cur].p [p2].br [k], j);

//--- потомок 1: P1 с поддеревом из P2
   CopyProg(gen [nxt].p [c1], gen [cur].p [p1]);
   Splice(gen [nxt].p [c1].br [k],
          gen [cur].p [p1].br [k], i, ie,
          gen [cur].p [p2].br [k], j, je);
   if(Depth(gen [nxt].p [c1].br [k]) > maxDepth)
      CopyBranch(gen [nxt].p [c1].br [k], gen [cur].p [p1].br [k]);   // [GP I]

//--- потомок 2: P2 с поддеревом из P1
   if(c2 < 0)
      return;

   CopyProg(gen [nxt].p [c2], gen [cur].p [p2]);
   Splice(gen [nxt].p [c2].br [k],
          gen [cur].p [p2].br [k], j, je,
          gen [cur].p [p1].br [k], i, ie);
   if(Depth(gen [nxt].p [c2].br [k]) > maxDepth)
      CopyBranch(gen [nxt].p [c2].br [k], gen [cur].p [p2].br [k]);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| новое поколение: round((1-pCross)*M) копий, остальное кроссовер  |
//+------------------------------------------------------------------+
void C_AO_GP::Breed()
  {
   int nxt = cur ^ 1;
   PrepRoulette();

   int nR  = (int)MathRound((1.0 - pCross) * popSize);
   int cnt = 0;

   for(; cnt < nR && cnt < popSize; cnt++)
      CopyProg(gen [nxt].p [cnt], gen [cur].p [Spin()]);

   while(cnt < popSize)
     {
      int  p1  = Spin();
      int  p2  = Spin();
      bool two = cnt + 1 < popSize;

      Crossover(p1, p2, cnt, two ? cnt + 1 : -1);
      cnt += two ? 2 : 1;
     }
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                            Moving                                |
//+------------------------------------------------------------------+
void C_AO_GP::Moving()
  {
   if(!revision)
      Ramped();                  // поколение 0
   else
     {
      Breed();                   // следующее поколение из текущего
      cur ^= 1;
     }

   long nodes = 0;
   for(int i = 0; i < popSize; i++)
      for(int k = 0; k < coords; k++)
        {
         a [i].c [k] = Wrap(BranchValue(gen [cur].p [i].br [k]), k);
         nodes += gen [cur].p [i].br [k].len;
        }

   avgSize = (double)nodes / ((double)popSize * coords);
  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                           Revision                               |
//| поколенческая замена без элитизма: сохраняется только лучший     |
//| за прогон (cB / fB) — это и есть результат по Коза               |
//+------------------------------------------------------------------+
void C_AO_GP::Revision()
  {
   for(int i = 0; i < popSize; i++)
     {
      gen [cur].p [i].f = a [i].f;

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
