//+——————————————————————————————————————————————————————————————————+
//|                                                      C_AO_BChimp |
//|                                  Copyright 2007-2026, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

#include "#C_AO.mqh"

//+------------------------------------------------------------------+
//| Chimp Optimization Algorithm — континуальное ядро (BChimp).      |
//|                                                                  |
//| Источник: Ayeche, Alti — Enhanced Binary Chimp Optimization      |
//| (Hum-Cent Intell Syst, 2023), поверх базового ChOA (Khishe,      |
//| Mosavi, 2020). Исходник — бинарный feature-selection (BChimp1/   |
//| BChimp2). Для C_AO бинаризация снята: остаётся континуальное     |
//| ядро ChOA, позиция = среднее четырёх кандидатов от вожаков       |
//| (ядро BChimp2 до сигмоид-трансфера).                             |
//|                                                                  |
//| Содержательная механика (сохранена полностью):                   |
//|   - четвёрка вожаков: attacker / barrier / chaser / driver —     |
//|     четыре ЛУЧШИХ позиции (элитарно, best-ever);                 |
//|   - четыре группы с РАЗНЫМИ законами динамических коэффициентов  |
//|     (корневой ^1/3 и кубический ^3), что разводит баланс         |
//|     exploration/exploitation по вожакам;                         |
//|   - хаотический множитель m[i] по Piecewise-карте, статический   |
//|     на агента (как в источнике — генерируется один раз).         |
//|                                                                  |
//| Расчёт на агента/координату для каждого вожака k=0..3:           |
//|     A_k = a·C1G_k·(2·rand − 1)   — СИММЕТРИЧЕН (E=0)             |
//|     C_k = C2G_k·rand                                             |
//|     D_k = | C_k·Lead_k − m_i·x_i |                               |
//|     X_k = Lead_k − A_k·D_k                                       |
//|   x_new = (X_0 + X_1 + X_2 + X_3) / 4                            |
//| где a = 2 − 2·(t/T) (спад 2→0), амплитуда a·C1G_k → 0.           |
//|                                                                  |
//| Важно: в исходнике A = a·(2·C1G·rand − 1) c E[A]=a·(C1G−1)≠0 —   |
//| смещение сносит континуальную популяцию в угол rangeMin (в       |
//| бинарном домене скрыто сигмоид-трансфером). Здесь A приведён     |
//| к симметричной канонической форме ChOA: A = f·(2r−1), f=a·C1G.   |
//|                                                                  |
//| Версия ЭТАЛОННАЯ: greedy-приёмки/откатов нет (в источнике их     |
//| нет), популяция движется свободно, элитарны только вожаки.       |
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| вожак стаи: координаты + фитнес                                  |
//+------------------------------------------------------------------+
struct S_Leader
 {
  double c           [];   // координаты
  double             f;      // фитнес

  void               Init(const int coords)
   {
    ArrayResize(c, coords);
    f = -DBL_MAX;
   }
 };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//| персональное состояние шимпанзе                                  |
//+------------------------------------------------------------------+
struct S_Chimp
 {
  double             m;      // хаотический множитель (Piecewise-карта), статичен
 };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
class C_AO_BChimp : public C_AO
 {
public:
                    ~C_AO_BChimp() {}
                     C_AO_BChimp()
   {
    ao_name = "BChimp";
    ao_desc = "Chimp Optimization Algorithm";
    ao_link = "https://www.mql5.com/ru/articles/23132";

    popSize = 30;

    ArrayResize(params, 1);
    params [0].name = "popSize";
    params [0].val  = popSize;
   }

  void               SetParams()
   {
    popSize = (int)params [0].val;
   }

  bool               Init(const double &rangeMinP  [],
                          const double &rangeMaxP  [],
                          const double &rangeStepP [],
                          const int     epochsP = 0);

  void               Moving();
  void               Revision();

private:
  //--- четыре вожака: 0=attacker, 1=barrier, 2=chaser, 3=driver
  S_Leader           lead [4];

  //--- персональное состояние агентов
  S_Chimp            chimp [];

  //--- счётчики итераций (для отношения t/T)
  int                epochs;
  int                epochsDen;
  int                epochNow;

  //--- попытка вставить агента в элитарную лесенку вожаков
  void               InsertLeader(const int i);
 };
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                              Init                                |
//+------------------------------------------------------------------+
bool C_AO_BChimp::Init(const double &rangeMinP  [],
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

//--- вожаки
  for(int k = 0; k < 4; k++)
    lead [k].Init(coords);

//--- персональное состояние
  ArrayResize(chimp, popSize);

//--- хаотический множитель m[i] — Piecewise-карта (index=6 в источнике),
//    P=0.4, x0=0.7, Value=1; значение берётся ДО обновления
  double P = 0.4;
  double x = 0.7;
  for(int i = 0; i < popSize; i++)
   {
    chimp [i].m = x;
    if(x >= 0.0 && x < P)
      x = x / P;
    else
      if(x >= P && x < 0.5)
        x = (x - P) / (0.5 - P);
      else
        if(x >= 0.5 && x < 1.0 - P)
          x = (1.0 - P - x) / (0.5 - P);
        else
          x = (1.0 - x) / P;
   }

//--- счётчики итераций
  epochs    = epochsP;
  epochsDen = (epochs > 1) ? epochs - 1 : 1;
  epochNow  = 0;

  return true;
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                            Moving                                |
//|                                                                  |
//|  1. первый прогон: cP -> c (внешний цикл оценит FF);             |
//|  2. далее: a и четыре группы коэффициентов на итерацию;          |
//|  3. для каждого агента/координаты — четыре кандидата от вожаков, |
//|     усреднение, проекция в границы (SeInDiSp).                   |
//+------------------------------------------------------------------+
void C_AO_BChimp::Moving()
 {
//--- первый прогон: cP -> c
  if(!revision)
   {
    for(int i = 0; i < popSize; i++)
      for(int c = 0; c < coords; c++)
        a [i].c [c] = u.SeInDiSp(a [i].cP [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    return;
   }

//--- отношение t/T (зажато в 1.0); a: спад 2 -> 0
  epochNow++;
  double tRatio = (double)epochNow / (double)epochsDen;
  if(tRatio > 1.0)
    tRatio = 1.0;

  double aCoef = 2.0 - 2.0 * tRatio;

//--- групповые динамические коэффициенты
//    u1 = (t/T)^(1/3) — корневой закон,  u3 = (t/T)^3 — кубический
  double u1 = MathPow(tRatio, 1.0 / 3.0);
  double u3 = tRatio * tRatio * tRatio;

  double C1G [4], C2G [4];
  C1G [0] = 1.95 - 2.0 * u1;
  C2G [0] = 2.0 * u1 + 0.5;   // группа 1 (attacker)
  C1G [1] = 1.95 - 2.0 * u1;
  C2G [1] = 2.0 * u3 + 0.5;   // группа 2 (barrier)
  C1G [2] = -2.0 * u3 + 2.5;
  C2G [2] = 2.0 * u1 + 0.5;   // группа 3 (chaser)
  C1G [3] = -2.0 * u3 + 2.5;
  C2G [3] = 2.0 * u3 + 0.5;   // группа 4 (driver)

//--- движение
  for(int i = 0; i < popSize; i++)
   {
    double mi = chimp [i].m;

    for(int c = 0; c < coords; c++)
     {
      double xi  = a [i].c [c];
      double sum = 0.0;

      for(int k = 0; k < 4; k++)
       {
        //--- A симметричен относительно 0: амплитуда aCoef*C1G[k] -> 0
        //    к концу прогона (E[A]=0, нет сноса в угол). C масштабирован
        //    к канону [0..~2.5]. Знак C1G[k] роли не играет — множится
        //    на симметричное (2r-1).
        double Ak = aCoef * C1G [k] * (2.0 * u.RNDprobab() - 1.0);
        double Ck = C2G [k] * u.RNDprobab();
        double Dk = MathAbs(Ck * lead [k].c [c] - mi * xi);
        sum += lead [k].c [c] - Ak * Dk;
       }

      double v = sum / 4.0;
      a [i].c [c] = u.SeInDiSp (v, rangeMin [c], rangeMax [c], rangeStep [c]);
     }
   }
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                           Revision                               |
//|                                                                  |
//|  Обновление личных/глобального лучших + элитарная актуализация   |
//|  четвёрки вожаков (best-ever). Откатов нет — версия эталонная.   |
//+------------------------------------------------------------------+
void C_AO_BChimp::Revision()
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

//--- элитарная актуализация вожаков (инициализация на первом вызове
//    и поддержание top-4 на последующих)
  for(int i = 0; i < popSize; i++)
    InsertLeader(i);

  if(!revision)
    revision = true;
 }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|        InsertLeader — вставка агента в лесенку вожаков           |
//|  Вожаки отсортированы по убыванию фитнеса (lead[0] — лучший).    |
//+------------------------------------------------------------------+
void C_AO_BChimp::InsertLeader(const int i)
 {
  double f = a [i].f;

  int p = -1;
  if(f > lead [0].f)
    p = 0;
  else
    if(f > lead [1].f)
      p = 1;
    else
      if(f > lead [2].f)
        p = 2;
      else
        if(f > lead [3].f)
          p = 3;

  if(p < 0)
    return;

//--- сдвиг вниз от хвоста к позиции вставки
  for(int s = 3; s > p; s--)
   {
    lead [s].f = lead [s - 1].f;
    ArrayCopy(lead [s].c, lead [s - 1].c, 0, 0, coords);
   }

//--- вставка
  lead [p].f = f;
  ArrayCopy(lead [p].c, a [i].c, 0, 0, coords);
 }
//+------------------------------------------------------------------+
