//+------------------------------------------------------------------+
//|                                  Copyright 2007-2026, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//+------------------------------------------------------------------+
#property copyright "Copyright 2007-2026, JQS aka Joo."
#property link      "https://www.mql5.com/ru/users/joo"
#property version   "1.00"
#property script_show_inputs

#include <Math\AOs\TestFunctions.mqh>
#include <Math\AOs\PopulationAO\#C_AO.mqh>
#include <Math\AOs\PopulationAO\#C_AO_enum.mqh>

//————————————————————————————————————————————————————————————————————
input string AOparam            = "----------------"; //AO parameters-----------
input E_AO   AOexactly_P        = NONE_AO;

input string TestStand_P        = "----------------"; //Test stand---------------
input int    NumbTestFuncRuns_P = 1e4;   // Количество вычислений функции
input int    NumberRepetTest_P  = 10;    // Количество повторных прогонов
//————————————————————————————————————————————————————————————————————


//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void OnStart()
 {
  C_AO *AO = SelectAO(AOexactly_P);

  if(AO == NULL)
   {
    Print("AO is not selected...");
    return;
   }

  Print(AO.GetName(), "|", AO.GetDesc(), "|", AO.GetParams());

//--- создаём 5 функций
#define NFUNCS 5
  C_Function *funcs [NFUNCS];
  funcs [0] = new C_Hilly();
  funcs [1] = new C_Forest();
  funcs [2] = new C_Megacity();
  funcs [3] = new C_Peaks();
  funcs [4] = new C_Skin();

//--- 10 координат: пара 0-1 = Hilly, 2-3 = Forest, 4-5 = Megacity,
//                   6-7 = Peaks, 8-9 = Skin
  int params = NFUNCS * 2;

  double rangeMin  [];
  double rangeMax  [];
  double rangeStep [];
  ArrayResize(rangeMin,  params);
  ArrayResize(rangeMax,  params);
  ArrayResize(rangeStep, params);

  for(int i = 0; i < NFUNCS; i++)
   {
    rangeMin  [i * 2]     = funcs [i].GetMinRangeX();
    rangeMax  [i * 2]     = funcs [i].GetMaxRangeX();
    rangeStep [i * 2]     = 0.0;

    rangeMin  [i * 2 + 1] = funcs [i].GetMinRangeY();
    rangeMax  [i * 2 + 1] = funcs [i].GetMaxRangeY();
    rangeStep [i * 2 + 1] = 0.0;
   }

//--- запуск тестов
  int    epochCount = NumbTestFuncRuns_P / (int)AO.params [0].val;
  double aveResult  = 0.0;

  Print("=============================");
  Print("Composite test: Hilly + Forest + Megacity + Peaks + Skin");
  Print("Coordinates: ", params, "; Epochs: ", epochCount,
        "; Repeats: ", NumberRepetTest_P);
  Print("=============================");

  for(int test = 0; test < NumberRepetTest_P; test++)
   {
    if(!AO.Init(rangeMin, rangeMax, rangeStep, epochCount))
      break;

    for(int epochCNT = 1; epochCNT <= epochCount && !IsStopped(); epochCNT++)
     {
      AO.Moving();

      //--- расчёт fitness для каждого агента
      for(int set = 0; set < ArraySize(AO.a); set++)
       {
        double sum = 0.0;
        for(int i = 0; i < NFUNCS; i++)
         {
          sum += funcs [i].Core(AO.a [set].c [i * 2],
                                AO.a [set].c [i * 2 + 1]);
         }
        AO.a [set].f = sum / NFUNCS;
       }

      AO.Revision();
     }

    Print("Run ", test + 1, "/", NumberRepetTest_P, ": ", AO.fB);
    aveResult += AO.fB;
   }

  aveResult /= (double)NumberRepetTest_P;

  Print("=============================");
  Print("Average result: ", DoubleToString(aveResult, 10),
        " (", DoubleToString(aveResult * 100.0, 2), "%)");
  Print("=============================");

//--- освобождение
  for(int i = 0; i < NFUNCS; i++)
    delete funcs [i];
  delete AO;
 }
//+------------------------------------------------------------------+
//+------------------------------------------------------------------+
