//+————————————————————————————————————————————————————————————————————————————+
//|                                                                       C_AO |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

#include "#C_AO.mqh"


#include "AO_BGA_Binary_Genetic_Algorithm.mqh"                                  //1 BGA
#include "AO_(P_O)ES_Evolution_Strategies.mqh"                                  //2 (P+O)ES
#include "AO_SDSm_StochasticDiffusionSearch.mqh"                                //3 SDSm
//4 ESG
//5 SIA
//6 DE
#include "AO_BSA_BirdSwarmAlgorithm.mqh"                                        //7 BSA
//8 HS
//9 SSG
//10 (PO)ES
#include "AO_WOA_WhaleOptimizationAlgorithm.mqh"                                //11 WOAm
//12 ACOm
//13 BFO-GA
//14 MEC
//15 IWO
//16 Micro-AIS
//17 COAm
//18 SDOm
//19 NMm
//20 FAm
//21 GSA
//22 BFO
//23 ABC
//24 BA
//25 SA
//26 IWDm
//27 PSO
//28 MA
//29 SFL
//30 FSS
//31 RND
#include "AO_GWO_GreyWolfOptimizer.mqh"                                         //32 GWO
//33 CSS
//34 EM







//#include "AO_ESG_Evolution_of_Social_Groups.mqh";
//#include "AO_DE_Differential_Evolution.mqh"
//#include "AO_(PO)ES_Evolution_Strategies.mqh"
//#include "AO_COAm_Cuckoo_Optimization_Algorithm.mqh"


//——————————————————————————————————————————————————————————————————————————————
enum E_AO
{
  AO_BGA,                                                                       //1 BGA
  AO_P_O_ES,                                                                    //2 (P+O)ES
  AO_SDSm,                                                                      //3 SDSm
  //4 ESG
  //5 SIA
  //6 DE
  AO_BSA,                                                                       //7 BSA
  //8 HS
  //9 SSG
  //10 (PO)ES
  AO_WOAm,                                                                      //11 WOAm
  //12 ACOm
  //13 BFO-GA
  //14 MEC
  //15 IWO
  //16 Micro-AIS
  //17 COAm
  //18 SDOm
  //19 NMm
  //20 FAm
  //21 GSA
  //22 BFO
  //23 ABC
  //24 BA
  //25 SA
  //26 IWDm
  //27 PSO
  //28 MA
  //29 SFL
  //30 FSS
  //31 RND
  AO_GWO,                                                                       //32 GWO
  //33 CSS
  //34 EM

  //AO_ESG,
  //AO_DE,

  //AO_PO_ES,
  //AO_COAm,

  AO_NONE
};
C_AO *SelectAO (E_AO a)
{
  C_AO *ao;
  switch (a)
  {
    case  AO_BGA    : ao = new C_AO_BGA    (); return (GetPointer (ao));        //1 BGA
    case  AO_P_O_ES : ao = new C_AO_P_O_ES (); return (GetPointer (ao));        //2 (P+O)ES
    case  AO_SDSm   : ao = new C_AO_SDSm   (); return (GetPointer (ao));        //3 SDSm
    case  AO_BSA    : ao = new C_AO_BSA    (); return (GetPointer (ao));        //7 BSA
    case  AO_WOAm   : ao = new C_AO_WOAm   (); return (GetPointer (ao));        //11 WOAm
    case  AO_GWO    : ao = new C_AO_GWO    (); return (GetPointer (ao));        //32 GWO

    default: ao = NULL; return NULL;
  }
}
//——————————————————————————————————————————————————————————————————————————————