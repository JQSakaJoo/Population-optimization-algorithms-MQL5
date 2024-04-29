//+————————————————————————————————————————————————————————————————————————————+
//|                                                                       C_AO |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

#include "#C_AO.mqh"

#include "AO_BGA_Binary_Genetic_Algorithm.mqh"
#include "AO_(P_O)ES_Evolution_Strategies.mqh"
#include "AO_SDSm_StochasticDiffusionSearch.mqh"
#include "AO_ESG_Evolution_of_Social_Groups.mqh";
#include "AO_SIA_SimulatedIsotropicAnnealing.mqh";
//#include "AO_DE_DifferentialEvolution.mqh";
#include "AO_BSA_BirdSwarmAlgorithm.mqh"
//HS
//SSG
//(PO)ES
#include "AO_BSO_BrainStormOptimization.mqh";
#include "AO_WOA_WhaleOptimizationAlgorithm.mqh"
//ACOm
#include "AO_TSEA_TurtleShellEvolutionAlgorithm.mqh"
//BFO-GA
//MEC
//IWO
//Micro-AIS
//COAm
//SDOm
//NMm
//FAm
//GSA
//BFO
//ABC
//BA
//SA
//IWDm
//PSO
#include "AO_Boids_BoidsAlgorithm.mqh";
//MA
//SFL
//FSS
//RND
#include "AO_GWO_GreyWolfOptimizer.mqh"
//CSS
//EM

//——————————————————————————————————————————————————————————————————————————————
enum E_AO
{
  NONE_AO,

  AO_BGA,
  AO_P_O_ES,
  AO_SDSm,
  AO_ESG,
  AO_SIA,
  //AO_DE,
  AO_BSA,
  //HS
  //SSG
  //(PO)ES
  AO_BSO,
  AO_WOAm,
  //ACOm
  AO_TSEA,
  //BFO-GA
  //MEC
  //IWO
  //Micro-AIS
  //COAm
  //SDOm
  //NMm
  //FAm
  //GSA
  //BFO
  //ABC
  //BA
  //SA
  //IWDm
  //PSO
  AO_Boids,
  //MA
  //SFL
  //FSS
  //RND
  AO_GWO,
  //CSS
  //EM
};
C_AO *SelectAO (E_AO a)
{
  C_AO *ao;
  switch (a)
  {
    case  AO_BGA    : ao = new C_AO_BGA    (); return (GetPointer (ao));
    case  AO_P_O_ES : ao = new C_AO_P_O_ES (); return (GetPointer (ao));
    case  AO_SDSm   : ao = new C_AO_SDSm   (); return (GetPointer (ao));
    case  AO_ESG    : ao = new C_AO_ESG    (); return (GetPointer (ao));
    case  AO_SIA    : ao = new C_AO_SIA    (); return (GetPointer (ao));
    //case  AO_DE     : ao = new C_AO_DE     (); return (GetPointer (ao));
    case  AO_BSA    : ao = new C_AO_BSA    (); return (GetPointer (ao));
    //HS
    //SSG
    //(PO)ES
    case  AO_BSO    : ao = new C_AO_BSO    (); return (GetPointer (ao));
    case  AO_WOAm   : ao = new C_AO_WOAm   (); return (GetPointer (ao));
    //ACOm
    case  AO_TSEA   : ao = new C_AO_TSEA   (); return (GetPointer (ao));
    //BFO-GA
    //MEC
    //IWO
    //Micro-AIS
    //COAm
    //SDOm
    //NMm
    //FAm
    //GSA
    //BFO
    //ABC
    //BA
    //SA
    //IWDm
    //PSO
    case  AO_Boids  : ao = new C_AO_Boids  (); return (GetPointer (ao));
    //MA
    //SFL
    //FSS
    //RND
    case  AO_GWO    : ao = new C_AO_GWO    (); return (GetPointer (ao));
    //CSS
    //EM

    default:
      ao = NULL; return NULL;
  }
}
//——————————————————————————————————————————————————————————————————————————————