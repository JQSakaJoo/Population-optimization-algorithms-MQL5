//+————————————————————————————————————————————————————————————————————————————+
//|                                                                       C_AO |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

#include "#C_AO.mqh"

#include "AO_ANS_AcrossNeighbourhoodSearch.mqh"
#include "AO_CLA_CodeLockAlgorithm.mqh"
#include "AO_AMOm_AnimalMigrationOptimization.mqh"
#include "AO_(P_O)ES_Evolution_Strategies.mqh"
#include "AO_CTA_CometTailAlgorithm.mqh"
#include "AO_SDSm_StochasticDiffusionSearch.mqh"
#include "AO_AAm_ArcheryAlgorithm.mqh"
#include "AO_ESG_Evolution_of_Social_Groups.mqh"
#include "AO_SIA_SimulatedIsotropicAnnealing.mqh"
#include "AO_ACS_ArtificialCooperativeSearch.mqh"
#include "AO_ASO_AnarchicSocietyOptimization.mqh"
#include "AO_AOSm_AtomicOrbitalSearch.mqh"
#include "AO_TSEA_TurtleShellEvolutionAlgorithm.mqh"
#include "AO_DE_DifferentialEvolution.mqh"
#include "AO_CRO_ChemicalReactionOptimisation.mqh"
#include "AO_BSA_BirdSwarmAlgorithm.mqh"
//HS
//SSG
#include "AO_BCOm_BacterialChemotaxisOptimization.mqh"
#include "AO_ABO_AfricanBuffaloOptimization.mqh"
#include "AO_(PO)ES_Evolution_Strategies.mqh"
#include "AO_TS_TabuSearch.mqh"
#include "AO_BSO_BrainStormOptimization.mqh";
#include "AO_WOA_WhaleOptimizationAlgorithm.mqh"
#include "AO_AEFA_ArtificialElectricFieldAlgorithm.mqh"
#include "AO_AEO_ArtificialEcosystemBasedOptimization.mqh"
//ACOm
//BFO-GA
#include "AO_SOA_SimpleOptimizationAlgorithm.mqh"
#include "AO_ABHA_ArtificialBeehiveAlgorithm.mqh"
#include "AO_ACMO_AtmosphereCloudsModelOptimization.mqh"
#include "AO_ADAMm_AdaptiveMomentEstimation.mqh"
#include "AO_ASHA_ArtificialShoweringAlgorithm.mqh"
#include "AO_ASBO_AdaptiveSocialBehaviorOptimization.mqh"
//MEC
//IWO
//Micro-AIS
#include "AO_COAm_CuckooOptimizationAlgorithm.mqh"
//SDOm
//NMm
//FAm
//GSA
//BFO
//ABC
//BA
#include "AO_AAA_ArtificialAlgaeAlgorithm.mqh"
//SA
//IWDm
//PSO
#include "AO_Boids_BoidsAlgorithm.mqh";
//MA
//SFL
//FSS
//RND
#include "AO_GWO_GreyWolfOptimizer.mqh"
#include "AO_AOA_ArithmeticOptimizationAlgorithm.mqh"
//CSS
//EM
#include "AO_BGA_Binary_Genetic_Algorithm.mqh"
#include "AO_RW_RandomWalk.mqh"

//——————————————————————————————————————————————————————————————————————————————
enum E_AO
{
  NONE_AO,

  AO_ANS,
  AO_CLA,
  AO_AMOm,
  AO_P_O_ES,
  AO_CTA,
  AO_SDSm,
  AO_AAm,
  AO_ESG,
  AO_SIA,
  AO_ACS,
  AO_ASO,
  AO_AOSm,
  AO_TSEA,
  AO_DE,
  AO_CRO,
  AO_BSA,
  //HS
  //SSG
  AO_BCOm,
  AO_ABO,
  AO_PO_ES,
  AO_TSm,
  AO_BSO,
  AO_WOAm,
  AO_AEFA,
  AO_AEO,
  //ACOm
  //BFO-GA
  AO_SOA,
  AO_ABHA,
  AO_ACMO,
  AO_ADAMm,
  AO_ASHA,
  AO_ASBO,
  //MEC
  //IWO
  //Micro-AIS
  AO_COAm,
  //SDOm
  //NMm
  //FAm
  //GSA
  //BFO
  //ABC
  //BA
  AO_AAA,
  //SA
  //IWDm
  //PSO
  AO_Boids,
  //MA
  //SFL
  //FSS
  //RND
  AO_GWO,
  AO_AOA,
  //CSS
  //EM
  AO_BGA,
  AO_RW,
};
C_AO *SelectAO (E_AO a)
{
  C_AO *ao;
  switch (a)
  {
    case  AO_ANS    : ao = new C_AO_ANS    (); return (GetPointer (ao));
    case  AO_CLA    : ao = new C_AO_CLA    (); return (GetPointer (ao));
    case  AO_AMOm   : ao = new C_AO_AMOm   (); return (GetPointer (ao));
    case  AO_P_O_ES : ao = new C_AO_P_O_ES (); return (GetPointer (ao));
    case  AO_CTA    : ao = new C_AO_CTA    (); return (GetPointer (ao));
    case  AO_SDSm   : ao = new C_AO_SDSm   (); return (GetPointer (ao));
    case  AO_AAm    : ao = new C_AO_AAm    (); return (GetPointer (ao));
    case  AO_ESG    : ao = new C_AO_ESG    (); return (GetPointer (ao));
    case  AO_SIA    : ao = new C_AO_SIA    (); return (GetPointer (ao));
    case  AO_ACS    : ao = new C_AO_ACS    (); return (GetPointer (ao));
    case  AO_ASO    : ao = new C_AO_ASO    (); return (GetPointer (ao));
    case  AO_AOSm   : ao = new C_AO_AOSm   (); return (GetPointer (ao));
    case  AO_TSEA   : ao = new C_AO_TSEA   (); return (GetPointer (ao));
    case  AO_DE     : ao = new C_AO_DE     (); return (GetPointer (ao));
    case  AO_CRO    : ao = new C_AO_CRO    (); return (GetPointer (ao));
    case  AO_BSA    : ao = new C_AO_BSA    (); return (GetPointer (ao));
    //HS
    //SSG
    case  AO_BCOm   : ao = new C_AO_BCOm   (); return (GetPointer (ao));
    case  AO_ABO    : ao = new C_AO_ABO    (); return (GetPointer (ao));
    case  AO_PO_ES  : ao = new C_AO_PO_ES  (); return (GetPointer (ao));
    case  AO_TSm    : ao = new C_AO_TSm    (); return (GetPointer (ao));
    case  AO_BSO    : ao = new C_AO_BSO    (); return (GetPointer (ao));
    case  AO_WOAm   : ao = new C_AO_WOAm   (); return (GetPointer (ao));
    case  AO_AEFA   : ao = new C_AO_AEFA   (); return (GetPointer (ao));
    case  AO_AEO    : ao = new C_AO_AEO    (); return (GetPointer (ao));
    //ACOm
    //BFO-GA
    case  AO_SOA    : ao = new C_AO_SOA    (); return (GetPointer (ao));
    case  AO_ABHA   : ao = new C_AO_ABHA   (); return (GetPointer (ao));
    case  AO_ACMO   : ao = new C_AO_ACMO   (); return (GetPointer (ao));
    case  AO_ADAMm  : ao = new C_AO_ADAMm  (); return (GetPointer (ao));
    case  AO_ASHA   : ao = new C_AO_ASHA   (); return (GetPointer (ao));
    case  AO_ASBO   : ao = new C_AO_ASBO   (); return (GetPointer (ao));
    //MEC
    //IWO
    //Micro-AIS
    case  AO_COAm   : ao = new C_AO_COAm   (); return (GetPointer (ao));
    //SDOm
    //NMm
    //FAm
    //GSA
    //BFO
    //ABC
    //BA
    case  AO_AAA    : ao = new C_AO_AAA    (); return (GetPointer (ao));
    //SA
    //IWDm
    //PSO
    case  AO_Boids  : ao = new C_AO_Boids  (); return (GetPointer (ao));
    //MA
    //SFL
    //FSS
    //RND
    case  AO_GWO    : ao = new C_AO_GWO    (); return (GetPointer (ao));
    case  AO_AOA    : ao = new C_AO_AOA    (); return (GetPointer (ao));
    //CSS
    //EM
    case  AO_BGA    : ao = new C_AO_BGA    (); return (GetPointer (ao));
    case  AO_RW     : ao = new C_AO_RW     (); return (GetPointer (ao));

    default:
      ao = NULL; return NULL;
  }
}
//——————————————————————————————————————————————————————————————————————————————