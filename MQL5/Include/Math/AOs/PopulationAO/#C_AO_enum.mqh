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
#include "AO_TETA_TimeEvolutionTravelAlgorithm.mqh"
#include "AO_SDSm_StochasticDiffusionSearch.mqh"
#include "AO_BOAm_BilliardsOptimizationAlgorithm.mqh"
#include "AO_AAm_ArcheryAlgorithm.mqh"
#include "AO_ESG_Evolution_of_Social_Groups.mqh"
#include "AO_SIA_SimulatedIsotropicAnnealing.mqh"
#include "AO_EOm_ExtremalOptimization.mqh"
#include "AO_BBO_BiogeographyBasedOptimization.mqh"
#include "AO_ACS_ArtificialCooperativeSearch.mqh"
#include "AO_DA_DialecticalAlgorithm.mqh"
#include "AO_BHAm_BlackHoleAlgorithm.mqh"
#include "AO_ASO_AnarchicSocietyOptimization.mqh"
#include "AO_RFO_RoyalFlushOptimization.mqh"
#include "AO_AOSm_AtomicOrbitalSearch.mqh"
#include "AO_TSEA_TurtleShellEvolutionAlgorithm.mqh"
#include "AO_BSA_BacktrackingSearchAlgorithm.mqh"
#include "AO_DE_DifferentialEvolution.mqh"
#include "AO_SRA_SuccessfulRestaurateurAlgorithm.mqh"
#include "AO_CRO_ChemicalReactionOptimisation.mqh"
#include "AO_BIO_BloodInheritanceOptimization.mqh"
#include "AO_BSA_BirdSwarmAlgorithm.mqh"
#include "AO_DEA_DolphinEcholocationAlgorithm.mqh"
//HS
//SSG
#include "AO_BCOm_BacterialChemotaxisOptimization.mqh"
#include "AO_ABO_AfricanBuffaloOptimization.mqh"
#include "AO_(PO)ES_Evolution_Strategies.mqh"
#include "AO_FBA_FractalBasedAlgorithm.mqh"
#include "AO_TS_TabuSearch.mqh"
#include "AO_BSO_BrainStormOptimization.mqh";
#include "AO_WOA_WhaleOptimizationAlgorithm.mqh"
#include "AO_AEFA_ArtificialElectricFieldAlgorithm.mqh"
#include "AO_AEO_ArtificialEcosystemBasedOptimization.mqh"
#include "AO_CAm_CamelAlgorithm.mqh"
//ACOm
#include "AO_CMAES_CovarianceMatrixAdaptationEvolutionStrategy.mqh"
//BFO-GA
#include "AO_SOA_SimpleOptimizationAlgorithm.mqh"
#include "AO_ABHA_ArtificialBeehiveAlgorithm.mqh"
#include "AO_CLA_L_CompetitiveLearningAlgorithm.mqh"
#include "AO_ACMO_AtmosphereCloudsModelOptimization.mqh"
#include "AO_ADAMm_AdaptiveMomentEstimation.mqh"
#include "AO_CoSO_CommunityOfScientistOptimization.mqh"
#include "AO_CGO_ChaosGameOptimization.mqh"
#include "AO_ATAm_ArtificialTribeAlgorithm.mqh"
#include "AO_A3_ArtificialAtomAlgorithm.mqh"
#include "AO_CROm_CoralReefsOptimization.mqh"
#include "AO_CFO_CentralForceOptimization.mqh"
#include "AO_ASHA_ArtificialShoweringAlgorithm.mqh"
#include "AO_ASBO_AdaptiveSocialBehaviorOptimization.mqh"
#include "AO_ES_EagleStrategy.mqh"
#include "AO_BRO_BattleRoyaleOptimizer.mqh"
//MEC
#include "AO_CSA_CircleSearchAlgorithm.mqh"
//IWO
//Micro-AIS
#include "AO_DOS_DeterministicOscillatorySearch.mqh"
#include "AO_EMA_ExchangeMarketAlgorithm.mqh"
#include "AO_COAm_CuckooOptimizationAlgorithm.mqh"
//SDOm
//NMm
#include "AO_COA_ChaosOptimizationAlgorithm.mqh"
#include "AO_BBBC_BigBangBigCrunch.mqh"
#include "AO_CPA_CyclicParthenogenesisAlgorithm.mqh"
//FAm
//GSA
//BFO
//ABC
//BA
#include "AO_AAA_ArtificialAlgaeAlgorithm.mqh"
//SA
//IWDm
//PSO
#include "AO_Boids_BoidsAlgorithm.mqh"
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

  AO_ANS,        //ANS     (Across Neighbourhood Search)
  AO_CLA,        //CLA     (Code Lock Algorithm, joo)
  AO_AMOm,       //AMOm    (Animal Migration Optimization, M)
  AO_P_O_ES,     //(P+O)ES (Evolution Strategies)
  AO_CTA,        //CTA     (Comet Tail Algorithm, joo)
  AO_TETA,       //TETA    (Time Evolution Travel Algorithm, joo)
  AO_SDSm,       //SDSm    (Stochastic Diffusion Search, M)
  AO_BOAm,       //BOAm    (Billiards Optimization Algorithm, M)
  AO_AAm,        //AAm     (Archery Algorithm, M)
  AO_ESG,        //ESG     (Evolution of Social Groups, joo)
  AO_SIA,        //SIA     (Simulated Isotropic Annealing, joo)
  AO_EOm,        //EOm     (Extremal Optimization M)
  AO_BBO,        //BBO     (Biogeography Based Optimization)
  AO_ACS,        //ACS     (Artificial Cooperative Search)
  AO_DA,         //DA      (Dialectical Algorithm)
  AO_BHAm,       //BHAm    (Black Hole Algorithm, M)
  AO_ASO,        //ASO     (Anarchic Society Optimization)
  AO_RFO,        //RFO     (Royal Flush Optimization)
  AO_AOSm,       //AOSm    (Atomic Orbital Search, M)
  AO_TSEA,       //TSEA    (Turtle Shell Evolution Algorithm, joo)
  AO_BSA_Backtr, //BSA     (Backtracking Search Algorithm)
  AO_DE,         //DE      (Differential Evolution)
  AO_SRA,        //SRA     (Successful Restaurateur Algorithm, joo)
  AO_CRO,        //CRO     (Chemical Reaction Optimisation)
  AO_BIO,        //BIO     (Blood Inheritance Optimization, joo)
  AO_BSA,        //BSA     (Bird Swarm Algorithm)
  AO_DEA,        //DEA     (Dolphin Echolocation Algorithm)
  //HS
  //SSG
  AO_BCOm,       //BCOm    (Bacterial Chemotaxis Optimization, M)
  AO_ABO,        //ABO     (African Buffalo Optimization)
  AO_PO_ES,      //(PO)ES  (Evolution_Strategies)
  AO_FBA,        //FBA     (Fractal Based Algorithm)
  AO_TSm,        //TS      (Tabu Search, M)
  AO_BSO,        //BSO     (Brain Storm Optimization)
  AO_WOAm,       //WOA     (Whale Optimization Algorithm)
  AO_AEFA,       //AEFA    (Artificial Electric Field Algorithm)
  AO_AEO,        //AEO     (Artificial Ecosystem Based Optimization)
  AO_CAm,        //CAm     (Camel Algorithm, M)
  //ACOm
  AO_CMAES,      //CMAES   (Covariance Matrix Adaptation Evolution Strategy)
  //BFO-GA
  AO_SOA,        //SOA     (Simple Optimization Algorithm)
  AO_ABHA,       //ABHA    (Artificial Beehive Algorithm)
  AO_CLA_l,      //CLA_l   (Competitive Learning Algorithm)
  AO_ACMO,       //ACMO    (Atmosphere Clouds Model Optimization)
  AO_ADAMm,      //ADAMm   (Adaptive Moment Estimation, M)
  CoSO,          //CoSO    (Community of Scientist Optimization)
  AO_CGO,        //CGO     (Chaos Game Optimization)
  AO_ATAm,       //ATAm    (Artificial Tribe Algorithm, M)
  AO_A3,         //A3      (Artificial Atom Algorithm)
  AO_CROm_coral, //CROm    (Coral Reefs Optimization, M)
  AO_CFO,        //CFO     (Central Force Optimization)
  AO_ASHA,       //ASHA    (Artificial Showering Algorithm)
  AO_ASBO,       //ASBO    (Adaptive Social Behavior Optimization)
  AO_ES,         //ES      (Eagle Strategy)
  AO_BRO,        //BRO     (Battle Royale Optimizer)
  //MEC
  AO_CSA,        //CSA     (Circle Search Algorithm)
  //IWO
  //Micro-AIS
  AO_DOS,        //DOS     (Deterministic Oscillatory Search)
  AO_EMA,        //EMA     (Exchange Market Algorithm)
  AO_COAm,       //COAm    (Cuckoo Optimization Algorithm, M)
  //SDOm
  //NMm
  AO_COA_chaos,  //COA     (Chaos Optimization Algorithm)
  AO_BBBC,       //BBBC    (BigBang Big Crunch)
  AO_CPA,        //CPA     (Cyclic Parthenogenesis Algorithm)
  //FAm
  //GSA
  //BFO
  //ABC
  //BA
  AO_AAA,        //AAA     (Artificial Algae Algorithm)
  //SA
  //IWDm
  //PSO
  AO_Boids,      //Boids   (Boids Algorithm)
  //MA
  //SFL
  //FSS
  //RND
  AO_GWO,        //GWO     (Grey Wolf Optimizer)
  AO_AOA,        //AOA     (Arithmetic Optimization Algorithm)
  //CSS
  //EM
  AO_BGA,        //BGA     (Binary Genetic Algorithm)
  AO_RW,         //RW      (Random Walk)
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
    case  AO_TETA   : ao = new C_AO_TETA   (); return (GetPointer (ao));
    case  AO_SDSm   : ao = new C_AO_SDSm   (); return (GetPointer (ao));
    case  AO_BOAm   : ao = new C_AO_BOAm   (); return (GetPointer (ao));
    case  AO_AAm    : ao = new C_AO_AAm    (); return (GetPointer (ao));
    case  AO_ESG    : ao = new C_AO_ESG    (); return (GetPointer (ao));
    case  AO_SIA    : ao = new C_AO_SIA    (); return (GetPointer (ao));
    case  AO_EOm    : ao = new C_AO_EOm    (); return (GetPointer (ao));
    case  AO_BBO    : ao = new C_AO_BBO    (); return (GetPointer (ao));
    case  AO_ACS    : ao = new C_AO_ACS    (); return (GetPointer (ao));
    case  AO_DA     : ao = new C_AO_DA     (); return (GetPointer (ao));
    case  AO_BHAm   : ao = new C_AO_BHAm   (); return (GetPointer (ao));
    case  AO_ASO    : ao = new C_AO_ASO    (); return (GetPointer (ao));
    case  AO_RFO    : ao = new C_AO_RFO    (); return (GetPointer (ao));
    case  AO_AOSm   : ao = new C_AO_AOSm   (); return (GetPointer (ao));
    case  AO_TSEA   : ao = new C_AO_TSEA   (); return (GetPointer (ao));
    case  AO_BSA_Backtr : ao = new C_AO_BSA_Backtracking (); return (GetPointer (ao));
    case  AO_DE     : ao = new C_AO_DE     (); return (GetPointer (ao));
    case  AO_SRA    : ao = new C_AO_SRA    (); return (GetPointer (ao));
    case  AO_CRO    : ao = new C_AO_CRO    (); return (GetPointer (ao));
    case  AO_BIO    : ao = new C_AO_BIO    (); return (GetPointer (ao));
    case  AO_BSA    : ao = new C_AO_BSA    (); return (GetPointer (ao));
    case  AO_DEA    : ao = new C_AO_DEA    (); return (GetPointer (ao));
    //HS
    //SSG
    case  AO_BCOm   : ao = new C_AO_BCOm   (); return (GetPointer (ao));
    case  AO_ABO    : ao = new C_AO_ABO    (); return (GetPointer (ao));
    case  AO_PO_ES  : ao = new C_AO_PO_ES  (); return (GetPointer (ao));
    case  AO_FBA    : ao = new C_AO_FBA    (); return (GetPointer (ao));
    case  AO_TSm    : ao = new C_AO_TSm    (); return (GetPointer (ao));
    case  AO_BSO    : ao = new C_AO_BSO    (); return (GetPointer (ao));
    case  AO_WOAm   : ao = new C_AO_WOAm   (); return (GetPointer (ao));
    case  AO_AEFA   : ao = new C_AO_AEFA   (); return (GetPointer (ao));
    case  AO_AEO    : ao = new C_AO_AEO    (); return (GetPointer (ao));
    case  AO_CAm    : ao = new C_AO_CAm    (); return (GetPointer (ao));
    //ACOm
    case  AO_CMAES  : ao = new C_AO_CMAES  (); return (GetPointer (ao));
    //BFO-GA
    case  AO_SOA    : ao = new C_AO_SOA    (); return (GetPointer (ao));
    case  AO_ABHA   : ao = new C_AO_ABHA   (); return (GetPointer (ao));
    case  AO_CLA_l  : ao = new C_AO_CLA_l  (); return (GetPointer (ao));
    case  AO_ACMO   : ao = new C_AO_ACMO   (); return (GetPointer (ao));
    case  AO_ADAMm  : ao = new C_AO_ADAMm  (); return (GetPointer (ao));
    case  AO_CGO    : ao = new C_AO_CGO    (); return (GetPointer (ao));
    case  AO_ATAm   : ao = new C_AO_ATAm   (); return (GetPointer (ao));
    case  AO_A3     : ao = new C_AO_A3     (); return (GetPointer (ao));
    case  AO_CROm_coral: ao = new C_AO_CROm (); return(GetPointer (ao));
    case  AO_CFO    : ao = new C_AO_CFO    (); return (GetPointer (ao));
    case  AO_ASHA   : ao = new C_AO_ASHA   (); return (GetPointer (ao));
    case  AO_ASBO   : ao = new C_AO_ASBO   (); return (GetPointer (ao));
    case  AO_ES     : ao = new C_AO_ES     (); return (GetPointer (ao));
    case  AO_BRO    : ao = new C_AO_BRO    (); return (GetPointer (ao));
    //MEC
    case  AO_CSA    : ao = new C_AO_CSA    (); return (GetPointer (ao));
    //IWO
    //Micro-AIS
    case  AO_DOS    : ao = new C_AO_DOS    (); return (GetPointer (ao));
    case  AO_EMA    : ao = new C_AO_EMA    (); return (GetPointer (ao));
    case  AO_COAm   : ao = new C_AO_COAm   (); return (GetPointer (ao));
    //SDOm
    //NMm
    case  AO_COA_chaos : ao = new C_AO_COA_chaos (); return (GetPointer (ao));
    case  AO_BBBC   : ao = new C_AO_BBBC   (); return (GetPointer (ao));
    case  AO_CPA    : ao = new C_AO_CPA    (); return (GetPointer (ao));
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