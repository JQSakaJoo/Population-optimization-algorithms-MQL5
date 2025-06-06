//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_BIO |
//|                                            Copyright 2007-2025, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/17246

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
class C_AO_BIO : public C_AO
{
  public: //--------------------------------------------------------------------
  C_AO_BIO ()
  {
    ao_name = "BIO";
    ao_desc = "Blood Inheritance Optimization (joo)";
    ao_link = "https://www.mql5.com/ru/articles/17246";

    popSize = 50; // размер популяции

    ArrayResize (params, 1);
    params [0].name = "popSize"; params [0].val = popSize;
  }

  void SetParams ()
  {
    popSize = (int)params [0].val;
  }

  bool Init (const double &rangeMinP  [],  // минимальные значения
             const double &rangeMaxP  [],  // максимальные значения
             const double &rangeStepP [],  // шаг изменения
             const int     epochsP = 0);   // количество эпох

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------


  private: //-------------------------------------------------------------------
  struct S_Papa
  {
      int bTypes [];
  };
  struct S_Mama
  {
      S_Papa pa [4];
  };
  S_Mama ma [4];

  S_AO_Agent p [];

  int  GetBloodType     (int ind);
  void GetBloodMutation (double &gene, int indGene, int bloodType);
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_BIO::Init (const double &rangeMinP  [],
                     const double &rangeMaxP  [],
                     const double &rangeStepP [],
                     const int     epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  ArrayResize (p, popSize * 2);
  for (int i = 0; i < popSize * 2; i++) p [i].Init (coords);

  //1-1
  ArrayResize (ma [0].pa [0].bTypes, 1);

  ma [0].pa [0].bTypes [0] = 1;

  //2-2
  ArrayResize (ma [1].pa [1].bTypes, 2);

  ma [1].pa [1].bTypes [0] = 1;
  ma [1].pa [1].bTypes [1] = 2;

  //3-3
  ArrayResize (ma [2].pa [2].bTypes, 2);

  ma [2].pa [2].bTypes [0] = 1;
  ma [2].pa [2].bTypes [1] = 3;

  //1-2; 2-1
  ArrayResize (ma [0].pa [1].bTypes, 2);
  ArrayResize (ma [1].pa [0].bTypes, 2);

  ma [0].pa [1].bTypes [0] = 1;
  ma [0].pa [1].bTypes [1] = 2;

  ma [1].pa [0].bTypes [0] = 1;
  ma [1].pa [0].bTypes [1] = 2;

  //1-3; 3-1
  ArrayResize (ma [0].pa [2].bTypes, 2);
  ArrayResize (ma [2].pa [0].bTypes, 2);

  ma [0].pa [2].bTypes [0] = 1;
  ma [0].pa [2].bTypes [1] = 3;

  ma [2].pa [0].bTypes [0] = 1;
  ma [2].pa [0].bTypes [1] = 3;

  //1-4; 4-1
  ArrayResize (ma [0].pa [3].bTypes, 2);
  ArrayResize (ma [3].pa [0].bTypes, 2);

  ma [0].pa [3].bTypes [0] = 2;
  ma [0].pa [3].bTypes [1] = 3;

  ma [3].pa [0].bTypes [0] = 2;
  ma [3].pa [0].bTypes [1] = 3;

  //2-3; 3-2
  ArrayResize (ma [1].pa [2].bTypes, 4);
  ArrayResize (ma [2].pa [1].bTypes, 4);

  ma [1].pa [2].bTypes [0] = 1;
  ma [1].pa [2].bTypes [1] = 2;
  ma [1].pa [2].bTypes [2] = 3;
  ma [1].pa [2].bTypes [3] = 4;

  ma [2].pa [1].bTypes [0] = 1;
  ma [2].pa [1].bTypes [1] = 2;
  ma [2].pa [1].bTypes [2] = 3;
  ma [2].pa [1].bTypes [3] = 4;

  //2-4; 4-2; 3-4; 4-3; 4-4
  ArrayResize (ma [1].pa [3].bTypes, 3);
  ArrayResize (ma [3].pa [1].bTypes, 3);
  ArrayResize (ma [2].pa [3].bTypes, 3);
  ArrayResize (ma [3].pa [2].bTypes, 3);
  ArrayResize (ma [3].pa [3].bTypes, 3);

  ma [1].pa [3].bTypes [0] = 2;
  ma [1].pa [3].bTypes [1] = 3;
  ma [1].pa [3].bTypes [2] = 4;

  ma [3].pa [1].bTypes [0] = 2;
  ma [3].pa [1].bTypes [1] = 3;
  ma [3].pa [1].bTypes [2] = 4;

  ma [2].pa [3].bTypes [0] = 2;
  ma [2].pa [3].bTypes [1] = 3;
  ma [2].pa [3].bTypes [2] = 4;

  ma [3].pa [2].bTypes [0] = 2;
  ma [3].pa [2].bTypes [1] = 3;
  ma [3].pa [2].bTypes [2] = 4;

  ma [3].pa [3].bTypes [0] = 2;
  ma [3].pa [3].bTypes [1] = 3;
  ma [3].pa [3].bTypes [2] = 4;

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BIO::Moving ()
{
  //----------------------------------------------------------------------------
  if (!revision)
  {
    for (int i = 0; i < popSize; i++)
    {
      for (int j = 0; j < coords; j++)
      {
        a [i].c [j] = u.RNDfromCI (rangeMin [j], rangeMax [j]);
        a [i].c [j] = u.SeInDiSp (a [i].c [j], rangeMin [j], rangeMax [j], rangeStep [j]);
      }
    }
    revision = true;
    return;
  }

  //----------------------------------------------------------------------------
  double rnd        = 0.0;
  int    papIND     = 0;
  int    mamIND     = 0;
  int    pBloodType = 0;
  int    mBloodType = 0;
  int    cBloodType = 0;
  int    bloodIND   = 0;

  for (int i = 0; i < popSize; i++)
  {
    rnd = u.RNDprobab ();
    rnd *= rnd;
    papIND = (int)u.Scale (rnd, 0.0, 1.0, 0, popSize - 1);

    rnd = u.RNDprobab ();
    rnd *= rnd;
    mamIND = (int)u.Scale (rnd, 0.0, 1.0, 0, popSize - 1);

    pBloodType = GetBloodType (papIND);
    mBloodType = GetBloodType (mamIND);

    for (int c = 0; c < coords; c++)
    {
      bloodIND   = MathRand () % ArraySize (ma [mBloodType - 1].pa [pBloodType - 1].bTypes);
      cBloodType = ma [mBloodType - 1].pa [pBloodType - 1].bTypes [bloodIND];

      if (cBloodType == 1) a [i].c [c] = cB [c];
      else
      {
        if (u.RNDbool () < 0.5) a [i].c [c] = p [papIND].c [c];
        else                    a [i].c [c] = p [mamIND].c [c];

        GetBloodMutation (a [i].c [c], c, cBloodType);
        a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
int C_AO_BIO::GetBloodType (int ind)
{
  if (ind % 4 == 0) return 1;
  if (ind % 4 == 1) return 2;
  if (ind % 4 == 2) return 3;
  if (ind % 4 == 3) return 4;

  return 1;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void  C_AO_BIO::GetBloodMutation (double &gene, int indGene, int bloodType)
{
  switch (bloodType)
  {
    case 2:
      gene = u.PowerDistribution (gene, rangeMin [indGene], rangeMax [indGene], 20);
      return;
    case 3:
      gene += (cB [indGene] - gene) * u.RNDprobab ();
      return;
    default:
    {
      gene = rangeMax [indGene] - (gene - rangeMin [indGene]);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BIO::Revision ()
{
  //----------------------------------------------------------------------------
  for (int i = 0; i < popSize; i++)
  {
    // Обновляем лучшее глобальное решение
    if (a [i].f > fB)
    {
      fB = a [i].f;
      ArrayCopy (cB, a [i].c, 0, 0, WHOLE_ARRAY);
    }
  }

  //----------------------------------------------------------------------------
  for (int i = 0; i < popSize; i++)
  {
    p [popSize + i] = a [i];
  }

  S_AO_Agent pT []; ArrayResize (pT, popSize * 2);
  u.Sorting (p, pT, popSize * 2);
}
//——————————————————————————————————————————————————————————————————————————————