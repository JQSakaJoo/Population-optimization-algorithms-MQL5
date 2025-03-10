//+————————————————————————————————————————————————————————————————————————————+
//|                                                                  C_AO_BOAm |
//|                                            Copyright 2007-2025, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/17325

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
class C_AO_BOAm : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_BOAm () { }
  C_AO_BOAm ()
  {
    ao_name = "BOAm";
    ao_desc = "Billiards Optimization Algorithm M";
    ao_link = "https://www.mql5.com/ru/articles/17325";

    popSize    = 50;  // число шаров (агентов)
    numPockets = 25;  // число карманов на бильярдном столе

    ArrayResize (params, 2);

    params [0].name = "popSize";    params [0].val = popSize;
    params [1].name = "numPockets"; params [1].val = numPockets;
  }

  void SetParams ()
  {
    popSize    = (int)params [0].val;
    numPockets = (int)params [1].val;
  }

  bool Init (const double &rangeMinP  [],  // минимальные значения
             const double &rangeMaxP  [],  // максимальные значения
             const double &rangeStepP [],  // шаг изменения
             const int     epochsP = 0);   // количество эпох

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  int numPockets;       // число карманов (лучших решений)

  private: //-------------------------------------------------------------------
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
//--- Инициализация
bool C_AO_BOAm::Init (const double &rangeMinP  [],
                      const double &rangeMaxP  [],
                      const double &rangeStepP [],
                      const int epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
//--- Основной шаг алгоритма
void C_AO_BOAm::Moving ()
{
  //----------------------------------------------------------------------------
  // Начальная инициализация
  if (!revision)
  {
    for (int p = 0; p < popSize; p++)
    {
      for (int c = 0; c < coords; c++)
      {
        a [p].c  [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        a [p].c  [c] = u.SeInDiSp (a [p].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
        a [p].cB [c] = a [p].c [c];  // Сохраняем начальную позицию
      }
    }

    revision = true;
    return;
  }

  //----------------------------------------------------------------------------
  for (int p = 0; p < popSize; p++)
  {
    for (int c = 0; c < coords; c++)
    {
      int pocketID = u.RNDminusOne (numPockets);

      //a [p].c [c] = a [p].cB [c] + u.RNDprobab () * (a [pocketID].cB [c] - u.RNDintInRange (1, 2) * a [p].cB [c]);
      a [p].c [c] = a [p].cB [c] + u.RNDprobab () * (a [pocketID].cB [c] - a [p].cB [c]) * u.RNDintInRange (1, 2);
      a [p].c [c] = u.SeInDiSp (a [p].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
//--- Обновление лучшего решения с учетом жадного выбора и вероятности принятия худших решений
void C_AO_BOAm::Revision ()
{
  int bestIND = -1;

  for (int i = 0; i < popSize; i++)
  {
    if (a [i].f > fB)
    {
      fB = a [i].f;
      bestIND = i;
    }

    if (a [i].f > a [i].fB)
    {
      a [i].fB = a [i].f;
      ArrayCopy (a [i].cB, a [i].c, 0, 0, WHOLE_ARRAY);
    }
  }

  if (bestIND != -1) ArrayCopy (cB, a [bestIND].c, 0, 0, WHOLE_ARRAY);

  S_AO_Agent aT []; ArrayResize (aT, popSize);
  u.Sorting_fB (a, aT, popSize);
}
//——————————————————————————————————————————————————————————————————————————————