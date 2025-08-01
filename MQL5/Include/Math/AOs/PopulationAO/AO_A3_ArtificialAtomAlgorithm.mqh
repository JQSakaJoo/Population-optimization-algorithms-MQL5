//+------------------------------------------------------------------+
//|                                  Copyright 2007-2025, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//+------------------------------------------------------------------+
#include "#C_AO.mqh"
//————————————————————————————————————————————————————————————————————
class C_AO_A3 : public C_AO
{
  public: //----------------------------------------------------------
  ~C_AO_A3 () { }
  C_AO_A3 ()
  {
    ao_name = "A3";
    ao_desc = "Artificial Atom Algorithm";
    ao_link = "https://www.mql5.com/ru/articles/18958";

    popSize      = 10;    // количество атомов (m)
    covalentRate = 0.1;   // коэффициент ковалентной связи (β)

    ArrayResize (params, 2);

    params [0].name = "popSize";      params [0].val = popSize;
    params [1].name = "covalentRate"; params [1].val = covalentRate;
  }

  void SetParams ()
  {
    popSize      = (int)params [0].val;
    covalentRate = params      [1].val;
  }

  bool Init (const double &rangeMinP  [],  // минимальные значения
             const double &rangeMaxP  [],  // максимальные значения
             const double &rangeStepP [],  // шаг изменения
             const int     epochsP = 0);   // количество эпох

  void Moving   ();
  void Revision ();

  //------------------------------------------------------------------
  double covalentRate;       // коэффициент ковалентной связи (β)

  private: //---------------------------------------------------------
  int    covalentCount;      // количество атомов для ковалентной связи
};
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Инициализация
bool C_AO_A3::Init (const double &rangeMinP  [],
                    const double &rangeMaxP  [],
                    const double &rangeStepP [],
                    const int     epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //------------------------------------------------------------------
  covalentCount = (int)MathFloor (popSize * covalentRate);
  if (covalentCount < 1) covalentCount = 1;
  if (covalentCount >= popSize) covalentCount = popSize - 1;

  return true;
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Основной шаг алгоритма
void C_AO_A3::Moving ()
{
  // Начальная инициализация популяции
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

  //------------------------------------------------------------------
  int    ind = 0;

  for (int i = 0; i < popSize; i++)
  {
    if (i <= covalentCount)
    {
      for (int c = 0; c < coords; c++)
      {
        if (u.RNDprobab () < covalentRate)
        {
          a [i].c [c] = a [i].c [c] + u.RNDprobab () * (cB [c] - a [i].c [c]) * covalentCount;//(1.0 - covalentCount);
        }
        else
        {
          a [i].c [c] = u.PowerDistribution (cB [c], rangeMin [c], rangeMax [c], 20);
        }
        a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }
    else
    {
      for (int c = 0; c < coords; c++)
      {
        ind = u.RNDintInRange (0, covalentCount);

        a [i].c [c] = a [i].c [c] + u.RNDprobab () * (a [ind].c [c] - a [i].cW [c]) * (1.0 - covalentCount);//covalentCount;
        a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
//--- Обновление лучшего и худшего решений
void C_AO_A3::Revision ()
{
  // Обновляем локальное худшее решение
  for (int i = 0; i < popSize; i++)
  {
    if (a [i].f < a [i].fW)
    {
      fW = a [i].f;
      ArrayCopy (a [i].cW, a [i].c, 0, 0, WHOLE_ARRAY);
    }
  }

  static S_AO_Agent aT []; ArrayResize (aT, popSize);
  u.Sorting (a, aT, popSize);

  if (a [0].f > fB)
  {
    fB = a [0].f;
    ArrayCopy (cB, a [0].c, 0, 0, WHOLE_ARRAY);
  }

  if (a [popSize - 1].f < fW)
  {
    fW = a [popSize - 1].f;
    ArrayCopy (cW, a [popSize - 1].c, 0, 0, WHOLE_ARRAY);
  }
}
//————————————————————————————————————————————————————————————————————