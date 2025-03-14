//+————————————————————————————————————————————————————————————————————————————+
//|                                                                  C_AO_TETA |
//|                                            Copyright 2007-2025, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/16963

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
// TETA - Time Evolution Travel Algorithm
// Алгоритм оптимизации, основанный на концепции перемещения между параллельными вселенными
// через изменение ключевых временых якорей (событий) в жизни
class C_AO_TETA : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_TETA () { }
  C_AO_TETA ()
  {
    ao_name = "TETA";
    ao_desc = "Time Evolution Travel Algorithm (joo)";
    ao_link = "https://www.mql5.com/ru/articles/16963";

    popSize = 50; // количество параллельных вселенных в популяции

    ArrayResize (params, 1);
    params [0].name = "popSize"; params [0].val = popSize;
  }

  void SetParams ()
  {
    popSize = (int)params [0].val;
  }

  bool Init (const double &rangeMinP  [],  // минимальные значения для якорей
             const double &rangeMaxP  [],  // максимальные значения для якорей
             const double &rangeStepP [],  // шаг изменения якорей
             const int     epochsP = 0);   // количество эпох поиска

  void Moving ();
  void Revision ();

  private: //-------------------------------------------------------------------
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_TETA::Init (const double &rangeMinP  [], // минимальные значения для якорей
                      const double &rangeMaxP  [], // максимальные значения для якорей
                      const double &rangeStepP [], // шаг изменения якорей
                      const int     epochsP = 0)   // количество эпох поиска
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_TETA::Moving ()
{
  //----------------------------------------------------------------------------
  if (!revision)
  {
    // Инициализация начальных значений якорей во всех параллельных вселенных
    for (int i = 0; i < popSize; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        a [i].c [c] = u.SeInDiSp  (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }

    revision = true;
    return;
  }

  //----------------------------------------------------------------------------
  double rnd  = 0.0;
  double val  = 0.0;
  int    pair = 0.0;

  for (int i = 0; i < popSize; i++)
  {
    for (int c = 0; c < coords; c++)
    {
      // Генерация вероятности, которая определяет как шанс выбора вселенной,
      // так и силу изменения якорей
      rnd  = u.RNDprobab ();
      rnd *= rnd;

      // Выбор вселенной для обмена опытом
      pair = (int)u.Scale (rnd, 0.0, 1.0, 0, popSize - 1);

      if (i != pair)
      {
        if (i < pair)
        {
          // Если текущая вселенная более благоприятна:
          // Небольшое изменение якоря (пропорционально rnd) для поиска лучшего баланса
          val = a [i].c [c] + (rnd) * (a [pair].cB [c] - a [i].cB [c]);
        }
        else
        {
          if (u.RNDprobab () > rnd)
          {
            // Если текущая вселенная менее благоприятна:
            // Значительное изменение якоря (пропорционально 1.0 - rnd)
            val = a [i].cB [c] + (1.0 - rnd) * (a [pair].cB [c] - a [i].cB [c]);
          }
          else
          {
            // Полное принятие конфигурации якоря из более успешной вселенной
            val = a [pair].cB [c];
          }
        }
      }
      else
      {
        // Локальная настройка якоря через гауссово распределение
        val = u.GaussDistribution (cB [c], rangeMin [c], rangeMax [c], 1);
      }

      a [i].c [c] = u.SeInDiSp  (val, rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_TETA::Revision ()
{
  for (int i = 0; i < popSize; i++)
  {
    // Обновление глобально лучшей конфигурации якорей
    if (a [i].f > fB)
    {
      fB = a [i].f;
      ArrayCopy (cB, a [i].c);
    }

    // Обновление лучшей известной конфигурации якорей для каждой вселенной
    if (a [i].f > a [i].fB)
    {
      a [i].fB = a [i].f;
      ArrayCopy (a [i].cB, a [i].c);
    }
  }

  // Сортировка вселенных по степени их благоприятности
  static S_AO_Agent aT []; ArrayResize (aT, popSize);
  u.Sorting_fB (a, aT, popSize);
}
//——————————————————————————————————————————————————————————————————————————————