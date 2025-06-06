//+————————————————————————————————————————————————————————————————————————————+
//|                                                                    C_AO_DA |
//|                                            Copyright 2007-2025, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/16999

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
// Класс реализующий диалектический алгоритм оптимизации
class C_AO_DA : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_DA () { }
  C_AO_DA ()
  {
    ao_name = "DA";
    ao_desc = "Dialectical Algorithm";
    ao_link = "https://www.mql5.com/ru/articles/16999";

    popSize = 50;       // размер популяции
    k1      = 40;       // спекулятивные мыслители
    k2      = 1;        // соседи

    ArrayResize (params, 3);
    params [0].name = "popSize"; params [0].val = popSize;
    params [1].name = "k1";      params [1].val = k1;
    params [2].name = "k2";      params [2].val = k2;
  }

  // Установка параметров алгоритма
  void SetParams ()
  {
    popSize = (int)params [0].val;
    k1      = (int)params [1].val;
    k2      = (int)params [2].val;
  }

  bool Init (const double &rangeMinP  [], // минимальный диапазон поиска
             const double &rangeMaxP  [], // максимальный диапазон поиска
             const double &rangeStepP [], // шаг поиска
             const int     epochsP = 0);  // количество эпох

  void Moving   ();    // Перемещение агентов в пространстве поиска
  void Revision ();    // Пересмотр и обновление лучших найденных решений

  //----------------------------------------------------------------------------
  int k1;       // количество спекулятивных мыслителей
  int k2;       // количество соседей для анализа

  private: //-------------------------------------------------------------------
  // Вычисление евклидового расстояния между двумя векторами
  double EuclideanDistance (const double &vec1 [], const double &vec2 [], const int dim)
  {
    double sum  = 0;
    double diff = 0.0;

    for (int i = 0; i < dim; i++)
    {
      diff = vec1 [i] - vec2 [i];
      sum += diff * diff;
    }
    return MathSqrt (sum);
  }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_DA::Init (const double &rangeMinP  [], // минимальный диапазон поиска
                    const double &rangeMaxP  [], // максимальный диапазон поиска
                    const double &rangeStepP [], // шаг поиска
                    const int     epochsP = 0)   // количество эпох
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
// Реализация движения агентов в пространстве поиска
void C_AO_DA::Moving ()
{
  //----------------------------------------------------------------------------
  // Начальная инициализация позиций агентов случайным образом
  if (!revision)
  {
    for (int i = 0; i < popSize; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }

    revision = true;
    return;
  }

  //----------------------------------------------------------------------------
  // Обновление позиции лучшего мыслителя
  for (int c = 0; c < coords; c++)
  {
    a [0].c [c] = a [0].cB [c] + u.RNDprobab () * (a [1].c [c] - a [0].c [c]);
    a [0].c [c] = u.SeInDiSp (a [0].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
  }

  //----------------------------------------------------------------------------
  double dist_next   = -DBL_MAX;  // максимальное расстояние до следующего соседа
  double dist_prev   = -DBL_MAX;  // максимальное расстояние до предыдущего соседа
  double dist        = 0.0;       // текущее расстояние
  int    antiNextIND = 0;         // индекс наиболее удаленного следующего соседа
  int    antiPrevIND = 0;         // индекс наиболее удаленного предыдущего соседа
  int    antiIND     = 0;         // выбранный индекс для обновления позиции

  // Обновление позиций спекулятивных мыслителей -------------------------------
  for (int i = k2; i < k1; i++)
  {
    // Поиск наиболее удаленного предыдущего соседа
    for (int j = 1; j <= k2; j++)
    {
      dist = EuclideanDistance (a [i].cB, a [i - j].cB, coords);
      if (dist > dist_prev)
      {
        dist_prev   = dist;
        antiPrevIND = i - j;
      }
    }

    // Поиск наиболее удаленного следующего соседа
    for (int j = 1; j <= k2; j++)
    {
      dist = EuclideanDistance (a [i].cB, a [i + j].cB, coords);
      if (dist > dist_next)
      {
        dist_next = dist;
        antiNextIND  = i + j;
      }
    }

    // Выбор наиболее удаленного соседа для обновления позиции
    if (dist_prev > dist_next) antiIND = antiPrevIND;
    else                       antiIND = antiNextIND;

    // Обновление координат спекулятивного мыслителя
    for (int c = 0; c < coords; c++)
    {
      a [i].c [c] = a [i].cB [c] + u.RNDbool () * (a [antiIND].c [c] - a [i].c [c]);
      //a [i].c [c] = a [i].cB [c] + u.RNDprobab () * (a [antiIND].c [c] - a [i].c [c]);
      a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }

  // Обновление позиций практических мыслителей --------------------------------
  for (int i = k1; i < popSize; i++)
  {
    // Случайный выбор двух спекулятивных мыслителей
    antiNextIND = u.RNDintInRange (0, k1 - 1);
    antiPrevIND = u.RNDintInRange (0, k1 - 1);

    if (antiNextIND == antiPrevIND) antiNextIND = antiPrevIND + 1;

    // Расчет расстояний до выбранных мыслителей
    dist_next = EuclideanDistance (a [i].cB, a [antiNextIND].cB, coords);
    dist_prev = EuclideanDistance (a [i].cB, a [antiPrevIND].cB, coords);

    // Выбор ближайшего мыслителя для обновления позиции
    if (dist_prev < dist_next) antiIND = antiPrevIND;
    else                       antiIND = antiNextIND;

    // Обновление координат практического мыслителя
    for (int c = 0; c < coords; c++)
    {
      a [i].c [c] = a [i].cB [c] + u.RNDprobab () * (a [antiIND].c [c] - a [i].c [c]);
      a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
// Пересмотр и обновление лучших найденных решений
void C_AO_DA::Revision ()
{
  int ind = -1;

  // Обновление лучших найденных решений для каждого агента
  for (int i = 0; i < popSize; i++)
  {
    // Обновление глобального лучшего решения
    if (a [i].f > fB)
    {
      fB = a [i].f;
      ind = i;
    }

    // Обновление личного лучшего решения агента
    if (a [i].f > a [i].fB)
    {
      a [i].fB = a [i].f;
      ArrayCopy (a [i].cB, a [i].c, 0, 0, WHOLE_ARRAY);
    }
  }

  // Обновление координат глобального лучшего решения
  if (ind != -1) ArrayCopy (cB, a [ind].c, 0, 0, WHOLE_ARRAY);

  // Сортировка агентов по их лучшим найденным решениям
  static S_AO_Agent aT []; ArrayResize (aT, popSize);
  u.Sorting_fB (a, aT, popSize);
}
//——————————————————————————————————————————————————————————————————————————————