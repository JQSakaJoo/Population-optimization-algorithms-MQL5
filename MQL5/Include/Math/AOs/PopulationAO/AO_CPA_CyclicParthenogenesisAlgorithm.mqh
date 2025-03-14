//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_CPA |
//|                                            Copyright 2007-2025, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/16877

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
// Класс реализующий Алгоритм Циклического Партеногенеза (CPA)
// Наследуется от базового класса оптимизации
class C_AO_CPA : public C_AO
{
  public:
  C_AO_CPA (void)
  {
    ao_name = "CPA";
    ao_desc = "Cyclic Parthenogenesis Algorithm";
    ao_link = "https://www.mql5.com/ru/articles/16877";

    popSize = 50;       // общий размер популяции Na

    Nc      = 10;       // количество колоний
    Fr      = 0.2;      // соотношение особей женского пола
    Pf      = 0.9;      // вероятность перелета между колониями
    alpha1  = 0.3;      // коэффициент масштабирования для партеногенеза
    alpha2  = 0.9;      // коэффициент масштабирования для спаривания

    ArrayResize (params, 6);

    // Установка параметров алгоритма
    params [0].name = "popSize";     params [0].val = popSize;
    params [1].name = "Nc";          params [1].val = Nc;
    params [2].name = "Fr";          params [2].val = Fr;
    params [3].name = "Pf";          params [3].val = Pf;
    params [4].name = "alpha1_init"; params [4].val = alpha1;
    params [5].name = "alpha2_init"; params [5].val = alpha2;
  }

  void SetParams ()
  {
    popSize = (int)params [0].val;

    Nc      = (int)params [1].val;
    Fr      = params      [2].val;
    Pf      = params      [3].val;
    alpha1  = params      [4].val;
    alpha2  = params      [5].val;
  }

  bool Init (const double &rangeMinP  [], // минимальный диапазон поиска
             const double &rangeMaxP  [], // максимальный диапазон поиска
             const double &rangeStepP [], // шаг поиска
             const int     epochsP = 0);  // количество эпох

  void Moving   ();         // функция перемещения особей
  void Revision ();         // функция пересмотра и обновления позиций

  //----------------------------------------------------------------------------
  int    Nc;                // количество колоний
  double Fr;                // соотношение особей женского пола
  double Pf;                // вероятность перелета между колониями

  private: //-------------------------------------------------------------------
  int    epochs;            // общее количество эпох
  int    epochNow;          // текущая эпоха
  int    Nm;                // количество особей в каждой колонии
  double alpha1;            // коэффициент масштабирования для партеногенеза
  double alpha2;            // коэффициент масштабирования для спаривания
  int    fNumber;           // количество особей женского пола в колонии
  int    mNumber;           // количество особей мужского пола в колонии

  S_AO_Agent aT [];         // временная колония для сортировки
  void SortFromTo (S_AO_Agent &p [], S_AO_Agent &pTemp [], int from, int count); // функция сортировки агентов
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
// Инициализация алгоритма с заданными параметрами поиска
bool C_AO_CPA::Init (const double &rangeMinP  [],
                     const double &rangeMaxP  [],
                     const double &rangeStepP [],
                     const int     epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  epochs   = epochsP;
  epochNow = 0;
  // Вычисление размера колонии и количества особей каждого пола
  Nm       = popSize / Nc;
  fNumber  = int(Nm * Fr); if (fNumber < 1) fNumber = 1;
  mNumber  = Nm - fNumber;

  ArrayResize (aT, Nm);

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
// Основная функция перемещения особей в пространстве поиска
void C_AO_CPA::Moving ()
{
  epochNow++;
  //----------------------------------------------------------------------------
  // Начальная случайная инициализация позиций, если это первая итерация
  if (!revision)
  {
    for (int i = 0; i < popSize; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        // Генерация случайной позиции в заданном диапазоне
        a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        a [i].c [c] = u.SeInDiSp  (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }

    revision = true;
    return;
  }

  //----------------------------------------------------------------------------
  // Вычисление коэффициента уменьшения силы поиска с течением времени
  double k    = (epochs - epochNow)/(double)epochs;
  int    ind  = 0;
  int    indF = 0;

  // Обработка каждой колонии
  for (int col = 0; col < Nc; col++)
  {
    // Обновление позиций особей женского пола (партеногенез)
    for (int f = 0; f < fNumber; f++)
    {
      ind = col * Nm + f;

      for (int c = 0; c < coords; c++)
      {
        // Партеногенетическое обновление позиции с использованием нормального распределения
        a [ind].c [c] = a [ind].cP [c] + alpha1 * k * u.GaussDistribution (0.0, -1.0, 1.0, 8) * (rangeMax [c] - rangeMin [c]);
      }
    }

    // Обновление позиций особей мужского пола (спаривание)
    for (int m = fNumber; m < Nm; m++)
    {
      ind = col * Nm + m;

      // Выбор случайной особи женского пола для спаривания
      indF = u.RNDintInRange (ind, col * Nm + fNumber - 1);

      for (int c = 0; c < coords; c++)
      {
        // Обновление позиции на основе выбранной особи женского пола
        a [ind].c [c] = a [ind].cP [c] + alpha2 * u.RNDprobab () * (a [indF].cP [c] - a [ind].cP [c]);
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
// Функция пересмотра позиций и обмена информацией между колониями
void C_AO_CPA::Revision ()
{
  // Поиск и обновление лучшего решения
  int ind = -1;

  for (int i = 0; i < popSize; i++)
  {
    if (a [i].f > fB)
    {
      fB = a [i].f;
      ind = i;
    }
  }

  if (ind != -1) ArrayCopy (cB, a [ind].c, 0, 0, WHOLE_ARRAY);

  //----------------------------------------------------------------------------
  // Сохранение текущих позиций
  for (int i = 0; i < popSize; i++)
  {
    ArrayCopy (a [i].cP, a [i].c, 0, 0, WHOLE_ARRAY);
  }

  // Сортировка особей в каждой колонии по значению целевой функции
  for (int col = 0; col < Nc; col++)
  {
    ind = col * Nm;
    SortFromTo (a, aT, ind, Nm);
  }

  // Механизм перелета (миграции) между колониями
  if (u.RNDprobab () < Pf)
  {
    int indCol_1 = 0;
    int indCol_2 = 0;

    // Выбор двух случайных различных колоний
    indCol_1 = u.RNDminusOne (Nc);
    do indCol_2 = u.RNDminusOne (Nc);
    while (indCol_1 == indCol_2);

    // Обеспечение, чтобы лучшее решение было в первой колонии
    if (a [indCol_1 * Nm].f < a [indCol_2 * Nm].f)
    {
      int temp = indCol_1;
      indCol_1 = indCol_2;
      indCol_2 = temp;
    }

    // Копирование лучшего решения в худшую колонию
    ArrayCopy (a [indCol_2 * Nm + Nm - 1].cP, a [indCol_1 * Nm].cP, 0, 0, WHOLE_ARRAY);

    // Пересортировка колонии после миграции
    SortFromTo (a, aT, indCol_2 * Nm, Nm);
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
// Вспомогательная функция сортировки агентов по значению целевой функции
void C_AO_CPA::SortFromTo (S_AO_Agent &p [], S_AO_Agent &pTemp [], int from, int count)
{
  int    cnt = 1;
  int    t0  = 0;
  double t1  = 0.0;
  int    ind [];
  double val [];
  ArrayResize (ind, count);
  ArrayResize (val, count);

  // Копирование значений для сортировки
  for (int i = 0; i < count; i++)
  {
    ind [i] = i + from;
    val [i] = p [i + from].f;
  }

  // Сортировка пузырьком по убыванию
  while (cnt > 0)
  {
    cnt = 0;
    for (int i = 0; i < count - 1; i++)
    {
      if (val [i] < val [i + 1])
      {
        // Обмен индексов
        t0 = ind [i + 1];
        ind [i + 1] = ind [i];
        ind [i] = t0;

        // Обмен значений
        t1 = val [i + 1];
        val [i + 1] = val [i];
        val [i] = t1;

        cnt++;
      }
    }
  }

  // Применение результатов сортировки
  for (int i = 0;    i < count; i++)        pTemp [i] = p [ind [i]];
  for (int i = from; i < from + count; i++) p     [i] = pTemp  [i - from];
}
//——————————————————————————————————————————————————————————————————————————————