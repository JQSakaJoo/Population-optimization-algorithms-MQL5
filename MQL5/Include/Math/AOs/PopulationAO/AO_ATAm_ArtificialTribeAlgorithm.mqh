//+————————————————————————————————————————————————————————————————————————————+
//|                                                                  C_AO_ATAm |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/16588

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
class C_AO_ATAm : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_ATAm () { }
  C_AO_ATAm ()
  {
    ao_name = "ATAm";
    ao_desc = "Artificial Tribe Algorithm M";
    ao_link = "https://www.mql5.com/ru/articles/16588";

    popSize      = 50;   // Размер популяции
    AT_criterion = 0.9;  // Критерий оценки текущей ситуации
    AT_w         = 0.8;  // Глобальный инерционный вес

    ArrayResize (params, 3);

    // Инициализация параметров
    params [0].name = "popSize";      params [0].val = popSize;
    params [1].name = "AT_criterion"; params [1].val = AT_criterion;
    params [2].name = "AT_w";         params [2].val = AT_w;
  }

  void SetParams () // Метод для установки параметров
  {
    popSize      = (int)params [0].val;
    AT_criterion = params     [1].val;
    AT_w         = params     [2].val;
  }

  bool Init (const double &rangeMinP  [], // Минимальный диапазон поиска
             const double &rangeMaxP  [], // Максимальный диапазон поиска
             const double &rangeStepP [], // Шаг поиска
             const int     epochsP = 0);  // Количество эпох

  void Moving   ();       // Метод перемещения
  void Revision ();       // Метод ревизии

  //----------------------------------------------------------------------------
  double AT_criterion;    // Критерий оценки текущей ситуации
  double AT_w;            // Глобальный инерционный вес

  private: //-------------------------------------------------------------------
  double prevBestFitness; //Предыдущее лучшее решение
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_ATAm::Init (const double &rangeMinP  [], // Минимальный диапазон поиска
                      const double &rangeMaxP  [], // Максимальный диапазон поиска
                      const double &rangeStepP [], // Шаг поиска
                      const int     epochsP = 0)   // Количество эпох
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false; // Инициализация стандартных параметров

  //----------------------------------------------------------------------------
  prevBestFitness = -DBL_MAX;
  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_ATAm::Moving ()
{
  // Начальное случайное позиционирование
  if (!revision) // Если еще не было ревизии
  {
    for (int i = 0; i < popSize; i++) // Для каждой частицы
    {
      for (int c = 0; c < coords; c++) // Для каждой координаты
      {
        a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);                             // Генерация случайной позиции
        a [i].c [c] = u.SeInDiSp  (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]); // Приведение к дискретным значениям
      }
    }

    revision = true; // Установка флага ревизии
    return;          // Выход из метода
  }

  //----------------------------------------------------------------------------
  // Проверка критерия существования
  double diff = (fB - prevBestFitness) / (fB - fW);

  double Xi   = 0.0;
  double Xi_1 = 0.0;
  double Yi   = 0.0;
  double Yi_1 = 0.0;
  double Xs   = 0.0;
  double Xg   = 0.0;
  int    p    = 0;
  double r1   = 0.0;
  double r2   = 0.0;

  if (diff > AT_criterion)
  {
    // Поведение распространения (хорошая ситуация)
    for (int i = 0; i < popSize; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        p  = u.RNDminusOne (popSize);
        r1 = u.RNDprobab ();

        Xi = a [i].cP [c];
        Yi = a [p].cP [c];

        Xi_1 = r1 * Xi + (1.0 - r1) * Yi;
        Yi_1 = r1 * Yi + (1.0 - r1) * Xi;

        a [i].c [c] = u.SeInDiSp  (Xi_1, rangeMin [c], rangeMax [c], rangeStep [c]);
        a [p].c [c] = u.SeInDiSp  (Yi_1, rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }
  }
  else
  {
    // Поведение миграции (плохая ситуация)
    for (int i = 0; i < popSize; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        /*
        r1 = u.RNDprobab ();
        r2 = u.RNDprobab ();

        Xi = a [i].cP [c];
        Xs = a [i].cB [c];
        Xg = cB [c];

        Xi_1 = Xi + r1 * (Xs - Xi) + AT_w * r2 * (Xg - Xi);

        a [i].c [c] = u.SeInDiSp (Xi_1, rangeMin [c], rangeMax [c], rangeStep [c]);
        */
        
        if (u.RNDprobab () < diff)
        {
          Xi_1 = u.GaussDistribution (cB [c], rangeMin [c], rangeMax [c], 1);
          a [i].c [c] = u.SeInDiSp (Xi_1, rangeMin [c], rangeMax [c], rangeStep [c]);
        }
        else
        {
          r1 = u.RNDprobab ();
          r2 = u.RNDprobab ();

          Xi = a [i].cP [c];
          Xs = a [i].cB [c];
          Xg = cB [c];

          Xi_1 = Xi + r1 * (Xs - Xi) + AT_w * r2 * (Xg - Xi);

          a [i].c [c] = u.SeInDiSp (Xi_1, rangeMin [c], rangeMax [c], rangeStep [c]);
        }
        
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_ATAm::Revision ()
{
  //----------------------------------------------------------------------------
  int    indB  = -1;                // Индекс лучшей частицы
  double tempB = fB;

  for (int i = 0; i < popSize; i++) // Для каждой частицы
  {
    if (a [i].f > fB)               // Если значение функции лучше, чем текущее лучшее
    {
      fB   = a [i].f;               // Обновление лучшего значения функции
      indB = i;                     // Сохранение индекса лучшей частицы
    }

    if (a [i].f < fW)               // Если значение функции хуже, чем текущее худшее
    {
      fW   = a [i].f;               // Обновление худшего значения функции
    }

    if (a [i].f > a [i].fB)
    {
      a [i].fB = a [i].f;
      ArrayCopy (a [i].cB, a [i].c, 0, 0, WHOLE_ARRAY);
    }

    ArrayCopy (a [i].cP, a [i].c, 0, 0, WHOLE_ARRAY);
  }

  if (indB != -1)
  {
    prevBestFitness = tempB;
    ArrayCopy (cB, a [indB].c, 0, 0, WHOLE_ARRAY); // Копирование координат лучшей частицы
  }
}
//——————————————————————————————————————————————————————————————————————————————
