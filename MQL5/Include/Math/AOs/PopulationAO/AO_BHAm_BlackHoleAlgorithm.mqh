//+————————————————————————————————————————————————————————————————————————————+
//|                                                                  C_AO_BHAm |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/16655

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
class C_AO_BHAm : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_BHAm () { }
  C_AO_BHAm ()
  {
    ao_name = "BHAm";
    ao_desc = "Black Hole Algorithm M";
    ao_link = "https://www.mql5.com/ru/articles/16655";

    popSize = 50;   // Размер популяции

    ArrayResize (params, 1);

    // Инициализация параметров
    params [0].name = "popSize"; params [0].val = popSize;
  }

  void SetParams () // Метод для установки параметров
  {
    popSize = (int)params [0].val;
  }

  bool Init (const double &rangeMinP  [], // Минимальный диапазон поиска
             const double &rangeMaxP  [], // Максимальный диапазон поиска
             const double &rangeStepP [], // Шаг поиска
             const int     epochsP = 0);  // Количество эпох

  void Moving   ();       // Метод перемещения
  void Revision ();       // Метод ревизии

  private: //-------------------------------------------------------------------
  int blackHoleIndex;    // Индекс черной дыры (лучшего решения)
};

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_BHAm::Init (const double &rangeMinP  [],
                      const double &rangeMaxP  [],
                      const double &rangeStepP [],
                      const int     epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  blackHoleIndex = 0; // Инициализация индекса черной дыры
  return true;
}
//——————————————————————————————————————————————————————————————————————————————
/*
//——————————————————————————————————————————————————————————————————————————————
void C_AO_BHA::Moving ()
{
  // Начальное случайное позиционирование при первом запуске
  if (!revision)
  {
    for (int i = 0; i < popSize; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        // Генерация случайной позиции в допустимом диапазоне
        a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        // Приведение к дискретным значениям согласно шагу
        a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }
    revision = true;
    return;
  }

  // Расчет суммы фитнес-значений для радиуса горизонта событий
  double sumFitness = 0.0;
  for (int i = 0; i < popSize; i++)
  {
    sumFitness += a [i].f;
  }

  // Расчет радиуса горизонта событий
  // R = fitBH / Σfiti
  double eventHorizonRadius = a [blackHoleIndex].f / sumFitness;

  // Обновление позиций звезд
  for (int i = 0; i < popSize; i++)
  {
    // Пропускаем черную дыру
    if (i == blackHoleIndex) continue;

    // Расчет расстояния до черной дыры
    double distance = 0.0;
    for (int c = 0; c < coords; c++)
    {
      double diff = a [blackHoleIndex].c [c] - a [i].c [c];
      distance += diff * diff;
    }
    distance = sqrt (distance);

    // Проверка пересечения горизонта событий
    if (distance < eventHorizonRadius)
    {
      // Звезда поглощена - создаем новую случайную звезду
      for (int c = 0; c < coords; c++)
      {
        a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }
    else
    {
      // Обновление позиции звезды по формуле:
      // Xi(t+1) = Xi(t) + rand × (XBH - Xi(t))
      for (int c = 0; c < coords; c++)
      {
        double rnd = u.RNDfromCI (0.0, 1.0);
        double newPosition = a [i].c [c] + rnd * (a [blackHoleIndex].c [c] - a [i].c [c]);

        // Проверка и коррекция границ
        newPosition = u.SeInDiSp (newPosition, rangeMin [c], rangeMax [c], rangeStep [c]);
        a [i].c [c] = newPosition;
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————
*/

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BHAm::Moving ()
{
  // Начальное случайное позиционирование при первом запуске
  if (!revision)
  {
    for (int i = 0; i < popSize; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        // Генерация случайной позиции в допустимом диапазоне
        a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        // Приведение к дискретным значениям согласно шагу
        a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }
    revision = true;
    return;
  }

  //----------------------------------------------------------------------------
  // Расчет среднего значения фитнес-значений для радиуса горизонта событий
  double aveFit = 0.0;
  double maxFit = fB;
  double minFit = a [0].f;

  for (int i = 0; i < popSize; i++)
  {
    aveFit += a [i].f;
    if (a [i].f < minFit) minFit = a [i].f;
  }
  aveFit /= popSize;

  // Расчет радиуса горизонта событий
  double eventHorizonRadius = (aveFit - minFit) / (maxFit - minFit);

  // Обновление позиций звезд
  for (int i = 0; i < popSize; i++)
  {
    // Пропускаем черную дыру
    if (i == blackHoleIndex) continue;

    for (int c = 0; c < coords; c++)
    {
      if (u.RNDprobab () < eventHorizonRadius * 0.01)
      {
        a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
      else
      {
        double rnd = u.RNDfromCI (0.0, 1.0);
        double newPosition = a [i].c [c] + rnd * (a [blackHoleIndex].c [c] - a [i].c [c]);

        // Проверка и коррекция границ
        newPosition = u.SeInDiSp (newPosition, rangeMin [c], rangeMax [c], rangeStep [c]);
        a [i].c [c] = newPosition;
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_BHAm::Revision ()
{
  // Поиск лучшего решения (черной дыры)
  for (int i = 0; i < popSize; i++)
  {
    if (a [i].f > fB)
    {
      fB = a [i].f;
      blackHoleIndex = i;
      ArrayCopy (cB, a [i].c, 0, 0, WHOLE_ARRAY);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————