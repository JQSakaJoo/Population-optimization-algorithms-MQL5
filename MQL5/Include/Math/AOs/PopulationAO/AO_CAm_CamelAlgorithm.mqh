//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_CAm |
//|                                            Copyright 2007-2025, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/18057

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
class C_AO_CAm : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_CAm () { }
  C_AO_CAm ()
  {
    ao_name = "CAm";
    ao_desc = "Camel Algorithm M";
    ao_link = "https://www.mql5.com/ru/articles/18057";

    popSize   = 50;     // размер популяции (camel caravan)
    Tmin      = 50;     // минимальная температура
    Tmax      = 100;    // максимальная температура
    omega     = 0.8;    // фактор нагрузки для supply
    dyingRate = 0.01;   // скорость "смерти" верблюдов
    alpha     = 0.9;    // параметр видимости для эффекта оазиса

    ArrayResize (params, 6);

    params [0].name = "popSize";    params [0].val = popSize;
    params [1].name = "Tmin";       params [1].val = Tmin;
    params [2].name = "Tmax";       params [2].val = Tmax;
    params [3].name = "omega";      params [3].val = omega;
    params [4].name = "dyingRate";  params [4].val = dyingRate;
    params [5].name = "alpha";      params [5].val = alpha;
  }

  void SetParams ()
  {
    popSize   = (int)params [0].val;
    Tmin      = params      [1].val;
    Tmax      = params      [2].val;
    omega     = params      [3].val;
    dyingRate = params      [4].val;
    alpha     = params      [5].val;
  }

  bool Init (const double &rangeMinP  [],  // минимальные значения
             const double &rangeMaxP  [],  // максимальные значения
             const double &rangeStepP [],  // шаг изменения
             const int     epochsP = 0);   // количество эпох

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  double Tmin;         // минимальная температура
  double Tmax;         // максимальная температура
  double omega;        // фактор нагрузки для supply
  double dyingRate;    // скорость "смерти" верблюдов
  double alpha;        // параметр видимости для эффекта оазиса

  private: //-------------------------------------------------------------------
  double temperature [];   // текущая температура для каждого верблюда
  double supply      [];   // текущий запас воды и пищи для каждого верблюда
  double endurance   [];   // текущая выносливость для каждого верблюда
  double initialSupply;    // начальный запас (обычно 1.0)
  double initialEndurance; // начальная выносливость (обычно 1.0)
  int    traveledSteps;    // количество пройденных шагов
  int    totalSteps;       // общее количество шагов

  // Вспомогательные методы
  void   InitializePopulation ();
  void   UpdateFactors        ();
  void   UpdatePositions      ();
  void   ApplyOasisEffect     ();
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_CAm::Init (const double &rangeMinP  [],  // минимальные значения
                     const double &rangeMaxP  [],  // максимальные значения
                     const double &rangeStepP [],  // шаг изменения
                     const int     epochsP = 0)    // количество эпох
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  // Инициализация массивов для каждого верблюда
  ArrayResize (temperature, popSize);
  ArrayResize (supply,      popSize);
  ArrayResize (endurance,   popSize);

  // Установка начальных значений
  initialSupply    = 1.0;
  initialEndurance = 1.0;
  traveledSteps    = 0;
  totalSteps       = epochsP;
  
  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Основной метод оптимизации                                                 |
//+----------------------------------------------------------------------------+
void C_AO_CAm::Moving ()
{
  // Первая итерация - инициализация начальной популяции
  if (!revision)
  {
    InitializePopulation ();
    revision = true;
    return;
  }

  // Увеличиваем счетчик пройденных шагов
  traveledSteps++;

  // Основной процесс оптимизации
  // 1. Обновляем факторы (температура, запас, выносливость)
  UpdateFactors ();

  // 2. Обновляем позиции верблюдов
  UpdatePositions ();

  // 3. Применяем эффект оазиса (обновление запасов и выносливости)
  ApplyOasisEffect ();

  // 4. Сохраняем состояние верблюдов
  for (int i = 0; i < popSize; i++) a [i].fP = a [i].f;
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Обновление лучшего решения                                                 |
//+----------------------------------------------------------------------------+
void C_AO_CAm::Revision ()
{
  // Поиск лучшего решения в текущей популяции
  for (int i = 0; i < popSize; i++)
  {
    // Обновление лучшего решения
    if (a [i].f > fB)
    {
      fB = a [i].f;
      ArrayCopy (cB, a [i].c, 0, 0, WHOLE_ARRAY);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Инициализация начальной популяции                                          |
//+----------------------------------------------------------------------------+
void C_AO_CAm::InitializePopulation ()
{
  // Инициализация начальной популяции равномерно по всему пространству
  for (int i = 0; i < popSize; i++)
  {
    for (int c = 0; c < coords; c++)
    {
      // Генерация случайных координат в допустимых пределах
      a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
      // Округление до ближайшего допустимого шага
      a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }

    // Инициализация факторов для каждого верблюда
    supply      [i] = initialSupply;
    endurance   [i] = initialEndurance;
  }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Обновление факторов (температура, запас, выносливость)                     |
//+----------------------------------------------------------------------------+
void C_AO_CAm::UpdateFactors ()
{
  double journeyRatio = (double)traveledSteps / (double)totalSteps;

  for (int i = 0; i < popSize; i++)
  {
    // Обновление температуры - случайное значение в диапазоне [Tmin, Tmax],
    // формула (1): Tnow = (Tmax - Tmin) * Rand(0,1) + Tmin
    temperature [i] = u.RNDfromCI (Tmin, Tmax);

    // Обновление запаса - уменьшается с течением времени,
    // формула (2): Snow = Spast * (1 - ω * Traveled steps / Total journey steps)
    supply [i] = supply [i] * (1.0 - omega * journeyRatio);

    // Обновление выносливости - зависит от температуры и времени,
    // формула (4): Enow = Epast * (1 - Tnow/Tmax) * (1 - Traveled steps / Total journey steps)
    double temperatureRatio = temperature [i] / Tmax;
    endurance [i] = endurance [i] * (1.0 - temperatureRatio) * (1.0 - journeyRatio);
  }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Обновление позиций верблюдов                                               |
//+----------------------------------------------------------------------------+
void C_AO_CAm::UpdatePositions ()
{
  for (int i = 0; i < popSize; i++)
  {
    /*
    // Проверка на "смерть" верблюда (quicksand, storm, etc.)
    if (u.RNDprobab () < dyingRate)
    {
      // Генерируем новую позицию случайным образом
      for (int c = 0; c < coords; c++)
      {
        a [i].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      }
      continue;
    }
    */

    // Обновление позиции-------------------------------------------------------
    double delta = u.RNDfromCI (-1.0, 1.0); // Фактор случайного блуждания

    // Обновляем каждую координату
    for (int c = 0; c < coords; c++)
    {
      /*
      // Применяем формулу обновления из статьи
      double enduranceFactor = (1.0 - endurance [i] / initialEndurance);
      double supplyFactor    = MathExp (1.0 - supply [i] / initialSupply);

      // Обновление позиции
      a [i].c [c] = a [i].c [c] + delta * enduranceFactor * supplyFactor * (cB [c] - a [i].c [c]);

      // Проверка на выход за границы и корректировка до допустимого значения
      a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
      */

      // Проверка на "смерть" верблюда (quicksand, storm, etc.)
      if (u.RNDprobab () < dyingRate)
      {
        // Генерируем новую позицию относительно координаты оазиса по нормальному распределению
        a [i].c [c] = u.GaussDistribution (cB [c], rangeMin [c], rangeMax [c], 8);
      }
      else
      {
        // Применяем формулу обновления из статьи
        double enduranceFactor = (1.0 - endurance [i] / initialEndurance);
        double supplyFactor    = MathExp (1.0 - supply [i] / initialSupply);

        // Обновление позиции
        a [i].c [c] = a [i].c [c] + delta * enduranceFactor * supplyFactor * (cB [c] - a [i].c [c]);
      }

      // Проверка на выход за границы и корректировка до допустимого значения
      a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//+----------------------------------------------------------------------------+
//| Применение эффекта оазиса (обновление запасов и выносливости)              |
//+----------------------------------------------------------------------------+
void C_AO_CAm::ApplyOasisEffect ()
{
  for (int i = 0; i < popSize; i++)
  {
    // Условие для обнаружения оазиса:
    // 1) Верблюд должен "видеть" оазис (случайная вероятность, зависящая от alpha)
    // 2) Текущее решение должно быть лучше, чем в предыдущей итерации

    if (u.RNDprobab () > (1.0 - alpha) && a [i].f > a [i].fP)
    {
      // Обнаружен оазис, пополняем запасы и выносливость
      supply    [i] = initialSupply;
      endurance [i] = initialEndurance;
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————
