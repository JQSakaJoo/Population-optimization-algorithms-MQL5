//+————————————————————————————————————————————————————————————————————————————+
//|                                                                  C_AO_AOSm |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/16315

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
struct S_Layer
{
    int    pc;  // счетчик частиц
    double BEk; // энергия связи
    double LEk; // минимальная энергия

    void Init ()
    {
      pc  = 0;
      BEk = 0.0;
      LEk = 0.0;
    }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
struct S_Atom
{
    S_Layer layers [];  // массив слоев атома

    void Init (int layersNumb)
    {
      ArrayResize (layers, layersNumb);
      for (int i = 0; i < layersNumb; i++) layers [i].Init ();
    }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
struct S_Electron
{
    int layerID [];  // массив идентификаторов слоев для электрона

    void Init (int coords)
    {
      ArrayResize (layerID, coords);
    }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
// Класс реализующий алгоритм атомарной оптимизации (Atomic Orbital Search)
class C_AO_AOSm : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_AOSm () { }
  C_AO_AOSm ()
  {
    ao_name = "AOSm";
    ao_desc = "Atomic Orbital Search M";
    ao_link = "https://www.mql5.com/ru/articles/16315";

    popSize         = 50;      // размер популяции

    maxLayers       = 10;      // максимальное количество слоев
    photonEmissions = 20;      // количество излучений фотонов
    PR              = 0.1;     // коэффициент перехода фотонов

    ArrayResize (params, 4);

    params [0].name = "popSize";         params [0].val = popSize;
    params [1].name = "maxLayers";       params [1].val = maxLayers;
    params [2].name = "photonEmissions"; params [2].val = photonEmissions;
    params [3].name = "photonRate";      params [3].val = PR;
  }

  // Установка параметров алгоритма
  void SetParams ()
  {
    popSize         = (int)params [0].val;
    maxLayers       = (int)params [1].val;
    photonEmissions = (int)params [2].val;
    PR              = params      [3].val;
  }

  // Инициализация алгоритма с заданными параметрами поиска
  bool Init (const double &rangeMinP  [], // минимальный диапазон поиска
             const double &rangeMaxP  [], // максимальный диапазон поиска
             const double &rangeStepP [], // шаг поиска
             const int     epochsP = 0);  // количество эпох

  void Moving   (); // Метод перемещения частиц
  void Revision (); // Метод пересмотра лучших решений

  //----------------------------------------------------------------------------
  int    maxLayers;       // максимальное количество слоев
  int    photonEmissions; // количество излучений фотонов
  double PR;              // коэффициент перехода фотонов

  private: //-------------------------------------------------------------------
  int        atomStage;        // текущая стадия атома
  int        currentLayers []; // текущее количество слоев для соответствующего атома (координаты)
  S_Atom     atoms     [];     // атомы, размер соответствует количеству координат
  S_Electron electrons [];     // электроны, соответствует размеру популяции

  // Получение орбитальной полосы для заданных параметров
  int  GetOrbitBandID          (int layersNumb, double min, double max, double center, double p);

  // Распределение частиц в пространстве поиска
  void DistributeParticles     ();

  // Генерация слоев в атомах
  void LayersGenerationInAtoms ();

  // Обновление идентификаторов электронов
  void UpdateElectronsIDs      ();

  // Расчет параметров слоев
  void CalcLayerParams         ();

  // Обновление положения электронов
  void UpdateElectrons         ();
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
// Инициализация алгоритма
bool C_AO_AOSm::Init (const double &rangeMinP  [],
                      const double &rangeMaxP  [],
                      const double &rangeStepP [],
                      const int epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  atomStage     = 0;

  ArrayResize (currentLayers, coords);

  ArrayResize (atoms, coords);
  for (int i = 0; i < coords; i++) atoms [i].Init (maxLayers);

  ArrayResize (electrons,   popSize);
  for (int i = 0; i < popSize; i++) ArrayResize (electrons [i].layerID, coords);

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
// Метод перемещения частиц
void C_AO_AOSm::Moving ()
{
  // Начальное случайное позиционирование
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
  // Выполнение соответствующего этапа алгоритма
  if (atomStage == 0)
  {
    DistributeParticles ();
  }
  else
  {
    LayersGenerationInAtoms ();
    UpdateElectronsIDs      ();
    CalcLayerParams         ();
    UpdateElectrons         ();
  }

  // Переход к следующей стадии
  atomStage++;
  if (atomStage > photonEmissions) atomStage = 0;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
// Определение орбитальной полосы для частицы
int C_AO_AOSm::GetOrbitBandID (int layersNumb, double min, double max, double center, double p)
{
  // Расчет ширины полос слева и справа от центра
  double leftWidth  = (center - min) / layersNumb;
  double rightWidth = (max - center) / layersNumb;

  // Определение номера полосы
  if (p < center)
  {
    // Объект находится слева от центра
    for (int i = 1; i <= layersNumb; i++)
    {
      if (p >= center - i * leftWidth) return i - 1;
    }
    return layersNumb - 1; // Если объект находится в крайней левой полосе
  }
  else
    if (p > center)
    {
      // Объект находится справа от центра
      for (int i = 1; i <= layersNumb; i++)
      {
        if (p <= center + i * rightWidth) return i - 1;
      }
      return layersNumb - 1; // Если объект находится в крайней правой полосе
    }
    else
    {
      // Объект находится в центре
      return 0;
    }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
// Распределение частиц в пространстве поиска
void C_AO_AOSm::DistributeParticles ()
{
  for (int i = 0; i < popSize; i++)
  {
    for (int c = 0; c < coords; c++)
    {
      // Используем гауссово распределение для позиционирования частиц
      a [i].c [c] = u.GaussDistribution (cB [c], rangeMin [c], rangeMax [c], 8);
      a [i].c [c] = u.SeInDiSp (a [i].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
// Генерация случайного количества слоев для каждого атома
void C_AO_AOSm::LayersGenerationInAtoms ()
{
  for (int i = 0; i < coords; i++)
  {
    currentLayers [i] = u.RNDintInRange (1, maxLayers);
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
// Обновление идентификаторов слоев для каждого электрона
void C_AO_AOSm::UpdateElectronsIDs ()
{
  for (int i = 0; i < popSize; i++)
  {
    for (int c = 0; c < coords; c++)
    {
      electrons [i].layerID [c] = GetOrbitBandID (currentLayers [c], rangeMin [c], rangeMax [c], cB [c], a [i].c [c]);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
// Расчет параметров для каждого слоя
void C_AO_AOSm::CalcLayerParams ()
{
  double energy;

  // Обработка каждой координаты (атома)
  for (int c = 0; c < coords; c++)
  {
    atoms [c].Init (maxLayers);

    // Обработка каждого слоя
    for (int L = 0; L < currentLayers [c]; L++)
    {
      energy = -DBL_MAX;

      // Обработка каждого электрона
      for (int e = 0; e < popSize; e++)
      {
        if (electrons [e].layerID [c] == L)
        {
          atoms [c].layers [L].pc++;
          atoms [c].layers [L].BEk += a [e].f;

          if (a [e].f > energy)
          {
            energy = a [e].f;
            atoms [c].layers [L].LEk = a [e].c [c];
          }
        }
      }

      // Расчет средних значений для слоя
      if (atoms [c].layers [L].pc != 0)
      {
        atoms [c].layers [L].BEk /= atoms [c].layers [L].pc;
      }
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
// Обновление положения электронов
void C_AO_AOSm::UpdateElectrons ()
{
  double α;      // коэффициент скорости
  double φ;      // вероятность перехода
  double newPos; // новая позиция
  double LE;     // лучшая энергия
  int    lID;    // идентификатор слоя

  // Обработка каждой частицы
  for (int p = 0; p < popSize; p++)
  {
    for (int c = 0; c < coords; c++)
    {
      φ = u.RNDprobab ();

      if (φ < PR)
      {
        // Случайный переход к центру
        newPos = cB [c];
      }
      else
      {
        lID = electrons [p].layerID [c];

        α = u.RNDfromCI (-1.0, 1.0);

        // Если текущая энергия частицы меньше средней энергии слоя
        if (a [p].f < atoms [c].layers [lID].BEk)
        {
          // Движение в сторону глобального оптимума----------------------------
          LE     = cB [c];

          newPos = a [p].cB [c]+ α * (LE - a [p].cB [c]);
        }
        else
        {
          // Движение в сторону локального оптимума слоя------------------------
          LE     = atoms [c].layers [lID].LEk;

          newPos = a [p].cB [c]+ α * (LE - a [p].cB [c]);
        }
      }

      // Ограничение новой позиции в пределах диапазона поиска с учетом шага
      a [p].c [c] = u.SeInDiSp (newPos, rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
// Метод пересмотра лучших решений
void C_AO_AOSm::Revision ()
{
  int bestIndex = -1;

  // Поиск лучшего решения в текущей итерации
  for (int i = 0; i < popSize; i++)
  {
    // Обновление глобального лучшего решения
    if (a [i].f > fB)
    {
      fB = a [i].f;
      bestIndex = i;
    }

    // Обновление персонального лучшего решения
    if (a [i].f > a [i].fB)
    {
      a [i].fB = a [i].f;
      ArrayCopy (a [i].cB, a [i].c, 0, 0, WHOLE_ARRAY);
    }
  }

  // Обновление лучших координат если найдено лучшее решение
  if (bestIndex != -1) ArrayCopy (cB, a [bestIndex].c, 0, 0, WHOLE_ARRAY);
}
//——————————————————————————————————————————————————————————————————————————————