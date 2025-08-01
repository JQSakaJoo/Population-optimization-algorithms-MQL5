//+——————————————————————————————————————————————————————————————————+
//|                                                       C_AO_CLA_l |
//|                                  Copyright 2007-2025, Andrey Dik |
//|                                https://www.mql5.com/ru/users/joo |
//———————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/18857

#include "#C_AO.mqh"

//————————————————————————————————————————————————————————————————————
class C_AO_CLA_l : public C_AO
{
  public: //----------------------------------------------------------
  ~C_AO_CLA_l () { }
  C_AO_CLA_l ()
  {
    ao_name = "CLA_L";
    ao_desc = "Competitive Learning Algorithm";
    ao_link = "https://www.mql5.com/ru/articles/18857";

    popSize        = 198;
    numClasses     = 3;
    beta           = 0.3;
    gamma          = 0.8;
    deltaIter      = 2;

    ArrayResize (params, 5);

    params [0].name = "popSize";     params [0].val = popSize;
    params [1].name = "numClasses";  params [1].val = numClasses;
    params [2].name = "beta";        params [2].val = beta;
    params [3].name = "gamma";       params [3].val = gamma;
    params [4].name = "deltaIter";   params [4].val = deltaIter;
  }

  void SetParams ()
  {
    popSize    = (int)params [0].val;
    numClasses = (int)params [1].val;
    beta       = params      [2].val;
    gamma      = params      [3].val;
    deltaIter  = (int)params [4].val;
  }

  bool Init (const double &rangeMinP  [],
             const double &rangeMaxP  [],
             const double &rangeStepP [],
             const int     epochsP = 0);

  void Moving   ();
  void Revision ();
  void Injection (const int popPos, const int coordPos, const double value) { }

  //------------------------------------------------------------------
  int    numClasses;
  double beta;
  double gamma;
  int    deltaIter;

  private: //---------------------------------------------------------
  int    currentIter;
  int    studentsPerClass;
  int    totalIters;

  // Структуры для классов
  S_AO_Agent teachers        [];
  double     classRanks      [];
  double     classTotalCosts [];

  // Временные массивы
  double     CL [];             // Среднее знание учителей

  // Вспомогательные методы
  void   UpdateTeachersAndCosts  ();
  void   UpdateClassRanks        ();
  void   UpdateStudentsKnowledge ();
};
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
bool C_AO_CLA_l::Init (const double &rangeMinP  [],
                       const double &rangeMaxP  [],
                       const double &rangeStepP [],
                       const int epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //------------------------------------------------------------------
  currentIter = 0;
  totalIters = epochsP;
  studentsPerClass = popSize / numClasses;

  if (studentsPerClass < 1)
  {
    Print ("Ошибка: слишком мало студентов на класс");
    return false;
  }

  // Корректировка размера популяции
  //spopSize = studentsPerClass * numClasses;
  ArrayResize (a, popSize);
  for (int i = 0; i < popSize; i++) a [i].Init (coords);

  // Инициализация структур классов
  ArrayResize (teachers, numClasses);
  ArrayResize (classRanks, numClasses);
  ArrayResize (classTotalCosts, numClasses);

  for (int i = 0; i < numClasses; i++)
  {
    teachers        [i].Init (coords);
    classRanks      [i] = 1.0;
    classTotalCosts [i] = -DBL_MAX;
  }

  // Временные массивы
  ArrayResize (CL, coords);

  return true;
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
void C_AO_CLA_l::Moving ()
{
  // Начальная инициализация популяции
  if (!revision)
  {
    for (int i = 0; i < popSize; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        double val = u.RNDfromCI (rangeMin [c], rangeMax [c]);
        a [i].c [c] = u.SeInDiSp (val, rangeMin [c], rangeMax [c], rangeStep [c]);
      }
    }
    revision = true;
    return;
  }

  currentIter++;

  // Обновление знаний студентов
  UpdateStudentsKnowledge ();
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
void C_AO_CLA_l::UpdateStudentsKnowledge ()
{
  // Расчет факторов обучения для текущей итерации
  double progress = (double)currentIter / (double)MathMax (totalIters, 100);

  // Вычисление среднего знания всех учителей (CL)
  ArrayInitialize (CL, 0.0);
  int validTeachers = 0;

  for (int k = 0; k < numClasses; k++)
  {
    // Проверяем, что учитель инициализирован
    if (teachers [k].f > -DBL_MAX)
    {
      for (int c = 0; c < coords; c++)
      {
        CL [c] += teachers [k].c [c];
      }
      validTeachers++;
    }
  }

  if (validTeachers > 0)
  {
    for (int c = 0; c < coords; c++)
    {
      CL [c] /= validTeachers;
    }
  }

  // Обновление каждого студента
  for (int i = 0; i < popSize; i++)
  {
    int classIdx = i / studentsPerClass;

    // Teaching Factor - уменьшается с прогрессом
    double TF = MathExp (-0.6 * progress * classRanks [classIdx]);
    TF = MathMax (0.1, MathMin (1.0, TF));

    // Confirmatory Factor
    double CF = MathExp (-0.5 * progress * classRanks [classIdx]);
    CF = MathMax (0.1, MathMin (1.0, CF));

    // Обновление позиции студента
    for (int c = 0; c < coords; c++)
    {
      double newPos = a [i].c [c];

      // a) Teacher Learning - учимся у учителя класса
      if (teachers [classIdx].f > -DBL_MAX)
      {
        newPos += TF * (teachers [classIdx].c [c] - a [i].c [c]);
      }

      // b) Personal Learning - учимся у своего лучшего решения
      if (currentIter > deltaIter && a [i].fB > -DBL_MAX)
      {
        double PF = u.RNDprobab ();
        newPos += PF * (a [i].cB [c] - a [i].c [c]);
      }

      // c) Confirmatory Learning - учимся у среднего всех учителей
      if (currentIter > deltaIter && validTeachers > 0)
      {
        double rnd = u.RNDprobab ();
        if (rnd >= gamma) // Участвует с вероятностью (1-gamma)
        {
          newPos += CF * (CL [c] - a [i].c [c]);
        }
      }

      // Применение границ
      a [i].c [c] = u.SeInDiSp (newPos, rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
void C_AO_CLA_l::Revision ()
{
  for (int i = 0; i < popSize; i++)
  {
    // Обновление персональных лучших позиций
    if (a [i].f > a [i].fB)
    {
      a [i].fB = a [i].f;
      ArrayCopy (a [i].cB, a [i].c);
    }
    // Обновление глобального лучшего
    if (a [i].f > fB)
    {
      fB = a [i].f;
      ArrayCopy (cB, a [i].c);
    }
  }

  // Обновление учителей и costs
  UpdateTeachersAndCosts ();

  // Обновление рангов классов
  UpdateClassRanks ();
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
void C_AO_CLA_l::UpdateTeachersAndCosts ()
{
  for (int k = 0; k < numClasses; k++)
  {
    int startIdx = k * studentsPerClass;
    int endIdx   = startIdx + studentsPerClass;

    double bestFitness   = -DBL_MAX;
    int    bestIdx       = startIdx;
    double sumFitness    = 0.0;
    int    validStudents = 0;

    // Находим лучшего студента (учителя) в классе
    for (int i = startIdx; i < endIdx; i++)
    {
      if (a [i].f > -DBL_MAX) // Проверка валидности
      {
        sumFitness += a [i].f;
        validStudents++;

        if (a [i].f > bestFitness)
        {
          bestFitness = a [i].f;
          bestIdx = i;
        }
      }
    }

    // Обновляем учителя
    if (validStudents > 0)
    {
      ArrayCopy (teachers [k].c, a [bestIdx].c);
      teachers [k].f = bestFitness;

      // Расчет total cost класса
      double avgFitness = sumFitness / validStudents;
      classTotalCosts [k] = bestFitness + beta * avgFitness;
    }
  }
}
//————————————————————————————————————————————————————————————————————

//————————————————————————————————————————————————————————————————————
void C_AO_CLA_l::UpdateClassRanks ()
{
  // Находим min и max costs среди валидных классов
  double minCost   = DBL_MAX;
  double maxCost   = -DBL_MAX;
  int validClasses = 0;

  for (int k = 0; k < numClasses; k++)
  {
    if (classTotalCosts [k] > -DBL_MAX)
    {
      if (classTotalCosts [k] < minCost) minCost = classTotalCosts [k];
      if (classTotalCosts [k] > maxCost) maxCost = classTotalCosts [k];
      validClasses++;
    }
  }

  if (validClasses == 0 || maxCost - minCost < 1e-10)
  {
    // Все классы одинаковые или нет валидных
    for (int k = 0; k < numClasses; k++) classRanks [k] = 1.0;
  }
  else
  {
    // Ранжирование: лучшие классы (высокий cost) получают низкий ранг
    for (int k = 0; k < numClasses; k++)
    {
      if (classTotalCosts [k] > -DBL_MAX)
      {
        // Нормализация от 0 до 1
        double normalized = (classTotalCosts [k] - minCost) / (maxCost - minCost);

        // Инверсия: лучшие получают ранг близкий к 1
        classRanks [k] = 1.0 + (1.0 - normalized) * (numClasses - 1.0);

        // Ограничение
        classRanks [k] = MathMax (1.0, MathMin ((double)numClasses, classRanks [k]));
      }
      else
      {
        classRanks [k] = numClasses; // Худший ранг для неинициализированных
      }
    }
  }
}
//————————————————————————————————————————————————————————————————————