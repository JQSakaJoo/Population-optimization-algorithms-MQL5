//+————————————————————————————————————————————————————————————————————————————+
//|                                                                   C_AO_RFO |
//|                                            Copyright 2007-2025, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/17063

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
// Структура для представления отдельной "руки"
struct S_RFO_Agent
{
    double card [];       // карты
    double f;             // значение функции пригодности ("ценность руки")
    int    cardRanks [];  // номера секторов ("ранги карт")

    void Init (int coords)
    {
      ArrayResize (cardRanks, coords);
      ArrayResize (card,      coords);
      f = -DBL_MAX;       // инициализация минимальным значением
    }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
class C_AO_RFO : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_RFO () { }
  C_AO_RFO ()
  {
    ao_name = "RFO";
    ao_desc = "Royal Flush Optimization (joo)";
    ao_link = "https://www.mql5.com/ru/articles/17063";

    popSize      = 50;      // размер "покерного стола" (популяции)
    deckSize     = 1000;    // количество "карт" в колоде (секторов)
    dealerBluff  = 0.03;    // вероятность "блефа" (мутации)

    ArrayResize (params, 3);

    params [0].name = "popSize";      params [0].val = popSize;
    params [1].name = "deckSize";     params [1].val = deckSize;
    params [2].name = "dealerBluff";  params [2].val = dealerBluff;
  }

  void SetParams ()
  {
    popSize     = (int)params [0].val;
    deckSize    = (int)params [1].val;
    dealerBluff = params      [2].val;
  }

  bool Init (const double &rangeMinP  [],  // минимальные значения
             const double &rangeMaxP  [],  // максимальные значения
             const double &rangeStepP [],  // шаг изменения
             const int     epochsP = 0);   // количество эпох

  void Moving   ();
  void Revision ();

  //----------------------------------------------------------------------------
  int    deckSize;         // количество секторов в измерении
  double dealerBluff;      // вероятность мутации

  S_RFO_Agent deck [];     // основная колода (популяция)
  S_RFO_Agent tempDeck []; // временная колода для сортировки
  S_RFO_Agent hand [];     // текущая рука (потомки)

  private: //-------------------------------------------------------------------
  int cutPoints;           // количество точек разрезания
  int tempCuts  [];        // временные индексы точек разрезания
  int finalCuts [];        // финальные индексы с учетом начала и конца

  void   Evolution ();     // основной процесс эволюции
  double DealCard (int rank, int suit);  // преобразование сектора в реальное значение
  void   ShuffleRanks (int &ranks []);   // мутация рангов
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_RFO::Init (const double &rangeMinP  [],
                     const double &rangeMaxP  [],
                     const double &rangeStepP [],
                     const int     epochsP = 0)
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  // Инициализация структур для хранения "рук" и "колод"
  ArrayResize (hand, popSize);
  for (int i = 0; i < popSize; i++) hand [i].Init (coords);

  ArrayResize (deck,     popSize * 2);
  ArrayResize (tempDeck, popSize * 2);
  for (int i = 0; i < popSize * 2; i++)
  {
    deck     [i].Init (coords);
    tempDeck [i].Init (coords);
  }

  // Инициализация массивов для точек разрезания
  cutPoints = 3;  // три точки разрезания для "трехкарточного" кроссовера
  ArrayResize (tempCuts,  cutPoints);
  ArrayResize (finalCuts, cutPoints + 2);

  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_RFO::Moving ()
{
  //----------------------------------------------------------------------------
  if (!revision)
  {
    // Инициализация начальной "раздачи"
    for (int i = 0; i < popSize; i++)
    {
      for (int c = 0; c < coords; c++)
      {
        hand [i].cardRanks [c] = u.RNDminusOne (deckSize);
        hand [i].card      [c] = DealCard (hand [i].cardRanks [c], c);
        hand [i].card      [c] = u.SeInDiSp (hand [i].card [c], rangeMin [c], rangeMax [c], rangeStep [c]);

        a [i].c [c] = hand [i].card [c];
      }
    }

    revision = true;
    return;
  }

  //----------------------------------------------------------------------------
  Evolution ();
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_RFO::Revision ()
{
  // Поиск лучшей "комбинации"
  int bestHand = -1;

  for (int i = 0; i < popSize; i++)
  {
    if (a [i].f > fB)
    {
      fB = a [i].f;
      bestHand = i;
    }
  }

  if (bestHand != -1) ArrayCopy (cB, a [bestHand].c, 0, 0, WHOLE_ARRAY);

  //----------------------------------------------------------------------------
  // Обновление значений пригодности
  for (int i = 0; i < popSize; i++)
  {
    hand [i].f = a [i].f;
  }

  // Добавление текущих рук в общую колоду
  for (int i = 0; i < popSize; i++)
  {
    deck [popSize + i] = hand [i];
  }

  // Сортировка колоды по ценности комбинаций
  u.Sorting (deck, tempDeck, popSize * 2);
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_RFO::Evolution ()
{
  // Для каждой позиции за столом
  for (int i = 0; i < popSize; i++)
  {
    // Выбор оппонента с учетом его рейтинга (квадрат вероятности для усиления селекции)
    double rnd = u.RNDprobab ();
    rnd *= rnd;
    int opponent = (int)u.Scale (rnd, 0.0, 1.0, 0, popSize - 1);

    // Копирование текущей руки
    hand [i] = deck [i];

    // Определение точек разрезания для обмена картами
    for (int j = 0; j < cutPoints; j++)
    {
      tempCuts [j] = u.RNDminusOne (coords);
    }

    // Сортировка точек разрезания и добавление границ
    ArraySort (tempCuts);
    ArrayCopy (finalCuts, tempCuts, 1, 0, WHOLE_ARRAY);
    finalCuts [0] = 0;
    finalCuts [cutPoints + 1] = coords - 1;

    // Случайный выбор начальной точки для обмена
    int startPoint = u.RNDbool ();

    // Обмен картами между руками
    for (int j = startPoint; j < cutPoints + 2; j += 2)
    {
      if (j < cutPoints + 1)
      {
        for (int len = finalCuts [j]; len < finalCuts [j + 1]; len++) hand [i].cardRanks [len] = deck [opponent].cardRanks [len];
      }
    }

    // Возможность "блефа" (мутации)
    ShuffleRanks (hand [i].cardRanks);

    // Преобразование рангов в реальные значения
    for (int c = 0; c < coords; c++)
    {
      hand [i].card [c] = DealCard (hand [i].cardRanks [c], c);
      hand [i].card [c] = u.SeInDiSp (hand [i].card [c], rangeMin [c], rangeMax [c], rangeStep [c]);

      a [i].c [c] = hand [i].card [c];
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
double C_AO_RFO::DealCard (int rank, int suit)
{
  // Преобразование ранга карты в реальное значение с случайным смещением внутри сектора
  double suitRange = (rangeMax [suit] - rangeMin [suit]) / deckSize;
  return rangeMin [suit] + (u.RNDprobab () + rank) * suitRange;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_RFO::ShuffleRanks (int &ranks [])
{
  // Перетасовка рангов (мутация)
  for (int i = 0; i < coords; i++)
  {
    if (u.RNDprobab () < dealerBluff) ranks [i] = (int)MathRand () % deckSize;
  }
}
//——————————————————————————————————————————————————————————————————————————————
