//+————————————————————————————————————————————————————————————————————————————+
//|                                                                    C_AO_RW |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//Article: https://www.mql5.com/ru/articles/

#include "#C_AO.mqh"

//——————————————————————————————————————————————————————————————————————————————
class C_FastRandomFloat
{
  private:
  uint state;

  // Быстрый целочисленный генератор
  uint xorshift ()
  {
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
  }

  public:
  void Init (uint seed)
  {
    state = seed;
  }

  // Получить случайное число в диапазоне [min, max]
  double nextFloat (const double min, const double max)
  {
    // Используем младшие 23 бита для мантиссы (достаточно для float)
    return min + (xorshift () & 0x7FFFFF) * (max - min) / 0x7FFFFF;
  }

  // Перегруженная версия для работы по умолчанию в диапазоне [0,1]
  double nextFloat ()
  {
    return (xorshift () & 0x7FFFFF) / (double)0x7FFFFF;
  }
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
class C_AO_RW : public C_AO
{
  public: //--------------------------------------------------------------------
  ~C_AO_RW () { }
  C_AO_RW ()
  {
    ao_name = "RW";
    ao_desc = "Random Walk";
    ao_link = "https://www.mql5.com/ru/articles/";

    popSize = 50;   //population size

    ArrayResize (params, 1);

    params [0].name = "popSize"; params [0].val = popSize;
  }

  void SetParams ()
  {
    popSize = (int)params [0].val;
  }

  bool Init (const double &rangeMinP  [], //minimum search range
             const double &rangeMaxP  [], //maximum search range
             const double &rangeStepP [], //step search
             const int     epochsP = 0);  //number of epochs

  void Moving   ();
  void Revision ();

  private: //-------------------------------------------------------------------
  C_FastRandomFloat rnd;
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
bool C_AO_RW::Init (const double &rangeMinP  [], //minimum search range
                    const double &rangeMaxP  [], //maximum search range
                    const double &rangeStepP [], //step search
                    const int     epochsP = 0)   //number of epochs
{
  if (!StandardInit (rangeMinP, rangeMaxP, rangeStepP)) return false;

  //----------------------------------------------------------------------------
  rnd.Init (GetTickCount ());
  return true;
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_RW::Moving ()
{
  for (int w = 0; w < popSize; w++)
  {
    for (int c = 0; c < coords; c++)
    {
      a [w].c [c] = u.RNDfromCI (rangeMin [c], rangeMax [c]);
      a [w].c [c] = u.SeInDiSp  (a [w].c [c], rangeMin [c], rangeMax [c], rangeStep [c]);
    }
  }
}
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
void C_AO_RW::Revision ()
{
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
}
//——————————————————————————————————————————————————————————————————————————————
