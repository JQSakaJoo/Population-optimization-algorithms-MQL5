//+————————————————————————————————————————————————————————————————————————————+
//|                                                             C_AO_Utilities |
//|                                            Copyright 2007-2024, Andrey Dik |
//|                                          https://www.mql5.com/ru/users/joo |
//—————————————————————————————————————————————————————————————————————————————+

//——————————————————————————————————————————————————————————————————————————————
class C_LCG
{
  public:

  void Init (ulong initialSeed)
  {
    seed = initialSeed;
  }

  // Быстрый линейный конгруэнтный генератор
  ulong Rand()
  {
    seed = (1664525 * seed + 1013904223) & 0xFFFFFFFF; // Модуль 2^32
    return seed;
  }

  // Генерация случайного числа в диапазоне [min, max]
  double RNDfromCI(double min, double max)
  {
    return min + ((max - min) * Rand() / 4294967296.0);
  }

  // Генерация случайного целого числа в диапазоне [min, max]
  ulong RNDintInRange(int min, int max)
  {
    return min + (Rand() % (max - min + 1));
  }

  // Генерация случайного булевого значения
  bool RNDbool()
  {
    return Rand() % 2 == 0;
  }

  // Генерация случайного числа в диапазоне [0, 1)
  double RNDprobab()
  {
    return (double)Rand() / 4294967296.0; // 2^32
  }

  // Генерация случайного числа в диапазоне [0; number - 1]
  ulong RNDminusOne(int number)
  {
    return Rand() % number;
  }

  private:
  ulong seed;  // Текущее значение seed
};
//——————————————————————————————————————————————————————————————————————————————

//——————————————————————————————————————————————————————————————————————————————
class C_AO_Utilities
{
  public: //--------------------------------------------------------------------
  double Scale                  (double In, double InMIN, double InMAX, double OutMIN, double OutMAX);
  double Scale                  (double In, double InMIN, double InMAX, double OutMIN, double OutMAX,  bool revers);
  double RNDfromCI              (double min, double max);
  int    RNDintInRange          (int min, int max);
  bool   RNDbool                ();
  double RNDprobab              ();
  int    RNDminusOne            (int number);
  double SeInDiSp               (double In, double InMin, double InMax, double Step);
  void   DecimalToGray          (ulong decimalNumber, char &array []);
  void   IntegerToBinary        (ulong number, char &array []);
  ulong  GrayToDecimal          (const char &grayCode [], int startInd, int endInd);
  ulong  BinaryToInteger        (const char &binaryStr [], const int startInd, const int endInd);
  ulong  GetMaxDecimalFromGray  (int digitsInGrayCode);
  double GaussDistribution      (const double In, const double outMin, const double outMax, const double sigma);
  double PowerDistribution      (const double In, const double outMin, const double outMax, const double p);
  double LevyFlightDistribution (double levisPower);
  double LognormalDistribution  (double center, double min_value, double max_value, double peakDisplCoeff = 0.2);
  C_LCG  lcg;

  //----------------------------------------------------------------------------
  template<typename T>
  void Sorting (T &p [], T &pTemp [], int size)
  {
    int    cnt = 1;
    int    t0  = 0;
    double t1  = 0.0;
    int    ind [];
    double val [];

    ArrayResize (ind, size);
    ArrayResize (val, size);

    for (int i = 0; i < size; i++)
    {
      ind [i] = i;
      val [i] = p [i].f;
    }

    while (cnt > 0)
    {
      cnt = 0;
      for (int i = 0; i < size - 1; i++)
      {
        if (val [i] < val [i + 1])
        {
          t0 = ind [i + 1];
          t1 = val [i + 1];
          ind [i + 1] = ind [i];
          val [i + 1] = val [i];
          ind [i] = t0;
          val [i] = t1;
          cnt++;
        }
      }
    }

    for (int u = 0; u < size; u++) pTemp [u] = p [ind [u]];
    for (int u = 0; u < size; u++) p [u] = pTemp [u];
  }

  //----------------------------------------------------------------------------
  template<typename T>
  void Sorting_fB (T &p [], T &pTemp [], int size)
  {
    int    cnt = 1;
    int    t0  = 0;
    double t1  = 0.0;
    int    ind [];
    double val [];

    ArrayResize (ind, size);
    ArrayResize (val, size);

    for (int i = 0; i < size; i++)
    {
      ind [i] = i;
      val [i] = p [i].fB;
    }

    while (cnt > 0)
    {
      cnt = 0;
      for (int i = 0; i < size - 1; i++)
      {
        if (val [i] < val [i + 1])
        {
          t0 = ind [i + 1];
          t1 = val [i + 1];
          ind [i + 1] = ind [i];
          val [i + 1] = val [i];
          ind [i] = t0;
          val [i] = t1;
          cnt++;
        }
      }
    }

    for (int u = 0; u < size; u++) pTemp [u] = p [ind [u]];
    for (int u = 0; u < size; u++) p [u] = pTemp [u];
  }

  //----------------------------------------------------------------------------
  struct S_Roulette
  {
      double start;
      double end;
  };
  S_Roulette roulette [];

  template<typename T>
  void PreCalcRoulette (T &agents [])
  {
    int aPopSize = ArraySize (agents);
    roulette [0].start = agents [0].f;
    roulette [0].end   = roulette [0].start + (agents [0].f - agents [aPopSize - 1].f);

    for (int s = 1; s < aPopSize; s++)
    {
      if (s != aPopSize - 1)
      {
        roulette [s].start = roulette [s - 1].end;
        roulette [s].end   = roulette [s].start + (agents [s].f - agents [aPopSize - 1].f);
      }
      else
      {
        roulette [s].start = roulette [s - 1].end;
        roulette [s].end   = roulette [s].start + (agents [s - 1].f - agents [s].f) * 0.1;
      }
    }
  }
  int  SpinRoulette (int aPopSize);
};
//——————————————————————————————————————————————————————————————————————————————

//------------------------------------------------------------------------------
double C_AO_Utilities :: Scale (double In, double InMIN, double InMAX, double OutMIN, double OutMAX)
{
  if (OutMIN == OutMAX) return (OutMIN);
  if (InMIN == InMAX) return (double((OutMIN + OutMAX) / 2.0));
  else
  {
    if (In < InMIN) return OutMIN;
    if (In > InMAX) return OutMAX;

    return (((In - InMIN) * (OutMAX - OutMIN) / (InMAX - InMIN)) + OutMIN);
  }
}

//------------------------------------------------------------------------------
double C_AO_Utilities :: Scale (double In, double InMIN, double InMAX, double OutMIN, double OutMAX,  bool revers)
{
  if (OutMIN == OutMAX) return (OutMIN);
  if (InMIN == InMAX) return (double((OutMIN + OutMAX) / 2.0));
  else
  {
    if (In < InMIN) return revers ? OutMAX : OutMIN;
    if (In > InMAX) return revers ? OutMIN : OutMAX;

    double res = (((In - InMIN) * (OutMAX - OutMIN) / (InMAX - InMIN)) + OutMIN);

    if (!revers) return res;
    else         return (OutMAX + OutMIN) - res;
  }
}

//------------------------------------------------------------------------------
double C_AO_Utilities ::RNDfromCI (double min, double max)
{
  if (min == max) return min;
  if (min > max)
  {
    double temp = min;
    min = max;
    max = temp;
  }
  return min + ((max - min) * rand () / 32767.0);
}

//------------------------------------------------------------------------------
int C_AO_Utilities :: RNDintInRange (int min, int max)
{
  if (min == max) return min;
  if (min > max)
  {
    int temp = min;
    min = max;
    max = temp;
  }
  return min + rand () % (max - min + 1);
}

//------------------------------------------------------------------------------
bool C_AO_Utilities :: RNDbool ()
{
  return rand () % 2 == 0;
}

//------------------------------------------------------------------------------
double C_AO_Utilities :: RNDprobab ()
{
  return (double)rand () / 32767;
}

//------------------------------------------------------------------------------
//generating a random number in the range [0; number - 1]
int C_AO_Utilities :: RNDminusOne (int number)
{
  return MathRand () % number;
}

//------------------------------------------------------------------------------
// Choice in discrete space
double C_AO_Utilities :: SeInDiSp (double In, double InMin, double InMax, double Step)
{
  if (In <= InMin) return (InMin);
  if (In >= InMax) return (InMax);
  if (Step == 0.0) return (In);
  else return (InMin + Step * (double)MathRound ((In - InMin) / Step));
}

//------------------------------------------------------------------------------
//Converting a decimal number to a Gray code
void C_AO_Utilities ::DecimalToGray (ulong decimalNumber, char &array [])
{
  ulong grayCode = decimalNumber ^ (decimalNumber >> 1);
  IntegerToBinary (grayCode, array);
}

//Converting a decimal number to a binary number
void C_AO_Utilities ::IntegerToBinary (ulong number, char &array [])
{
  ArrayResize (array, 0);
  ulong temp;
  int cnt = 0;
  while (number > 0)
  {
    ArrayResize (array, cnt + 1);
    temp = number % 2;
    array [cnt] = (char)temp;
    number = number / 2;
    cnt++;
  }

  ArrayReverse (array, 0, WHOLE_ARRAY);
}

//Converting from Gray's code to a decimal number
ulong C_AO_Utilities ::GrayToDecimal (const char &grayCode [], int startInd, int endInd)
{
  ulong grayCodeS = BinaryToInteger (grayCode, startInd, endInd);
  ulong result = grayCodeS;

  while ((grayCodeS >>= 1) > 0)
  {
    result ^= grayCodeS;
  }
  return result;
}

//Converting a binary string to a decimal number
ulong C_AO_Utilities ::BinaryToInteger (const char &binaryStr [], const int startInd, const int endInd)
{
  ulong result = 0;
  if (startInd == endInd) return 0;

  for (int i = startInd; i <= endInd; i++)
  {
    result = (result << 1) + binaryStr [i];
  }
  return result;
}

//Calculation of the maximum possible ulong number using the Gray code for a given number of bits
ulong C_AO_Utilities ::GetMaxDecimalFromGray (int digitsInGrayCode)
{
  ulong maxValue = 1;

  for (int i = 1; i < digitsInGrayCode; i++)
  {
    maxValue <<= 1;
    maxValue |= 1;
  }

  return maxValue;
}

//------------------------------------------------------------------------------
double C_AO_Utilities :: GaussDistribution (const double In, const double outMin, const double outMax, const double sigma)
{
  double logN = 0.0;
  double u1   = RNDfromCI (0.0, 1.0);
  double u2   = RNDfromCI (0.0, 1.0);

  logN = u1 <= 0.0 ? 0.000000000000001 : u1;

  double z0 = sqrt (-2 * log (logN)) * cos (2 * M_PI * u2);

  double sigmaN = sigma > 8.583864105157389 ? 8.583864105157389 : sigma;

  if (z0 >=  sigmaN) z0 = RNDfromCI (0.0,     sigmaN);
  if (z0 <= -sigmaN) z0 = RNDfromCI (-sigmaN, 0.0);

  if (z0 >= 0.0) z0 =  Scale (z0,        0.0, sigmaN, 0.0, outMax - In, false);
  else           z0 = -Scale (fabs (z0), 0.0, sigmaN, 0.0, In - outMin, false);

  return In + z0;
}

//------------------------------------------------------------------------------
double C_AO_Utilities :: PowerDistribution (const double In, const double outMin, const double outMax, const double p)
{
  double rnd = RNDfromCI (-1.0, 1.0);
  double r   = pow (fabs (rnd), p);

  if (rnd >= 0.0) return In + Scale (r, 0.0, 1.0, 0.0, outMax - In, false);
  else            return In - Scale (r, 0.0, 1.0, 0.0, In - outMin, false);
}

//------------------------------------------------------------------------------
//A distribution function close to the Levy Flight distribution.
//The function generates numbers in the range [0.0;1.0], with the distribution shifted to 0.0.
double C_AO_Utilities :: LevyFlightDistribution (double levisPower)
{
  double min = pow (20, -levisPower); //calculate the minimum possible value
  double r = RNDfromCI (1.0, 20);     //generating a number in the range [1; 20]

  r = pow (r, -levisPower);           //we raise the number r to a power
  r = (r - min) / (1 - min);          //we scale the resulting number to [0; 1]

  return r;
}

//------------------------------------------------------------------------------
//The lognormal distribution of the species:  min|------P---C---P------|max
double C_AO_Utilities :: LognormalDistribution (double center, double min_value, double max_value, double peakDisplCoeff = 0.2)
{
  // Проверка правой границы
  if (min_value >= max_value)
  {
    return max_value;
  }

  // Проверка левой границы
  if (max_value <= min_value)
  {
    return min_value;
  }

  // Генерация случайного числа от 0 до 1
  double random = MathRand () / 32767.0;

  // Коррекция центра, если он выходит за границы
  if (center < min_value)
  {
    center = min_value;
    random = 1;
  }

  if (center > max_value)
  {
    center = max_value;
    random = 0;
  }

  // Расчет положения пиков
  double peak_left  = center - (center - min_value) * peakDisplCoeff;
  double peak_right = center + (max_value - center) * peakDisplCoeff;

  double result = 0.0;

  if (random < 0.5) // Левая часть распределения
  {
    // Расчет параметров для левой части
    double diff_center_peak = MathMax (center - peak_left, DBL_EPSILON);
    double diff_center_min  = MathMax (center - min_value, DBL_EPSILON);

    double mu_left = MathLog (diff_center_peak);
    double sigma_left = MathSqrt (2.0 * MathLog (MathMax (diff_center_min / diff_center_peak, DBL_EPSILON)) / 9.0);

    // Генерация случайных чисел для метода Бокса-Мюллера
    double u1 = MathRand () / 32767.0;
    double u2 = MathRand () / 32767.0;

    // Защита от нулевых значений
    u1 = MathMax (u1, DBL_EPSILON);

    // Применение метода Бокса-Мюллера
    double z = MathSqrt (-2.0 * MathLog (u1)) * MathCos (2.0 * M_PI * u2);

    // Расчет результата для левой части
    result = center - MathExp (mu_left + sigma_left * z);
  }
  else // Правая часть распределения
  {
    // Расчет параметров для правой части
    double diff_peak_center = MathMax (peak_right - center, DBL_EPSILON);
    double diff_max_center  = MathMax (max_value - center,  DBL_EPSILON);

    double mu_right    = MathLog  (diff_peak_center);
    double sigma_right = MathSqrt (2.0 * MathLog (MathMax (diff_max_center / diff_peak_center, DBL_EPSILON)) / 9.0);

    // Генерация случайных чисел для метода Бокса-Мюллера
    double u1 = MathRand () / 32767.0;
    double u2 = MathRand () / 32767.0;

    // Защита от нулевых значений
    u1 = MathMax (u1, DBL_EPSILON);

    // Применение метода Бокса-Мюллера
    double z = MathSqrt (-2.0 * MathLog (u1)) * MathCos (2.0 * M_PI * u2);

    // Расчет результата для правой части
    result = center + MathExp (mu_right + sigma_right * z);
  }

  // Проверка и коррекция результата, если он выходит за границы
  if (result < min_value || result > max_value) return RNDfromCI (min_value, max_value);

  return result;
}