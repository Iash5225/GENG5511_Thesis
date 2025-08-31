function Debye3=Debye3(Xvalue) 
%
% This function calculates the Debye function of order 3, defined as
%
% Debye3(x) = 3 * [Integral {0 to x} t^3/(exp(t)-1) dt] / (x ^ 3)
%
%      The code uses Chebyshev series whose coefficients are given to 20 decimal places.
%
%   MACHINE-DEPENDENT PARAMETERS:
%
%      NTERMS - INTEGER - The no. of elements of the array ADEB3.
%                         The recommended value is such that
%                             ABS(ADEB3(NTERMS)) < EPS/100,
%                         subject to 1 <= NTERMS <= 18
%
%      XLOW - DOUBLE PRECISION - The value below which
%                    DEBYE3 = 1 - 3x/8 + x*x/20 to machine precision.
%                    The recommended value is SQRT(8 * EPSNEG)
%
%      XUPPER - DOUBLE PRECISION - The value above which
%               DEBYE3 = (18*zeta(4)/x^3) - 3*exp(-x)(x^3+3x^2+6x+6)/x^3.
%                      The recommended value is -Log(2 * EPS)
%
%      XLIM1 - DOUBLE PRECISION - The value above which DEBYE3 = 18*zeta(4)/x^3
%                     The recommended value is -Log(XMIN)
%
%      XLIM2 - DOUBLE PRECISION - The value above which DEBYE3 = 0.0 to machine
%                     precision. The recommended value is CUBE ROOT(19/XMIN)
%

ZERO = 0;
PT375 = 0.375;
HALF = 0.5;
ONE = 1;
THREE = 3;
FOUR = 4;
SIX = 6;
SEVP5 = 7.5;
EIGHT = 8;
TWENTY = 20;
ONEHUN = 100;
DEBINF = 5.13299112734217E-02;
ADEB3(1) = 2.70773706832744;
ADEB3(2) = 0.340068135211092;
ADEB3(3) = -1.29451501844409E-02;
ADEB3(4) = 7.96375538017382E-04;
ADEB3(5) = -5.46360009590824E-05;
ADEB3(6) = 3.92430195988049E-06;
ADEB3(7) = -2.8940328235386E-07;
ADEB3(8) = 2.173176139625E-08;
ADEB3(9) = -1.65420999498E-09;
ADEB3(10) = 1.2727961892E-10;
ADEB3(11) = -9.87963459E-12;
ADEB3(12) = 7.725074E-13;
ADEB3(13) = -6.077972E-14;
ADEB3(14) = 4.80759E-15;
ADEB3(15) = -3.8204E-16;
ADEB3(16) = 3.048E-17;
ADEB3(17) = -2.44E-18;
ADEB3(18) = 2E-19;
ADEB3(19) = -2E-20;
%
      X = Xvalue;
%
%   Error test
%
      if (X < ZERO)
         Debye3 = ZERO;
         return
      end 
%
%   Compute the machine-dependent constants.
%   Machine dependent constants in VBA:
%     Smallest allowed negative number -2.2251E-308
%     Smallest allowed positive number 2.2251E-308
%     Largest allowed positive number 9.99999999999999E+307
%     Largest allowed negative number -9.99999999999999E+307
%     Largest allowed positive number via formula 1.7976931348623158e+308
%     Largest allowed negative number via formula -1.7976931348623158e+308
%
% Set T = smallest postive number (machine dependent)
      T = 2.2251E-308;
      XLIM1 = -log(T);
      XK = ONE/THREE;
      XKI = (ONE/DEBINF)^XK;
      RK = T^XK;
      XLIM2 = XKI/RK;
% Set T = smallest relative spacing (machine depedendent)
      T = 1.11E-16;
      XLOW = sqrt(T*EIGHT);
      XUPPER = -log(T + T);
      T = T/ONEHUN;
      for NTERMS = 19:-1:1 
         if (abs(ADEB3(NTERMS)) > T) 
             if (X <= FOUR) 
                 if (X < XLOW) 
                     Debye3 = ((X - SEVP5)*X + TWENTY)/TWENTY;
                 else
                     T = ((X.*X./EIGHT) - HALF) - HALF;
                     Debye3 = Cheval(NTERMS, ADEB3, T) - PT375.*X;
                 end
             else%   Code for x > 4.0
                 if (X > XLIM2) 
                     Debye3 = ZERO;
                 else
                     Debye3 = ONE./(DEBINF.*X.*X.*X);
                     if (X < XLIM1) 
                         EXPMX = exp(-X);
                         if (X > XUPPER) 
                             SUM = (((X + THREE).*X + SIX).*X + SIX)./(X.*X.*X);
                         else
                             SUM = ZERO;
                             RK = floor(XLIM1/X);
                             NEXP = floor(RK);
                             XK = RK.*X;
                             for i = NEXP:-1:1
                                 XKI = ONE/XK;
                                 T = (((SIX.*XKI + SIX).*XKI + THREE).*XKI + ONE)./RK;
                                 SUM = SUM.*EXPMX + T;
                                 RK = RK - ONE;
                                 XK = XK - X;
                             end
                         end
                         Debye3 = Debye3 - THREE.*SUM.*EXPMX;
                     end
                 end
             end
             break
         end
         
      end
%
%   Code for x <= 4.0
%
  
end 