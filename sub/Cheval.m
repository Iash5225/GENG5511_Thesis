function Cheval=Cheval(N, a , T)
%
%   This function evaluates a Chebyshev series, using the Clenshaw method with
%   Reinsch modification, as analysed in the paper by Oliver.
%
%   INPUT PARAMETERS
%
%       N - INTEGER - The no. of terms in the sequence
%
%       A - DOUBLE PRECISION ARRAY, dimension 0 to N - The coefficients of
%           the Chebyshev series
%
%       T - DOUBLE PRECISION - The value at which the series is to be evaluated
%

ZERO = 0;
HALF = 0.5;
TEST = 0.6;
TWO = 2;
U1 = ZERO;
%
%   If ABS ( T )  < 0.6 use the standard Clenshaw method
%
if (abs(T) < TEST) 
         U0 = ZERO;
         Tt = T + T;
         for i = N:-1:1
            U2 = U1;
            U1 = U0;
            U0 = Tt * U1 + a(i) - U2;
         end
         Cheval = (U0 - U2) / TWO;
else
%
%   If ABS ( T )  > =  0.6 use the Reinsch modification
%
         d1 = ZERO;
%
%   T > =  0.6 code
%
         if (T > ZERO) 
            Tt = (T - HALF) - HALF;
            Tt = Tt + Tt;
            for i = N:-1:1
               d2 = d1;
               U2 = U1;
               d1 = Tt * U2 + a(i) + d2;
               U1 = d1 + U2;
            end
            Cheval = (d1 + d2) / TWO;
         else
%
%   T < =  -0.6 code
%
            Tt = (T + HALF) + HALF;
            Tt = Tt + Tt;
            for i = N:-1:1
               d2 = d1;
               U2 = U1;
               d1 = Tt.*U2 + a(i) - d2;
               U1 = d1 - U2;
            end
            Cheval = (d1 - d2) / TWO;
         end
end
end 