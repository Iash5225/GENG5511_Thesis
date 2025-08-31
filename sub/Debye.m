function Debye=Debye(X)
%*************************************************************************
% Function to compute Cv/R for one Debye distribution
%*************************************************************************

lim = 9.99999999999999E+307;
lim = log(sqrt(lim)) - 1;
if (X < lim)
  Debye = 4 * Debye3(X) - 3 * X / (exp(X) - 1);
else
  Debye = 4* Debye3(X) + 3 * X;
end 
end