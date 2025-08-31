function Einstein=Einstein(X)
%*************************************************************************
%Function to compute Cv/R for one harmonic mode
%*************************************************************************

lim = 9.99999999999999E+307;
lim = log(sqrt(lim)) - 1;
if (X < lim) 
  Einstein = X.^2.*exp(X)./(exp(X) - 1).^2;
else
  Einstein = 0;
end 
end 
