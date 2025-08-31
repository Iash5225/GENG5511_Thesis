function poly_fit2=poly_fit2(b,T,Tt,pt)%
%     B1 = 0.3072;
%     B2 = 0.0007;
%     D1 = 0.01714;
%     D2 = -0.00025;



poly_fit2 = pt .* exp(Tt ./ T .* (b(1) .* (1 - T ./ Tt) .^ 1 + b(2) .* (1 - T ./ Tt) .^ 1.5 + b(3) .* (1 - T ./ Tt) .^ 5));
end