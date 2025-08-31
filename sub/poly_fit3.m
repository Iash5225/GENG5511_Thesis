function poly_fit3=poly_fit3(b,T,Tt,pt)%
%     B1 = 0.3072;
%     B2 = 0.0007;
%     D1 = 0.01714;
%     D2 = -0.00025;



poly_fit3 = pt .* (1 + b(1).*(T./Tt - 1)^b(2)+ b(3).*(T./Tt - 1)^b(4));
end 