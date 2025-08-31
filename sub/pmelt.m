function p_melt = pmelt(T)
    % Constants
    pt = 0.068891;  % MPa
    Tt = 83.806;    % K
    b1 = 1506.5415;
    b2 = 1.73136;
    b3 = 4677.1597;
    b4 = 0.9849295;
    
    % Calculate the melting pressure
    p_melt = pt .* (1 + b1 .* (T ./ Tt - 1).^b2 + b3 .* (T ./ Tt - 1).^b4);
end