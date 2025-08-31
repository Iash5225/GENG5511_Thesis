function p_sub = psub(T)
    % Constants
    pt = 0.068891;  % MPa
    Tt = 83.806;    % K
%     d4 = -9.231;
%     d5 = -4.954;
%     d6 = 7.043;
%     
    d4 = -10.763;
    d5 = -1.526;
    d6 = -0.4245;
    % Calculate the sublimation pressure
    p_sub = pt .* exp((Tt ./ T) .* (d4 .* (1 - T ./ Tt)^1 + d5 .* (1 - T ./ Tt).^1.5 + d6 .* (1 - T ./ Tt)^5));
end