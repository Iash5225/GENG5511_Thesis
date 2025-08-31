function props = computeThermoPropsTV(T, V, parameters)
    % Function to compute thermodynamic properties based on temperature (T)
    % and pressure (p), using given model parameters.
    %
    % T: Temperature (K)
    % p: Pressure (MPa)
    % parameters: Model parameters (array)
    
    % Constants
    R = 8.31451; % Universal gas constant (J K^-1 mol^-1)
    v00 = parameters(1); % Reference molar volume (cm3/mol)
    a1 = parameters(2);
    a2 = parameters(3);
    a3 = parameters(4);
    a(1) = 3;
    % Initialize arrays and constants
    Th = parameters(10:15); % Characteristic temperatures
    g = parameters(16:21); % Gruneisen parameters
    q = parameters(22:27); % Exponent for volume dependence of Gamma
    aa = parameters(28);
    bb = parameters(29);
    cc = parameters(30);
    s = parameters(31);
    % Solve for volume (v) at given pressure and temperature
%     eps = 1e-10 * p;
%     v = 26 * exp(-p / 25000); % Initial guess for root finder
    v = V;
%     for j = 1:100
        [pcalc, dpdv, Theta, Gamma, lnz] = evaluateArgon(T, v, v00, a, a1, a2, a3, Th, g, q, aa, bb, cc);
    
    % Compute other properties
    z = T / Th(1);
    fz = z^4 / (1 + bb * z^2);
    fz1 = 2 * z^3 * (2 + bb * z^2) / (1 + bb * z^2)^2;
    fz2 = 2 * z^2 * (6 + 3 * bb * z^2 + bb^2 * z^4) / (1 + bb * z^2)^3;
    fv = exp(cc * (v - v00) / v00);

    % Compute specific heat capacities, compressibilities, and other properties
    cv = a(1) * R * (Debye3(Theta(1) / T) - (Theta(1) / T) * Debye3D(Theta(1) / T));
    dpdt = Gamma(1) * cv;
    cv = cv - (aa * R * T / Th(1)) * fz2 * fv;
    dpdt = (dpdt / v) - aa * (cc / v00) * R * fz1 * fv;
    cp = cv - T * (dpdt^2) / dpdv;
    KappaT = -1 / (v * dpdv);
    KappaS = -cv / (cp * v * dpdv);
    Alpha = -(dpdt / dpdv) / v;

    % Compute thermodynamic energies
    U = v00 * (0.5 * a1 * lnz^2 + a2 * (lnz^3) / 3 + 0.25 * a3 * (lnz^4));
    U = U + a(1) * R * T * Debye3(Theta(1) / T);
     U = U + aa * R * Th(1) * (fz - z * fz1) * fv;
    S = a(1) * R * ((4 / 3) * Debye3(Theta(1) / T) - log(1 - exp(-Theta(1) / T)));
    S = S - aa * R * fz1 * fv;
    
    % Thermodynamic potentials
    A = U - T * S;  % Helmholtz energy
    H = U + pcalc * v;  % Enthalpy
    G = H - T * S;  % Gibbs energy

    % Thermal Grüneisen parameter
    Gruneisen = Alpha * v / (cv * KappaT);

    % Store all properties in a matrix similar to `benzeneprops`
    props = zeros(1, 12);
    props(1) = pcalc;
    props(2) = KappaT;
    props(3) = KappaS;
    props(4) = Alpha;
    props(5) = cp;
    props(6) = cv;
    props(7) = Gruneisen;
    props(8) = U;
    props(9) = S;
    props(10) = A;
    props(11) = H;
    props(12) = G;
end

function [pcalc, dpdv,Theta,Gamma,lnz] = evaluateArgon(T, v, v00, a,a1, a2, a3, Th, g, q, aa, bb, cc)
    % Function to evaluate pressure (pcalc) and its derivative (dpdv) at a given T and v

    R = 8.31451;  % Universal gas constant (J K^-1 mol^-1)
    
    % Initialize variables
    Gamma = zeros(1, 1);
    Theta = zeros(1, 1);

    % Calculate Gamma and Theta for j = 1
    for j = 1:1
        Gamma(j) = g(j) * (v / v00) ^ q(j);
        if abs(q(j)) > 1e-10
            Theta(j) = Th(j) * exp((g(j) / q(j)) * (1 - (v / v00) ^ q(j)));
        else
            Theta(j) = Th(j) * (v00 / v) ^ g(j);
        end
    end

    % Calculate pressure components (pj) for j = 1
    z = v00 / v;
    lnz = log(z);
    
    pj_0 = z * (a1 * lnz + a2 * lnz^2 + a3 * lnz^3);
    pj_1 = a(1) * (R * T * Gamma(1) / v) * Debye3(Theta(1) / T);
    
    % Total pressure
    pcalc = pj_0 + pj_1;

    % Compute dp/dv for j = 1
    dpj_0 = -(z / v) * (a1 * (1 + lnz) + a2 * lnz * (2 + lnz) + a3 * lnz^2 * (3 + lnz));
    dpj_1 = (pj_1 / v) * (q(1) - 1) - a(1) * R * Theta(1) * (Gamma(1) / v)^2 * Debye3D(Theta(1) / T);
    
    % Total dp/dv
    dpdv = dpj_0 + dpj_1;

    % Apply corrections to pcalc and dpdv
    z = T / Th(1);
    exp_factor = exp(cc * (v - v00) / v00);
    correction = aa * (cc / v00) * R * Th(1) * (z^4 / (1 + bb * z^2)) * exp_factor;
    
    pcalc = pcalc - correction;
    dpdv = dpdv - aa * (cc / v00)^2 * R * Th(1) * (z^4 / (1 + bb * z^2)) * exp_factor;
end