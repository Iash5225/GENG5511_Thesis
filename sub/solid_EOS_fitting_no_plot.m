 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;path(path,[pwd,'\..\..\SUB']);
%fontsize = 20; markersize = 11; linewidth = 0.8;
fontsize = 14; markersize = 7; linewidth = 0.9;
mymarker = {'s','d','x','o','^','>','<','+','*','o','x','s','d','^','v','>','<','+','*','o','x','s','d','^','v','>','<','+','*','o','x','s','d','^','v','>','<','+','*','o','x','s','d','^','v','>','<','+','*','o','x','s','d','^','v','>','<'};
mycolor = [  0 0 0; 112 48 160;192 0 0; 1 175 146;222 110 38;   0 0 255; 150 150 150;      95 58 91;
                72 113 57;  27 71 116;  222 110 38;  139 44 42;
                0 200 0;   255 0 240; 
                92 103 177;  71 30 118;  100 200 0;  239 144 42;
                120 100 255;   55 200 80;200  20  150;   25  105 88; 88 10 198;
                100 55 241; 254 120 62; 165 158 171; 224 21 138; 155 100 8; 84 184 93;
                193 233 41; 250 199 64;200 175 41;127 217 16;0  0  0;   255  0 0;   0 0 255;     95 58 91;
                72 113 57;  27 71 116;  222 110 38;  139 44 42;
                0 200 0;   255 0 240; 
                92 103 177;  71 30 118;  100 200 0;  239 144 42;
                120 100 255;   55 200 80;200  20  150;   25  105 88; 88 10 198;
                100 55 241; 254 120 62; 165 158 171; 224 21 138; 155 100 8; 84 184 93;
                193 233 41; 250 199 64;200 175 41;127 217 16;]/256;%note colour needs to be adjusted manually
%preparation: read txt files
Lplot = 0;% 0 not plot fitted results, 1 plot fitted results
[Year_Vm_sub,Author_Vm_sub,T_Vm_sub,Vm_sub] = textread('../evaluation data/cell volume/Vm_sublimation.txt','%s%s%f%f','headerlines',2);
for ii = 1:length(Vm_sub)
    p_Vm_sub(ii) = psub(T_Vm_sub(ii));
end
[Year_Vm_melt,Author_Vm_melt,T_Vm_melt,Vm_melt] = textread('../evaluation data/cell volume/Vm_melting.txt','%s%s%f%f','headerlines',2);
for ii = 1:length(Vm_melt)
    p_Vm_melt(ii) = pmelt(T_Vm_melt(ii));
end
[Year_Vm_highp,Author_Vm_highp,T_Vm_highp,Vm_highp,p_Vm_highp] = textread('../evaluation data/cell volume/Vm_high_pressure.txt','%s%s%f%f%f','headerlines',2);
[Year_cp_sub,Author_cp_sub,T_cp_sub,cp_sub] = textread('../evaluation data/heat capacity/cp_sublimation.txt','%s%s%f%f','headerlines',2);
for ii = 1:length(cp_sub)
    p_cp_sub(ii) = psub(T_cp_sub(ii));
end
[Year_alpha_sub,Author_alpha_sub,T_alpha_sub,alpha_sub] = textread('../evaluation data/thermal expansion/alpha_sublimation.txt','%s%s%f%f','headerlines',2);
for ii = 1:length(alpha_sub)
    p_alpha_sub(ii) = psub(T_alpha_sub(ii));
end
[Year_BetaT_sub,Author_BetaT_sub,T_BetaT_sub,BetaT_sub] = textread('../evaluation data/bulk modulus/BetaT_sublimation.txt','%s%s%f%f','headerlines',2);
for ii = 1:length(BetaT_sub)
    p_BetaT_sub(ii) = psub(T_BetaT_sub(ii));
end
[Year_BetaS_sub,Author_BetaS_sub,T_BetaS_sub,BetaS_sub] = textread('../evaluation data/bulk modulus/BetaS_sublimation.txt','%s%s%f%f','headerlines',2);
for ii = 1:length(BetaS_sub)
    p_BetaS_sub(ii) = psub(T_BetaS_sub(ii));
end
[Year_sub,Author_sub,T_sub,p_sub,G_fluid_sub,V_fluid_sub] = textread('../evaluation data/sublimation/sublimation_for_fitting.txt','%s%s%f%f%f%f','headerlines',2);
[Year_melt,Author_melt,T_melt,p_melt,G_fluid_melt,V_fluid_melt] = textread('../evaluation data/melting/melting_for_fitting.txt','%s%s%f%f%f%f','headerlines',2);
[Year_H_sub,Author_H_sub,T_H_sub,delta_H_sub,H_fluid_sub] = textread('../evaluation data/enthalpy/enthalpy of sublimation for fitting.txt','%s%s%f%f%f','headerlines',2);
for ii = 1:length(T_H_sub)
    p_H_sub(ii) = psub(T_H_sub(ii));
end

[Year_H_melt,Author_H_melt,T_H_melt,delta_H_melt,H_fluid_melt] = textread('../evaluation data/enthalpy/enthalpy of fusion for fitting.txt','%s%s%f%f%f','headerlines',2);
for ii = 1:length(T_H_melt)
    p_H_melt(ii) = pmelt(T_H_melt(ii));
end

params_init = [22.501,2614.9,6518,10,0,0,0,0,0,88.8,0,0,0,0,0,2.659,0,0,0,0,0,0.2987,0,0,0,0,0,0.0501,0.4482,3.2797,129.2201];
lb = [22.501,0,0,0,0,0,0,0,0,20,0,0,0,0,0,-5,0,0,0,0,0,-10,0,0,0,0,0,0,0,0,0];
ub = [22.501,10000,10000,10000,0,0,0,0,0,300,0,0,0,0,0,5,0,0,0,0,0,10,0,0,0,0,0,100,100,100,1000];

% Optimization options
% options = optimset('Display', 'iter', 'TolFun', 1e-25, 'MaxIter', 2000000);
options = optimset('Display', 'iter', 'TolFun', 1e-25, 'MaxIter', 20000000, 'OutputFcn', @outfun);



% Perform the optimization using fmincon
[params_fit, fval] = fmincon(@(params) combined_cost_function(params, T_Vm_sub, p_Vm_sub, Vm_sub, ...
                           T_Vm_melt, p_Vm_melt, Vm_melt, ...
                           T_Vm_highp, p_Vm_highp, Vm_highp, ...
                           T_cp_sub, p_cp_sub, cp_sub, ...
                           T_alpha_sub, p_alpha_sub, alpha_sub, ...
                           T_BetaT_sub, p_BetaT_sub, BetaT_sub, ...
                           T_BetaS_sub, p_BetaS_sub, BetaS_sub, ...
                           T_sub, p_sub, G_fluid_sub, V_fluid_sub, ...
                           T_melt, p_melt, G_fluid_melt,V_fluid_melt,T_H_sub, p_H_sub, delta_H_sub, H_fluid_sub, ...
                           T_H_melt, p_H_melt, delta_H_melt, H_fluid_melt), ...
                           params_init, [], [], [], [], lb, ub, [], options);

%need to: plot for each property
%for data verses fitted curve
%check the extrapolation behaviour
% Compute fitted values using the optimized parameters


function stop = outfun(x, optimValues, state)
    stop = false;
    switch state
        case 'iter'
            % Calculate and display individual deviations at each iteration
            [~, deviations] = combined_cost_function(x, T_Vm_sub, p_Vm_sub, Vm_sub, ...
                            T_Vm_melt, p_Vm_melt, Vm_melt, ...
                            T_Vm_highp, p_Vm_highp, Vm_highp, ...
                            T_cp_sub, p_cp_sub, cp_sub, ...
                            T_alpha_sub, p_alpha_sub, alpha_sub, ...
                            T_BetaT_sub, p_BetaT_sub, BetaT_sub, ...
                            T_BetaS_sub, p_BetaS_sub, BetaS_sub, ...
                            T_sub, p_sub, G_fluid_sub, V_fluid_sub, ...
                            T_melt, p_melt, G_fluid_melt, V_fluid_melt, T_H_sub, p_H_sub, delta_H_sub, H_fluid_sub, ...
                            T_H_melt, p_H_melt, delta_H_melt, H_fluid_melt);

            disp(['Current total_deviation: ', num2str(optimValues.fval)]);
            disp(['Vm_sub deviation: ', num2str(deviations.Vm_sub)]);
            disp(['Vm_melt deviation: ', num2str(deviations.Vm_melt)]);
            % Display other deviations as needed
        otherwise
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%deviation function
function [total_deviation, deviations] = combined_cost_function(params,T_Vm_sub, p_Vm_sub,Vm_sub,T_Vm_melt, p_Vm_melt,Vm_melt,T_Vm_highp, p_Vm_highp,Vm_highp,T_cp_sub, p_cp_sub,cp_sub,T_alpha_sub,p_alpha_sub,alpha_sub,T_BetaT_sub, p_BetaT_sub,BetaT_sub,T_BetaS_sub, p_BetaS_sub,BetaS_sub,T_sub, p_sub,G_fluid_sub,V_fluid_sub,T_melt,p_melt,G_fluid_melt,V_fluid_melt,T_H_sub,p_H_sub,delta_H_sub,H_fluid_sub,T_H_melt,p_H_melt,delta_H_melt,H_fluid_melt)
    % Compute individual property deviations
St_REFPROP = 131.1500;%J/mol/K
Ht_REFPROP = 1689.10;%J/molProps(K1038,K1039,12,$B$7:$B$36)+K1038*O1041
Tt = 83.8058; pt = 0.068891;%unit: K and MPa
deltaS_triple =  params(31) - St_REFPROP;
props_Triple = computeThermoProps(Tt,pt, params);
Ht_fitted = props_Triple(12) + Tt*params(31);
deltaH_triple = Ht_fitted - Ht_REFPROP;
    for ii = 1:length(T_Vm_sub)
        props = computeThermoProps(T_Vm_sub(ii), p_Vm_sub(ii), params);
        Vm_sub_deviation(ii) = 100*(Vm_sub(ii) - props(1)) ./ Vm_sub(ii);
    end
    Vm_sub_dev = sqrt(sumsqr(Vm_sub_deviation)/length(Vm_sub_deviation));
    for ii = 1:length(T_Vm_melt)
        props = computeThermoProps(T_Vm_melt(ii), p_Vm_melt(ii), params);
        Vm_melt_deviation(ii) = 100*(Vm_melt(ii) - props(1)) ./ Vm_melt(ii);       
    end
    Vm_melt_dev = sqrt(sumsqr(Vm_melt_deviation)/length(Vm_melt_deviation));
    for ii = 1:length(T_Vm_highp)
        props = computeThermoProps(T_Vm_highp(ii), p_Vm_highp(ii), params);
        Vm_highp_deviation(ii) = 100*(Vm_highp(ii) - props(1)) ./ Vm_highp(ii);       
    end
    Vm_highp_dev = sqrt(sumsqr(Vm_highp_deviation)/length(Vm_highp_deviation));
    for ii = 1:length(T_cp_sub)
        props = computeThermoProps(T_cp_sub(ii), p_cp_sub(ii), params);
        cp_sub_deviation(ii) = 100*(cp_sub(ii) - props(5)) ./ cp_sub(ii);
    end    
    cp_sub_dev = sqrt(sumsqr(cp_sub_deviation)/length(cp_sub_deviation));
    for ii = 1:length(T_alpha_sub)
        props = computeThermoProps(T_alpha_sub(ii), p_alpha_sub(ii), params);
        alpha_sub_deviation(ii) = 100*(alpha_sub(ii) - props(4)) ./ alpha_sub(ii);
    end        
    alpha_sub_dev = sqrt(sumsqr(alpha_sub_deviation)/length(alpha_sub_deviation));    
    for ii = 1:length(T_BetaT_sub)
        props = computeThermoProps(T_BetaT_sub(ii), p_BetaT_sub(ii), params);
        BetaT_sub_deviation(ii) = real(100*(BetaT_sub(ii) - props(2)) ./ BetaT_sub(ii));
    end  
    BetaT_sub_dev = sqrt(sumsqr(BetaT_sub_deviation)/length(BetaT_sub_deviation));    
    for ii = 1:length(T_BetaS_sub)
        props = computeThermoProps(T_BetaS_sub(ii), p_BetaS_sub(ii), params);
        BetaS_sub_deviation(ii) = 100*(BetaS_sub(ii) - props(3)) ./ BetaS_sub(ii);
    end      
    BetaS_sub_dev = sqrt(sumsqr(BetaS_sub_deviation)/length(BetaS_sub_deviation));    
%     T_H_sub,p_H_sub,H_sub,H_fluid_sub
    for ii = 1:length(T_H_sub)
        H_solid_sub(ii) = H_fluid_sub(ii)-1000*delta_H_sub(ii);
        props = computeThermoProps(T_H_sub(ii), p_H_sub(ii), params);
        H_solid_sub_fitted(ii) = props(11) - deltaH_triple;
        H_solid_sub_deviation(ii) = 100*(H_solid_sub(ii) - H_solid_sub_fitted(ii)) / H_solid_sub_fitted(ii) ;
    end
    H_solid_sub_dev = sqrt(sumsqr(H_solid_sub_deviation)/length(H_solid_sub_deviation));    
    for ii = 1:length(T_H_melt)
        H_solid_melt(ii) = H_fluid_melt(ii)-1000*delta_H_melt(ii);
        props = computeThermoProps(T_H_melt(ii), p_H_melt(ii), params);
        H_solid_melt_fitted(ii) = props(11) - deltaH_triple;
        H_solid_melt_deviation(ii) = 100*(H_solid_melt(ii) - H_solid_melt_fitted(ii)) / H_solid_melt_fitted(ii) ;
    end    
    H_solid_melt_dev = sqrt(sumsqr(H_solid_melt_deviation)/length(H_solid_melt_deviation));    
    
    for ii = 1:length(p_sub)
        props = computeThermoProps(T_sub(ii), p_sub(ii)/10^6, params);
        delta_G_sub(ii) = G_fluid_sub(ii)-props(12) + deltaH_triple - T_sub(ii) * deltaS_triple;
        p_fitted_sub(ii) = p_sub(ii) - delta_G_sub(ii)/(V_fluid_sub(ii) - props(1))*10^6;%unit Pa
        p_sub_deviation(ii) = 100 * (p_sub(ii) - p_fitted_sub(ii))/p_sub(ii);
    end  
    
    p_sub_dev = sqrt(sumsqr(p_sub_deviation)/length(p_sub_deviation));    
    
    for ii = 1:length(p_melt)
        props = computeThermoProps(T_melt(ii), p_melt(ii), params);
        delta_G_melt(ii) = G_fluid_melt(ii)-props(12) + deltaH_triple - T_melt(ii) * deltaS_triple;
        p_fitted_melt(ii) = p_melt(ii) - delta_G_melt(ii)/(V_fluid_melt(ii) - props(1));%unit MPa
        p_melt_deviation(ii) = 100 * (p_melt(ii) - p_fitted_melt(ii))/p_melt(ii);
    end     
    p_melt_dev = sqrt(sumsqr(p_melt_deviation)/length(p_melt_deviation));  
    
    % Store individual deviations in a struct for easy access
    deviations.Vm_sub = Vm_sub_dev;
    deviations.Vm_melt = Vm_melt_dev;
    deviations.Vm_highp = Vm_highp_dev;
    deviations.cp_sub = cp_sub_dev;
    deviations.alpha_sub = alpha_sub_dev;
    deviations.BetaT_sub = BetaT_sub_dev;
    deviations.BetaS_sub = BetaS_sub_dev;
    deviations.H_solid_sub = H_solid_sub_dev;
    deviations.H_solid_melt = H_solid_melt_dev;
    deviations.p_sub = p_sub_dev;
    deviations.p_melt = p_melt_dev;
    
    total_deviation = Vm_sub_dev + Vm_melt_dev + Vm_highp_dev + cp_sub_dev + alpha_sub_dev + BetaT_sub_dev + BetaS_sub_dev + H_solid_sub_dev + H_solid_melt_dev + p_sub_dev + p_melt_dev;
 
end