%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;path(path,[pwd,'\..\..\SUB']);
%fontsize = 20; markersize = 11; linewidth = 0.8;
fontsize = 13.5; markersize = 7; linewidth = 0.9;
mymarker = {'s','d','x','o','^','>','<','+','*','v','pentagram','h','s','d','^','v','>','<','+','*','o','x','s','d','^','v','>','<','+','*','o','x','s','d','^','v','>','<','+','*','o','x','s','d','^','v','>','<','+','*','o','x','s','d','^','v','>','<'};
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
[Year_Vm_sub,Author_Vm_sub,T_Vm_sub,Vm_sub] = textread('../evaluation data/cell volume/Vm_sublimation.txt','%s%s%f%f','headerlines',2);
for ii = 1:length(Vm_sub)
    p_Vm_sub(ii) = psub(T_Vm_sub(ii));
end
[Year_Vm_sub_all,Author_Vm_sub_all,T_Vm_sub_all,Vm_sub_all] = textread('../evaluation data/cell volume/Vm_sublimation_all_data.txt','%s%s%f%f','headerlines',2);
for ii = 1:length(Vm_sub_all)
    p_Vm_sub_all(ii) = psub(T_Vm_sub_all(ii));
end
[Year_Vm_melt,Author_Vm_melt,T_Vm_melt,Vm_melt] = textread('../evaluation data/cell volume/Vm_melting.txt','%s%s%f%f','headerlines',2);
for ii = 1:length(Vm_melt)
    p_Vm_melt(ii) = pmelt(T_Vm_melt(ii));
end
[Year_Vm_melt_all,Author_Vm_melt_all,T_Vm_melt_all,Vm_melt_all] = textread('../evaluation data/cell volume/Vm_melting_all_data.txt','%s%s%f%f','headerlines',2);
for ii = 1:length(Vm_melt_all)
    p_Vm_melt_all(ii) = pmelt(T_Vm_melt_all(ii));
end
[Year_Vm_highp,Author_Vm_highp,T_Vm_highp,Vm_highp,p_Vm_highp] = textread('../evaluation data/cell volume/Vm_high_pressure.txt','%s%s%f%f%f','headerlines',2);
[Year_Vm_highp_all,Author_Vm_highp_all,T_Vm_highp_all,Vm_highp_all,p_Vm_highp_all] = textread('../evaluation data/cell volume/Vm_high_pressure_all_data.txt','%s%s%f%f%f','headerlines',2);
[Year_cp_sub,Author_cp_sub,T_cp_sub,cp_sub] = textread('../evaluation data/heat capacity/cp_sublimation.txt','%s%s%f%f','headerlines',2);
for ii = 1:length(cp_sub)
    p_cp_sub(ii) = psub(T_cp_sub(ii));
end
[Year_cp_sub_all,Author_cp_sub_all,T_cp_sub_all,cp_sub_all] = textread('../evaluation data/heat capacity/cp_sublimation_all_data.txt','%s%s%f%f','headerlines',2);
for ii = 1:length(cp_sub_all)
    p_cp_sub_all(ii) = psub(T_cp_sub_all(ii));
end
[Year_alpha_sub,Author_alpha_sub,T_alpha_sub,alpha_sub] = textread('../evaluation data/thermal expansion/alpha_sublimation.txt','%s%s%f%f','headerlines',2);
for ii = 1:length(alpha_sub)
    p_alpha_sub(ii) = psub(T_alpha_sub(ii));
end
[Year_alpha_sub_all,Author_alpha_sub_all,T_alpha_sub_all,alpha_sub_all] = textread('../evaluation data/thermal expansion/alpha_sublimation_all_data.txt','%s%s%f%f','headerlines',2);
for ii = 1:length(alpha_sub_all)
    p_alpha_sub_all(ii) = psub(T_alpha_sub_all(ii));
end
[Year_BetaT_sub,Author_BetaT_sub,T_BetaT_sub,BetaT_sub] = textread('../evaluation data/bulk modulus/BetaT_sublimation.txt','%s%s%f%f','headerlines',2);
BetaT_sub = 1./BetaT_sub;
for ii = 1:length(BetaT_sub)
    p_BetaT_sub(ii) = psub(T_BetaT_sub(ii));
end
[Year_BetaT_sub_all,Author_BetaT_sub_all,T_BetaT_sub_all,BetaT_sub_all] = textread('../evaluation data/bulk modulus/BetaT_sublimation_all_data.txt','%s%s%f%f','headerlines',2);
BetaT_sub_all = 1./BetaT_sub_all;
for ii = 1:length(BetaT_sub_all)
    p_BetaT_sub_all(ii) = psub(T_BetaT_sub_all(ii));
end
[Year_BetaS_sub,Author_BetaS_sub,T_BetaS_sub,BetaS_sub] = textread('../evaluation data/bulk modulus/BetaS_sublimation.txt','%s%s%f%f','headerlines',2);
BetaS_sub = 1./BetaS_sub;
for ii = 1:length(BetaS_sub)
    p_BetaS_sub(ii) = psub(T_BetaS_sub(ii));
end
[Year_BetaS_sub_all,Author_BetaS_sub_all,T_BetaS_sub_all,BetaS_sub_all] = textread('../evaluation data/bulk modulus/BetaS_sublimation_all_data.txt','%s%s%f%f','headerlines',2);
BetaS_sub_all = 1./BetaS_sub_all;
for ii = 1:length(BetaS_sub_all)
    p_BetaS_sub_all(ii) = psub(T_BetaS_sub_all(ii));
end
[Year_sub,Author_sub,T_sub,p_sub,G_fluid_sub,V_fluid_sub] = textread('../evaluation data/sublimation/sublimation_for_fitting.txt','%s%s%f%f%f%f','headerlines',2);
[Year_sub_all,Author_sub_all,T_sub_all,p_sub_all,G_fluid_sub_all,V_fluid_sub_all] = textread('../evaluation data/sublimation/sublimation_all_data.txt','%s%s%f%f%f%f','headerlines',2);
[Year_KT_highp_all,Author_KT_highp_all,T_KT_highp_all,KT_highp_all,p_KT_highp_all] = textread('../evaluation data/bulk modulus/BetaT_high_pressure_all_data.txt','%s%s%f%f%f','headerlines',2);
KT_highp_all = 1./KT_highp_all;
[Year_melt,Author_melt,T_melt,p_melt,G_fluid_melt,V_fluid_melt] = textread('../evaluation data/melting/melting_for_fitting.txt','%s%s%f%f%f%f','headerlines',2);
[Year_melt_all,Author_melt_all,T_melt_all,p_melt_all,G_fluid_melt_all,V_fluid_melt_all] = textread('../evaluation data/melting/melting_all_data.txt','%s%s%f%f%f%f','headerlines',2);

[Year_H_sub,Author_H_sub,T_H_sub,delta_H_sub,H_fluid_sub] = textread('../evaluation data/enthalpy/enthalpy of sublimation for fitting.txt','%s%s%f%f%f','headerlines',2);
for ii = 1:length(T_H_sub)
    p_H_sub(ii) = psub(T_H_sub(ii));
end

% 
[Year_H_melt,Author_H_melt,T_H_melt,delta_H_melt,H_fluid_melt] = textread('../evaluation data/enthalpy/enthalpy of fusion for fitting.txt','%s%s%f%f%f','headerlines',2);
for ii = 1:length(T_H_melt)
    p_H_melt(ii) = pmelt(T_H_melt(ii));
end
[Year_H_melt_all,Author_H_melt_all,T_H_melt_all,delta_H_melt_all,H_fluid_melt_all] = textread('../evaluation data/enthalpy/enthalpy of fusion all data.txt','%s%s%f%f%f','headerlines',2);
for ii = 1:length(T_H_melt_all)
    p_H_melt_all(ii) = pmelt(T_H_melt_all(ii));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choice = 1;%Vm along sublimation
% choice = 2;%Vm along melting
% choice = 3;%cp along sublimation
% choice = 4;%alpha along sublimation
% choice = 5;%KT
% choice = 6;%KS
% choice = 7;%enthalpy sublimation
% choice = 8;%enthalpy melting
choice = 9;%sublimation
% choice = 10;%melting
% choice = 11;%high pressure volume
% choice = 12;%KT high pressure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 1

ndata_Vm_sub=length(T_Vm_sub);xxcount_Vm_sub=1; Nxstart_Vm_sub = 1;yy_Vm_sub=Year_Vm_sub(1);AA_Vm_sub = Author_Vm_sub(1);
if ndata_Vm_sub == 1
    Nxend_Vm_sub(xxcount_Vm_sub) = 1;
else
    for ii = 2:ndata_Vm_sub
        if ~strcmp(yy_Vm_sub,string(Year_Vm_sub(ii)))||~strcmp(AA_Vm_sub,Author_Vm_sub(ii))
            Nxend_Vm_sub(xxcount_Vm_sub) = ii-1;
            xxcount_Vm_sub = xxcount_Vm_sub + 1;
            Nxstart_Vm_sub(xxcount_Vm_sub) = ii;
            yy_Vm_sub=Year_Vm_sub(ii);AA_Vm_sub=Author_Vm_sub(ii);
        end
        if ii == ndata_Vm_sub
            Nxend_Vm_sub(xxcount_Vm_sub) = ndata_Vm_sub;
        end
    end
end

    for ii = 1:length(T_Vm_sub)
        props_Vm_sub = computeThermoPropsForResult(T_Vm_sub(ii), p_Vm_sub(ii));
        Vm_sub_deviation(ii) = 100*(Vm_sub(ii) - props_Vm_sub(1)) ./ Vm_sub(ii);
    end


kk=1;
for ii=1:xxcount_Vm_sub
    [Tmax_Vm_sub(kk,:) Tmin_Vm_sub(kk,:)]=maxmin(T_Vm_sub(Nxstart_Vm_sub(ii):Nxend_Vm_sub(ii)));
    if Tmin_Vm_sub(kk,:)<83.4
        Tmin_Vm_sub(kk,:) = round (Tmin_Vm_sub(kk,:),0);
    end
    if Tmax_Vm_sub(kk,:)<83.4
        Tmax_Vm_sub(kk,:) = round (Tmax_Vm_sub(kk,:),0);
    end
        T_range_Vm_sub = strcat(string(Tmin_Vm_sub),"-",string(Tmax_Vm_sub));
     N_Vm_sub(kk,:)=Nxend_Vm_sub(ii)-Nxstart_Vm_sub(ii)+1;
     RMS_Vm_sub(kk,:)=sqrt(sumsqr(Vm_sub_deviation(Nxstart_Vm_sub(ii):Nxend_Vm_sub(ii)))/N_Vm_sub(kk,:));
     AAD_Vm_sub(kk,:)=mean(abs(Vm_sub_deviation(Nxstart_Vm_sub(ii):Nxend_Vm_sub(ii))));
     kk=kk+1;
     
end    
    if ~isfolder('../../paper writing/table/cell volume')
        mkdir('../../paper writing/table/cell volume')
    end
    fid_name = fopen('../../paper writing/table/cell volume/sublimation.txt','w');
    fprintf(fid_name,'    Year              Author                T/K            N         RMS        AAD\n');
    for ii = 1:xxcount_Vm_sub
        fprintf(fid_name,'%8s%20s%20s%12d%12.2f%12.2f\n',string(Year_Vm_sub(Nxstart_Vm_sub(ii))),string(Author_Vm_sub(Nxstart_Vm_sub(ii))),T_range_Vm_sub(ii),N_Vm_sub(ii),RMS_Vm_sub(ii),AAD_Vm_sub(ii));
    end    
    
ndata_Vm_sub_all=length(T_Vm_sub_all);xxcount_Vm_sub_all=1; Nxstart_Vm_sub_all = 1;yy_Vm_sub_all=Year_Vm_sub_all(1);AA_Vm_sub_all = Author_Vm_sub_all(1);
if ndata_Vm_sub_all == 1
    Nxend_Vm_sub_all(xxcount_Vm_sub_all) = 1;
else
    for ii = 2:ndata_Vm_sub_all
        if ~strcmp(yy_Vm_sub_all,string(Year_Vm_sub_all(ii)))||~strcmp(AA_Vm_sub_all,Author_Vm_sub_all(ii))
            Nxend_Vm_sub_all(xxcount_Vm_sub_all) = ii-1;
            xxcount_Vm_sub_all = xxcount_Vm_sub_all + 1;
            Nxstart_Vm_sub_all(xxcount_Vm_sub_all) = ii;
            yy_Vm_sub_all=Year_Vm_sub_all(ii);AA_Vm_sub_all=Author_Vm_sub_all(ii);
        end
        if ii == ndata_Vm_sub_all
            Nxend_Vm_sub_all(xxcount_Vm_sub_all) = ndata_Vm_sub_all;
        end
    end
end

    for ii = 1:length(T_Vm_sub_all)
        props_Vm_sub_all = computeThermoPropsForResult(T_Vm_sub_all(ii), p_Vm_sub_all(ii));
        Vm_sub_deviation_all(ii) = 100*(Vm_sub_all(ii) - props_Vm_sub_all(1)) ./ Vm_sub_all(ii);
    end
 kk = 1;
for ii=1:xxcount_Vm_sub_all
    [Tmax_Vm_sub_all(kk,:) Tmin_Vm_sub_all(kk,:)]=maxmin(T_Vm_sub_all(Nxstart_Vm_sub_all(ii):Nxend_Vm_sub_all(ii)));
    if Tmin_Vm_sub_all(kk,:)<83.4
        Tmin_Vm_sub_all(kk,:) = round (Tmin_Vm_sub_all(kk,:),0);
    end
    if Tmax_Vm_sub_all(kk,:)<83.4
        Tmax_Vm_sub_all(kk,:) = round (Tmax_Vm_sub_all(kk,:),0);
    end
            T_range_Vm_sub_all = strcat(string(Tmin_Vm_sub_all),"-",string(Tmax_Vm_sub_all));
     N_Vm_sub_all(kk,:)=Nxend_Vm_sub_all(ii)-Nxstart_Vm_sub_all(ii)+1;
     RMS_Vm_sub_all(kk,:)=sqrt(sumsqr(Vm_sub_deviation_all(Nxstart_Vm_sub_all(ii):Nxend_Vm_sub_all(ii)))/N_Vm_sub_all(kk,:));
     AAD_Vm_sub_all(kk,:)=mean(abs(Vm_sub_deviation_all(Nxstart_Vm_sub_all(ii):Nxend_Vm_sub_all(ii))));
     kk=kk+1;
     
end        
    
     if ~isfolder('../../paper writing/table/cell volume')
        mkdir('../../paper writing/table/cell volume')
    end
    fid_name = fopen('../../paper writing/table/cell volume/sublimation all data.txt','w');
    fprintf(fid_name,'    Year              Author                T/K            N         RMS        AAD\n');
    for ii = 1:xxcount_Vm_sub_all
        fprintf(fid_name,'%8s%20s%20s%12d%12.2f%12.2f\n',string(Year_Vm_sub_all(Nxstart_Vm_sub_all(ii))),string(Author_Vm_sub_all(Nxstart_Vm_sub_all(ii))),T_range_Vm_sub_all(ii),N_Vm_sub_all(ii),RMS_Vm_sub_all(ii),AAD_Vm_sub_all(ii));
    end     
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 2

ndata_Vm_melt=length(T_Vm_melt);xxcount_Vm_melt=1; Nxstart_Vm_melt = 1;yy_Vm_melt=Year_Vm_melt(1);AA_Vm_melt = Author_Vm_melt(1);
if ndata_Vm_melt == 1
    Nxend_Vm_melt(xxcount_Vm_melt) = 1;
else
    for ii = 2:ndata_Vm_melt
        if ~strcmp(yy_Vm_melt,string(Year_Vm_melt(ii)))||~strcmp(AA_Vm_melt,Author_Vm_melt(ii))
            Nxend_Vm_melt(xxcount_Vm_melt) = ii-1;
            xxcount_Vm_melt = xxcount_Vm_melt + 1;
            Nxstart_Vm_melt(xxcount_Vm_melt) = ii;
            yy_Vm_melt=Year_Vm_melt(ii);AA_Vm_melt=Author_Vm_melt(ii);
        end
        if ii == ndata_Vm_melt
            Nxend_Vm_melt(xxcount_Vm_melt) = ndata_Vm_melt;
        end
    end
end

    for ii = 1:length(T_Vm_melt)
        props_Vm_melt = computeThermoPropsForResult(T_Vm_melt(ii), p_Vm_melt(ii));
        Vm_melt_deviation(ii) = 100*(Vm_melt(ii) - props_Vm_melt(1)) ./ Vm_melt(ii);
    end


kk=1;
for ii=1:xxcount_Vm_melt
    [Tmax_Vm_melt(kk,:) Tmin_Vm_melt(kk,:)]=maxmin(T_Vm_melt(Nxstart_Vm_melt(ii):Nxend_Vm_melt(ii)));
    if Tmin_Vm_melt(kk,:)>83.99
        Tmin_Vm_melt(kk,:) = round (Tmin_Vm_melt(kk,:),0);
    end
    if Tmax_Vm_melt(kk,:)>83.99
        Tmax_Vm_melt(kk,:) = round (Tmax_Vm_melt(kk,:),0);
    end
        T_range_Vm_melt = strcat(string(Tmin_Vm_melt),"-",string(Tmax_Vm_melt));
     N_Vm_melt(kk,:)=Nxend_Vm_melt(ii)-Nxstart_Vm_melt(ii)+1;
     RMS_Vm_melt(kk,:)=sqrt(sumsqr(Vm_melt_deviation(Nxstart_Vm_melt(ii):Nxend_Vm_melt(ii)))/N_Vm_melt(kk,:));
     AAD_Vm_melt(kk,:)=mean(abs(Vm_melt_deviation(Nxstart_Vm_melt(ii):Nxend_Vm_melt(ii))));
     kk=kk+1;
     
end    
    if ~isfolder('../../paper writing/table/cell volume')
        mkdir('../../paper writing/table/cell volume')
    end
    fid_name = fopen('../../paper writing/table/cell volume/melting.txt','w');
    fprintf(fid_name,'    Year              Author                T/K            N         RMS        AAD\n');
    for ii = 1:xxcount_Vm_melt
        fprintf(fid_name,'%8s%20s%20s%12d%12.2f%12.2f\n',string(Year_Vm_melt(Nxstart_Vm_melt(ii))),string(Author_Vm_melt(Nxstart_Vm_melt(ii))),T_range_Vm_melt(ii),N_Vm_melt(ii),RMS_Vm_melt(ii),AAD_Vm_melt(ii));
    end    
    
ndata_Vm_melt_all=length(T_Vm_melt_all);xxcount_Vm_melt_all=1; Nxstart_Vm_melt_all = 1;yy_Vm_melt_all=Year_Vm_melt_all(1);AA_Vm_melt_all = Author_Vm_melt_all(1);
if ndata_Vm_melt_all == 1
    Nxend_Vm_melt_all(xxcount_Vm_melt_all) = 1;
else
    for ii = 2:ndata_Vm_melt_all
        if ~strcmp(yy_Vm_melt_all,string(Year_Vm_melt_all(ii)))||~strcmp(AA_Vm_melt_all,Author_Vm_melt_all(ii))
            Nxend_Vm_melt_all(xxcount_Vm_melt_all) = ii-1;
            xxcount_Vm_melt_all = xxcount_Vm_melt_all + 1;
            Nxstart_Vm_melt_all(xxcount_Vm_melt_all) = ii;
            yy_Vm_melt_all=Year_Vm_melt_all(ii);AA_Vm_melt_all=Author_Vm_melt_all(ii);
        end
        if ii == ndata_Vm_melt_all
            Nxend_Vm_melt_all(xxcount_Vm_melt_all) = ndata_Vm_melt_all;
        end
    end
end

    for ii = 1:length(T_Vm_melt_all)
        props_Vm_melt_all = computeThermoPropsForResult(T_Vm_melt_all(ii), p_Vm_melt_all(ii));
        Vm_melt_deviation_all(ii) = 100*(Vm_melt_all(ii) - props_Vm_melt_all(1)) ./ Vm_melt_all(ii);
    end
 kk = 1;
for ii=1:xxcount_Vm_melt_all
    [Tmax_Vm_melt_all(kk,:) Tmin_Vm_melt_all(kk,:)]=maxmin(T_Vm_melt_all(Nxstart_Vm_melt_all(ii):Nxend_Vm_melt_all(ii)));
    if Tmin_Vm_melt_all(kk,:)>83.99
        Tmin_Vm_melt_all(kk,:) = round (Tmin_Vm_melt_all(kk,:),0);
    end
    if Tmax_Vm_melt_all(kk,:)>83.99
        Tmax_Vm_melt_all(kk,:) = round (Tmax_Vm_melt_all(kk,:),0);
    end
            T_range_Vm_melt_all = strcat(string(Tmin_Vm_melt_all),"-",string(Tmax_Vm_melt_all));
     N_Vm_melt_all(kk,:)=Nxend_Vm_melt_all(ii)-Nxstart_Vm_melt_all(ii)+1;
     RMS_Vm_melt_all(kk,:)=sqrt(sumsqr(Vm_melt_deviation_all(Nxstart_Vm_melt_all(ii):Nxend_Vm_melt_all(ii)))/N_Vm_melt_all(kk,:));
     AAD_Vm_melt_all(kk,:)=mean(abs(Vm_melt_deviation_all(Nxstart_Vm_melt_all(ii):Nxend_Vm_melt_all(ii))));
     kk=kk+1;
     
end        
    
     if ~isfolder('../../paper writing/table/cell volume')
        mkdir('../../paper writing/table/cell volume')
    end
    fid_name = fopen('../../paper writing/table/cell volume/melting all data.txt','w');
    fprintf(fid_name,'    Year              Author                T/K            N         RMS        AAD\n');
    for ii = 1:xxcount_Vm_melt_all
        fprintf(fid_name,'%8s%20s%20s%12d%12.2f%12.2f\n',string(Year_Vm_melt_all(Nxstart_Vm_melt_all(ii))),string(Author_Vm_melt_all(Nxstart_Vm_melt_all(ii))),T_range_Vm_melt_all(ii),N_Vm_melt_all(ii),RMS_Vm_melt_all(ii),AAD_Vm_melt_all(ii));
    end     
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 3

ndata_cp_sub=length(T_cp_sub);xxcount_cp_sub=1; Nxstart_cp_sub = 1;yy_cp_sub=Year_cp_sub(1);AA_cp_sub = Author_cp_sub(1);
if ndata_cp_sub == 1
    Nxend_cp_sub(xxcount_cp_sub) = 1;
else
    for ii = 2:ndata_cp_sub
        if ~strcmp(yy_cp_sub,string(Year_cp_sub(ii)))||~strcmp(AA_cp_sub,Author_cp_sub(ii))
            Nxend_cp_sub(xxcount_cp_sub) = ii-1;
            xxcount_cp_sub = xxcount_cp_sub + 1;
            Nxstart_cp_sub(xxcount_cp_sub) = ii;
            yy_cp_sub=Year_cp_sub(ii);AA_cp_sub=Author_cp_sub(ii);
        end
        if ii == ndata_cp_sub
            Nxend_cp_sub(xxcount_cp_sub) = ndata_cp_sub;
        end
    end
end

    for ii = 1:length(T_cp_sub)
        props_cp_sub = computeThermoPropsForResult(T_cp_sub(ii), p_cp_sub(ii));
        cp_sub_deviation(ii) = 100*(cp_sub(ii) - props_cp_sub(5)) ./ cp_sub(ii);
    end


kk=1;
for ii=1:xxcount_cp_sub
    [Tmax_cp_sub(kk,:) Tmin_cp_sub(kk,:)]=maxmin(T_cp_sub(Nxstart_cp_sub(ii):Nxend_cp_sub(ii)));
    if Tmin_cp_sub(kk,:)<83.4
        Tmin_cp_sub(kk,:) = round (Tmin_cp_sub(kk,:),0);
    end
    if Tmax_cp_sub(kk,:)<83.4
        Tmax_cp_sub(kk,:) = round (Tmax_cp_sub(kk,:),0);
    end
        T_range_cp_sub = strcat(string(Tmin_cp_sub),"-",string(Tmax_cp_sub));
     N_cp_sub(kk,:)=Nxend_cp_sub(ii)-Nxstart_cp_sub(ii)+1;
     RMS_cp_sub(kk,:)=sqrt(sumsqr(cp_sub_deviation(Nxstart_cp_sub(ii):Nxend_cp_sub(ii)))/N_cp_sub(kk,:));
     AAD_cp_sub(kk,:)=mean(abs(cp_sub_deviation(Nxstart_cp_sub(ii):Nxend_cp_sub(ii))));
     kk=kk+1;
     
end    
    if ~isfolder('../../paper writing/table/heat capacity')
        mkdir('../../paper writing/table/heat capacity')
    end
    fid_name = fopen('../../paper writing/table/heat capacity/sublimation.txt','w');
    fprintf(fid_name,'    Year              Author                T/K            N         RMS        AAD\n');
    for ii = 1:xxcount_cp_sub
        fprintf(fid_name,'%8s%20s%20s%12d%12.2f%12.2f\n',string(Year_cp_sub(Nxstart_cp_sub(ii))),string(Author_cp_sub(Nxstart_cp_sub(ii))),T_range_cp_sub(ii),N_cp_sub(ii),RMS_cp_sub(ii),AAD_cp_sub(ii));
    end    
    
ndata_cp_sub_all=length(T_cp_sub_all);xxcount_cp_sub_all=1; Nxstart_cp_sub_all = 1;yy_cp_sub_all=Year_cp_sub_all(1);AA_cp_sub_all = Author_cp_sub_all(1);
if ndata_cp_sub_all == 1
    Nxend_cp_sub_all(xxcount_cp_sub_all) = 1;
else
    for ii = 2:ndata_cp_sub_all
        if ~strcmp(yy_cp_sub_all,string(Year_cp_sub_all(ii)))||~strcmp(AA_cp_sub_all,Author_cp_sub_all(ii))
            Nxend_cp_sub_all(xxcount_cp_sub_all) = ii-1;
            xxcount_cp_sub_all = xxcount_cp_sub_all + 1;
            Nxstart_cp_sub_all(xxcount_cp_sub_all) = ii;
            yy_cp_sub_all=Year_cp_sub_all(ii);AA_cp_sub_all=Author_cp_sub_all(ii);
        end
        if ii == ndata_cp_sub_all
            Nxend_cp_sub_all(xxcount_cp_sub_all) = ndata_cp_sub_all;
        end
    end
end

    for ii = 1:length(T_cp_sub_all)
        props_cp_sub_all = computeThermoPropsForResult(T_cp_sub_all(ii), p_cp_sub_all(ii));

        cp_sub_deviation_all(ii) = 100*(cp_sub_all(ii) - props_cp_sub_all(5)) ./ cp_sub_all(ii);
        if T_cp_sub_all(ii) == 0
            cp_sub_deviation_all(ii) = 0;
        end        
    end
 kk = 1;
for ii=1:xxcount_cp_sub_all
    [Tmax_cp_sub_all(kk,:) Tmin_cp_sub_all(kk,:)]=maxmin(T_cp_sub_all(Nxstart_cp_sub_all(ii):Nxend_cp_sub_all(ii)));
    if Tmin_cp_sub_all(kk,:)<83.4
        Tmin_cp_sub_all(kk,:) = round (Tmin_cp_sub_all(kk,:),0);
    end
    if Tmax_cp_sub_all(kk,:)<83.4
        Tmax_cp_sub_all(kk,:) = round (Tmax_cp_sub_all(kk,:),0);
    end
            T_range_cp_sub_all = strcat(string(Tmin_cp_sub_all),"-",string(Tmax_cp_sub_all));
     N_cp_sub_all(kk,:)=Nxend_cp_sub_all(ii)-Nxstart_cp_sub_all(ii)+1;
     RMS_cp_sub_all(kk,:)=sqrt(sumsqr(cp_sub_deviation_all(Nxstart_cp_sub_all(ii):Nxend_cp_sub_all(ii)))/N_cp_sub_all(kk,:));
     AAD_cp_sub_all(kk,:)=mean(abs(cp_sub_deviation_all(Nxstart_cp_sub_all(ii):Nxend_cp_sub_all(ii))));
     kk=kk+1;
     
end        
    
     if ~isfolder('../../paper writing/table/heat capacity')
        mkdir('../../paper writing/table/heat capacity')
    end
    fid_name = fopen('../../paper writing/table/heat capacity/sublimation all data.txt','w');
    fprintf(fid_name,'    Year              Author                T/K            N         RMS        AAD\n');
    for ii = 1:xxcount_cp_sub_all
        fprintf(fid_name,'%8s%20s%20s%12d%12.2f%12.2f\n',string(Year_cp_sub_all(Nxstart_cp_sub_all(ii))),string(Author_cp_sub_all(Nxstart_cp_sub_all(ii))),T_range_cp_sub_all(ii),N_cp_sub_all(ii),RMS_cp_sub_all(ii),AAD_cp_sub_all(ii));
    end     
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 4

ndata_alpha_sub=length(T_alpha_sub);xxcount_alpha_sub=1; Nxstart_alpha_sub = 1;yy_alpha_sub=Year_alpha_sub(1);AA_alpha_sub = Author_alpha_sub(1);
if ndata_alpha_sub == 1
    Nxend_alpha_sub(xxcount_alpha_sub) = 1;
else
    for ii = 2:ndata_alpha_sub
        if ~strcmp(yy_alpha_sub,string(Year_alpha_sub(ii)))||~strcmp(AA_alpha_sub,Author_alpha_sub(ii))
            Nxend_alpha_sub(xxcount_alpha_sub) = ii-1;
            xxcount_alpha_sub = xxcount_alpha_sub + 1;
            Nxstart_alpha_sub(xxcount_alpha_sub) = ii;
            yy_alpha_sub=Year_alpha_sub(ii);AA_alpha_sub=Author_alpha_sub(ii);
        end
        if ii == ndata_alpha_sub
            Nxend_alpha_sub(xxcount_alpha_sub) = ndata_alpha_sub;
        end
    end
end

    for ii = 1:length(T_alpha_sub)
        props_alpha_sub = computeThermoPropsForResult(T_alpha_sub(ii), p_alpha_sub(ii));
        alpha_sub_deviation(ii) = 100*(alpha_sub(ii) - props_alpha_sub(4)) ./ alpha_sub(ii);
    end


kk=1;
for ii=1:xxcount_alpha_sub
    [Tmax_alpha_sub(kk,:) Tmin_alpha_sub(kk,:)]=maxmin(T_alpha_sub(Nxstart_alpha_sub(ii):Nxend_alpha_sub(ii)));
    if Tmin_alpha_sub(kk,:)<83.4
        Tmin_alpha_sub(kk,:) = round (Tmin_alpha_sub(kk,:),0);
    end
    if Tmax_alpha_sub(kk,:)<83.4
        Tmax_alpha_sub(kk,:) = round (Tmax_alpha_sub(kk,:),0);
    end
        T_range_alpha_sub = strcat(string(Tmin_alpha_sub),"-",string(Tmax_alpha_sub));
     N_alpha_sub(kk,:)=Nxend_alpha_sub(ii)-Nxstart_alpha_sub(ii)+1;
     RMS_alpha_sub(kk,:)=sqrt(sumsqr(alpha_sub_deviation(Nxstart_alpha_sub(ii):Nxend_alpha_sub(ii)))/N_alpha_sub(kk,:));
     AAD_alpha_sub(kk,:)=mean(abs(alpha_sub_deviation(Nxstart_alpha_sub(ii):Nxend_alpha_sub(ii))));
     kk=kk+1;
     
end    
    if ~isfolder('../../paper writing/table/thermal expansivity')
        mkdir('../../paper writing/table/thermal expansivity')
    end
    fid_name = fopen('../../paper writing/table/thermal expansivity/sublimation.txt','w');
    fprintf(fid_name,'    Year              Author                T/K            N         RMS        AAD\n');
    for ii = 1:xxcount_alpha_sub
        fprintf(fid_name,'%8s%20s%20s%12d%12.2f%12.2f\n',string(Year_alpha_sub(Nxstart_alpha_sub(ii))),string(Author_alpha_sub(Nxstart_alpha_sub(ii))),T_range_alpha_sub(ii),N_alpha_sub(ii),RMS_alpha_sub(ii),AAD_alpha_sub(ii));
    end    
    
ndata_alpha_sub_all=length(T_alpha_sub_all);xxcount_alpha_sub_all=1; Nxstart_alpha_sub_all = 1;yy_alpha_sub_all=Year_alpha_sub_all(1);AA_alpha_sub_all = Author_alpha_sub_all(1);
if ndata_alpha_sub_all == 1
    Nxend_alpha_sub_all(xxcount_alpha_sub_all) = 1;
else
    for ii = 2:ndata_alpha_sub_all
        if ~strcmp(yy_alpha_sub_all,string(Year_alpha_sub_all(ii)))||~strcmp(AA_alpha_sub_all,Author_alpha_sub_all(ii))
            Nxend_alpha_sub_all(xxcount_alpha_sub_all) = ii-1;
            xxcount_alpha_sub_all = xxcount_alpha_sub_all + 1;
            Nxstart_alpha_sub_all(xxcount_alpha_sub_all) = ii;
            yy_alpha_sub_all=Year_alpha_sub_all(ii);AA_alpha_sub_all=Author_alpha_sub_all(ii);
        end
        if ii == ndata_alpha_sub_all
            Nxend_alpha_sub_all(xxcount_alpha_sub_all) = ndata_alpha_sub_all;
        end
    end
end

    for ii = 1:length(T_alpha_sub_all)
        props_alpha_sub_all = computeThermoPropsForResult(T_alpha_sub_all(ii), p_alpha_sub_all(ii));

        alpha_sub_deviation_all(ii) = real(100*(alpha_sub_all(ii) - props_alpha_sub_all(4)) ./ alpha_sub_all(ii));
        if T_alpha_sub_all(ii) == 0
            alpha_sub_deviation_all(ii) = 0;
        end        
    end
 kk = 1;
for ii=1:xxcount_alpha_sub_all
    [Tmax_alpha_sub_all(kk,:) Tmin_alpha_sub_all(kk,:)]=maxmin(T_alpha_sub_all(Nxstart_alpha_sub_all(ii):Nxend_alpha_sub_all(ii)));
    if Tmin_alpha_sub_all(kk,:)<83.4
        Tmin_alpha_sub_all(kk,:) = round (Tmin_alpha_sub_all(kk,:),0);
    end
    if Tmax_alpha_sub_all(kk,:)<83.4
        Tmax_alpha_sub_all(kk,:) = round (Tmax_alpha_sub_all(kk,:),0);
    end
            T_range_alpha_sub_all = strcat(string(Tmin_alpha_sub_all),"-",string(Tmax_alpha_sub_all));
     N_alpha_sub_all(kk,:)=Nxend_alpha_sub_all(ii)-Nxstart_alpha_sub_all(ii)+1;
     RMS_alpha_sub_all(kk,:)=sqrt(sumsqr(alpha_sub_deviation_all(Nxstart_alpha_sub_all(ii):Nxend_alpha_sub_all(ii)))/N_alpha_sub_all(kk,:));
     AAD_alpha_sub_all(kk,:)=mean(abs(alpha_sub_deviation_all(Nxstart_alpha_sub_all(ii):Nxend_alpha_sub_all(ii))));
     kk=kk+1;
     
end        
    
     if ~isfolder('../../paper writing/table/thermal expansivity')
        mkdir('../../paper writing/table/thermal expansivity')
    end
    fid_name = fopen('../../paper writing/table/thermal expansivity/sublimation all data.txt','w');
    fprintf(fid_name,'    Year              Author                T/K            N         RMS        AAD\n');
    for ii = 1:xxcount_alpha_sub_all
        fprintf(fid_name,'%8s%20s%20s%12d%12.2f%12.2f\n',string(Year_alpha_sub_all(Nxstart_alpha_sub_all(ii))),string(Author_alpha_sub_all(Nxstart_alpha_sub_all(ii))),T_range_alpha_sub_all(ii),N_alpha_sub_all(ii),RMS_alpha_sub_all(ii),AAD_alpha_sub_all(ii));
    end     
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 5

ndata_BetaT_sub=length(T_BetaT_sub);xxcount_BetaT_sub=1; Nxstart_BetaT_sub = 1;yy_BetaT_sub=Year_BetaT_sub(1);AA_BetaT_sub = Author_BetaT_sub(1);
if ndata_BetaT_sub == 1
    Nxend_BetaT_sub(xxcount_BetaT_sub) = 1;
else
    for ii = 2:ndata_BetaT_sub
        if ~strcmp(yy_BetaT_sub,string(Year_BetaT_sub(ii)))||~strcmp(AA_BetaT_sub,Author_BetaT_sub(ii))
            Nxend_BetaT_sub(xxcount_BetaT_sub) = ii-1;
            xxcount_BetaT_sub = xxcount_BetaT_sub + 1;
            Nxstart_BetaT_sub(xxcount_BetaT_sub) = ii;
            yy_BetaT_sub=Year_BetaT_sub(ii);AA_BetaT_sub=Author_BetaT_sub(ii);
        end
        if ii == ndata_BetaT_sub
            Nxend_BetaT_sub(xxcount_BetaT_sub) = ndata_BetaT_sub;
        end
    end
end

    for ii = 1:length(T_BetaT_sub)
        props_BetaT_sub = computeThermoPropsForResult(T_BetaT_sub(ii), p_BetaT_sub(ii));
        BetaT_sub_deviation(ii) = 100*(BetaT_sub(ii) - 1./props_BetaT_sub(2)) ./ BetaT_sub(ii);
    end


kk=1;
for ii=1:xxcount_BetaT_sub
    [Tmax_BetaT_sub(kk,:) Tmin_BetaT_sub(kk,:)]=maxmin(T_BetaT_sub(Nxstart_BetaT_sub(ii):Nxend_BetaT_sub(ii)));
    if Tmin_BetaT_sub(kk,:)<83.4
        Tmin_BetaT_sub(kk,:) = round (Tmin_BetaT_sub(kk,:),0);
    end
    if Tmax_BetaT_sub(kk,:)<83.4
        Tmax_BetaT_sub(kk,:) = round (Tmax_BetaT_sub(kk,:),0);
    end
        T_range_BetaT_sub = strcat(string(Tmin_BetaT_sub),"-",string(Tmax_BetaT_sub));
     N_BetaT_sub(kk,:)=Nxend_BetaT_sub(ii)-Nxstart_BetaT_sub(ii)+1;
     RMS_BetaT_sub(kk,:)=sqrt(sumsqr(BetaT_sub_deviation(Nxstart_BetaT_sub(ii):Nxend_BetaT_sub(ii)))/N_BetaT_sub(kk,:));
     AAD_BetaT_sub(kk,:)=mean(abs(BetaT_sub_deviation(Nxstart_BetaT_sub(ii):Nxend_BetaT_sub(ii))));
     kk=kk+1;
     
end    
    if ~isfolder('../../paper writing/table/bulk modulus')
        mkdir('../../paper writing/table/bulk modulus')
    end
    fid_name = fopen('../../paper writing/table/bulk modulus/BetaT sublimation.txt','w');
    fprintf(fid_name,'    Year              Author                T/K            N         RMS        AAD\n');
    for ii = 1:xxcount_BetaT_sub
        fprintf(fid_name,'%8s%20s%20s%12d%12.2f%12.2f\n',string(Year_BetaT_sub(Nxstart_BetaT_sub(ii))),string(Author_BetaT_sub(Nxstart_BetaT_sub(ii))),T_range_BetaT_sub(ii),N_BetaT_sub(ii),RMS_BetaT_sub(ii),AAD_BetaT_sub(ii));
    end    
    
ndata_BetaT_sub_all=length(T_BetaT_sub_all);xxcount_BetaT_sub_all=1; Nxstart_BetaT_sub_all = 1;yy_BetaT_sub_all=Year_BetaT_sub_all(1);AA_BetaT_sub_all = Author_BetaT_sub_all(1);
if ndata_BetaT_sub_all == 1
    Nxend_BetaT_sub_all(xxcount_BetaT_sub_all) = 1;
else
    for ii = 2:ndata_BetaT_sub_all
        if ~strcmp(yy_BetaT_sub_all,string(Year_BetaT_sub_all(ii)))||~strcmp(AA_BetaT_sub_all,Author_BetaT_sub_all(ii))
            Nxend_BetaT_sub_all(xxcount_BetaT_sub_all) = ii-1;
            xxcount_BetaT_sub_all = xxcount_BetaT_sub_all + 1;
            Nxstart_BetaT_sub_all(xxcount_BetaT_sub_all) = ii;
            yy_BetaT_sub_all=Year_BetaT_sub_all(ii);AA_BetaT_sub_all=Author_BetaT_sub_all(ii);
        end
        if ii == ndata_BetaT_sub_all
            Nxend_BetaT_sub_all(xxcount_BetaT_sub_all) = ndata_BetaT_sub_all;
        end
    end
end

    for ii = 1:length(T_BetaT_sub_all)
        props_BetaT_sub_all = computeThermoPropsForResult(T_BetaT_sub_all(ii), p_BetaT_sub_all(ii));

        BetaT_sub_deviation_all(ii) = real(100*(BetaT_sub_all(ii) - 1./props_BetaT_sub_all(2)) ./ BetaT_sub_all(ii));
        if T_BetaT_sub_all(ii) == 0
            BetaT_sub_deviation_all(ii) = 0;
        end        
    end
 kk = 1;
for ii=1:xxcount_BetaT_sub_all
    [Tmax_BetaT_sub_all(kk,:) Tmin_BetaT_sub_all(kk,:)]=maxmin(T_BetaT_sub_all(Nxstart_BetaT_sub_all(ii):Nxend_BetaT_sub_all(ii)));
    if Tmin_BetaT_sub_all(kk,:)<83.4
        Tmin_BetaT_sub_all(kk,:) = round (Tmin_BetaT_sub_all(kk,:),0);
    end
    if Tmax_BetaT_sub_all(kk,:)<83.4
        Tmax_BetaT_sub_all(kk,:) = round (Tmax_BetaT_sub_all(kk,:),0);
    end
            T_range_BetaT_sub_all = strcat(string(Tmin_BetaT_sub_all),"-",string(Tmax_BetaT_sub_all));
     N_BetaT_sub_all(kk,:)=Nxend_BetaT_sub_all(ii)-Nxstart_BetaT_sub_all(ii)+1;
     RMS_BetaT_sub_all(kk,:)=sqrt(sumsqr(BetaT_sub_deviation_all(Nxstart_BetaT_sub_all(ii):Nxend_BetaT_sub_all(ii)))/N_BetaT_sub_all(kk,:));
     AAD_BetaT_sub_all(kk,:)=mean(abs(BetaT_sub_deviation_all(Nxstart_BetaT_sub_all(ii):Nxend_BetaT_sub_all(ii))));
     kk=kk+1;
     
end        
    
     if ~isfolder('../../paper writing/table/bulk modulus')
        mkdir('../../paper writing/table/bulk modulus')
    end
    fid_name = fopen('../../paper writing/table/bulk modulus/BetaT sublimation all data.txt','w');
    fprintf(fid_name,'    Year              Author                T/K            N         RMS        AAD\n');
    for ii = 1:xxcount_BetaT_sub_all
        fprintf(fid_name,'%8s%20s%20s%12d%12.2f%12.2f\n',string(Year_BetaT_sub_all(Nxstart_BetaT_sub_all(ii))),string(Author_BetaT_sub_all(Nxstart_BetaT_sub_all(ii))),T_range_BetaT_sub_all(ii),N_BetaT_sub_all(ii),RMS_BetaT_sub_all(ii),AAD_BetaT_sub_all(ii));
    end     
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 6

ndata_BetaS_sub=length(T_BetaS_sub);xxcount_BetaS_sub=1; Nxstart_BetaS_sub = 1;yy_BetaS_sub=Year_BetaS_sub(1);AA_BetaS_sub = Author_BetaS_sub(1);
if ndata_BetaS_sub == 1
    Nxend_BetaS_sub(xxcount_BetaS_sub) = 1;
else
    for ii = 2:ndata_BetaS_sub
        if ~strcmp(yy_BetaS_sub,string(Year_BetaS_sub(ii)))||~strcmp(AA_BetaS_sub,Author_BetaS_sub(ii))
            Nxend_BetaS_sub(xxcount_BetaS_sub) = ii-1;
            xxcount_BetaS_sub = xxcount_BetaS_sub + 1;
            Nxstart_BetaS_sub(xxcount_BetaS_sub) = ii;
            yy_BetaS_sub=Year_BetaS_sub(ii);AA_BetaS_sub=Author_BetaS_sub(ii);
        end
        if ii == ndata_BetaS_sub
            Nxend_BetaS_sub(xxcount_BetaS_sub) = ndata_BetaS_sub;
        end
    end
end

    for ii = 1:length(T_BetaS_sub)
        props_BetaS_sub = computeThermoPropsForResult(T_BetaS_sub(ii), p_BetaS_sub(ii));
        BetaS_sub_deviation(ii) = 100*(BetaS_sub(ii) - 1./props_BetaS_sub(3)) ./ BetaS_sub(ii);
    end


kk=1;
for ii=1:xxcount_BetaS_sub
    [Tmax_BetaS_sub(kk,:) Tmin_BetaS_sub(kk,:)]=maxmin(T_BetaS_sub(Nxstart_BetaS_sub(ii):Nxend_BetaS_sub(ii)));
    if Tmin_BetaS_sub(kk,:)<83.4
        Tmin_BetaS_sub(kk,:) = round (Tmin_BetaS_sub(kk,:),0);
    end
    if Tmax_BetaS_sub(kk,:)<83.4
        Tmax_BetaS_sub(kk,:) = round (Tmax_BetaS_sub(kk,:),0);
    end
        T_range_BetaS_sub = strcat(string(Tmin_BetaS_sub),"-",string(Tmax_BetaS_sub));
     N_BetaS_sub(kk,:)=Nxend_BetaS_sub(ii)-Nxstart_BetaS_sub(ii)+1;
     RMS_BetaS_sub(kk,:)=sqrt(sumsqr(BetaS_sub_deviation(Nxstart_BetaS_sub(ii):Nxend_BetaS_sub(ii)))/N_BetaS_sub(kk,:));
     AAD_BetaS_sub(kk,:)=mean(abs(BetaS_sub_deviation(Nxstart_BetaS_sub(ii):Nxend_BetaS_sub(ii))));
     kk=kk+1;
     
end    
    if ~isfolder('../../paper writing/table/bulk modulus')
        mkdir('../../paper writing/table/bulk modulus')
    end
    fid_name = fopen('../../paper writing/table/bulk modulus/BetaS sublimation.txt','w');
    fprintf(fid_name,'    Year              Author                T/K            N         RMS        AAD\n');
    for ii = 1:xxcount_BetaS_sub
        fprintf(fid_name,'%8s%20s%20s%12d%12.2f%12.2f\n',string(Year_BetaS_sub(Nxstart_BetaS_sub(ii))),string(Author_BetaS_sub(Nxstart_BetaS_sub(ii))),T_range_BetaS_sub(ii),N_BetaS_sub(ii),RMS_BetaS_sub(ii),AAD_BetaS_sub(ii));
    end    
    
ndata_BetaS_sub_all=length(T_BetaS_sub_all);xxcount_BetaS_sub_all=1; Nxstart_BetaS_sub_all = 1;yy_BetaS_sub_all=Year_BetaS_sub_all(1);AA_BetaS_sub_all = Author_BetaS_sub_all(1);
if ndata_BetaS_sub_all == 1
    Nxend_BetaS_sub_all(xxcount_BetaS_sub_all) = 1;
else
    for ii = 2:ndata_BetaS_sub_all
        if ~strcmp(yy_BetaS_sub_all,string(Year_BetaS_sub_all(ii)))||~strcmp(AA_BetaS_sub_all,Author_BetaS_sub_all(ii))
            Nxend_BetaS_sub_all(xxcount_BetaS_sub_all) = ii-1;
            xxcount_BetaS_sub_all = xxcount_BetaS_sub_all + 1;
            Nxstart_BetaS_sub_all(xxcount_BetaS_sub_all) = ii;
            yy_BetaS_sub_all=Year_BetaS_sub_all(ii);AA_BetaS_sub_all=Author_BetaS_sub_all(ii);
        end
        if ii == ndata_BetaS_sub_all
            Nxend_BetaS_sub_all(xxcount_BetaS_sub_all) = ndata_BetaS_sub_all;
        end
    end
end

    for ii = 1:length(T_BetaS_sub_all)
        props_BetaS_sub_all = computeThermoPropsForResult(T_BetaS_sub_all(ii), p_BetaS_sub_all(ii));

        BetaS_sub_deviation_all(ii) = real(100*(BetaS_sub_all(ii) - 1./props_BetaS_sub_all(3)) ./ BetaS_sub_all(ii));
        if T_BetaS_sub_all(ii) == 0
            BetaS_sub_deviation_all(ii) = 0;
        end        
    end
 kk = 1;
for ii=1:xxcount_BetaS_sub_all
    [Tmax_BetaS_sub_all(kk,:) Tmin_BetaS_sub_all(kk,:)]=maxmin(T_BetaS_sub_all(Nxstart_BetaS_sub_all(ii):Nxend_BetaS_sub_all(ii)));
    if Tmin_BetaS_sub_all(kk,:)<83.4
        Tmin_BetaS_sub_all(kk,:) = round (Tmin_BetaS_sub_all(kk,:),0);
    end
    if Tmax_BetaS_sub_all(kk,:)<83.4
        Tmax_BetaS_sub_all(kk,:) = round (Tmax_BetaS_sub_all(kk,:),0);
    end
            T_range_BetaS_sub_all = strcat(string(Tmin_BetaS_sub_all),"-",string(Tmax_BetaS_sub_all));
     N_BetaS_sub_all(kk,:)=Nxend_BetaS_sub_all(ii)-Nxstart_BetaS_sub_all(ii)+1;
     RMS_BetaS_sub_all(kk,:)=sqrt(sumsqr(BetaS_sub_deviation_all(Nxstart_BetaS_sub_all(ii):Nxend_BetaS_sub_all(ii)))/N_BetaS_sub_all(kk,:));
     AAD_BetaS_sub_all(kk,:)=mean(abs(BetaS_sub_deviation_all(Nxstart_BetaS_sub_all(ii):Nxend_BetaS_sub_all(ii))));
     kk=kk+1;
     
end        
    
     if ~isfolder('../../paper writing/table/bulk modulus')
        mkdir('../../paper writing/table/bulk modulus')
    end
    fid_name = fopen('../../paper writing/table/bulk modulus/BetaS sublimation all data.txt','w');
    fprintf(fid_name,'    Year              Author                T/K            N         RMS        AAD\n');
    for ii = 1:xxcount_BetaS_sub_all
        fprintf(fid_name,'%8s%20s%20s%12d%12.2f%12.2f\n',string(Year_BetaS_sub_all(Nxstart_BetaS_sub_all(ii))),string(Author_BetaS_sub_all(Nxstart_BetaS_sub_all(ii))),T_range_BetaS_sub_all(ii),N_BetaS_sub_all(ii),RMS_BetaS_sub_all(ii),AAD_BetaS_sub_all(ii));
    end     
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 7
    ndata_H_sub=length(T_H_sub);xxcount_H_sub=1; Nxstart_H_sub = 1;yy_H_sub=Year_H_sub(1);AA_H_sub = Author_H_sub(1);
    if ndata_H_sub == 1
        Nxend_H_sub(xxcount_H_sub) = 1;
    else
        for ii = 2:ndata_H_sub
            if ~strcmp(yy_H_sub,string(Year_H_sub(ii)))||~strcmp(AA_H_sub,Author_H_sub(ii))
                Nxend_H_sub(xxcount_H_sub) = ii-1;
                xxcount_H_sub = xxcount_H_sub + 1;
                Nxstart_H_sub(xxcount_H_sub) = ii;
                yy_H_sub=Year_H_sub(ii);AA_H_sub=Author_H_sub(ii);
            end
            if ii == ndata_H_sub
                Nxend_H_sub(xxcount_H_sub) = ndata_H_sub;
            end
        end
    end   

    St_REFPROP = 131.1500;%J/mol/K
    Ht_REFPROP = 1689.10;%J/molProps(K1038,K1039,12,$B$7:$B$36)+K1038*O1041
    Tt = 83.8058; pt = 0.068891;%unit: K and MPa
    deltaS_triple =  130.37 - St_REFPROP;
    props_Triple = computeThermoPropsForResult(Tt,pt);
    Ht_fitted = props_Triple(12) + Tt*130.37;
    deltaH_triple = Ht_fitted - Ht_REFPROP;

    for ii = 1:length(T_H_sub)
            H_solid_sub(ii) = H_fluid_sub(ii)-1000*delta_H_sub(ii);
            fitted_props_H_sub = computeThermoPropsForResult(T_H_sub(ii), p_H_sub(ii));

            H_solid_sub_fitted(ii) = fitted_props_H_sub(11) - deltaH_triple;
            H_solid_sub_deviation(ii) = 100*(H_solid_sub(ii) - H_solid_sub_fitted(ii)) / H_solid_sub(ii) ;   

    end 
kk=1;
for ii=1:xxcount_H_sub
    [Tmax_H_sub(kk,:) Tmin_H_sub(kk,:)]=maxmin(T_H_sub(Nxstart_H_sub(ii):Nxend_H_sub(ii)));
    if Tmin_H_sub(kk,:)<83.4
        Tmin_H_sub(kk,:) = round (Tmin_H_sub(kk,:),0);
    end
    if Tmax_H_sub(kk,:)<83.4
        Tmax_H_sub(kk,:) = round (Tmax_H_sub(kk,:),0);
    end
        T_range_H_sub = strcat(string(Tmin_H_sub),"-",string(Tmax_H_sub));
     N_H_sub(kk,:)=Nxend_H_sub(ii)-Nxstart_H_sub(ii)+1;
     RMS_H_sub(kk,:)=sqrt(sumsqr(H_solid_sub_deviation(Nxstart_H_sub(ii):Nxend_H_sub(ii)))/N_H_sub(kk,:));
     AAD_H_sub(kk,:)=mean(abs(H_solid_sub_deviation(Nxstart_H_sub(ii):Nxend_H_sub(ii))));
     kk=kk+1;
     
end    
    if ~isfolder('../../paper writing/table/enthalpy')
        mkdir('../../paper writing/table/enthalpy')
    end
    fid_name = fopen('../../paper writing/table/enthalpy/H sublimation.txt','w');
    fprintf(fid_name,'    Year              Author                T/K            N         RMS        AAD\n');
    for ii = 1:xxcount_H_sub
        fprintf(fid_name,'%8s%20s%20s%12d%12.2f%12.2f\n',string(Year_H_sub(Nxstart_H_sub(ii))),string(Author_H_sub(Nxstart_H_sub(ii))),T_range_H_sub(ii),N_H_sub(ii),RMS_H_sub(ii),AAD_H_sub(ii));
    end        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 8
    ndata_H_melt=length(T_H_melt);xxcount_H_melt=1; Nxstart_H_melt = 1;yy_H_melt=Year_H_melt(1);AA_H_melt = Author_H_melt(1);
    if ndata_H_melt == 1
        Nxend_H_melt(xxcount_H_melt) = 1;
    else
        for ii = 2:ndata_H_melt
            if ~strcmp(yy_H_melt,string(Year_H_melt(ii)))||~strcmp(AA_H_melt,Author_H_melt(ii))
                Nxend_H_melt(xxcount_H_melt) = ii-1;
                xxcount_H_melt = xxcount_H_melt + 1;
                Nxstart_H_melt(xxcount_H_melt) = ii;
                yy_H_melt=Year_H_melt(ii);AA_H_melt=Author_H_melt(ii);
            end
            if ii == ndata_H_melt
                Nxend_H_melt(xxcount_H_melt) = ndata_H_melt;
            end
        end
    end   

    St_REFPROP = 131.1500;%J/mol/K
    Ht_REFPROP = 1689.10;%J/molProps(K1038,K1039,12,$B$7:$B$36)+K1038*O1041
    Tt = 83.8058; pt = 0.068891;%unit: K and MPa
    deltaS_triple =  130.37 - St_REFPROP;
    props_Triple = computeThermoPropsForResult(Tt,pt);
    Ht_fitted = props_Triple(12) + Tt*130.37;
    deltaH_triple = Ht_fitted - Ht_REFPROP;

    for ii = 1:length(T_H_melt)
            H_solid_melt(ii) = H_fluid_melt(ii)-1000*delta_H_melt(ii);
            fitted_props_H_melt = computeThermoPropsForResult(T_H_melt(ii), p_H_melt(ii));

            H_solid_melt_fitted(ii) = fitted_props_H_melt(11) - deltaH_triple;
            H_solid_melt_deviation(ii) = 100*(H_solid_melt(ii) - H_solid_melt_fitted(ii)) / H_solid_melt(ii) ;   

    end 
kk=1;
for ii=1:xxcount_H_melt
    [Tmax_H_melt(kk,:) Tmin_H_melt(kk,:)]=maxmin(T_H_melt(Nxstart_H_melt(ii):Nxend_H_melt(ii)));
    if Tmin_H_melt(kk,:)>83.99
        Tmin_H_melt(kk,:) = round (Tmin_H_melt(kk,:),0);
    end
    if Tmax_H_melt(kk,:)>83.99
        Tmax_H_melt(kk,:) = round (Tmax_H_melt(kk,:),0);
    end
        T_range_H_melt = strcat(string(Tmin_H_melt),"-",string(Tmax_H_melt));
     N_H_melt(kk,:)=Nxend_H_melt(ii)-Nxstart_H_melt(ii)+1;
     RMS_H_melt(kk,:)=sqrt(sumsqr(H_solid_melt_deviation(Nxstart_H_melt(ii):Nxend_H_melt(ii)))/N_H_melt(kk,:));
     AAD_H_melt(kk,:)=mean(abs(H_solid_melt_deviation(Nxstart_H_melt(ii):Nxend_H_melt(ii))));
     kk=kk+1;
     
end    
    if ~isfolder('../../paper writing/table/enthalpy')
        mkdir('../../paper writing/table/enthalpy')
    end
    fid_name = fopen('../../paper writing/table/enthalpy/H melting.txt','w');
    fprintf(fid_name,'    Year              Author                T/K            N         RMS        AAD\n');
    for ii = 1:xxcount_H_melt
        fprintf(fid_name,'%8s%20s%20s%12d%12.2f%12.2f\n',string(Year_H_melt(Nxstart_H_melt(ii))),string(Author_H_melt(Nxstart_H_melt(ii))),T_range_H_melt(ii),N_H_melt(ii),RMS_H_melt(ii),AAD_H_melt(ii));
    end     
    ndata_H_melt_all=length(T_H_melt_all);xxcount_H_melt_all=1; Nxstart_H_melt_all = 1;yy_H_melt_all=Year_H_melt_all(1);AA_H_melt_all = Author_H_melt_all(1);
    if ndata_H_melt_all == 1
        Nxend_H_melt_all(xxcount_H_melt_all) = 1;
    else
        for ii = 2:ndata_H_melt_all
            if ~strcmp(yy_H_melt_all,string(Year_H_melt_all(ii)))||~strcmp(AA_H_melt_all,Author_H_melt_all(ii))
                Nxend_H_melt_all(xxcount_H_melt_all) = ii-1;
                xxcount_H_melt_all = xxcount_H_melt_all + 1;
                Nxstart_H_melt_all(xxcount_H_melt_all) = ii;
                yy_H_melt_all=Year_H_melt(ii);AA_H_melt_all=Author_H_melt_all(ii);
            end
            if ii == ndata_H_melt_all
                Nxend_H_melt_all(xxcount_H_melt_all) = ndata_H_melt_all;
            end
        end
    end   

    St_REFPROP = 131.1500;%J/mol/K
    Ht_REFPROP = 1689.10;%J/molProps(K1038,K1039,12,$B$7:$B$36)+K1038*O1041
    Tt = 83.8058; pt = 0.068891;%unit: K and MPa
    deltaS_triple =  130.37 - St_REFPROP;
    props_Triple = computeThermoPropsForResult(Tt,pt);
    Ht_fitted = props_Triple(12) + Tt*130.37;
    deltaH_triple = Ht_fitted - Ht_REFPROP;

    for ii = 1:length(T_H_melt_all)
            H_solid_melt_all(ii) = H_fluid_melt_all(ii)-1000*delta_H_melt_all(ii);
            fitted_props_H_melt_all = computeThermoPropsForResult(T_H_melt_all(ii), p_H_melt_all(ii));

            H_solid_melt_fitted_all(ii) = fitted_props_H_melt_all(11) - deltaH_triple;
            H_solid_melt_deviation_all(ii) = 100*(H_solid_melt_all(ii) - H_solid_melt_fitted_all(ii)) / H_solid_melt_all(ii) ;   

    end 
kk=1;
for ii=1:xxcount_H_melt_all
    [Tmax_H_melt_all(kk,:) Tmin_H_melt_all(kk,:)]=maxmin(T_H_melt_all(Nxstart_H_melt_all(ii):Nxend_H_melt_all(ii)));
    if Tmin_H_melt_all(kk,:)>83.99
        Tmin_H_melt_all(kk,:) = round (Tmin_H_melt_all(kk,:),0);
    end
    if Tmax_H_melt_all(kk,:)>83.99
        Tmax_H_melt_all(kk,:) = round (Tmax_H_melt_all(kk,:),0);
    end
        T_range_H_melt_all = strcat(string(Tmin_H_melt_all),"-",string(Tmax_H_melt_all));
     N_H_melt_all(kk,:)=Nxend_H_melt_all(ii)-Nxstart_H_melt_all(ii)+1;
     RMS_H_melt_all(kk,:)=sqrt(sumsqr(H_solid_melt_deviation_all(Nxstart_H_melt_all(ii):Nxend_H_melt_all(ii)))/N_H_melt_all(kk,:));
     AAD_H_melt_all(kk,:)=mean(abs(H_solid_melt_deviation_all(Nxstart_H_melt_all(ii):Nxend_H_melt_all(ii))));
     kk=kk+1;
     
end    
    if ~isfolder('../../paper writing/table/enthalpy')
        mkdir('../../paper writing/table/enthalpy')
    end
    fid_name = fopen('../../paper writing/table/enthalpy/H melting all.txt','w');
    fprintf(fid_name,'    Year              Author                T/K            N         RMS        AAD\n');
    for ii = 1:xxcount_H_melt_all
        fprintf(fid_name,'%8s%20s%20s%12d%12.2f%12.2f\n',string(Year_H_melt_all(Nxstart_H_melt_all(ii))),string(Author_H_melt_all(Nxstart_H_melt_all(ii))),T_range_H_melt_all(ii),N_H_melt_all(ii),RMS_H_melt_all(ii),AAD_H_melt_all(ii));
    end       
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 9
    ndata_sub=length(T_sub);xxcount_sub=1; Nxstart_sub = 1;yy_sub=Year_sub(1);AA_sub = Author_sub(1);
    if ndata_sub == 1
        Nxend_sub(xxcount_sub) = 1;
    else
        for ii = 2:ndata_sub
            if ~strcmp(yy_sub,string(Year_sub(ii)))||~strcmp(AA_sub,Author_sub(ii))
                Nxend_sub(xxcount_sub) = ii-1;
                xxcount_sub = xxcount_sub + 1;
                Nxstart_sub(xxcount_sub) = ii;
                yy_sub=Year_sub(ii);AA_sub=Author_sub(ii);
            end
            if ii == ndata_sub
                Nxend_sub(xxcount_sub) = ndata_sub;
            end
        end
    end    


    St_REFPROP = 131.1500;%J/mol/K
    Ht_REFPROP = 1689.10;%J/molProps(K1038,K1039,12,$B$7:$B$36)+K1038*O1041
    Tt = 83.8058; pt = 0.068891;%unit: K and MPa
    deltaS_triple =  130.37 - St_REFPROP;
    props_Triple = computeThermoPropsForResult(Tt,pt);
    Ht_fitted = props_Triple(12) + Tt*130.37;
    deltaH_triple = Ht_fitted - Ht_REFPROP;

    for ii = 1:length(T_sub)

        fitted_props_sub = computeThermoPropsForResult(T_sub(ii), p_sub(ii)/10^6);
        delta_G_sub(ii,:) = G_fluid_sub(ii) - fitted_props_sub(12) + deltaH_triple - T_sub(ii) * deltaS_triple;
        p_fitted_sub(ii) = p_sub(ii) - delta_G_sub(ii)/(V_fluid_sub(ii) - fitted_props_sub(1)) * 10^6; % unit Pa
        p_sub_deviation(ii) = 100 * (p_sub(ii) - p_fitted_sub(ii)) / p_sub(ii);

    end    
kk=1;

for ii=1:xxcount_sub
    [Tmax_sub(kk,:) Tmin_sub(kk,:)]=maxmin(T_sub(Nxstart_sub(ii):Nxend_sub(ii)));
    if Tmin_sub(kk,:)<83.4
        Tmin_sub(kk,:) = round (Tmin_sub(kk,:),0);
    end
    if Tmax_sub(kk,:)<83.4
        Tmax_sub(kk,:) = round (Tmax_sub(kk,:),0);
    end
        T_range_sub = strcat(string(Tmin_sub),"-",string(Tmax_sub));
     N_sub(kk,:)=Nxend_sub(ii)-Nxstart_sub(ii)+1;
     RMS_sub(kk,:)=sqrt(sumsqr(p_sub_deviation(Nxstart_sub(ii):Nxend_sub(ii)))/N_sub(kk,:));
     AAD_sub(kk,:)=mean(abs(p_sub_deviation(Nxstart_sub(ii):Nxend_sub(ii))));
     kk=kk+1;
     
end    
    if ~isfolder('../../paper writing/table/sublimation')
        mkdir('../../paper writing/table/sublimation')
    end
    fid_name = fopen('../../paper writing/table/sublimation/sublimation.txt','w');
    fprintf(fid_name,'    Year              Author                T/K            N         RMS        AAD\n');
    for ii = 1:xxcount_sub
        fprintf(fid_name,'%8s%20s%20s%12d%12.2f%12.2f\n',string(Year_sub(Nxstart_sub(ii))),string(Author_sub(Nxstart_sub(ii))),T_range_sub(ii),N_sub(ii),RMS_sub(ii),AAD_sub(ii));
    end         
    
    ndata_sub_all=length(T_sub_all);xxcount_sub_all=1; Nxstart_sub_all = 1;yy_sub_all=Year_sub_all(1);AA_sub_all = Author_sub_all(1);
    if ndata_sub_all == 1
        Nxend_sub_all(xxcount_sub_all) = 1;
    else
        for ii = 2:ndata_sub_all
            if ~strcmp(yy_sub_all,string(Year_sub_all(ii)))||~strcmp(AA_sub_all,Author_sub_all(ii))
                Nxend_sub_all(xxcount_sub_all) = ii-1;
                xxcount_sub_all = xxcount_sub_all + 1;
                Nxstart_sub_all(xxcount_sub_all) = ii;
                yy_sub_all=Year_sub_all(ii);AA_sub_all=Author_sub_all(ii);
            end
            if ii == ndata_sub_all
                Nxend_sub_all(xxcount_sub_all) = ndata_sub_all;
            end
        end
    end    


    St_REFPROP = 131.1500;%J/mol/K
    Ht_REFPROP = 1689.10;%J/molProps(K1038,K1039,12,$B$7:$B$36)+K1038*O1041
    Tt = 83.8058; pt = 0.068891;%unit: K and MPa
    deltaS_triple =  130.37 - St_REFPROP;
    props_Triple = computeThermoPropsForResult(Tt,pt);
    Ht_fitted = props_Triple(12) + Tt*130.37;
    deltaH_triple = Ht_fitted - Ht_REFPROP;

    for ii = 1:length(T_sub_all)

        fitted_props_sub_all = computeThermoPropsForResult(T_sub_all(ii), p_sub_all(ii)/10^6);
        delta_G_sub_all(ii,:) = G_fluid_sub_all(ii) - fitted_props_sub_all(12) + deltaH_triple - T_sub_all(ii) * deltaS_triple;
        p_fitted_sub_all(ii) = p_sub_all(ii) - delta_G_sub_all(ii)/(V_fluid_sub_all(ii) - fitted_props_sub_all(1)) * 10^6; % unit Pa
        p_sub_deviation_all(ii) = 100 * (p_sub_all(ii) - p_fitted_sub_all(ii)) / p_sub_all(ii);

    end    
kk=1;

for ii=1:xxcount_sub_all
    [Tmax_sub_all(kk,:) Tmin_sub_all(kk,:)]=maxmin(T_sub_all(Nxstart_sub_all(ii):Nxend_sub_all(ii)));
    if Tmin_sub_all(kk,:)<83.4
        Tmin_sub_all(kk,:) = round (Tmin_sub_all(kk,:),0);
    end
    if Tmax_sub_all(kk,:)<83.4
        Tmax_sub_all(kk,:) = round (Tmax_sub_all(kk,:),0);
    end
        T_range_sub_all = strcat(string(Tmin_sub_all),"-",string(Tmax_sub_all));
     N_sub_all(kk,:)=Nxend_sub_all(ii)-Nxstart_sub_all(ii)+1;
     RMS_sub_all(kk,:)=sqrt(sumsqr(p_sub_deviation_all(Nxstart_sub_all(ii):Nxend_sub_all(ii)))/N_sub_all(kk,:));
     AAD_sub_all(kk,:)=mean(abs(p_sub_deviation_all(Nxstart_sub_all(ii):Nxend_sub_all(ii))));
     kk=kk+1;
     
end    
    if ~isfolder('../../paper writing/table/sublimation')
        mkdir('../../paper writing/table/sublimation')
    end
    fid_name = fopen('../../paper writing/table/sublimation/sublimation_all.txt','w');
    fprintf(fid_name,'    Year              Author                T/K            N         RMS        AAD\n');
    for ii = 1:xxcount_sub_all
        fprintf(fid_name,'%8s%20s%20s%12d%12.2f%12.2f\n',string(Year_sub_all(Nxstart_sub_all(ii))),string(Author_sub_all(Nxstart_sub_all(ii))),T_range_sub_all(ii),N_sub_all(ii),RMS_sub_all(ii),AAD_sub_all(ii));
    end         
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 10
    ndata_melt=length(T_melt);xxcount_melt=1; Nxstart_melt = 1;yy_melt=Year_melt(1);AA_melt = Author_melt(1);
    if ndata_melt == 1
        Nxend_melt(xxcount_melt) = 1;
    else
        for ii = 2:ndata_melt
            if ~strcmp(yy_melt,string(Year_melt(ii)))||~strcmp(AA_melt,Author_melt(ii))
                Nxend_melt(xxcount_melt) = ii-1;
                xxcount_melt = xxcount_melt + 1;
                Nxstart_melt(xxcount_melt) = ii;
                yy_melt=Year_melt(ii);AA_melt=Author_melt(ii);
            end
            if ii == ndata_melt
                Nxend_melt(xxcount_melt) = ndata_melt;
            end
        end
    end    


    St_REFPROP = 131.1500;%J/mol/K
    Ht_REFPROP = 1689.10;%J/molProps(K1038,K1039,12,$B$7:$B$36)+K1038*O1041
    Tt = 83.8058; pt = 0.068891;%unit: K and MPa
    deltaS_triple =  130.37 - St_REFPROP;
    props_Triple = computeThermoPropsForResult(Tt,pt);
    Ht_fitted = props_Triple(12) + Tt*130.37;
    deltaH_triple = Ht_fitted - Ht_REFPROP;

    for ii = 1:length(T_melt)

        fitted_props_melt = computeThermoPropsForResult(T_melt(ii), p_melt(ii));
        delta_G_melt(ii,:) = G_fluid_melt(ii) - fitted_props_melt(12) + deltaH_triple - T_melt(ii) * deltaS_triple;
        p_fitted_melt(ii) = p_melt(ii) - delta_G_melt(ii)/(V_fluid_melt(ii) - fitted_props_melt(1)) ; % unit MPa
        p_melt_deviation(ii) = 100 * (p_melt(ii) - p_fitted_melt(ii)) / p_melt(ii);

    end    
kk=1;

for ii=1:xxcount_melt
    [Tmax_melt(kk,:) Tmin_melt(kk,:)]=maxmin(T_melt(Nxstart_melt(ii):Nxend_melt(ii)));
    if Tmin_melt(kk,:)>83.99
        Tmin_melt(kk,:) = round (Tmin_melt(kk,:),0);
    end
    if Tmax_melt(kk,:)>83.99
        Tmax_melt(kk,:) = round (Tmax_melt(kk,:),0);
    end
        T_range_melt = strcat(string(Tmin_melt),"-",string(Tmax_melt));
     N_melt(kk,:)=Nxend_melt(ii)-Nxstart_melt(ii)+1;
     RMS_melt(kk,:)=sqrt(sumsqr(p_melt_deviation(Nxstart_melt(ii):Nxend_melt(ii)))/N_melt(kk,:));
     AAD_melt(kk,:)=mean(abs(p_melt_deviation(Nxstart_melt(ii):Nxend_melt(ii))));
     kk=kk+1;
     
end    
    if ~isfolder('../../paper writing/table/melting')
        mkdir('../../paper writing/table/melting')
    end
    fid_name = fopen('../../paper writing/table/melting/melting.txt','w');
    fprintf(fid_name,'    Year              Author                T/K            N         RMS        AAD\n');
    for ii = 1:xxcount_melt
        fprintf(fid_name,'%8s%20s%20s%12d%12.2f%12.2f\n',string(Year_melt(Nxstart_melt(ii))),string(Author_melt(Nxstart_melt(ii))),T_range_melt(ii),N_melt(ii),RMS_melt(ii),AAD_melt(ii));
    end         
    
    ndata_melt_all=length(T_melt_all);xxcount_melt_all=1; Nxstart_melt_all = 1;yy_melt_all=Year_melt_all(1);AA_melt_all = Author_melt_all(1);
    if ndata_melt_all == 1
        Nxend_melt_all(xxcount_melt_all) = 1;
    else
        for ii = 2:ndata_melt_all
            if ~strcmp(yy_melt_all,string(Year_melt_all(ii)))||~strcmp(AA_melt_all,Author_melt_all(ii))
                Nxend_melt_all(xxcount_melt_all) = ii-1;
                xxcount_melt_all = xxcount_melt_all + 1;
                Nxstart_melt_all(xxcount_melt_all) = ii;
                yy_melt_all=Year_melt_all(ii);AA_melt_all=Author_melt_all(ii);
            end
            if ii == ndata_melt_all
                Nxend_melt_all(xxcount_melt_all) = ndata_melt_all;
            end
        end
    end    


    St_REFPROP = 131.1500;%J/mol/K
    Ht_REFPROP = 1689.10;%J/molProps(K1038,K1039,12,$B$7:$B$36)+K1038*O1041
    Tt = 83.8058; pt = 0.068891;%unit: K and MPa
    deltaS_triple =  130.37 - St_REFPROP;
    props_Triple = computeThermoPropsForResult(Tt,pt);
    Ht_fitted = props_Triple(12) + Tt*130.37;
    deltaH_triple = Ht_fitted - Ht_REFPROP;

    for ii = 1:length(T_melt_all)

        fitted_props_melt_all = computeThermoPropsForResult(T_melt_all(ii), p_melt_all(ii));
        delta_G_melt_all(ii,:) = G_fluid_melt_all(ii) - fitted_props_melt_all(12) + deltaH_triple - T_melt_all(ii) * deltaS_triple;
        p_fitted_melt_all(ii) = p_melt_all(ii) - delta_G_melt_all(ii)/(V_fluid_melt_all(ii) - fitted_props_melt_all(1)) ; % unit MPa
        p_melt_deviation_all(ii) = 100 * (p_melt_all(ii) - p_fitted_melt_all(ii)) / p_melt_all(ii);

    end    
kk=1;

for ii=1:xxcount_melt_all
    [Tmax_melt_all(kk,:) Tmin_melt_all(kk,:)]=maxmin(T_melt_all(Nxstart_melt_all(ii):Nxend_melt_all(ii)));
    if Tmin_melt_all(kk,:)>83.99
        Tmin_melt_all(kk,:) = round (Tmin_melt_all(kk,:),0);
    end
    if Tmax_melt_all(kk,:)>83.99
        Tmax_melt_all(kk,:) = round (Tmax_melt_all(kk,:),0);
    end
        T_range_melt_all = strcat(string(Tmin_melt_all),"-",string(Tmax_melt_all));
     N_melt_all(kk,:)=Nxend_melt_all(ii)-Nxstart_melt_all(ii)+1;
     RMS_melt_all(kk,:)=sqrt(sumsqr(p_melt_deviation_all(Nxstart_melt_all(ii):Nxend_melt_all(ii)))/N_melt_all(kk,:));
     AAD_melt_all(kk,:)=mean(abs(p_melt_deviation_all(Nxstart_melt_all(ii):Nxend_melt_all(ii))));
     kk=kk+1;
     
end    
    if ~isfolder('../../paper writing/table/melting')
        mkdir('../../paper writing/table/melting')
    end
    fid_name = fopen('../../paper writing/table/melting/melting_all.txt','w');
    fprintf(fid_name,'    Year              Author                T/K            N         RMS        AAD\n');
    for ii = 1:xxcount_melt_all
        fprintf(fid_name,'%8s%20s%20s%12d%12.2f%12.2f\n',string(Year_melt_all(Nxstart_melt_all(ii))),string(Author_melt_all(Nxstart_melt_all(ii))),T_range_melt_all(ii),N_melt_all(ii),RMS_melt_all(ii),AAD_melt_all(ii));
    end         
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 11
ndata_Vm_highp=length(T_Vm_highp);xxcount_Vm_highp=1; Nxstart_Vm_highp = 1;AA_Vm_highp=Author_Vm_highp(1);YY_Vm_highp=Year_Vm_highp(1);
for ii = 2:ndata_Vm_highp         
        if ~strcmp(AA_Vm_highp,Author_Vm_highp(ii))||~strcmp(YY_Vm_highp,Year_Vm_highp(ii))
            Nxend_Vm_highp(xxcount_Vm_highp) = ii-1;
            xxcount_Vm_highp = xxcount_Vm_highp + 1;
            Nxstart_Vm_highp(xxcount_Vm_highp) = ii;
            AA_Vm_highp=Author_Vm_highp(ii);YY_Vm_highp=Year_Vm_highp(ii);
        end
        if ii == ndata_Vm_highp
            Nxend_Vm_highp(xxcount_Vm_highp) = ndata_Vm_highp;
        end
end

for ii = 1:length(T_Vm_highp)
            props_Vm_highp = computeThermoPropsForResult(T_Vm_highp(ii), p_Vm_highp(ii));
            Vm_highp_deviation(ii,:) = 100*(Vm_highp(ii) - props_Vm_highp(1)) ./ Vm_highp(ii);    
end
for jj=1:xxcount_Vm_highp
    ndata_Vm_highp=length(T_Vm_highp(Nxstart_Vm_highp(jj):Nxend_Vm_highp(jj)));yycount_Vm_highp=1; Nystart_Vm_highp(jj,yycount_Vm_highp) = Nxstart_Vm_highp(jj);TT_Vm_highp=T_Vm_highp(Nxstart_Vm_highp(jj));
    for ii = 2:ndata_Vm_highp   

            %if ~strcmp(AA,Author(Nxstart(jj)+ii-1))&&~strcmp(YY,Year(Nxstart(jj)+ii-1))
            if abs(T_Vm_highp(Nxstart_Vm_highp(jj)+ii-1) - TT_Vm_highp)> 0.05
                Nyend_Vm_highp(jj,yycount_Vm_highp) = Nxstart_Vm_highp(jj)+ii-2;
                yycount_Vm_highp = yycount_Vm_highp + 1;
                Nystart_Vm_highp(jj,yycount_Vm_highp) = Nxstart_Vm_highp(jj)+ii-1;
                TT_Vm_highp=T_Vm_highp(Nxstart_Vm_highp(jj)+ii-1);
            end
            if ii == ndata_Vm_highp
                Nyend_Vm_highp(jj,yycount_Vm_highp) = ndata_Vm_highp + Nxstart_Vm_highp(jj)-1;
            end
    end     
end
[a1_Vm_highp a2_Vm_highp]=size(Nystart_Vm_highp);


    
 
 for ii = 1:a1_Vm_highp
    for jj=1:a2_Vm_highp
        Vm_highp_dev = [];
        if Nystart_Vm_highp(ii,jj)~=0
        ll1_Vm_highp(ii,jj)=string(strcat(Author_Vm_highp(Nystart_Vm_highp(ii,jj)),Year_Vm_highp(Nystart_Vm_highp(ii,jj)),' - ',{' '},char(string(T_Vm_highp(Nystart_Vm_highp(ii,jj)))),{' '},'K'));            
        for  kkk = 1:length(p_Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)))
            props_Vm_highp = computeThermoPropsForResult(T_Vm_highp(Nystart_Vm_highp(ii,jj)+kkk-1), p_Vm_highp(Nystart_Vm_highp(ii,jj)+kkk-1));
            Vm_highp_dev(kkk,:) = 100*(Vm_highp(Nystart_Vm_highp(ii,jj)+kkk-1) - props_Vm_highp(1)) ./ Vm_highp(Nystart_Vm_highp(ii,jj)+kkk-1);
        end
        

        end
    end    
 end

kk=1;
for ii=1:xxcount_Vm_highp
    [Tmax_Vm_highp(kk,:) Tmin_Vm_highp(kk,:)]=maxmin(T_Vm_highp(Nxstart_Vm_highp(ii):Nxend_Vm_highp(ii)));
    [pmax_Vm_highp(kk,:) pmin_Vm_highp(kk,:)]=maxmin(p_Vm_highp(Nxstart_Vm_highp(ii):Nxend_Vm_highp(ii)));
%     if Tmin_Vm_highp(kk,:)<83.4
        Tmin_Vm_highp(kk,:) = round (Tmin_Vm_highp(kk,:),0);
        pmin_Vm_highp(kk,:) = round (pmin_Vm_highp(kk,:),0);
%     end
%     if Tmax_Vm_highp(kk,:)<83.4
        Tmax_Vm_highp(kk,:) = round (Tmax_Vm_highp(kk,:),0);
        pmax_Vm_highp(kk,:) = round (pmax_Vm_highp(kk,:),0);
%     end
        T_range_Vm_highp = strcat(string(Tmin_Vm_highp),"-",string(Tmax_Vm_highp));
        p_range_Vm_highp = strcat(string(pmin_Vm_highp),"-",string(pmax_Vm_highp));
     N_Vm_highp(kk,:)=Nxend_Vm_highp(ii)-Nxstart_Vm_highp(ii)+1;
     RMS_Vm_highp(kk,:)=sqrt(sumsqr(Vm_highp_deviation(Nxstart_Vm_highp(ii):Nxend_Vm_highp(ii)))/N_Vm_highp(kk,:));
     AAD_Vm_highp(kk,:)=mean(abs(Vm_highp_deviation(Nxstart_Vm_highp(ii):Nxend_Vm_highp(ii))));
     kk=kk+1;
     
end    
    if ~isfolder('../../paper writing/table/cell volume/high pressure')
        mkdir('../../paper writing/table/cell volume/high pressure')
    end
    fid_name = fopen('../../paper writing/table/cell volume/high pressure/high pressure.txt','w');
    fprintf(fid_name,'    Year              Author                T/K               p/MPa           N         RMS         AAD\n');
    for ii = 1:xxcount_Vm_highp
        fprintf(fid_name,'%8s%20s%20s%20s%12d%12.2f%12.2f\n',string(Year_Vm_highp(Nxstart_Vm_highp(ii))),string(Author_Vm_highp(Nxstart_Vm_highp(ii))),T_range_Vm_highp(ii),p_range_Vm_highp(ii),N_Vm_highp(ii),RMS_Vm_highp(ii),AAD_Vm_highp(ii));
    end     
ndata_Vm_highp_all=length(T_Vm_highp_all);xxcount_Vm_highp_all=1; Nxstart_Vm_highp_all = 1;AA_Vm_highp_all=Author_Vm_highp_all(1);YY_Vm_highp_all=Year_Vm_highp_all(1);
for ii = 2:ndata_Vm_highp_all         
        if ~strcmp(AA_Vm_highp_all,Author_Vm_highp_all(ii))||~strcmp(YY_Vm_highp_all,Year_Vm_highp_all(ii))
            Nxend_Vm_highp_all(xxcount_Vm_highp_all) = ii-1;
            xxcount_Vm_highp_all = xxcount_Vm_highp_all + 1;
            Nxstart_Vm_highp_all(xxcount_Vm_highp_all) = ii;
            AA_Vm_highp_all=Author_Vm_highp_all(ii);YY_Vm_highp_all=Year_Vm_highp_all(ii);
        end
        if ii == ndata_Vm_highp_all
            Nxend_Vm_highp_all(xxcount_Vm_highp_all) = ndata_Vm_highp_all;
        end
end

for ii = 1:length(T_Vm_highp_all)
            props_Vm_highp_all = computeThermoPropsForResult(T_Vm_highp_all(ii), p_Vm_highp_all(ii));
            Vm_highp_all_deviation(ii,:) = 100*(Vm_highp_all(ii) - props_Vm_highp_all(1)) ./ Vm_highp_all(ii);    
end
for jj=1:xxcount_Vm_highp_all
    ndata_Vm_highp_all=length(T_Vm_highp_all(Nxstart_Vm_highp_all(jj):Nxend_Vm_highp_all(jj)));yycount_Vm_highp_all=1; Nystart_Vm_highp_all(jj,yycount_Vm_highp_all) = Nxstart_Vm_highp_all(jj);TT_Vm_highp_all=T_Vm_highp_all(Nxstart_Vm_highp_all(jj));
    for ii = 2:ndata_Vm_highp_all   

            %if ~strcmp(AA,Author(Nxstart(jj)+ii-1))&&~strcmp(YY,Year(Nxstart(jj)+ii-1))
            if abs(T_Vm_highp_all(Nxstart_Vm_highp_all(jj)+ii-1) - TT_Vm_highp_all)> 0.05
                Nyend_Vm_highp_all(jj,yycount_Vm_highp_all) = Nxstart_Vm_highp_all(jj)+ii-2;
                yycount_Vm_highp_all = yycount_Vm_highp_all + 1;
                Nystart_Vm_highp_all(jj,yycount_Vm_highp_all) = Nxstart_Vm_highp_all(jj)+ii-1;
                TT_Vm_highp_all=T_Vm_highp_all(Nxstart_Vm_highp_all(jj)+ii-1);
            end
            if ii == ndata_Vm_highp_all
                Nyend_Vm_highp_all(jj,yycount_Vm_highp_all) = ndata_Vm_highp_all + Nxstart_Vm_highp_all(jj)-1;
            end
    end     
end
[a1_Vm_highp_all a2_Vm_highp_all]=size(Nystart_Vm_highp_all);


    
 
 for ii = 1:a1_Vm_highp_all
    for jj=1:a2_Vm_highp_all
        Vm_highp_all_dev = [];
        if Nystart_Vm_highp_all(ii,jj)~=0
        ll1_Vm_highp_all(ii,jj)=string(strcat(Author_Vm_highp_all(Nystart_Vm_highp_all(ii,jj)),Year_Vm_highp_all(Nystart_Vm_highp_all(ii,jj)),' - ',{' '},char(string(T_Vm_highp_all(Nystart_Vm_highp_all(ii,jj)))),{' '},'K'));            
        for  kkk = 1:length(p_Vm_highp_all(Nystart_Vm_highp_all(ii,jj):Nyend_Vm_highp_all(ii,jj)))
            props_Vm_highp_all = computeThermoPropsForResult(T_Vm_highp_all(Nystart_Vm_highp_all(ii,jj)+kkk-1), p_Vm_highp_all(Nystart_Vm_highp_all(ii,jj)+kkk-1));
            Vm_highp_all_dev(kkk,:) = 100*(Vm_highp_all(Nystart_Vm_highp_all(ii,jj)+kkk-1) - props_Vm_highp_all(1)) ./ Vm_highp_all(Nystart_Vm_highp_all(ii,jj)+kkk-1);
        end
        

        end
    end    
 end

kk=1;
for ii=1:xxcount_Vm_highp_all
    [Tmax_Vm_highp_all(kk,:) Tmin_Vm_highp_all(kk,:)]=maxmin(T_Vm_highp_all(Nxstart_Vm_highp_all(ii):Nxend_Vm_highp_all(ii)));
    [pmax_Vm_highp_all(kk,:) pmin_Vm_highp_all(kk,:)]=maxmin(p_Vm_highp_all(Nxstart_Vm_highp_all(ii):Nxend_Vm_highp_all(ii)));
%     if Tmin_Vm_highp_all(kk,:)<83.4
        Tmin_Vm_highp_all(kk,:) = round (Tmin_Vm_highp_all(kk,:),0);
        pmin_Vm_highp_all(kk,:) = round (pmin_Vm_highp_all(kk,:),0);
%     end
%     if Tmax_Vm_highp_all(kk,:)<83.4
        Tmax_Vm_highp_all(kk,:) = round (Tmax_Vm_highp_all(kk,:),0);
        pmax_Vm_highp_all(kk,:) = round (pmax_Vm_highp_all(kk,:),0);
%     end
        T_range_Vm_highp_all = strcat(string(Tmin_Vm_highp_all),"-",string(Tmax_Vm_highp_all));
        p_range_Vm_highp_all = strcat(string(pmin_Vm_highp_all),"-",string(pmax_Vm_highp_all));
     N_Vm_highp_all(kk,:)=Nxend_Vm_highp_all(ii)-Nxstart_Vm_highp_all(ii)+1;
     RMS_Vm_highp_all(kk,:)=sqrt(sumsqr(Vm_highp_all_deviation(Nxstart_Vm_highp_all(ii):Nxend_Vm_highp_all(ii)))/N_Vm_highp_all(kk,:));
     AAD_Vm_highp_all(kk,:)=mean(abs(Vm_highp_all_deviation(Nxstart_Vm_highp_all(ii):Nxend_Vm_highp_all(ii))));
     kk=kk+1;
     
end    
    if ~isfolder('../../paper writing/table/cell volume/high pressure')
        mkdir('../../paper writing/table/cell volume/high pressure')
    end
    fid_name = fopen('../../paper writing/table/cell volume/high pressure/high pressure_all.txt','w');
    fprintf(fid_name,'    Year              Author                T/K               p/MPa           N         RMS         AAD\n');
    for ii = 1:xxcount_Vm_highp_all
        fprintf(fid_name,'%8s%20s%20s%20s%12d%12.2f%12.2f\n',string(Year_Vm_highp_all(Nxstart_Vm_highp_all(ii))),string(Author_Vm_highp_all(Nxstart_Vm_highp_all(ii))),T_range_Vm_highp_all(ii),p_range_Vm_highp_all(ii),N_Vm_highp_all(ii),RMS_Vm_highp_all(ii),AAD_Vm_highp_all(ii));
    end    
 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 12
ndata_KT_highp_all=length(T_KT_highp_all);xxcount_KT_highp_all=1; Nxstart_KT_highp_all = 1;AA_KT_highp_all=Author_KT_highp_all(1);YY_KT_highp_all=Year_KT_highp_all(1);
for ii = 2:ndata_KT_highp_all         
        if ~strcmp(AA_KT_highp_all,Author_KT_highp_all(ii))||~strcmp(YY_KT_highp_all,Year_KT_highp_all(ii))
            Nxend_KT_highp_all(xxcount_KT_highp_all) = ii-1;
            xxcount_KT_highp_all = xxcount_KT_highp_all + 1;
            Nxstart_KT_highp_all(xxcount_KT_highp_all) = ii;
            AA_KT_highp_all=Author_KT_highp_all(ii);YY_KT_highp_all=Year_KT_highp_all(ii);
        end
        if ii == ndata_KT_highp_all
            Nxend_KT_highp_all(xxcount_KT_highp_all) = ndata_KT_highp_all;
        end
end

for ii = 1:length(T_KT_highp_all)
            props_KT_highp_all = computeThermoPropsForResult(T_KT_highp_all(ii), p_KT_highp_all(ii));
            KT_highp_all_deviation(ii,:) = 100*(KT_highp_all(ii) - 1./props_KT_highp_all(2)) ./ KT_highp_all(ii);    
end
for jj=1:xxcount_KT_highp_all
    ndata_KT_highp_all=length(T_KT_highp_all(Nxstart_KT_highp_all(jj):Nxend_KT_highp_all(jj)));yycount_KT_highp_all=1; Nystart_KT_highp_all(jj,yycount_KT_highp_all) = Nxstart_KT_highp_all(jj);TT_KT_highp_all=T_KT_highp_all(Nxstart_KT_highp_all(jj));
    for ii = 2:ndata_KT_highp_all   

            %if ~strcmp(AA,Author(Nxstart(jj)+ii-1))&&~strcmp(YY,Year(Nxstart(jj)+ii-1))
            if abs(T_KT_highp_all(Nxstart_KT_highp_all(jj)+ii-1) - TT_KT_highp_all)> 0.05
                Nyend_KT_highp_all(jj,yycount_KT_highp_all) = Nxstart_KT_highp_all(jj)+ii-2;
                yycount_KT_highp_all = yycount_KT_highp_all + 1;
                Nystart_KT_highp_all(jj,yycount_KT_highp_all) = Nxstart_KT_highp_all(jj)+ii-1;
                TT_KT_highp_all=T_KT_highp_all(Nxstart_KT_highp_all(jj)+ii-1);
            end
            if ii == ndata_KT_highp_all
                Nyend_KT_highp_all(jj,yycount_KT_highp_all) = ndata_KT_highp_all + Nxstart_KT_highp_all(jj)-1;
            end
    end     
end
[a1_KT_highp_all a2_KT_highp_all]=size(Nystart_KT_highp_all);


    
 
 for ii = 1:a1_KT_highp_all
    for jj=1:a2_KT_highp_all
        KT_highp_all_dev = [];
        if Nystart_KT_highp_all(ii,jj)~=0
        ll1_KT_highp_all(ii,jj)=string(strcat(Author_KT_highp_all(Nystart_KT_highp_all(ii,jj)),Year_KT_highp_all(Nystart_KT_highp_all(ii,jj)),' - ',{' '},char(string(T_KT_highp_all(Nystart_KT_highp_all(ii,jj)))),{' '},'K'));            
        for  kkk = 1:length(p_KT_highp_all(Nystart_KT_highp_all(ii,jj):Nyend_KT_highp_all(ii,jj)))
            props_KT_highp_all = computeThermoPropsForResult(T_KT_highp_all(Nystart_KT_highp_all(ii,jj)+kkk-1), p_KT_highp_all(Nystart_KT_highp_all(ii,jj)+kkk-1));
            KT_highp_all_dev(kkk,:) = 100*(KT_highp_all(Nystart_KT_highp_all(ii,jj)+kkk-1) - 1./props_KT_highp_all(2)) ./ KT_highp_all(Nystart_KT_highp_all(ii,jj)+kkk-1);
        end
        

        end
    end    
 end

kk=1;
for ii=1:xxcount_KT_highp_all
    [Tmax_KT_highp_all(kk,:) Tmin_KT_highp_all(kk,:)]=maxmin(T_KT_highp_all(Nxstart_KT_highp_all(ii):Nxend_KT_highp_all(ii)));
    [pmax_KT_highp_all(kk,:) pmin_KT_highp_all(kk,:)]=maxmin(p_KT_highp_all(Nxstart_KT_highp_all(ii):Nxend_KT_highp_all(ii)));
%     if Tmin_KT_highp_all(kk,:)<83.4
        Tmin_KT_highp_all(kk,:) = round (Tmin_KT_highp_all(kk,:),0);
        pmin_KT_highp_all(kk,:) = round (pmin_KT_highp_all(kk,:),0);
%     end
%     if Tmax_KT_highp_all(kk,:)<83.4
        Tmax_KT_highp_all(kk,:) = round (Tmax_KT_highp_all(kk,:),0);
        pmax_KT_highp_all(kk,:) = round (pmax_KT_highp_all(kk,:),0);
%     end
        T_range_KT_highp_all = strcat(string(Tmin_KT_highp_all),"-",string(Tmax_KT_highp_all));
        p_range_KT_highp_all = strcat(string(pmin_KT_highp_all),"-",string(pmax_KT_highp_all));
     N_KT_highp_all(kk,:)=Nxend_KT_highp_all(ii)-Nxstart_KT_highp_all(ii)+1;
     RMS_KT_highp_all(kk,:)=sqrt(sumsqr(KT_highp_all_deviation(Nxstart_KT_highp_all(ii):Nxend_KT_highp_all(ii)))/N_KT_highp_all(kk,:));
     AAD_KT_highp_all(kk,:)=mean(abs(KT_highp_all_deviation(Nxstart_KT_highp_all(ii):Nxend_KT_highp_all(ii))));
     kk=kk+1;
     
end    
    if ~isfolder('../../paper writing/table/bulk modulus/high pressure')
        mkdir('../../paper writing/table/bulk modulus/high pressure')
    end
    fid_name = fopen('../../paper writing/table/bulk modulus/high pressure/KT_all.txt','w');
    fprintf(fid_name,'    Year              Author                T/K               p/MPa           N         RMS         AAD\n');
    for ii = 1:xxcount_KT_highp_all
        fprintf(fid_name,'%8s%20s%20s%20s%12d%12.2f%12.2f\n',string(Year_KT_highp_all(Nxstart_KT_highp_all(ii))),string(Author_KT_highp_all(Nxstart_KT_highp_all(ii))),T_range_KT_highp_all(ii),p_range_KT_highp_all(ii),N_KT_highp_all(ii),RMS_KT_highp_all(ii),AAD_KT_highp_all(ii));
    end    
 
end