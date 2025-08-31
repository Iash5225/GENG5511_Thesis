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
Lplot = 1;% 0 not plot fitted results, 1 plot fitted results
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

% params_init = [22.501,2614.9,6518,10,0,0,0,0,0,88.8,0,0,0,0,0,2.659,0,0,0,0,0,0.2987,0,0,0,0,0,0.0501,0.4482,3.2797,129.2201];
%  params_init = [22.5010000000000,2666.57712657095,7182.10495276709,10.9738793460463,0,0,0,0,0,88.4187937382078,0,0,0,0,0,2.64968221299643,0,0,0,0,0,0.0654976996469354,0,0,0,0,0,0.00315255988663863,0.751332648200701,2.85818383962980,130.802847769017];
% params_init = [22.5010000000000,2666.57712657095,7182.10495276709,10.9738793460463,0,0,0,0,0,88.4187937382078,0,0,0,0,0,2.64968221299643,0,0,0,0,0,0.0654976996469354,0,0,0,0,0,0.00315255988663863,0.751332648200701,2.85818383962980,130.802847769017];
% params_init = [22.5010000000000,2613.36513197677,6578.91823295083,0.0479456358591460,0,0,0,0,0,88.7597567375365,0,0,0,0,0,2.70590227942827,0,0,0,0,0,0.215325374971741,0,0,0,0,0,0.000108354627003047,1.22941489874094,40.1401450079614,130.844240762090];
% params_init = [22.5010000000000,2652.75143379056,6764.67416880259,13.3270930005519,0,0,0,0,0,86.0144705897861,0,0,0,0,0,2.77291040952701,0,0,0,0,0,0.218180728566795,0,0,0,0,0,0.0392401424864426,1.29209316748706,7.53419859973110,130.831775423514];
% params_init = [22.5010000000000,2650.66771721636,6785.17253936522,7.56883396146445,0,0,0,0,0,85.8764565484837,0,0,0,0,0,2.77251504982636,0,0,0,0,0,0.159028250445040,0,0,0,0,1.49011611938477e-08,0.0514363814775774,1.18772044364822,5.89338184883926,130.560654696078];
% params_init = [22.5550000000000,2652.61707352136,6789.85675106941,7.03167243246522,0,0,0,0,0,85.8713527567015,0,0,0,0,0,2.77455541271994,0,0,0,0,0,0.168126634286510,0,0,0,0,0,0.0531435506621265,1.22677902896779,5.91290901831193,130.559391749957];
% params_init = [22.5550000000000,2652.74895279117,6951.21407167918,0.420340270319679,0,0,0,0,0,85.7898103422512,0,0,0,0,0,2.70081913270137,0,0,0,0,0,0.109595017680243,0,0,0,0,0,0.0368775749277790,0.828530956043237,4.88844291487550,130.423272215059];
% params_init = [22.5550000000000,2652.74895279117,6951.21407167918,0.420340270336971,0,0,0,0,0,85.7898103422497,0,0,0,0,0,2.70081913276557,0,0,0,0,0,0.109595017560591,0,0,0,0,0,0.0368775749292855,0.828530955386335,4.88844291460752,130.423272215062];
% params_init = [22.5550000000000,2672.75850009876,6949.65942735616,37.4426690934005,0,0,0,0,0,85.7934248179273,0,0,0,0,0,2.72643751288280,0,0,0,0,0,0.170743661090858,0,0,0,0,0,0.0422480934716084,0.914355989686453,5.15001804301763,130.381434124793];
% params_init = [22.555,2672.86748983869,6949.65115605911,37.4430581602660,0,0,0,0,0,85.7773190654182,0,0,0,0,0,2.75779856999379,0,0,0,0,0,0.181208946004173,0,0,0,0,0,0.0559743462286786,1.15134309868676,5.35640975026643,130.287044493036];
% params_init = [22.555,2672.88384644784,6949.64990867093,37.4435047856296,0,0,0,0,0,85.8645580253154,0,0,0,0,0,2.78453949246262,0,0,0,0,0,0.150127659969475,0,0,0,0,0,0.0639733915576612,1.28270820832069,5.45226468940159,130.163550381449];
% params_init = [22.555,2672.88280896774,6949.65158502347,37.4170137655152,0,0,0,0,0,86.0285735331876,0,0,0,0,0,2.74947423008234,0,0,0,0,0,-0.0371313605412209,0,0,0,0,0,0.0271439066613556,0.968498322095311,6.78918293945055,130.435543778892];
% params_init = [22.555,2672.88280879896,6949.65158504018,37.4170137631036,0,0,0,0,0,86.0285746658828,0,0,0,0,0,2.74945243680585,0,0,0,0,0,-0.0370929666095603,0,0,0,0,0,0.0271175324818954,0.968497787558168,6.78918208039884,130.435544702695];
% params_init = [22.5550000000000,2639.21343385169,6954.15546240242,37.0996433482041,0,0,0,0,0,86.1433815863236,0,0,0,0,0,2.67697380297446,0,0,0,0,0,0.182791166501241,0,0,0,0,0,0.0254337093410858,0.686606591910504,6.09971263788264,130.690631220896];
% params_init = [22.555,2672.88280890923,6949.65158502979,37.4170137647429,0,0,0,0,0,86.0285742450053,0,0,0,0,0,2.74947011234005,0,0,0,0,0,-0.0371198364030973,0,0,0,0,0,0.0254558949658529,0.968498351398830,6.78918260312645,130.435543948130];
% params_init = [22.555,2672.88279572213,6949.65158646638,37.4170136246946,0,0,0,0,0,86.0287778158475,0,0,0,0,0,2.74736916230632,0,0,0,0,0,-0.0346183798671545,0,0,0,0,0,0.0254580754821122,0.968430956603284,6.78913884013286,130.435638681147];
% params_init = [22.555,2636.70610603933,6953.82558914833,37.1358081092509,0,0,0,0,0,86.1450070955185,0,0,0,0,0,2.67518976903739,0,0,0,0,0,0.190408996510545,0,0,0,0,0,0.0259302768075426,0.777320455870285,6.47592602341366,130.774166555759];
% params_init = [22.5550000000000,2637.33617409871,6953.87087622058,37.1403958142949,0,0,0,0,0,86.1479475125384,0,0,0,0,0,2.67836449970903,0,0,0,0,0,0.143633992903819,0,0,0,0,0,0.0236914330853278,1.03523886979435,7.87434341364638,130.849056986844];
% params_init = [22.5550000000000,2638.87967596836,6956.74590423197,37.2264289017128,0,0,0,0,0,86.5542220285932,0,0,0,0,0,2.68942313224628,0,0,0,0,0,0.0396970532290793,0,0,0,0,0,0.0198332211569660,0.542711028545775,5.88110414300863,130.339058033350];
%  params_init = [22.555,2638.87960550061,6956.74623655026,37.2264149734175,0,0,0,0,0,86.5553796191804,0,0,0,0,0,2.68923687171096,0,0,0,0,0,0.0385707185622832,0,0,0,0,0,0.0197503798713480,0.540748767571423,5.87629360011326,130.338171011153];
% params_init = [22.555,2638.33369709270,7299.74583312913,9.78588995212791,0,0,0,0,0,86.2779986081257,0,0,0,0,0,2.64054689502414,0,0,0,0,0,-0.121021067558611,0,0,0,0,0,0.00840619006753843,0.281416096516770,7.11745257258451,130.377125864504];
% params_init = [22.555,2656.52642414977,7298.23149507131,10.2141048204410,0,0,0,0,0,86.4362106997872,0,0,0,0,0,2.67962223740954,0,0,0,0,0,0.00242239049165781,0,0,0,0,0,0.0128218007844121,0.388169436431343,7.85061199420585,130.371937461416];
params_init = [22.5550000000000,2656.52642414977,7298.23149507131,10.2141048204410,0,0,0,0,0,86.4362106997872,0,0,0,0,0,2.67962223740954,0,0,0,0,0,0.00242239049165781,0,0,0,0,0,0.0128218007844121,0.388169436431343,7.85061199420585,130.371937461416];
lb = [22.555,0,0,0,0,0,0,0,0,20,0,0,0,0,0,-5,0,0,0,0,0,-10,0,0,0,0,0,0,0,0,0];
ub = [22.555,10000,10000,10000,0,0,0,0,0,300,0,0,0,0,0,5,0,0,0,0,0,10,0,0,0,0,0,100,100,100,1000];

% Optimization options
% options = optimset('Display', 'iter', 'TolFun', 1e-25, 'MaxIter', 2000000);
options = optimoptions('fmincon', 'Display', 'iter', 'TolFun', 1e-25, 'MaxIterations', 20000000, ...
                       'MaxFunctionEvaluations', 1e+06, 'OutputFcn', @outfun);

% Perform the optimization using fmincon
[params_fit, fval] = fmincon(@(params) combined_cost_function(params, T_Vm_sub, p_Vm_sub, Vm_sub, ...
                           T_Vm_melt, p_Vm_melt, Vm_melt, ...
                           T_Vm_highp, p_Vm_highp, Vm_highp, ...
                           T_cp_sub, p_cp_sub, cp_sub, ...
                           T_alpha_sub, p_alpha_sub, alpha_sub, ...
                           T_BetaT_sub, p_BetaT_sub, BetaT_sub, ...
                           T_BetaS_sub, p_BetaS_sub, BetaS_sub, ...
                           T_sub, p_sub,Year_sub, G_fluid_sub, V_fluid_sub, ...
                           T_melt, p_melt, G_fluid_melt, V_fluid_melt, T_H_sub, p_H_sub, delta_H_sub, H_fluid_sub, ...
                           T_H_melt, p_H_melt, delta_H_melt, H_fluid_melt), ...
                           params_init, [], [], [], [], lb, ub, [], ...
                           optimset('Display', 'iter', 'TolFun', 1e-25, 'MaxIter', 20000000, ...
                                    'OutputFcn', @(x, optimValues, state) outfun(x, optimValues, state, ...
                                    T_Vm_sub, p_Vm_sub, Vm_sub, ...
                                    T_Vm_melt, p_Vm_melt, Vm_melt, ...
                                    T_Vm_highp, p_Vm_highp, Vm_highp, ...
                                    T_cp_sub, p_cp_sub, cp_sub, ...
                                    T_alpha_sub, p_alpha_sub, alpha_sub, ...
                                    T_BetaT_sub, p_BetaT_sub, BetaT_sub, ...
                                    T_BetaS_sub, p_BetaS_sub, BetaS_sub, ...
                                    T_sub, p_sub, Year_sub,G_fluid_sub, V_fluid_sub, ...
                                    T_melt, p_melt, G_fluid_melt, V_fluid_melt, T_H_sub, p_H_sub, delta_H_sub, H_fluid_sub, ...
                                    T_H_melt, p_H_melt, delta_H_melt, H_fluid_melt)));




%need to: plot for each property
%for data verses fitted curve
%check the extrapolation behaviour
% Compute fitted values using the optimized parameters
if Lplot == 1
T_plot_calc = 1:400;
for ii = 1:length(T_plot_calc)
    if T_plot_calc(ii)<83.806
        p_plot_calc(ii) = psub(T_plot_calc(ii));
    else
        p_plot_calc(ii) = pmelt(T_plot_calc(ii));
    end

    fitted_props(ii,:) = computeThermoProps(T_plot_calc(ii), p_plot_calc(ii), params_fit);
end
% Example: Compare experimental and fitted molar volume for sublimation

for ii = 1:length(Vm_sub)
    p_Vm_sub(ii) = psub(T_Vm_sub(ii));
end

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
figure(1);
for ii = 1:xxcount_Vm_sub
        plot(T_Vm_sub(Nxstart_Vm_sub(ii):Nxend_Vm_sub(ii)),Vm_sub(Nxstart_Vm_sub(ii):Nxend_Vm_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_Vm_sub(Nxstart_Vm_sub(ii)),"-",Author_Vm_sub(Nxstart_Vm_sub(ii)))));
        hold on    
end
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
for ii = 1:xxcount_Vm_melt
        plot(T_Vm_melt(Nxstart_Vm_melt(ii):Nxend_Vm_melt(ii)),Vm_melt(Nxstart_Vm_melt(ii):Nxend_Vm_melt(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_Vm_melt(Nxstart_Vm_melt(ii)),"-",Author_Vm_melt(Nxstart_Vm_melt(ii)))));
        hold on    
end
hold on;
plot(T_plot_calc, fitted_props(:,1), 'r-', 'LineWidth', 1.5);  % Fitted values
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
    ylabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/cell volume')
        mkdir('../solid data fitting/figure/cell volume')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/cell volume/Vm.png')
    
figure(2)%deviation plot Vm_sub

    for ii = 1:length(T_Vm_sub)
        props_Vm_sub = computeThermoProps(T_Vm_sub(ii), p_Vm_sub(ii),  params_fit);
        Vm_sub_deviation(ii) = 100*(Vm_sub(ii) - props_Vm_sub(1)) ./ Vm_sub(ii);
    end

for ii = 1:xxcount_Vm_sub
        plot(T_Vm_sub(Nxstart_Vm_sub(ii):Nxend_Vm_sub(ii)),Vm_sub_deviation(Nxstart_Vm_sub(ii):Nxend_Vm_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_Vm_sub(Nxstart_Vm_sub(ii)),"-",Author_Vm_sub(Nxstart_Vm_sub(ii)))));
        hold on    
end
            c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
            xlim([0 90]);
        ytickformat('%.1f');    
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
    ylabel(['100\cdot',char(8201),'(\itV\rm_m_,_e_x_p',char(hex2dec('2212')),'\it V\rm_m_,_c_a_l_c) /\it V\rm_m_,_c_a_l_c'],'fontsize',fontsize,'linewidth',linewidth)  
     set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/cell volume')
        mkdir('../solid data fitting/figure/cell volume')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/cell volume/Vm sublimation deviation.png')   
    
figure(3)%deviation plot Vm_melt

    for ii = 1:length(T_Vm_melt)
        props_Vm_melt = computeThermoProps(T_Vm_melt(ii), p_Vm_melt(ii),  params_fit);
        Vm_melt_deviation(ii) = 100*(Vm_melt(ii) - props_Vm_melt(1)) ./ Vm_melt(ii);
    end

for ii = 1:xxcount_Vm_melt
        plot(T_Vm_melt(Nxstart_Vm_melt(ii):Nxend_Vm_melt(ii)),Vm_melt_deviation(Nxstart_Vm_melt(ii):Nxend_Vm_melt(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_Vm_melt(Nxstart_Vm_melt(ii)),"-",Author_Vm_melt(Nxstart_Vm_melt(ii)))));
        hold on    
end
        c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
        xlim([80 400]);
        ytickformat('%.1f');        
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
    ylabel(['100\cdot',char(8201),'(\itV\rm_m_,_e_x_p',char(hex2dec('2212')),'\it V\rm_m_,_c_a_l_c) /\it V\rm_m_,_c_a_l_c'],'fontsize',fontsize,'linewidth',linewidth)  
     set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/cell volume')
        mkdir('../solid data fitting/figure/cell volume')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/cell volume/Vm melting deviation.png') 
    
    

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

figure(4)
%  my_dashline=[0 0 0 0;4 1 1 2;2 0.5 0.5 1;6 6 6 6;3 2 3 2;1 1 1 1;7 1 7 1;13 3 13 3;1 4 1 4];
for ii = 1:a1_Vm_highp
    for jj=1:a2_Vm_highp
        if Nystart_Vm_highp(ii,jj)~=0
            ll1_Vm_highp(ii,jj)=string(strcat(Author_Vm_highp(Nystart_Vm_highp(ii,jj)),Year_Vm_highp(Nystart_Vm_highp(ii,jj)),' - ',{' '},char(string(T_Vm_highp(Nystart_Vm_highp(ii,jj)))),{' '},'K'));            
            plot(p_Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)),Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(jj,:),'linewidth',linewidth,'DisplayName',ll1_Vm_highp(ii,jj)); 
        hold on;

            T_plot_Vm_highp = mean(T_Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)));
            pp=1;
            p_plot_Vm_highp=[];
            Vm_p_plot_Vm_highp=[];  
            p_plot_Vm_highp_min=max(min(p_Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)))-2,0);
            p_plot_Vm_highp_max=max(p_Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)))+30;     
            for jjk=p_plot_Vm_highp_min:1:p_plot_Vm_highp_max
                p_plot_Vm_highp(pp,:)=jjk;
                props_Vm_highp(pp,:)=computeThermoProps(T_plot_Vm_highp, p_plot_Vm_highp(pp,:),  params_fit);
    %             benzeneprops(T_plot,p_plot_Vm_highp(pp,:));
                 Vm_p_plot_Vm_highp(pp,:)=props_Vm_highp(pp,1);
                pp=pp+1;
            end        
                plot(p_plot_Vm_highp,Vm_p_plot_Vm_highp,'color',mycolor(jj,:),'linewidth',linewidth);
        end
    end    
end  
        ytickformat('%.1f');  
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth); 
    ylabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,10,800,750])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/cell volume/high pressure/')
        mkdir('../solid data fitting/figure/cell volume/high pressure/')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/cell volume/high pressure/high pressure.png')  
    
    
    figure(5)
 for ii = 1:a1_Vm_highp
    for jj=1:a2_Vm_highp
        Vm_highp_dev = [];
        if Nystart_Vm_highp(ii,jj)~=0
        ll1_Vm_highp(ii,jj)=string(strcat(Author_Vm_highp(Nystart_Vm_highp(ii,jj)),Year_Vm_highp(Nystart_Vm_highp(ii,jj)),' - ',{' '},char(string(T_Vm_highp(Nystart_Vm_highp(ii,jj)))),{' '},'K'));            
        for  kkk = 1:length(p_Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)))
            props_Vm_highp = computeThermoProps(T_Vm_highp(Nystart_Vm_highp(ii,jj)+kkk-1), p_Vm_highp(Nystart_Vm_highp(ii,jj)+kkk-1),  params_fit);
            Vm_highp_dev(kkk,:) = 100*(Vm_highp(Nystart_Vm_highp(ii,jj)+kkk-1) - props_Vm_highp(1)) ./ Vm_highp(Nystart_Vm_highp(ii,jj)+kkk-1);
        end
        plot(p_Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)),Vm_highp_dev,...
        mymarker{ii},'markersize',markersize,'color',mycolor(jj,:),'linewidth',linewidth,'DisplayName',ll1_Vm_highp(ii,jj)); 
    hold on;
        end
    end    
 end
        c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
    ylabel(['100\cdot',char(8201),'(\itV\rm_m_,_e_x_p',char(hex2dec('2212')),'\it V\rm_m_,_c_a_l_c) /\it V\rm_m_,_c_a_l_c'],'fontsize',fontsize,'linewidth',linewidth)  
     set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/cell volume')
        mkdir('../solid data fitting/figure/cell volume')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/cell volume/high pressure/high pressure_dev.png')  


for ii = 1:length(cp_sub)
    p_cp_sub(ii) = psub(T_cp_sub(ii));
end

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
figure(6);
for ii = 1:xxcount_cp_sub
        plot(T_cp_sub(Nxstart_cp_sub(ii):Nxend_cp_sub(ii)),cp_sub(Nxstart_cp_sub(ii):Nxend_cp_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_cp_sub(Nxstart_cp_sub(ii)),"-",Author_cp_sub(Nxstart_cp_sub(ii)))));
        hold on    
end

plot(T_plot_calc, fitted_props(:,5), 'r-', 'LineWidth', 1.5);  % Fitted values
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
                    xlim([0 85]);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
    ylabel('\it c_p\rm / (J K^-^1 mol^-^1)','fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/heat capacity')
        mkdir('../solid data fitting/figure/heat capacity')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/heat capacity/cp.png')
    
figure(7)%deviation plot cp_sub

    for ii = 1:length(T_cp_sub)
        props_cp_sub = computeThermoProps(T_cp_sub(ii), p_cp_sub(ii),  params_fit);
        cp_sub_deviation(ii) = 100*(cp_sub(ii) - props_cp_sub(5)) ./ cp_sub(ii);
    end

for ii = 1:xxcount_cp_sub
        plot(T_cp_sub(Nxstart_cp_sub(ii):Nxend_cp_sub(ii)),cp_sub_deviation(Nxstart_cp_sub(ii):Nxend_cp_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_cp_sub(Nxstart_cp_sub(ii)),"-",Author_cp_sub(Nxstart_cp_sub(ii)))));
        hold on    
end
            c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
            xlim([0 90]);
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
    ylabel(['100\cdot',char(8201),'(\itc_p\rm_,_e_x_p',char(hex2dec('2212')),'\it c_p\rm_,_c_a_l_c) /\it c_p\rm_,_c_a_l_c'],'fontsize',fontsize,'linewidth',linewidth)  
     set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/heat capacity')
        mkdir('../solid data fitting/figure/heat capacity')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/heat capacity/cp sublimation deviation.png')   
    

for ii = 1:length(alpha_sub)
    p_alpha_sub(ii) = psub(T_alpha_sub(ii));
end

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
figure(8);
for ii = 1:xxcount_alpha_sub
        plot(T_alpha_sub(Nxstart_alpha_sub(ii):Nxend_alpha_sub(ii)),alpha_sub(Nxstart_alpha_sub(ii):Nxend_alpha_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_alpha_sub(Nxstart_alpha_sub(ii)),"-",Author_alpha_sub(Nxstart_alpha_sub(ii)))));
        hold on    
end

plot(T_plot_calc, fitted_props(:,4), 'r-', 'LineWidth', 1.5);  % Fitted values
xlim([0 84]);
        ytickformat('%.1f'); 
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
    ylabel('\it \alpha \rm\cdot10^4\rm / K^-^1','fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/thermal expansivity')
        mkdir('../solid data fitting/figure/thermal expansivity')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/thermal expansivity/alpha.png')
    
figure(9)%deviation plot alpha_sub

    for ii = 1:length(T_alpha_sub)
        props_alpha_sub = computeThermoProps(T_alpha_sub(ii), p_alpha_sub(ii),  params_fit);
        alpha_sub_deviation(ii) = 100*(alpha_sub(ii) - props_alpha_sub(4)) ./ alpha_sub(ii);
    end

for ii = 1:xxcount_alpha_sub
        plot(T_alpha_sub(Nxstart_alpha_sub(ii):Nxend_alpha_sub(ii)),alpha_sub_deviation(Nxstart_alpha_sub(ii):Nxend_alpha_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_alpha_sub(Nxstart_alpha_sub(ii)),"-",Author_alpha_sub(Nxstart_alpha_sub(ii)))));
        hold on    
end
            c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
            xlim([0 90]);
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
    ylabel(['100\cdot',char(8201),'(\it\alpha\rm_e_x_p',char(hex2dec('2212')),'\it \alpha\rm_c_a_l_c) /\it \alpha\rm_c_a_l_c'],'fontsize',fontsize,'linewidth',linewidth)  
     set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/thermal expansivity')
        mkdir('../solid data fitting/figure/thermal expansivity')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/thermal expansivity/alpha sublimation deviation.png')       

for ii = 1:length(BetaT_sub)
    p_BetaT_sub(ii) = psub(T_BetaT_sub(ii));
end

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
figure(10);
for ii = 1:xxcount_BetaT_sub
        plot(T_BetaT_sub(Nxstart_BetaT_sub(ii):Nxend_BetaT_sub(ii)),BetaT_sub(Nxstart_BetaT_sub(ii):Nxend_BetaT_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_BetaT_sub(Nxstart_BetaT_sub(ii)),"-",Author_BetaT_sub(Nxstart_BetaT_sub(ii)))));
        hold on    
end

plot(T_plot_calc, fitted_props(:,2), 'r-', 'LineWidth', 1.5);  % Fitted values
xlim([0 85]);
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
ylabel('\it \beta_T \rm\cdot10^4\rm / MPa^-^1','fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/bulk modulus')
        mkdir('../solid data fitting/figure/bulk modulus')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/bulk modulus/BetaT.png')
    
figure(11)%deviation plot BetaT_sub

    for ii = 1:length(T_BetaT_sub)
        props_BetaT_sub = computeThermoProps(T_BetaT_sub(ii), p_BetaT_sub(ii),  params_fit);
        BetaT_sub_deviation(ii) = 100*(BetaT_sub(ii) - props_BetaT_sub(2)) ./ BetaT_sub(ii);
    end

for ii = 1:xxcount_BetaT_sub
        plot(T_BetaT_sub(Nxstart_BetaT_sub(ii):Nxend_BetaT_sub(ii)),BetaT_sub_deviation(Nxstart_BetaT_sub(ii):Nxend_BetaT_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_BetaT_sub(Nxstart_BetaT_sub(ii)),"-",Author_BetaT_sub(Nxstart_BetaT_sub(ii)))));
        hold on    
end
            c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
            xlim([0 90]);
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
    ylabel(['100\cdot',char(8201),'(\it\beta_T\rm_,_e_x_p',char(hex2dec('2212')),'\it \beta_T\rm_,_c_a_l_c) /\it \beta_T\rm_,_c_a_l_c'],'fontsize',fontsize,'linewidth',linewidth)  
     set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/bulk modulus')
        mkdir('../solid data fitting/figure/bulk modulus')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/bulk modulus/BetaT sublimation deviation.png') 

for ii = 1:length(BetaS_sub)
    p_BetaS_sub(ii) = psub(T_BetaS_sub(ii));
end

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
figure(12);
for ii = 1:xxcount_BetaS_sub
        plot(T_BetaS_sub(Nxstart_BetaS_sub(ii):Nxend_BetaS_sub(ii)),BetaS_sub(Nxstart_BetaS_sub(ii):Nxend_BetaS_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_BetaS_sub(Nxstart_BetaS_sub(ii)),"-",Author_BetaS_sub(Nxstart_BetaS_sub(ii)))));
        hold on    
end

plot(T_plot_calc, fitted_props(:,3), 'r-', 'LineWidth', 1.5);  % Fitted values
xlim([0 85]);
        ytickformat('%.1f'); 
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
ylabel('\it \beta_S \rm\cdot10^4\rm / MPa^-^1','fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/bulk modulus')
        mkdir('../solid data fitting/figure/bulk modulus')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/bulk modulus/BetaS.png')
    
figure(13)%deviation plot BetaS_sub

    for ii = 1:length(T_BetaS_sub)
        props_BetaS_sub = computeThermoProps(T_BetaS_sub(ii), p_BetaS_sub(ii),  params_fit);
        BetaS_sub_deviation(ii) = 100*(BetaS_sub(ii) - props_BetaS_sub(3)) ./ BetaS_sub(ii);
    end

for ii = 1:xxcount_BetaS_sub
        plot(T_BetaS_sub(Nxstart_BetaS_sub(ii):Nxend_BetaS_sub(ii)),BetaS_sub_deviation(Nxstart_BetaS_sub(ii):Nxend_BetaS_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_BetaS_sub(Nxstart_BetaS_sub(ii)),"-",Author_BetaS_sub(Nxstart_BetaS_sub(ii)))));
        hold on    
end
            c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
            xlim([0 90]);
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
    ylabel(['100\cdot',char(8201),'(\it\beta_S\rm_,_e_x_p',char(hex2dec('2212')),'\it \beta_S\rm_,_c_a_l_c) /\it \beta_S\rm_,_c_a_l_c'],'fontsize',fontsize,'linewidth',linewidth)  
     set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/bulk modulus')
        mkdir('../solid data fitting/figure/bulk modulus')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/bulk modulus/BetaS sublimation deviation.png') 
   
% 
% for ii = 1:length(BetaS_sub)
%     p_sub(ii) = psub(T_sub(ii));
% end

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
figure(14);
for ii = 1:xxcount_sub
        plot(T_sub(Nxstart_sub(ii):Nxend_sub(ii)),p_sub(Nxstart_sub(ii):Nxend_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_sub(Nxstart_sub(ii)),"-",Author_sub(Nxstart_sub(ii)))));
        hold on    
end

St_REFPROP = 131.1500;%J/mol/K
Ht_REFPROP = 1689.10;%J/molProps(K1038,K1039,12,$B$7:$B$36)+K1038*O1041
Tt = 83.8058; pt = 0.068891;%unit: K and MPa
deltaS_triple =  params_fit(31) - St_REFPROP;
props_Triple = computeThermoProps(Tt,pt, params_fit);
Ht_fitted = props_Triple(12) + Tt*params_fit(31);
deltaH_triple = Ht_fitted - Ht_REFPROP;

for ii = 1:length(T_sub)
    
    fitted_props_sub = computeThermoProps(T_sub(ii), p_sub(ii)/10^6, params_fit);
    delta_G_sub(ii,:) = G_fluid_sub(ii) - fitted_props_sub(12) + deltaH_triple - T_sub(ii) * deltaS_triple;
    p_fitted_sub(ii) = p_sub(ii) - delta_G_sub(ii)/(V_fluid_sub(ii) - fitted_props_sub(1)) * 10^6; % unit Pa
    p_sub_deviation(ii) = 100 * (p_sub(ii) - p_fitted_sub(ii)) / p_sub(ii);

end

% Sort T_sub and p_fitted_sub in ascending order
[T_sub_sorted, sortIdx] = sort(T_sub);
p_fitted_sub_sorted = p_fitted_sub(sortIdx);

% Plot the sorted values
plot(T_sub_sorted, p_fitted_sub_sorted, 'r-', 'LineWidth', 1.5); % Fitted values
set(gca, 'YScale', 'log');  % Set y-axis to logarithmic scale
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
ylabel('\it p \rm / Pa','fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/sublimation')
        mkdir('../solid data fitting/figure/sublimation')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/sublimation/sublimation.png')
    
figure(15)%deviation plot sublimation

for ii = 1:xxcount_sub
        plot(T_sub(Nxstart_sub(ii):Nxend_sub(ii)),p_sub_deviation(Nxstart_sub(ii):Nxend_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_sub(Nxstart_sub(ii)),"-",Author_sub(Nxstart_sub(ii)))));
        hold on    
end
            c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
            xlim([0 90]);
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
    ylabel(['100\cdot',char(8201),'(\itp\rm_,_e_x_p',char(hex2dec('2212')),'\it p\rm_,_c_a_l_c) /\it p\rm_,_c_a_l_c'],'fontsize',fontsize,'linewidth',linewidth)  
     set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/sublimation')
        mkdir('../solid data fitting/figure/sublimation')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/sublimation/Sublimation deviation.png')  
    
% 
% for ii = 1:length(BetaS_sub)
%     p_sub(ii) = psub(T_sub(ii));
% end

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
figure(16);
for ii = 1:xxcount_melt
        plot(T_melt(Nxstart_melt(ii):Nxend_melt(ii)),p_melt(Nxstart_melt(ii):Nxend_melt(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_melt(Nxstart_melt(ii)),"-",Author_melt(Nxstart_melt(ii)))));
        hold on    
end

St_REFPROP = 131.1500;%J/mol/K
Ht_REFPROP = 1689.10;%J/molProps(K1038,K1039,12,$B$7:$B$36)+K1038*O1041
Tt = 83.8058; pt = 0.068891;%unit: K and MPa
deltaS_triple =  params_fit(31) - St_REFPROP;
props_Triple = computeThermoProps(Tt,pt, params_fit);
Ht_fitted = props_Triple(12) + Tt*params_fit(31);
deltaH_triple = Ht_fitted - Ht_REFPROP;

for ii = 1:length(T_melt)
    
    fitted_props_melt = computeThermoProps(T_melt(ii), p_melt(ii), params_fit);
    delta_G_melt(ii,:) = G_fluid_melt(ii) - fitted_props_melt(12) + deltaH_triple - T_melt(ii) * deltaS_triple;
    p_fitted_melt(ii) = p_melt(ii) - delta_G_melt(ii)/(V_fluid_melt(ii) - fitted_props_melt(1)) ; % unit MPa
    p_melt_deviation(ii) = 100 * (p_melt(ii) - p_fitted_melt(ii)) / p_melt(ii);

end

% Sort T_melt and p_fitted_melt in ascending order
[T_melt_sorted, sortIdx] = sort(T_melt);
p_fitted_melt_sorted = p_fitted_melt(sortIdx);

% Plot the sorted values
plot(T_melt_sorted, p_fitted_melt_sorted, 'r-', 'LineWidth', 1.5); % Fitted values

        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
ylabel('\it p \rm / MPa','fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/melting')
        mkdir('../solid data fitting/figure/melting')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/melting/melting.png')
    
figure(17)%deviation plot melting

for ii = 1:xxcount_melt
        plot(T_melt(Nxstart_melt(ii):Nxend_melt(ii)),p_melt_deviation(Nxstart_melt(ii):Nxend_melt(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_melt(Nxstart_melt(ii)),"-",Author_melt(Nxstart_melt(ii)))));
        hold on    
end
            c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
            xlim([80 600]);
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
    ylabel(['100\cdot',char(8201),'(\itp\rm,_e_x_p',char(hex2dec('2212')),'\it p\rm_,_c_a_l_c) /\it p\rm_,_c_a_l_c'],'fontsize',fontsize,'linewidth',linewidth)  
     set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/melting')
        mkdir('../solid data fitting/figure/melting')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/melting/melting deviation.png')       
    

for ii = 1:length(BetaS_sub)
    p_sub(ii) = psub(T_sub(ii));
end

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
figure(18);

ytickformat('%.2f'); 
St_REFPROP = 131.1500;%J/mol/K
Ht_REFPROP = 1689.10;%J/molProps(K1038,K1039,12,$B$7:$B$36)+K1038*O1041
Tt = 83.8058; pt = 0.068891;%unit: K and MPa
deltaS_triple =  params_fit(31) - St_REFPROP;
props_Triple = computeThermoProps(Tt,pt, params_fit);
Ht_fitted = props_Triple(12) + Tt*params_fit(31);
deltaH_triple = Ht_fitted - Ht_REFPROP;

for ii = 1:length(T_H_sub)
        H_solid_sub(ii) = H_fluid_sub(ii)-1000*delta_H_sub(ii);
        fitted_props_H_sub = computeThermoProps(T_H_sub(ii), p_H_sub(ii), params_fit);

        H_solid_sub_fitted(ii) = fitted_props_H_sub(11) - deltaH_triple;
        H_solid_sub_deviation(ii) = 100*(H_solid_sub(ii) - H_solid_sub_fitted(ii)) / H_solid_sub_fitted(ii) ;   

end
for ii = 1:xxcount_H_sub
        plot(T_H_sub(Nxstart_H_sub(ii):Nxend_H_sub(ii)),H_solid_sub(Nxstart_H_sub(ii):Nxend_H_sub(ii))/1000,...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_H_sub(Nxstart_H_sub(ii)),"-",Author_H_sub(Nxstart_H_sub(ii)))));
        hold on    
end
% Sort T_H_sub and p_fitted_H_sub in ascending order
[T_H_sub_sorted, sortIdx] = sort(T_H_sub);
H_solid_sub_fitted_sorted = H_solid_sub_fitted(sortIdx);

% Plot the sorted values
plot(T_H_sub_sorted, H_solid_sub_fitted_sorted/1000, 'r-', 'LineWidth', 1.5); % Fitted values

        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
ylabel(['\it \DeltaH \rm/ (kJ mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/enthalpy')
        mkdir('../solid data fitting/figure/enthalpy')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/enthalpy/H_sub.png')
    
figure(19)%deviation plot enthalpy

for ii = 1:xxcount_H_sub
        plot(T_H_sub(Nxstart_H_sub(ii):Nxend_H_sub(ii)),H_solid_sub_deviation(Nxstart_H_sub(ii):Nxend_H_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_H_sub(Nxstart_H_sub(ii)),"-",Author_H_sub(Nxstart_H_sub(ii)))));
        hold on    
end
            c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
%             xlim([0 80]);
ytickformat('%.1f'); 
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
    ylabel(['100\cdot',char(8201),'(\it\DeltaH\rm_e_x_p',char(hex2dec('2212')),'\it \DeltaH\rm_c_a_l_c) /\it\DeltaH\rm_c_a_l_c'],'fontsize',fontsize,'linewidth',linewidth)  
     set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/enthalpy')
        mkdir('../solid data fitting/figure/enthalpy')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/enthalpy/H_sub_dev.png')       
 % 
% for ii = 1:length(BetaS_sub)
%     p_sub(ii) = psub(T_sub(ii));
% end

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
figure(20);


St_REFPROP = 131.1500;%J/mol/K
Ht_REFPROP = 1689.10;%J/molProps(K1038,K1039,12,$B$7:$B$36)+K1038*O1041
Tt = 83.8058; pt = 0.068891;%unit: K and MPa
deltaS_triple = params_fit(31) - St_REFPROP;
props_Triple = computeThermoProps(Tt,pt, params_fit);
Ht_fitted = props_Triple(12) + Tt*params_fit(31);
deltaH_triple = Ht_fitted - Ht_REFPROP;

for ii = 1:length(T_H_melt)
        H_solid_melt(ii) = H_fluid_melt(ii)-1000*delta_H_melt(ii);
        fitted_props_H_melt = computeThermoProps(T_H_melt(ii), p_H_melt(ii), params_fit);

        H_solid_melt_fitted(ii) = fitted_props_H_melt(11) - deltaH_triple;
        H_solid_melt_deviation(ii) = 100*(H_solid_melt(ii) - H_solid_melt_fitted(ii)) / H_solid_melt_fitted(ii) ;   

end
for ii = 1:xxcount_H_melt
        plot(T_H_melt(Nxstart_H_melt(ii):Nxend_H_melt(ii)),H_solid_melt(Nxstart_H_melt(ii):Nxend_H_melt(ii))/1000,...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_H_melt(Nxstart_H_melt(ii)),"-",Author_H_melt(Nxstart_H_melt(ii)))));
        hold on    
end
% Sort T_H_melt and p_fitted_H_melt in ascending order
[T_H_melt_sorted, sortIdx] = sort(T_H_melt);
H_solid_melt_fitted_sorted = H_solid_melt_fitted(sortIdx);

% Plot the sorted values
plot(T_H_melt_sorted, H_solid_melt_fitted_sorted/1000, 'r-', 'LineWidth', 1.5); % Fitted values

        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
ylabel(['\it \DeltaH \rm/ (kJ mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/enthalpy')
        mkdir('../solid data fitting/figure/enthalpy')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/enthalpy/H_melt.png')
    
figure(21)%deviation plot enthalpy

for ii = 1:xxcount_H_melt
        plot(T_H_melt(Nxstart_H_melt(ii):Nxend_H_melt(ii)),H_solid_melt_deviation(Nxstart_H_melt(ii):Nxend_H_melt(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_H_melt(Nxstart_H_melt(ii)),"-",Author_H_melt(Nxstart_H_melt(ii)))));
        hold on    
end
            c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
            xlim([80 600]);
            ytickformat('%.1f'); 
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
    ylabel(['100\cdot',char(8201),'(\it\DeltaH\rm_e_x_p',char(hex2dec('2212')),'\it \DeltaH\rm_c_a_l_c) /\it\DeltaH\rm_c_a_l_c'],'fontsize',fontsize,'linewidth',linewidth)  
     set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/enthalpy')
        mkdir('../solid data fitting/figure/enthalpy')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/enthalpy/H_melt_dev.png')
    
    figure(22)
    T1 = 0.0001;
    V_T1 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T1)
        props_T1 = computeThermoPropsTV(T1, V_T1(pp), params_fit);
        if props_T1(1)>0
            p_T1_plot(kk) = props_T1(1);
            V_T1_plot(kk) = V_T1(pp);
            kk = kk + 1;
        end
    end
 plot(V_T1_plot,p_T1_plot,...
            'color',mycolor(1,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 0 K");
  hold on  
     T2 = 50;
    V_T2 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T2)
        props_T2 = computeThermoPropsTV(T2, V_T2(pp), params_fit);
        if props_T2(1)>0
            p_T2_plot(kk) = props_T2(1);
            V_T2_plot(kk) = V_T2(pp);
            kk = kk + 1;
        end
    end
 plot(V_T2_plot,p_T2_plot,...
            '--','color',mycolor(2,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 50 K"); 
     hold on   
     
     T3 = 83.806;
    V_T3 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T3)
        props_T3 = computeThermoPropsTV(T3, V_T3(pp), params_fit);
        if props_T3(1)>0
            p_T3_plot(kk) = props_T3(1);
            V_T3_plot(kk) = V_T3(pp);
            kk = kk + 1;
        end
    end
 plot(V_T3_plot,p_T3_plot,...
            '-.','color',mycolor(3,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 83.806 K");  
        hold on
        
      T4 = 150;
    V_T4 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T4)
        props_T4 = computeThermoPropsTV(T4, V_T4(pp), params_fit);
        props_T4_check = computeThermoProps(T4, pmelt(T4), params_fit);          
        if props_T4(1)>0&&V_T4(pp)<props_T4_check(1)
            p_T4_plot(kk) = props_T4(1);
            V_T4_plot(kk) = V_T4(pp);
            kk = kk + 1;
        end
    end
 plot(V_T4_plot,p_T4_plot,...
            ':','color',mycolor(4,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 150 K");  
     hold on   
     T5 = 300;
    V_T5 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T5)
        props_T5 = computeThermoPropsTV(T5, V_T5(pp), params_fit);  
        props_T5_check = computeThermoProps(T5, pmelt(T5), params_fit);  
        if props_T5(1)>0&&V_T5(pp)<props_T5_check(1)
            p_T5_plot(kk) = props_T5(1);
            
            V_T5_plot(kk) = V_T5(pp);
            kk = kk + 1;
        end
    end

 dashline(V_T5_plot,p_T5_plot,4,1,2,1,'color',mycolor(5,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 300 K");
 
    
 T6 = [0.0001	2	4	6	8	10	12	14	16	18	20	22	24	26	28	30	32	34	36	38	40	42	44	46	48	50	52	54	56	58	60	62	64	66	68	70	72	74	76	78	80	82	83.806];
  for ii = 1:length(T6)

          p6(ii) = psub(T6(ii));
      props_T6 = computeThermoProps(T6(ii), p6(ii), params_fit);
      V_T6_plot(ii) = props_T6(1);
  end
 dashline(V_T6_plot,p6,6,6,6,6,'color',mycolor(6,:),'linewidth',linewidth,'DisplayName',"sublimation");   
 T7 = [83.806	84	86	88	90	92	94	96	98	100	102	104	106	108	110	112	114	116	118	120	122	124	126	128	130	132	134	136	138	140	142	144	146	148	150	152	154	156	158	160	162	164	166	168	170	172	174	176	178	180	182	184	186	188	190	192	194	196	198	200	202	204	206	208	210	212	214	216	218	220	222	224	226	228	230	232	234	236	238	240	242	244	246	248	250	252	254	256	258	260	262	264	266	268	270	272	274	276	278	280	282	284	286	288	290	292	294	296	298	300];
   for ii = 1:length(T7)

          p7(ii) = pmelt(T7(ii));
      props_T7 = computeThermoProps(T7(ii), p7(ii), params_fit);
      V_T7_plot(ii) = props_T7(1);
  end
 dashline(V_T7_plot,p7,7,1,7,1,'color',mycolor(7,:),'linewidth',linewidth,'DisplayName',"melting");
         legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);

    xlabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    ylabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth); 
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/evaluation')
        mkdir('../solid data fitting/figure/evaluation')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/evaluation/pV.png')
    
    figure(23)
    T1 = 0.0001;
    V_T1 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T1)
        props_T1 = computeThermoPropsTV(T1, V_T1(pp), params_fit);
        if props_T1(1)>0
            BetaS_T1_plot(kk) = props_T1(3);
            V_T1_plot(kk) = V_T1(pp);
            kk = kk + 1;
        end
    end
 plot(V_T1_plot,BetaS_T1_plot*10^4,...
            'color',mycolor(1,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 0 K");
  hold on  
     T2 = 50;
    V_T2 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T2)
        props_T2 = computeThermoPropsTV(T2, V_T2(pp), params_fit);
        if props_T2(1)>0
            BetaS_T2_plot(kk) = props_T2(3);
            V_T2_plot(kk) = V_T2(pp);
            kk = kk + 1;
        end
    end
 plot(V_T2_plot,BetaS_T2_plot*10^4,...
            '--','color',mycolor(2,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 50 K"); 
     hold on   
     
     T3 = 83.806;
    V_T3 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T3)
        props_T3 = computeThermoPropsTV(T3, V_T3(pp), params_fit);
        if props_T3(1)>0
            BetaS_T3_plot(kk) = props_T3(3);
            V_T3_plot(kk) = V_T3(pp);
            kk = kk + 1;
        end
    end
 plot(V_T3_plot,BetaS_T3_plot*10^4,...
            '-.','color',mycolor(3,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 83.806 K");  
        hold on
        
      T4 = 150;
    V_T4 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T4)
        props_T4 = computeThermoPropsTV(T4, V_T4(pp), params_fit);
        props_T4_check = computeThermoProps(T4, pmelt(T4), params_fit);          
        if props_T4(1)>0&&V_T4(pp)<props_T4_check(1)
            BetaS_T4_plot(kk) = props_T4(3);
            V_T4_plot(kk) = V_T4(pp);
            kk = kk + 1;
        end
    end
 plot(V_T4_plot,BetaS_T4_plot*10^4,...
            ':','color',mycolor(4,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 150 K");  
     hold on   
     T5 = 300;
    V_T5 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T5)
        props_T5 = computeThermoPropsTV(T5, V_T5(pp), params_fit);  
        props_T5_check = computeThermoProps(T5, pmelt(T5), params_fit);  
        if props_T5(1)>0&&V_T5(pp)<props_T5_check(1)
            BetaS_T5_plot(kk) = props_T5(3);
            
            V_T5_plot(kk) = V_T5(pp);
            kk = kk + 1;
        end
    end

 dashline(V_T5_plot,BetaS_T5_plot*10^4,4,1,2,1,'color',mycolor(5,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 300 K");
 
    V_T7_plot = [];
 T6 = [0.0001	2	4	6	8	10	12	14	16	18	20	22	24	26	28	30	32	34	36	38	40	42	44	46	48	50	52	54	56	58	60	62	64	66	68	70	72	74	76	78	80	82	83.806];
  for ii = 1:length(T6)

          p6(ii) = psub(T6(ii));
      props_T6 = computeThermoProps(T6(ii), p6(ii), params_fit);
      BetaS_T6_plot(ii) = props_T6(3);
      V_T6_plot(ii) = props_T6(1);
  end
 dashline(V_T6_plot,BetaS_T6_plot*10^4,6,6,6,6,'color',mycolor(6,:),'linewidth',linewidth,'DisplayName',"sublimation");   
 T7 = [83.806	84	86	88	90	92	94	96	98	100	102	104	106	108	110	112	114	116	118	120	122	124	126	128	130	132	134	136	138	140	142	144	146	148	150	152	154	156	158	160	162	164	166	168	170	172	174	176	178	180	182	184	186	188	190	192	194	196	198	200	202	204	206	208	210	212	214	216	218	220	222	224	226	228	230	232	234	236	238	240	242	244	246	248	250	252	254	256	258	260	262	264	266	268	270	272	274	276	278	280	282	284	286	288	290	292	294	296	298	300];
   for ii = 1:length(T7)

          p7(ii) = pmelt(T7(ii));
      props_T7 = computeThermoProps(T7(ii), p7(ii), params_fit);
       BetaS_T7_plot(ii) = props_T7(3);
      V_T7_plot(ii) = props_T7(1);
  end
 dashline(V_T7_plot,BetaS_T7_plot*10^4,7,1,7,1,'color',mycolor(7,:),'linewidth',linewidth,'DisplayName',"melting");
         legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);

    xlabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
ylabel(['\it \beta_S \rm\cdot10^4\rm / MPa^',char(hex2dec('2212')),'^1'],'fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/evaluation')
        mkdir('../solid data fitting/figure/evaluation')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/evaluation/KS.png')  
    
    figure(24)   
     T1 = 0.0001;
    V_T1 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T1)
        props_T1 = computeThermoPropsTV(T1, V_T1(pp), params_fit);
        if props_T1(1)>0
            alpha_T1_plot(kk) = props_T1(4);
            V_T1_plot(kk) = V_T1(pp);
            kk = kk + 1;
        end
    end
 plot(V_T1_plot,alpha_T1_plot*10^4,...
            'color',mycolor(1,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 0 K");
  hold on  
     T2 = 50;
    V_T2 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T2)
        props_T2 = computeThermoPropsTV(T2, V_T2(pp), params_fit);
        if props_T2(1)>0
            alpha_T2_plot(kk) = props_T2(4);
            V_T2_plot(kk) = V_T2(pp);
            kk = kk + 1;
        end
    end
 plot(V_T2_plot,alpha_T2_plot*10^4,...
            '--','color',mycolor(2,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 50 K"); 
     hold on   
     
     T3 = 83.806;
    V_T3 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T3)
        props_T3 = computeThermoPropsTV(T3, V_T3(pp), params_fit);
        if props_T3(1)>0
            alpha_T3_plot(kk) = props_T3(4);
            V_T3_plot(kk) = V_T3(pp);
            kk = kk + 1;
        end
    end
 plot(V_T3_plot,alpha_T3_plot*10^4,...
            '-.','color',mycolor(3,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 83.806 K");  
        hold on
        
      T4 = 150;
    V_T4 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T4)
        props_T4 = computeThermoPropsTV(T4, V_T4(pp), params_fit);
        props_T4_check = computeThermoProps(T4, pmelt(T4), params_fit);          
        if props_T4(1)>0&&V_T4(pp)<props_T4_check(1)
            alpha_T4_plot(kk) = props_T4(4);
            V_T4_plot(kk) = V_T4(pp);
            kk = kk + 1;
        end
    end
 plot(V_T4_plot,alpha_T4_plot*10^4,...
            ':','color',mycolor(4,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 150 K");  
     hold on   
     T5 = 300;
    V_T5 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T5)
        props_T5 = computeThermoPropsTV(T5, V_T5(pp), params_fit);  
        props_T5_check = computeThermoProps(T5, pmelt(T5), params_fit);  
        if props_T5(1)>0&&V_T5(pp)<props_T5_check(1)
            alpha_T5_plot(kk) = props_T5(4);
            
            V_T5_plot(kk) = V_T5(pp);
            kk = kk + 1;
        end
    end

 dashline(V_T5_plot,alpha_T5_plot*10^4,4,1,2,1,'color',mycolor(5,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 300 K");
 
    
 T6 = [0.0001	2	4	6	8	10	12	14	16	18	20	22	24	26	28	30	32	34	36	38	40	42	44	46	48	50	52	54	56	58	60	62	64	66	68	70	72	74	76	78	80	82	83.806];
  for ii = 1:length(T6)

          p6(ii) = psub(T6(ii));
      props_T6 = computeThermoProps(T6(ii), p6(ii), params_fit);
      alpha_T6_plot(ii) = props_T6(4);
      V_T6_plot(ii) = props_T6(1);
  end
 dashline(V_T6_plot,alpha_T6_plot*10^4,6,6,6,6,'color',mycolor(6,:),'linewidth',linewidth,'DisplayName',"sublimation");   
 T7 = [83.806	84	86	88	90	92	94	96	98	100	102	104	106	108	110	112	114	116	118	120	122	124	126	128	130	132	134	136	138	140	142	144	146	148	150	152	154	156	158	160	162	164	166	168	170	172	174	176	178	180	182	184	186	188	190	192	194	196	198	200	202	204	206	208	210	212	214	216	218	220	222	224	226	228	230	232	234	236	238	240	242	244	246	248	250	252	254	256	258	260	262	264	266	268	270	272	274	276	278	280	282	284	286	288	290	292	294	296	298	300];
   for ii = 1:length(T7)

          p7(ii) = pmelt(T7(ii));
      props_T7 = computeThermoProps(T7(ii), p7(ii), params_fit);
       alpha_T7_plot(ii) = props_T7(4);
      V_T7_plot(ii) = props_T7(1);
  end
 dashline(V_T7_plot,alpha_T7_plot*10^4,7,1,7,1,'color',mycolor(7,:),'linewidth',linewidth,'DisplayName',"melting");
         legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);

    xlabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
ylabel(['\it \alpha \rm\cdot10^4\rm / K^',char(hex2dec('2212')),'^1'],'fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/evaluation')
        mkdir('../solid data fitting/figure/evaluation')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/evaluation/alpha.png')   
    
    
     figure(25)
    T1 = 0.0001;
    V_T1 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T1)
        props_T1 = computeThermoPropsTV(T1, V_T1(pp), params_fit);
        if props_T1(1)>0
            Beta_T1_plot(kk) = props_T1(4)/props_T1(2);
            V_T1_plot(kk) = V_T1(pp);
            kk = kk + 1;
        end
    end
 plot(V_T1_plot,Beta_T1_plot,...
            'color',mycolor(1,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 0 K");
  hold on  
     T2 = 50;
    V_T2 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T2)
        props_T2 = computeThermoPropsTV(T2, V_T2(pp), params_fit);
        if props_T2(1)>0
            Beta_T2_plot(kk) = props_T2(4)/props_T2(2);
            V_T2_plot(kk) = V_T2(pp);
            kk = kk + 1;
        end
    end
 plot(V_T2_plot,Beta_T2_plot,...
            '--','color',mycolor(2,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 50 K"); 
     hold on   
     
     T3 = 83.806;
    V_T3 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T3)
        props_T3 = computeThermoPropsTV(T3, V_T3(pp), params_fit);
        if props_T3(1)>0
            Beta_T3_plot(kk) = props_T3(4)/props_T3(2);
            V_T3_plot(kk) = V_T3(pp);
            kk = kk + 1;
        end
    end
 plot(V_T3_plot,Beta_T3_plot,...
            '-.','color',mycolor(3,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 83.806 K");  
        hold on
        
      T4 = 150;
    V_T4 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T4)
        props_T4 = computeThermoPropsTV(T4, V_T4(pp), params_fit);
        props_T4_check = computeThermoProps(T4, pmelt(T4), params_fit);          
        if props_T4(1)>0&&V_T4(pp)<props_T4_check(1)
            Beta_T4_plot(kk) = props_T4(4)/props_T4(2);
            V_T4_plot(kk) = V_T4(pp);
            kk = kk + 1;
        end
    end
 plot(V_T4_plot,Beta_T4_plot,...
            ':','color',mycolor(4,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 150 K");  
     hold on   
     T5 = 300;
    V_T5 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T5)
        props_T5 = computeThermoPropsTV(T5, V_T5(pp), params_fit);  
        props_T5_check = computeThermoProps(T5, pmelt(T5), params_fit);  
        if props_T5(1)>0&&V_T5(pp)<props_T5_check(1)
            Beta_T5_plot(kk) = props_T5(4)/props_T5(2);
            
            V_T5_plot(kk) = V_T5(pp);
            kk = kk + 1;
        end
    end

 dashline(V_T5_plot,Beta_T5_plot,4,1,2,1,'color',mycolor(5,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 300 K");
 
    
 T6 = [0.0001	2	4	6	8	10	12	14	16	18	20	22	24	26	28	30	32	34	36	38	40	42	44	46	48	50	52	54	56	58	60	62	64	66	68	70	72	74	76	78	80	82	83.806];
  for ii = 1:length(T6)

          p6(ii) = psub(T6(ii));
      props_T6 = computeThermoProps(T6(ii), p6(ii), params_fit);
      Beta_T6_plot(ii) = props_T6(4)/props_T6(2);
      V_T6_plot(ii) = props_T6(1);
  end
 dashline(V_T6_plot,Beta_T6_plot,6,6,6,6,'color',mycolor(6,:),'linewidth',linewidth,'DisplayName',"sublimation");   
 T7 = [83.806	84	86	88	90	92	94	96	98	100	102	104	106	108	110	112	114	116	118	120	122	124	126	128	130	132	134	136	138	140	142	144	146	148	150	152	154	156	158	160	162	164	166	168	170	172	174	176	178	180	182	184	186	188	190	192	194	196	198	200	202	204	206	208	210	212	214	216	218	220	222	224	226	228	230	232	234	236	238	240	242	244	246	248	250	252	254	256	258	260	262	264	266	268	270	272	274	276	278	280	282	284	286	288	290	292	294	296	298	300];
   for ii = 1:length(T7)

          p7(ii) = pmelt(T7(ii));
      props_T7 = computeThermoProps(T7(ii), p7(ii), params_fit);
       Beta_T7_plot(ii) = props_T7(4)/props_T7(2);
      V_T7_plot(ii) = props_T7(1);
  end
 dashline(V_T7_plot,Beta_T7_plot,7,1,7,1,'color',mycolor(7,:),'linewidth',linewidth,'DisplayName',"melting");
         legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);

    xlabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
ylabel(['\it \beta \rm / MPa^',char(hex2dec('2212')),'^1'],'fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/evaluation')
        mkdir('../solid data fitting/figure/evaluation')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/evaluation/thermal pressure coefficient.png')  
    
    
      figure(26) 
      R = 8.31451;
 T6 = [0.0001	2	4	6	8	10	12	14	16	18	20	22	24	26	28	30	32	34	36	38	40	42	44	46	48	50	52	54	56	58	60	62	64	66	68	70	72	74	76	78	80	82	83.806];
  for ii = 1:length(T6)

          p6(ii) = psub(T6(ii));
      props_T6 = computeThermoProps(T6(ii), p6(ii), params_fit);
      cp_T6_plot_R(ii) = props_T6(5)/R;
      cv_T6_plot_R(ii) = props_T6(6)/R;
      V_T6_plot(ii) = props_T6(1);
  end
 plot(T6,cp_T6_plot_R,'color',mycolor(1,:),'linewidth',linewidth,'DisplayName',"cp-sublimation");  
 hold on
  plot(T6,cv_T6_plot_R,'color',mycolor(2,:),'linewidth',linewidth,'DisplayName',"cv-sublimation");  
 hold on
 T7 = [83.806	84	86	88	90	92	94	96	98	100	102	104	106	108	110	112	114	116	118	120	122	124	126	128	130	132	134	136	138	140	142	144	146	148	150	152	154	156	158	160	162	164	166	168	170	172	174	176	178	180	182	184	186	188	190	192	194	196	198	200	202	204	206	208	210	212	214	216	218	220	222	224	226	228	230	232	234	236	238	240	242	244	246	248	250	252	254	256	258	260	262	264	266	268	270	272	274	276	278	280	282	284	286	288	290	292	294	296	298	300	302	304	306	308	310	312	314	316	318	320	322	324	326	328	330	332	334	336	338	340	342	344	346	348	350	352	354	356	358	360	362	364	366	368	370	372	374	376	378	380	382	384	386	388	390	392	394	396	398	400	402	404	406	408	410	412	414	416	418	420	422	424	426	428	430	432	434	436	438	440	442	444	446	448	450	452	454	456	458	460	462	464	466	468	470	472	474	476	478	480	482	484	486	488	490	492	494	496	498	500	502	504	506	508	510	512	514	516	518	520	522	524	526	528	530	532	534	536	538	540	542	544	546	548	550	552	554	556	558	560	562	564	566	568	570	572	574	576	578	580	582	584	586	588	590	592	594	596	598	600	602	604	606	608	610	612	614	616	618	620	622	624	626	628	630	632	634	636	638	640	642	644	646	648	650	652	654	656	658	660	662	664	666	668	670	672	674	676	678	680	682	684	686	688	690	692	694	696	698	700	702	704	706	708	710	712	714	716	718	720	722	724	726	728	730	732	734	736	738	740	742	744	746	748	750	752	754	756	758	760	762	764	766	768	770	772	774	776	778	780	782	784	786	788	790	792	794	796	798	800];
   for ii = 1:length(T7)

          p7(ii) = pmelt(T7(ii));
      props_T7 = computeThermoProps(T7(ii), p7(ii), params_fit);
      cp_T7_plot_R(ii) = props_T7(5)/R;
      cv_T7_plot_R(ii) = props_T7(6)/R;
      V_T7_plot(ii) = props_T7(1);
  end
 plot(T7,cp_T7_plot_R,'color',mycolor(1,:),'linewidth',linewidth,'DisplayName',"cp-sublimation");
 hold on 
plot(T7,cv_T7_plot_R,'color',mycolor(2,:),'linewidth',linewidth,'DisplayName',"cv-sublimation");
 hold on  
 c(1)=plot([83.806,83.806],[0,10],'k--','linewidth',linewidth-0.2);
 ylim([0 5]);
 ytickformat('%.1f'); 
         legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);

    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
ylabel('\it C\rm_m/\itR','fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/evaluation')
        mkdir('../solid data fitting/figure/evaluation')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/evaluation/heat capacity.png')
    
    figure(27)
    T1 = 0.00001;
    V_T1 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T1)
        props_T1 = computeThermoPropsTV(T1, V_T1(pp), params_fit);
        if props_T1(1)>0
            Gamma_T1_plot(kk) = props_T1(7);
            V_T1_plot(kk) = V_T1(pp);
            kk = kk + 1;
        end
    end
 plot(V_T1_plot,Gamma_T1_plot,...
            'color',mycolor(1,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 0 K");
  hold on  
     T2 = 50;
    V_T2 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T2)
        props_T2 = computeThermoPropsTV(T2, V_T2(pp), params_fit);
        if props_T2(1)>0
            Gamma_T2_plot(kk) = props_T2(7);
            V_T2_plot(kk) = V_T2(pp);
            kk = kk + 1;
        end
    end
 plot(V_T2_plot,Gamma_T2_plot,...
            '--','color',mycolor(2,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 50 K"); 
     hold on   
     
     T3 = 83.806;
    V_T3 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T3)
        props_T3 = computeThermoPropsTV(T3, V_T3(pp), params_fit);
        if props_T3(1)>0
            Gamma_T3_plot(kk) = props_T3(7);
            V_T3_plot(kk) = V_T3(pp);
            kk = kk + 1;
        end
    end
 plot(V_T3_plot,Gamma_T3_plot,...
            '-.','color',mycolor(3,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 83.806 K");  
        hold on
        
      T4 = 150;
    V_T4 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T4)
        props_T4 = computeThermoPropsTV(T4, V_T4(pp), params_fit);
        props_T4_check = computeThermoProps(T4, pmelt(T4), params_fit);          
        if props_T4(1)>0&&V_T4(pp)<props_T4_check(1)
            Gamma_T4_plot(kk) = props_T4(7);
            V_T4_plot(kk) = V_T4(pp);
            kk = kk + 1;
        end
    end
 plot(V_T4_plot,Gamma_T4_plot,...
            ':','color',mycolor(4,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 150 K");  
     hold on   
     T5 = 300;
    V_T5 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T5)
        props_T5 = computeThermoPropsTV(T5, V_T5(pp), params_fit);  
        props_T5_check = computeThermoProps(T5, pmelt(T5), params_fit);  
        if props_T5(1)>0&&V_T5(pp)<props_T5_check(1)
            Gamma_T5_plot(kk) = props_T5(7);
            
            V_T5_plot(kk) = V_T5(pp);
            kk = kk + 1;
        end
    end

 dashline(V_T5_plot,Gamma_T5_plot,4,1,2,1,'color',mycolor(5,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 300 K");
 
    
 T6 = [0.0001	2	4	6	8	10	12	14	16	18	20	22	24	26	28	30	32	34	36	38	40	42	44	46	48	50	52	54	56	58	60	62	64	66	68	70	72	74	76	78	80	82	83.806];
  for ii = 1:length(T6)

          p6(ii) = psub(T6(ii));
      props_T6 = computeThermoProps(T6(ii), p6(ii), params_fit);
      Gamma_T6_plot(ii) = props_T6(7);
      V_T6_plot(ii) = props_T6(1);
  end
 dashline(V_T6_plot,Gamma_T6_plot,6,6,6,6,'color',mycolor(6,:),'linewidth',linewidth,'DisplayName',"sublimation");   
 T7 = [83.806	84	86	88	90	92	94	96	98	100	102	104	106	108	110	112	114	116	118	120	122	124	126	128	130	132	134	136	138	140	142	144	146	148	150	152	154	156	158	160	162	164	166	168	170	172	174	176	178	180	182	184	186	188	190	192	194	196	198	200	202	204	206	208	210	212	214	216	218	220	222	224	226	228	230	232	234	236	238	240	242	244	246	248	250	252	254	256	258	260	262	264	266	268	270	272	274	276	278	280	282	284	286	288	290	292	294	296	298	300];
V_T7_plot=[];
 for ii = 1:length(T7)

          p7(ii) = pmelt(T7(ii));
      props_T7 = computeThermoProps(T7(ii), p7(ii), params_fit);
       Gamma_T7_plot(ii) = props_T7(7);
      V_T7_plot(ii) = props_T7(1);
  end
 dashline(V_T7_plot,Gamma_T7_plot,7,1,7,1,'color',mycolor(7,:),'linewidth',linewidth,'DisplayName',"melting");
         legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);

    xlabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
ylabel('\it \gamma','fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/evaluation')
        mkdir('../solid data fitting/figure/evaluation')
    end
    print(gcf,'-dtiff','-r300','../solid data fitting/figure/evaluation/gamma.png')    
    
 end 
 
% Define the modified output function
function stop = outfun(x, optimValues, state, T_Vm_sub, p_Vm_sub, Vm_sub, ...
                       T_Vm_melt, p_Vm_melt, Vm_melt, ...
                       T_Vm_highp, p_Vm_highp, Vm_highp, ...
                       T_cp_sub, p_cp_sub, cp_sub, ...
                       T_alpha_sub, p_alpha_sub, alpha_sub, ...
                       T_BetaT_sub, p_BetaT_sub, BetaT_sub, ...
                       T_BetaS_sub, p_BetaS_sub, BetaS_sub, ...
                       T_sub, p_sub, Year_sub,G_fluid_sub, V_fluid_sub, ...
                       T_melt, p_melt, G_fluid_melt, V_fluid_melt, T_H_sub, p_H_sub, delta_H_sub, H_fluid_sub, ...
                       T_H_melt, p_H_melt, delta_H_melt, H_fluid_melt)
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
                            T_sub, p_sub, Year_sub,G_fluid_sub, V_fluid_sub, ...
                            T_melt, p_melt, G_fluid_melt, V_fluid_melt, T_H_sub, p_H_sub, delta_H_sub, H_fluid_sub, ...
                            T_H_melt, p_H_melt, delta_H_melt, H_fluid_melt);

            disp(['Current total_deviation: ', num2str(optimValues.fval)]);
            disp(['Vm_sub deviation: ', num2str(deviations.Vm_sub)]);
            disp(['Vm_melt deviation: ', num2str(deviations.Vm_melt)]);
            disp(['Vm_highp deviation: ', num2str(deviations.Vm_highp)]); 
            disp(['cp deviation: ', num2str(deviations.cp_sub)]);           
            disp(['alpha_sub deviation: ', num2str(deviations.alpha_sub)]);
            disp(['BetaT_sub deviation: ', num2str(deviations.BetaT_sub)]);
            disp(['BetaS_sub deviation: ', num2str(deviations.BetaS_sub)]);
            disp(['H_solid_sub deviation: ', num2str(deviations.H_solid_sub)]);
            disp(['H_solid_melt deviation: ', num2str(deviations.H_solid_melt)]);
            disp(['p_sub deviation: ', num2str(deviations.p_sub)]);
            disp(['p_melt deviation: ', num2str(deviations.p_melt)]); 
            disp(['Gamma T deviation: ', num2str(deviations.Gamma_T)]);             
            % Display other deviations as needed
        otherwise
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%deviation function
function [total_deviation, deviations] = combined_cost_function(params,T_Vm_sub, p_Vm_sub,Vm_sub,T_Vm_melt, p_Vm_melt,Vm_melt,T_Vm_highp, p_Vm_highp,Vm_highp,T_cp_sub, p_cp_sub,cp_sub,T_alpha_sub,p_alpha_sub,alpha_sub,T_BetaT_sub, p_BetaT_sub,BetaT_sub,T_BetaS_sub, p_BetaS_sub,BetaS_sub,T_sub, p_sub,Year_sub,G_fluid_sub,V_fluid_sub,T_melt,p_melt,G_fluid_melt,V_fluid_melt,T_H_sub,p_H_sub,delta_H_sub,H_fluid_sub,T_H_melt,p_H_melt,delta_H_melt,H_fluid_melt)
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
        if T_cp_sub(ii)<12
            cp_sub_deviation(ii) = 100*(cp_sub(ii) - props(5)) ./ cp_sub(ii);
        else
            cp_sub_deviation(ii) = 700*(cp_sub(ii) - props(5)) ./ cp_sub(ii);
        end
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
%         if T_sub(ii) < 40&& strcmp(Year_sub{ii}, '2018')
%             p_sub_deviation(ii) = 26 * p_sub_deviation(ii);
%         end
    end  
    
    p_sub_dev = sqrt(sumsqr(p_sub_deviation)/length(p_sub_deviation));    
    
    for ii = 1:length(p_melt)
        props = computeThermoProps(T_melt(ii), p_melt(ii), params);
        delta_G_melt(ii) = G_fluid_melt(ii)-props(12) + deltaH_triple - T_melt(ii) * deltaS_triple;
        p_fitted_melt(ii) = p_melt(ii) - delta_G_melt(ii)/(V_fluid_melt(ii) - props(1));%unit MPa
        p_melt_deviation(ii) = 100 * (p_melt(ii) - p_fitted_melt(ii))/p_melt(ii);
        if T_melt(ii)>500
            p_melt_deviation(ii) = p_melt_deviation(ii) * 5;
        end        
    end     
    p_melt_dev = sqrt(sumsqr(p_melt_deviation)/length(p_melt_deviation));  
    
T6 = [0.0001	2	4	6	8	10	12	14	16	18	20	22	24	26	28	30	32	34	36	38	40	42	44	46	48	50	52	54	56	58	60	62	64	66	68	70	72	74	76	78	80	82	83.806];
  for ii = 1:length(T6)
      p6(ii) = psub(T6(ii));
      props_T6 = computeThermoProps(T6(ii), p6(ii), params);
      Gamma_T6_fit(ii) = props_T6(7); 
      Vm_T6_fit(ii) = props_T6(1);    
      if ii>1
          slope(ii-1) = (Gamma_T6_fit(ii) - Gamma_T6_fit(ii-1))/(Vm_T6_fit(ii) - Vm_T6_fit(ii-1));
          if slope(ii-1)>0
              Gamma_T6_fit_dev(ii) = 45 + abs(200 * (Gamma_T6_fit(ii)-mean(Gamma_T6_fit))/mean(Gamma_T6_fit));
          else
              Gamma_T6_fit_dev(ii) = 50 * (Gamma_T6_fit(ii)-mean(Gamma_T6_fit))/mean(Gamma_T6_fit);
          end
      end
  end
    Gamma_T6_dev = sqrt(sumsqr(Gamma_T6_fit_dev)/length(Gamma_T6_fit_dev));
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
    deviations.Gamma_T = Gamma_T6_dev;
    
    total_deviation = Vm_sub_dev * 55 + Vm_melt_dev * 35 + Vm_highp_dev + cp_sub_dev * 30 + alpha_sub_dev  * 50 + BetaT_sub_dev * 20 + BetaS_sub_dev * 30 + H_solid_sub_dev * 25 + H_solid_melt_dev *2 + p_sub_dev * 2 + p_melt_dev *5 + Gamma_T6_dev*3.5;
 
end