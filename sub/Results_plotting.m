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
[Year_Vm_melt,Author_Vm_melt,T_Vm_melt,Vm_melt] = textread('../evaluation data/cell volume/Vm_melting.txt','%s%s%f%f','headerlines',2);
for ii = 1:length(Vm_melt)
    p_Vm_melt(ii) = pmelt(T_Vm_melt(ii));
end
[Year_Vm_highp,Author_Vm_highp,T_Vm_highp,Vm_highp,p_Vm_highp] = textread('../evaluation data/cell volume/Vm_high_pressure_1.txt','%s%s%f%f%f','headerlines',2);
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
[T_sub_fluid,p_sub_fluid,V_fluid_sub_fluid,G_fluid_sub_fluid] = textread('../evaluation data/sublimation/sublimation_fluid_property.txt','%f%f%f%f','headerlines',1);
[T_melt_fluid,p_melt_fluid,V_fluid_melt_fluid,G_fluid_melt_fluid] = textread('../evaluation data/melting/melting_fluid_property.txt','%f%f%f%f','headerlines',1);

[Year_H_melt,Author_H_melt,T_H_melt,delta_H_melt,H_fluid_melt] = textread('../evaluation data/enthalpy/enthalpy of fusion for fitting.txt','%s%s%f%f%f','headerlines',2);
for ii = 1:length(T_H_melt)
    p_H_melt(ii) = pmelt(T_H_melt(ii));
end
my_dashline=[0 0 0 0;4 1 1 2;2 0.5 0.5 1;3 2 3 2;1 1 1 1;7 1 7 1;6 6 6 6;1 4 1 4;13 3 13 3;];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choice = 1;%Vm along sublimation
% choice = 2;%Vm along melting
% choice = 3;%heat capacity along sublimation
% choice = 4;%thermal expansion
%  choice = 5;%isothermal bulk modulus
% choice = 6;%isentropic bulk modulus
% choice = 7;%enthalpy sublimation
% choice = 8;%enthalpy melting
choice = 9;%sublimation
% choice = 10;%melting
% choice = 11;%extrapolation cp plot
% choice = 12;%extrapolation gamma
% choice = 13;%extrapolation four plots 
% choice = 14;%volume high pressure
% choice = 14.5;%volume high pressure legend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 1
    [Year_Vm_sub,Author_Vm_sub,T_Vm_sub,Vm_sub] = textread('../evaluation data/cell volume/Vm_sublimation.txt','%s%s%f%f','headerlines',2);
    for ii = 1:length(Vm_sub)
        p_Vm_sub(ii) = psub(T_Vm_sub(ii));
    end
T_plot_calc = 1:83.806;
for ii = 1:length(T_plot_calc)
    if T_plot_calc(ii)<83.806
        p_plot_calc(ii) = psub(T_plot_calc(ii));
    else
        p_plot_calc(ii) = pmelt(T_plot_calc(ii));
    end

    fitted_props(ii,:) = computeThermoPropsForResult(T_plot_calc(ii), p_plot_calc(ii));
end
% Example: Compare experimental and fitted molar volume for sublimation


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
    ytickformat('%.1f');
                xlim([0 90]);
       c(1)= plot(T_plot_calc, fitted_props(:,1),...
            'color',mycolor(1,:),'linewidth',linewidth);
      annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
      
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
    ylabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[500,200,500,340])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../../paper writing/figure/cell volume')
        mkdir('../../paper writing/figure/cell volume')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/cell volume/sublimation.png')
    
figure(2)%deviation plot Vm_sub

    for ii = 1:length(T_Vm_sub)
        props_Vm_sub = computeThermoPropsForResult(T_Vm_sub(ii), p_Vm_sub(ii));
        Vm_sub_deviation(ii) = 100*(Vm_sub(ii) - props_Vm_sub(1)) ./ Vm_sub(ii);
    end

for ii = 1:xxcount_Vm_sub
        plot(T_Vm_sub(Nxstart_Vm_sub(ii):Nxend_Vm_sub(ii)),Vm_sub_deviation(Nxstart_Vm_sub(ii):Nxend_Vm_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_Vm_sub(Nxstart_Vm_sub(ii)),"-",Author_Vm_sub(Nxstart_Vm_sub(ii)))));
        hold on    
end
    annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfb','FontSize',fontsize,'FontName','Times');

            c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
            xlim([0 90]);
        ytickformat('%.1f');    
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
    ylabel(['100\cdot',char(8201),'(\itV\rm_m_,_e_x_p',char(hex2dec('2212')),'\it V\rm_m_,_c_a_l_c) /\it V\rm_m_,_e_x_p'],'fontsize',fontsize,'linewidth',linewidth)  
     set(gcf,'position',[500,200,500,340])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')

    print(gcf,'-dtiff','-r300','../../paper writing/figure/cell volume/sublimation deviation.png') 
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 2
    
T_plot_calc = 83.806:400;
for ii = 1:length(T_plot_calc)
    if T_plot_calc(ii)<83.806
        p_plot_calc(ii) = psub(T_plot_calc(ii));
    else
        p_plot_calc(ii) = pmelt(T_plot_calc(ii));
    end

    fitted_props(ii,:) = computeThermoPropsForResult(T_plot_calc(ii), p_plot_calc(ii));
end
% Example: Compare experimental and fitted molar volume for sublimation

figure(1);
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
       c(1)= plot(T_plot_calc, fitted_props(:,1),...
            'color',mycolor(1,:),'linewidth',linewidth);
      annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
        xlim([80 400]);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
    ylabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    set(gcf,'position',[500,200,500,340])
%     set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../../paper writing/figure/cell volume')
        mkdir('../../paper writing/figure/cell volume')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/cell volume/Vm melting.png')
  
    
figure(2)%deviation plot Vm_melt

    for ii = 1:length(T_Vm_melt)
        props_Vm_melt = computeThermoPropsForResult(T_Vm_melt(ii), p_Vm_melt(ii));
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
              annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfb','FontSize',fontsize,'FontName','Times');
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
    ylabel(['100\cdot',char(8201),'(\itV\rm_m_,_e_x_p',char(hex2dec('2212')),'\it V\rm_m_,_c_a_l_c) /\it V\rm_m_,_e_x_p'],'fontsize',fontsize,'linewidth',linewidth)  
%      set(gcf,'position',[300,100,650,540])
    set(gcf,'position',[500,200,500,340])    
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/cell volume')
        mkdir('../solid data fitting/figure/cell volume')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/cell volume/Vm melting deviation.png') 
    
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 3
    [Year_cp_sub,Author_cp_sub,T_cp_sub,cp_sub] = textread('../evaluation data/heat capacity/cp_sublimation.txt','%s%s%f%f','headerlines',2);
    for ii = 1:length(cp_sub)
        p_cp_sub(ii) = psub(T_cp_sub(ii));
    end
T_plot_calc = 1:83.806;
for ii = 1:length(T_plot_calc)
    if T_plot_calc(ii)<83.806
        p_plot_calc(ii) = psub(T_plot_calc(ii));
    else
        p_plot_calc(ii) = pmelt(T_plot_calc(ii));
    end

    fitted_props(ii,:) = computeThermoPropsForResult(T_plot_calc(ii), p_plot_calc(ii));
end
% Example: Compare experimental and fitted molar volume for sublimation


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
figure(1);
for ii = 1:xxcount_cp_sub
        plot(T_cp_sub(Nxstart_cp_sub(ii):Nxend_cp_sub(ii)),cp_sub(Nxstart_cp_sub(ii):Nxend_cp_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_cp_sub(Nxstart_cp_sub(ii)),"-",Author_cp_sub(Nxstart_cp_sub(ii)))));
        hold on    
end
%     ytickformat('%.1f');
                xlim([0 90]);
       c(1)= plot(T_plot_calc, fitted_props(:,5),...
            'color',mycolor(1,:),'linewidth',linewidth);
      annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
      
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
    ylabel(['\it c_p\rm / (J K^',char(hex2dec('2212')),'^1 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
%     ylabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[500,200,500,340])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../../paper writing/figure/heat capacity')
        mkdir('../../paper writing/figure/heat capacity')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/heat capacity/sublimation.png')
    
figure(2)%deviation plot cp_sub

    for ii = 1:length(T_cp_sub)
        props_cp_sub = computeThermoPropsForResult(T_cp_sub(ii), p_cp_sub(ii));
        cp_sub_deviation(ii) = 100*(cp_sub(ii) - props_cp_sub(5)) ./ cp_sub(ii);
    end

for ii = 1:xxcount_cp_sub
        plot(T_cp_sub(Nxstart_cp_sub(ii):Nxend_cp_sub(ii)),cp_sub_deviation(Nxstart_cp_sub(ii):Nxend_cp_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_cp_sub(Nxstart_cp_sub(ii)),"-",Author_cp_sub(Nxstart_cp_sub(ii)))));
        hold on    
end
    annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfb','FontSize',fontsize,'FontName','Times');

            c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
            xlim([0 90]);
%         ytickformat('%.1f');    
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
ylabel(['100\cdot (\itc_p\rm_,_e_x_p',char(hex2dec('2212')),' \itc_p\rm_,_c_a_l_c) / \itc_p\rm_,_e_x_p'],'fontsize',fontsize,'linewidth',linewidth);
set(gcf,'position',[500,200,500,340])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')

    print(gcf,'-dtiff','-r300','../../paper writing/figure/heat capacity/sublimation deviation.png') 
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 4
    [Year_alpha_sub,Author_alpha_sub,T_alpha_sub,alpha_sub] = textread('../evaluation data/thermal expansion/alpha_sublimation.txt','%s%s%f%f','headerlines',2);
    for ii = 1:length(alpha_sub)
        p_alpha_sub(ii) = psub(T_alpha_sub(ii));
    end
T_plot_calc = 1:83.806;
for ii = 1:length(T_plot_calc)
    if T_plot_calc(ii)<83.806
        p_plot_calc(ii) = psub(T_plot_calc(ii));
    else
        p_plot_calc(ii) = pmelt(T_plot_calc(ii));
    end

    fitted_props(ii,:) = computeThermoPropsForResult(T_plot_calc(ii), p_plot_calc(ii));
end
% Example: Compare experimental and fitted molar volume for sublimation


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
figure(1);
for ii = 1:xxcount_alpha_sub
        plot(T_alpha_sub(Nxstart_alpha_sub(ii):Nxend_alpha_sub(ii)),alpha_sub(Nxstart_alpha_sub(ii):Nxend_alpha_sub(ii))*10^4,...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_alpha_sub(Nxstart_alpha_sub(ii)),"-",Author_alpha_sub(Nxstart_alpha_sub(ii)))));
        hold on    
end
%     ytickformat('%.1f');
                xlim([0 90]);
       c(1)= plot(T_plot_calc, fitted_props(:,4)*10^4,...
            'color',mycolor(1,:),'linewidth',linewidth);
      annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
      
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
    ylabel('\it \alpha \rm\cdot10^4\rm / K^-^1','fontsize',fontsize,'linewidth',linewidth);
%     ylabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[500,200,500,340])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../../paper writing/figure/thermal expansion')
        mkdir('../../paper writing/figure/thermal expansion')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/thermal expansion/sublimation.png')
    
figure(2)%deviation plot alpha_sub

    for ii = 1:length(T_alpha_sub)
        props_alpha_sub = computeThermoPropsForResult(T_alpha_sub(ii), p_alpha_sub(ii));
        alpha_sub_deviation(ii) = 100*(alpha_sub(ii) - props_alpha_sub(4)) ./ alpha_sub(ii);
    end

for ii = 1:xxcount_alpha_sub
        plot(T_alpha_sub(Nxstart_alpha_sub(ii):Nxend_alpha_sub(ii)),alpha_sub_deviation(Nxstart_alpha_sub(ii):Nxend_alpha_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_alpha_sub(Nxstart_alpha_sub(ii)),"-",Author_alpha_sub(Nxstart_alpha_sub(ii)))));
        hold on    
end
    annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfb','FontSize',fontsize,'FontName','Times');

            c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
            xlim([0 90]);
%         ytickformat('%.1f');    
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
ylabel(['100\cdot (\it\alpha\rm_e_x_p',char(hex2dec('2212')),' \it\alpha\rm_c_a_l_c) / \it\alpha\rm_e_x_p'],'fontsize',fontsize,'linewidth',linewidth);
set(gcf,'position',[500,200,500,340])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')

    print(gcf,'-dtiff','-r300','../../paper writing/figure/thermal expansion/sublimation deviation.png') 
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 5
    [Year_KT_sub,Author_KT_sub,T_KT_sub,KT_sub] = textread('../evaluation data/bulk modulus/BetaT_sublimation.txt','%s%s%f%f','headerlines',2);
    KT_sub = 1./KT_sub;
    for ii = 1:length(KT_sub)
        p_KT_sub(ii) = psub(T_KT_sub(ii));
    end
T_plot_calc = 1:83.806;
for ii = 1:length(T_plot_calc)
    if T_plot_calc(ii)<83.806
        p_plot_calc(ii) = psub(T_plot_calc(ii));
    else
        p_plot_calc(ii) = pmelt(T_plot_calc(ii));
    end

    fitted_props(ii,:) = computeThermoPropsForResult(T_plot_calc(ii), p_plot_calc(ii));
end
% Example: Compare experimental and fitted molar volume for sublimation


ndata_KT_sub=length(T_KT_sub);xxcount_KT_sub=1; Nxstart_KT_sub = 1;yy_KT_sub=Year_KT_sub(1);AA_KT_sub = Author_KT_sub(1);
if ndata_KT_sub == 1
    Nxend_KT_sub(xxcount_KT_sub) = 1;
else
    for ii = 2:ndata_KT_sub
        if ~strcmp(yy_KT_sub,string(Year_KT_sub(ii)))||~strcmp(AA_KT_sub,Author_KT_sub(ii))
            Nxend_KT_sub(xxcount_KT_sub) = ii-1;
            xxcount_KT_sub = xxcount_KT_sub + 1;
            Nxstart_KT_sub(xxcount_KT_sub) = ii;
            yy_KT_sub=Year_KT_sub(ii);AA_KT_sub=Author_KT_sub(ii);
        end
        if ii == ndata_KT_sub
            Nxend_KT_sub(xxcount_KT_sub) = ndata_KT_sub;
        end
    end
end
figure(1);
for ii = 1:xxcount_KT_sub
        plot(T_KT_sub(Nxstart_KT_sub(ii):Nxend_KT_sub(ii)),KT_sub(Nxstart_KT_sub(ii):Nxend_KT_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_KT_sub(Nxstart_KT_sub(ii)),"-",Author_KT_sub(Nxstart_KT_sub(ii)))));
        hold on    
end
%     ytickformat('%.1f');
                xlim([0 90]);
       c(1)= plot(T_plot_calc, 1./fitted_props(:,2),...
            'color',mycolor(1,:),'linewidth',linewidth);
      annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
      
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
ylabel('\it K_T \rm / MPa','fontsize',fontsize,'linewidth',linewidth);
%     ylabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[500,200,500,340])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../../paper writing/figure/bulk modulus')
        mkdir('../../paper writing/figure/bulk modulus')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/bulk modulus/KT_sublimation.png')
    
figure(2)%deviation plot KT_sub

    for ii = 1:length(T_KT_sub)
        props_KT_sub = computeThermoPropsForResult(T_KT_sub(ii), p_KT_sub(ii));
        KT_sub_deviation(ii) = 100*(KT_sub(ii) - 1./props_KT_sub(2)) ./ KT_sub(ii);
    end

for ii = 1:xxcount_KT_sub
        plot(T_KT_sub(Nxstart_KT_sub(ii):Nxend_KT_sub(ii)),KT_sub_deviation(Nxstart_KT_sub(ii):Nxend_KT_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_KT_sub(Nxstart_KT_sub(ii)),"-",Author_KT_sub(Nxstart_KT_sub(ii)))));
        hold on    
end
    annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfb','FontSize',fontsize,'FontName','Times');

            c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
            xlim([0 90]);
%         ytickformat('%.1f');    
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
    ylabel(['100\cdot',char(8201),'(\itK_T\rm_,_e_x_p',char(hex2dec('2212')),'\it K_T\rm_,_c_a_l_c) /\it K_T\rm_,_e_x_p'],'fontsize',fontsize,'linewidth',linewidth)  
set(gcf,'position',[500,200,500,340])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')

    print(gcf,'-dtiff','-r300','../../paper writing/figure/bulk modulus/KT_sublimation deviation.png') 
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 6
    [Year_KS_sub,Author_KS_sub,T_KS_sub,KS_sub] = textread('../evaluation data/bulk modulus/BetaS_sublimation.txt','%s%s%f%f','headerlines',2);
    KS_sub = 1./KS_sub;
    for ii = 1:length(KS_sub)
        p_KS_sub(ii) = psub(T_KS_sub(ii));
    end
T_plot_calc = 1:83.806;
for ii = 1:length(T_plot_calc)
    if T_plot_calc(ii)<83.806
        p_plot_calc(ii) = psub(T_plot_calc(ii));
    else
        p_plot_calc(ii) = pmelt(T_plot_calc(ii));
    end

    fitted_props(ii,:) = computeThermoPropsForResult(T_plot_calc(ii), p_plot_calc(ii));
end
% Example: Compare experimental and fitted molar volume for sublimation


ndata_KS_sub=length(T_KS_sub);xxcount_KS_sub=1; Nxstart_KS_sub = 1;yy_KS_sub=Year_KS_sub(1);AA_KS_sub = Author_KS_sub(1);
if ndata_KS_sub == 1
    Nxend_KS_sub(xxcount_KS_sub) = 1;
else
    for ii = 2:ndata_KS_sub
        if ~strcmp(yy_KS_sub,string(Year_KS_sub(ii)))||~strcmp(AA_KS_sub,Author_KS_sub(ii))
            Nxend_KS_sub(xxcount_KS_sub) = ii-1;
            xxcount_KS_sub = xxcount_KS_sub + 1;
            Nxstart_KS_sub(xxcount_KS_sub) = ii;
            yy_KS_sub=Year_KS_sub(ii);AA_KS_sub=Author_KS_sub(ii);
        end
        if ii == ndata_KS_sub
            Nxend_KS_sub(xxcount_KS_sub) = ndata_KS_sub;
        end
    end
end
figure(1);
for ii = 1:xxcount_KS_sub
        plot(T_KS_sub(Nxstart_KS_sub(ii):Nxend_KS_sub(ii)),KS_sub(Nxstart_KS_sub(ii):Nxend_KS_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_KS_sub(Nxstart_KS_sub(ii)),"-",Author_KS_sub(Nxstart_KS_sub(ii)))));
        hold on    
end
%     ytickformat('%.1f');
                xlim([0 90]);
       c(1)= plot(T_plot_calc, 1./fitted_props(:,3),...
            'color',mycolor(1,:),'linewidth',linewidth);
      annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
      
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
ylabel('\it K_S \rm / MPa','fontsize',fontsize,'linewidth',linewidth);
%     ylabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[500,200,500,340])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../../paper writing/figure/bulk modulus')
        mkdir('../../paper writing/figure/bulk modulus')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/bulk modulus/KS_sublimation.png')
    
figure(2)%deviation plot KS_sub

    for ii = 1:length(T_KS_sub)
        props_KS_sub = computeThermoPropsForResult(T_KS_sub(ii), p_KS_sub(ii));
        KS_sub_deviation(ii) = 100*(KS_sub(ii) - 1./props_KS_sub(3)) ./ KS_sub(ii);
    end

for ii = 1:xxcount_KS_sub
        plot(T_KS_sub(Nxstart_KS_sub(ii):Nxend_KS_sub(ii)),KS_sub_deviation(Nxstart_KS_sub(ii):Nxend_KS_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_KS_sub(Nxstart_KS_sub(ii)),"-",Author_KS_sub(Nxstart_KS_sub(ii)))));
        hold on    
end
    annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfb','FontSize',fontsize,'FontName','Times');

            c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
            xlim([0 90]);
%         ytickformat('%.1f');    
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
    ylabel(['100\cdot',char(8201),'(\itK_S\rm_,_e_x_p',char(hex2dec('2212')),'\it K_S\rm_,_c_a_l_c) /\it K_S\rm_,_e_x_p'],'fontsize',fontsize,'linewidth',linewidth)  
set(gcf,'position',[500,200,500,340])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')

    print(gcf,'-dtiff','-r300','../../paper writing/figure/bulk modulus/KS_sublimation deviation.png') 
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 7
    figure(1);
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
    for ii = 1:xxcount_H_sub
            plot(T_H_sub(Nxstart_H_sub(ii):Nxend_H_sub(ii)),H_solid_sub(Nxstart_H_sub(ii):Nxend_H_sub(ii))/1000,...
                mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_H_sub(Nxstart_H_sub(ii)),"-",Author_H_sub(Nxstart_H_sub(ii)))));
            hold on    
    end
    % Sort T_H_sub and p_fitted_H_sub in ascending order
    [T_H_sub_sorted, sortIdx] = sort(T_H_sub);
    H_solid_sub_fitted_sorted = H_solid_sub_fitted(sortIdx);

    % Plot the sorted values
    plot(T_H_sub_sorted, H_solid_sub_fitted_sorted/1000,...
            'color',mycolor(1,:),'linewidth',linewidth);
      annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    ytickformat('%.1f'); 
%             legend show
%             legend('boxoff');
%             legend ('Location','bestoutside','NumColumns',1);
        xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
    ylabel(['\it \DeltaH \rm/ (kJ mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
        % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
        % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[500,200,500,340])
        set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
        if ~isfolder('../../paper writing/figure/enthalpy')
            mkdir('../../paper writing/figure/enthalpy')
        end
        print(gcf,'-dtiff','-r300','../../paper writing/figure/enthalpy/H_sub.png')

    figure(2)%deviation plot enthalpy

    for ii = 1:xxcount_H_sub
            plot(T_H_sub(Nxstart_H_sub(ii):Nxend_H_sub(ii)),H_solid_sub_deviation(Nxstart_H_sub(ii):Nxend_H_sub(ii)),...
                mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_H_sub(Nxstart_H_sub(ii)),"-",Author_H_sub(Nxstart_H_sub(ii)))));
            hold on    
    end
                c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
                xlim([74 84]);
    ytickformat('%.2f'); 
%             legend show
%             legend('boxoff');
%             legend ('Location','bestoutside','NumColumns',1);
        xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
        ylabel(['100\cdot',char(8201),'(\it\DeltaH\rm_e_x_p',char(hex2dec('2212')),'\it \DeltaH\rm_c_a_l_c) /\it\DeltaH\rm_e_x_p'],'fontsize',fontsize,'linewidth',linewidth)  
        annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfb','FontSize',fontsize,'FontName','Times');
        set(gcf,'position',[500,200,500,340])
        set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
        if ~isfolder('../../paper writing/figure/enthalpy')
            mkdir('../../paper writing/figure/enthalpy')
        end
        print(gcf,'-dtiff','-r300','../../paper writing/figure/enthalpy/H_sub_dev.png')       
     % 
    % for ii = 1:length(BetaS_sub)
    %     p_sub(ii) = psub(T_sub(ii));
    % end

 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 8
    figure(1);
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
    plot(T_H_melt_sorted, H_solid_melt_fitted_sorted/1000,...
            'color',mycolor(1,:),'linewidth',linewidth);
      annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
%     ytickformat('%.1f'); 
    xlim([80 380]);
    ylim([-10 32]);
%             legend show
%             legend('boxoff');
%             legend ('Location','bestoutside','NumColumns',1);
        xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
    ylabel(['\it \DeltaH \rm/ (kJ mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
        % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
        % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[500,200,500,340])
        set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
        if ~isfolder('../../paper writing/figure/enthalpy')
            mkdir('../../paper writing/figure/enthalpy')
        end
        print(gcf,'-dtiff','-r300','../../paper writing/figure/enthalpy/H_melt.png')

    figure(2)%deviation plot enthalpy

    for ii = 1:xxcount_H_melt
            plot(T_H_melt(Nxstart_H_melt(ii):Nxend_H_melt(ii)),H_solid_melt_deviation(Nxstart_H_melt(ii):Nxend_H_melt(ii)),...
                mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_H_melt(Nxstart_H_melt(ii)),"-",Author_H_melt(Nxstart_H_melt(ii)))));
            hold on    
    end
                c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
                    xlim([80 380]);
%                 xlim([74 84]);
%     ytickformat('%.2f'); 
%             legend show
%             legend('boxoff');
%             legend ('Location','bestoutside','NumColumns',1);
        xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
        ylabel(['100\cdot',char(8201),'(\it\DeltaH\rm_e_x_p',char(hex2dec('2212')),'\it \DeltaH\rm_c_a_l_c) /\it\DeltaH\rm_e_x_p'],'fontsize',fontsize,'linewidth',linewidth)  
        annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfb','FontSize',fontsize,'FontName','Times');
        set(gcf,'position',[500,200,500,340])
        set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
        if ~isfolder('../../paper writing/figure/enthalpy')
            mkdir('../../paper writing/figure/enthalpy')
        end
        print(gcf,'-dtiff','-r300','../../paper writing/figure/enthalpy/H_melt_dev.png')       
     % 
    % for ii = 1:length(BetaS_melt)
    %     p_melt(ii) = psub(T_melt(ii));
    % end

 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    figure(1);
for ii = 1:xxcount_sub
        plot(T_sub(Nxstart_sub(ii):Nxend_sub(ii)),p_sub(Nxstart_sub(ii):Nxend_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_sub(Nxstart_sub(ii)),"-",Author_sub(Nxstart_sub(ii)))));
        hold on    
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

% Sort T_sub and p_fitted_sub in ascending order
% [T_sub_sorted, sortIdx] = sort(T_sub);
% p_fitted_sub_sorted = p_fitted_sub(sortIdx);
% 
% % Plot the sorted values
% plot(T_sub_sorted, p_fitted_sub_sorted,...
%             'color',mycolor(1,:),'linewidth',linewidth);
for ii = 1:length(T_sub_fluid)
    p_aux(ii) = psub(T_sub_fluid(ii));
    fitted_props_sub_plot = computeThermoPropsForResult(T_sub_fluid(ii), p_sub_fluid(ii));
        delta_G_sub_plot(ii,:) = G_fluid_sub_fluid(ii) - fitted_props_sub_plot(12) + deltaH_triple - T_sub_fluid(ii) * deltaS_triple;
    p_fitted_sub_plot(ii) = p_sub_fluid(ii) - delta_G_sub_plot(ii)/(V_fluid_sub_fluid(ii)*10^6 - fitted_props_sub_plot(1)) ; % unit Pa
    p_aux_dev(ii) = 100 * (p_aux(ii) - p_fitted_sub_plot(ii)) / p_aux(ii);
   p_fitted_sub_plot(ii) = p_fitted_sub_plot(ii) * 10^6;
end
plot(T_sub_fluid, p_fitted_sub_plot,...
            'color',mycolor(1,:),'linewidth',linewidth);
set(gca, 'YScale', 'log');  % Set y-axis to logarithmic scale
ylim([1.5*10^-7 10^5]);
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
        annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
ylabel('\it p \rm / Pa','fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
        set(gcf,'position',[500,200,500,340])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
        if ~isfolder('../../paper writing/figure/sublimation')
            mkdir('../../paper writing/figure/sublimation')
        end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/sublimation/sublimation.png')
    
figure(2)%deviation plot sublimation

for ii = 1:xxcount_sub
        plot(T_sub(Nxstart_sub(ii):Nxend_sub(ii)),p_sub_deviation(Nxstart_sub(ii):Nxend_sub(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_sub(Nxstart_sub(ii)),"-",Author_sub(Nxstart_sub(ii)))));
        hold on    
end
    plot(T_sub_fluid,p_aux_dev,'--','color','r');
    hold on
            c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
            xlim([40 90]);
            ylim([-17 10]);

%     xlim([250 275])
%     ylim([-1.02 1.02])
%     ytickformat('%.2f')
%     xticks([250 270])            
            
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
        annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfb','FontSize',fontsize,'FontName','Times');
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
    ylabel(['100\cdot',char(8201),'(\itp\rm_e_x_p',char(hex2dec('2212')),'\it p\rm_c_a_l_c) /\it p\rm_e_x_p'],'fontsize',fontsize,'linewidth',linewidth)  
        set(gcf,'position',[500,200,500,340])
     set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
       
    axes('Position',[.565 .27 .18 .18])
    box on
%     plot(T_aux,aux_dev,'--','color','r');
    hold on
    plot([0,55],[0,0],'color','k','linewidth',linewidth-0.1)
    ylim([-150 -20]);
    xlim([20 60]);
    xticks([30 50])
    yticks([-120 -80 -40])    
 for ii = 1:xxcount_sub
        plot(T_sub(Nxstart_sub(ii):Nxend_sub(ii)),p_sub_deviation(Nxstart_sub(ii):Nxend_sub(ii)),...
            mymarker{ii},'markersize',markersize-2,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_sub(Nxstart_sub(ii)),"-",Author_sub(Nxstart_sub(ii)))));
        hold on    
end   
        plot(T_sub_fluid,p_aux_dev,'--','color','r');
    hold on                
        
    set(gca,'fontsize',fontsize-2,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../solid data fitting/figure/sublimation')
        mkdir('../solid data fitting/figure/sublimation')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/sublimation/Sublimation deviation.png')  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    figure(1);
for ii = 1:xxcount_melt
        plot(T_melt(Nxstart_melt(ii):Nxend_melt(ii)),p_melt(Nxstart_melt(ii):Nxend_melt(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_melt(Nxstart_melt(ii)),"-",Author_melt(Nxstart_melt(ii)))));
        hold on    
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
    p_fitted_melt(ii) = p_melt(ii) - delta_G_melt(ii)/(V_fluid_melt(ii) - fitted_props_melt(1)); % unit Pa
    p_melt_deviation(ii) = 100 * (p_melt(ii) - p_fitted_melt(ii)) / p_melt(ii);

end

% Sort T_melt and p_fitted_melt in ascending order
% [T_melt_sorted, sortIdx] = sort(T_melt);
% p_fitted_melt_sorted = p_fitted_melt(sortIdx);
% 
% % Plot the sorted values
% plot(T_melt_sorted, p_fitted_melt_sorted,...
%             'color',mycolor(1,:),'linewidth',linewidth);
for ii = 1:length(T_melt_fluid)
    p_aux(ii) = real(pmelt(T_melt_fluid(ii)));
    fitted_props_melt_plot = computeThermoPropsForResult(T_melt_fluid(ii), p_melt_fluid(ii));
        delta_G_melt_plot(ii,:) = G_fluid_melt_fluid(ii) - fitted_props_melt_plot(12) + deltaH_triple - T_melt_fluid(ii) * deltaS_triple;
    p_fitted_melt_plot(ii) = p_melt_fluid(ii) - delta_G_melt_plot(ii)/(V_fluid_melt_fluid(ii) - fitted_props_melt_plot(1)) ; % unit Pa
    p_aux_dev(ii) = 100 * (p_aux(ii) - p_fitted_melt_plot(ii)) / p_aux(ii);
%    p_fitted_melt_plot(ii) = p_fitted_melt_plot(ii) * 10^6;
end
ylim([0.1 10^4]);
plot(T_melt_fluid, p_fitted_melt_plot,...
            'color',mycolor(1,:),'linewidth',linewidth);
set(gca, 'YScale', 'log');  % Set y-axis to logarithmic scale
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
        annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
ylabel('\it p \rm / MPa','fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
        set(gcf,'position',[500,200,500,340])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
        if ~isfolder('../../paper writing/figure/melting')
            mkdir('../../paper writing/figure/melting')
        end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/melting/melting.png')
    
figure(2)%deviation plot melting

for ii = 1:xxcount_melt
        plot(T_melt(Nxstart_melt(ii):Nxend_melt(ii)),p_melt_deviation(Nxstart_melt(ii):Nxend_melt(ii)),...
            mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_melt(Nxstart_melt(ii)),"-",Author_melt(Nxstart_melt(ii)))));
        hold on    
end
    plot(T_melt_fluid,p_aux_dev,'--','color','r');
    hold on
            c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
%             xlim([40 90]);
%             ylim([-17 10]);

    xlim([0 800])
    ylim([-10 20])
%     ytickformat('%.2f')
%     xticks([250 270])            
            
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
        annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfb','FontSize',fontsize,'FontName','Times');
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);
    ylabel(['100\cdot',char(8201),'(\itp\rm_e_x_p',char(hex2dec('2212')),'\it p\rm_c_a_l_c) /\it p\rm_e_x_p'],'fontsize',fontsize,'linewidth',linewidth)  
        set(gcf,'position',[500,200,500,340])
     set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
       
%     axes('Position',[.565 .27 .18 .18])
%     box on
% %     plot(T_aux,aux_dev,'--','color','r');
%     hold on
%     plot([0,55],[0,0],'color','k','linewidth',linewidth-0.1)
% %     ylim([-150 -20]);
% %     xlim([20 60]);
%     xticks([30 50])
%     yticks([-120 -80 -40])    
%  for ii = 1:xxcount_melt
%         plot(T_melt(Nxstart_melt(ii):Nxend_melt(ii)),p_melt_deviation(Nxstart_melt(ii):Nxend_melt(ii)),...
%             mymarker{ii},'markersize',markersize-2,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',string(strcat(Year_melt(Nxstart_melt(ii)),"-",Author_melt(Nxstart_melt(ii)))));
%         hold on    
% end   
%         plot(T_melt_fluid,p_aux_dev,'--','color','r');
%     hold on                
        
%     set(gca,'fontsize',fontsize-2,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
%     if ~isfolder('../solid data fitting/figure/melting')
%         mkdir('../solid data fitting/figure/melting')
%     end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/melting/Sublimation deviation.png')  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if choice == 11
      figure(1) 
      R = 8.31451;
 T6 = [0.0001	2	4	6	8	10	12	14	16	18	20	22	24	26	28	30	32	34	36	38	40	42	44	46	48	50	52	54	56	58	60	62	64	66	68	70	72	74	76	78	80	82	83.806];
  for ii = 1:length(T6)

          p6(ii) = psub(T6(ii));
      props_T6 = computeThermoPropsForResult(T6(ii), p6(ii));
      cp_T6_plot_R(ii) = props_T6(5)/R;
      cv_T6_plot_R(ii) = props_T6(6)/R;
      V_T6_plot(ii) = props_T6(1);
  end
 plot(T6,cp_T6_plot_R,'color',mycolor(1,:),'linewidth',linewidth,'DisplayName',"cp-sublimation");  
 hold on
  plot(T6,cv_T6_plot_R, '--','color',mycolor(2,:),'linewidth',linewidth,'DisplayName',"cv-sublimation");  
 hold on
 T7 = [83.806	84	86	88	90	92	94	96	98	100	102	104	106	108	110	112	114	116	118	120	122	124	126	128	130	132	134	136	138	140	142	144	146	148	150	152	154	156	158	160	162	164	166	168	170	172	174	176	178	180	182	184	186	188	190	192	194	196	198	200	202	204	206	208	210	212	214	216	218	220	222	224	226	228	230	232	234	236	238	240	242	244	246	248	250	252	254	256	258	260	262	264	266	268	270	272	274	276	278	280	282	284	286	288	290	292	294	296	298	300	302	304	306	308	310	312	314	316	318	320	322	324	326	328	330	332	334	336	338	340	342	344	346	348	350	352	354	356	358	360	362	364	366	368	370	372	374	376	378	380	382	384	386	388	390	392	394	396	398	400	402	404	406	408	410	412	414	416	418	420	422	424	426	428	430	432	434	436	438	440	442	444	446	448	450	452	454	456	458	460	462	464	466	468	470	472	474	476	478	480	482	484	486	488	490	492	494	496	498	500	502	504	506	508	510	512	514	516	518	520	522	524	526	528	530	532	534	536	538	540	542	544	546	548	550	552	554	556	558	560	562	564	566	568	570	572	574	576	578	580	582	584	586	588	590	592	594	596	598	600	602	604	606	608	610	612	614	616	618	620	622	624	626	628	630	632	634	636	638	640	642	644	646	648	650	652	654	656	658	660	662	664	666	668	670	672	674	676	678	680	682	684	686	688	690	692	694	696	698	700	702	704	706	708	710	712	714	716	718	720	722	724	726	728	730	732	734	736	738	740	742	744	746	748	750	752	754	756	758	760	762	764	766	768	770	772	774	776	778	780	782	784	786	788	790	792	794	796	798	800];
   for ii = 1:length(T7)

          p7(ii) = pmelt(T7(ii));
      props_T7 = computeThermoPropsForResult(T7(ii), p7(ii));
      cp_T7_plot_R(ii) = props_T7(5)/R;
      cv_T7_plot_R(ii) = props_T7(6)/R;
      V_T7_plot(ii) = props_T7(1);
  end
 plot(T7,cp_T7_plot_R,'color',mycolor(1,:),'linewidth',linewidth,'DisplayName',"cp-sublimation");
 hold on 
plot(T7, cv_T7_plot_R, '--', 'color', mycolor(2,:), 'linewidth', linewidth, 'DisplayName', "cv-sublimation");
 hold on  
% Blue-green color can be set using RGB values (e.g., [0, 0.5, 0.5])
% c(1) = plot([83.806, 83.806], [0, 10], '-.', 'color',[1 175 146]./256, 'linewidth', linewidth);
dashline([83.806, 83.806], [0, 10],4,1,1,2,'color',	[1 175 146]./256,'linewidth',linewidth);
 ylim([0 5]);
%  ytickformat('%.1f'); 
%          legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);

    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
ylabel('\it C\rm_m/\itR','fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
        set(gcf,'position',[500,200,500,340])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../../paper writing/figure/extrapolation/')
        mkdir('../../paper writing/figure/extrapolation/')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/extrapolation/heat capacity.png')    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 12
    T1 = 0.00001;
    V_T1 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T1)
        props_T1 = computeThermoPropsTVForResult(T1, V_T1(pp));
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
        props_T2 = computeThermoPropsTVForResult(T2, V_T2(pp));
        if props_T2(1)>0
            Gamma_T2_plot(kk) = props_T2(7);
            V_T2_plot(kk) = V_T2(pp);
            kk = kk + 1;
        end
    end
 plot(V_T2_plot,Gamma_T2_plot,...
            '--','color',mycolor(2,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 50 K"); 
% dashline(V_T2_plot,Gamma_T2_plot,6,3,4,2,'color',mycolor(2,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 50 K"); 
     hold on   
     
     T3 = 83.806;
    V_T3 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T3)
        props_T3 = computeThermoPropsTVForResult(T3, V_T3(pp));
        if props_T3(1)>0
            Gamma_T3_plot(kk) = props_T3(7);
            V_T3_plot(kk) = V_T3(pp);
            kk = kk + 1;
        end
    end
 plot(V_T3_plot,Gamma_T3_plot,...
            '-.','color',mycolor(3,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 83.806 K");  
%  dashline(V_T3_plot,Gamma_T3_plot,70,7,70,7,'color',mycolor(3,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 83.806 K"); 
        hold on
        
      T4 = 150;
    V_T4 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T4)
        props_T4 = computeThermoPropsTVForResult(T4, V_T4(pp));
        props_T4_check = computeThermoPropsForResult(T4, pmelt(T4));          
        if props_T4(1)>0&&V_T4(pp)<props_T4_check(1)
            Gamma_T4_plot(kk) = props_T4(7);
            V_T4_plot(kk) = V_T4(pp);
            kk = kk + 1;
        end
    end
 plot(V_T4_plot,Gamma_T4_plot,...
            ':','color',mycolor(4,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 150 K");
%   dashline(V_T4_plot,Gamma_T4_plot,2, 2, 1, 1,'color',mycolor(4,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 150 K");
       
     hold on   
     T5 = 300;
    V_T5 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T5)
        props_T5 = computeThermoPropsTVForResult(T5, V_T5(pp));  
        props_T5_check = computeThermoPropsForResult(T5, pmelt(T5));  
        if props_T5(1)>0&&V_T5(pp)<props_T5_check(1)
            Gamma_T5_plot(kk) = props_T5(7);
            
            V_T5_plot(kk) = V_T5(pp);
            kk = kk + 1;
        end
    end

%  dashline(V_T5_plot,Gamma_T5_plot,1,1,1,1,'color',mycolor(5,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 300 K");
  plot(V_T5_plot,Gamma_T5_plot,...
            '-','color',mycolor(5,:),'linewidth',linewidth*2,'DisplayName',"\itT\rm = 300 K");
    
 T6 = [0.0001	2	4	6	8	10	12	14	16	18	20	22	24	26	28	30	32	34	36	38	40	42	44	46	48	50	52	54	56	58	60	62	64	66	68	70	72	74	76	78	80	82	83.806];
  for ii = 1:length(T6)

          p6(ii) = psub(T6(ii));
      props_T6 = computeThermoPropsForResult(T6(ii), p6(ii));
      Gamma_T6_plot(ii) = props_T6(7);
      V_T6_plot(ii) = props_T6(1);
  end
%  dashline(V_T6_plot,Gamma_T6_plot,6,6,6,6,'color',mycolor(6,:),'linewidth',linewidth,'DisplayName',"sublimation");   
  plot(V_T6_plot,Gamma_T6_plot,...
            '--','color',mycolor(6,:),'linewidth',linewidth*2,'DisplayName',"sublimation");
 T7 = [83.806	84	86	88	90	92	94	96	98	100	102	104	106	108	110	112	114	116	118	120	122	124	126	128	130	132	134	136	138	140	142	144	146	148	150	152	154	156	158	160	162	164	166	168	170	172	174	176	178	180	182	184	186	188	190	192	194	196	198	200	202	204	206	208	210	212	214	216	218	220	222	224	226	228	230	232	234	236	238	240	242	244	246	248	250	252	254	256	258	260	262	264	266	268	270	272	274	276	278	280	282	284	286	288	290	292	294	296	298	300];
V_T7_plot=[];
 for ii = 1:length(T7)

          p7(ii) = pmelt(T7(ii));
      props_T7 = computeThermoPropsForResult(T7(ii), p7(ii));
       Gamma_T7_plot(ii) = props_T7(7);
      V_T7_plot(ii) = props_T7(1);
  end
%  dashline(V_T7_plot,Gamma_T7_plot,7,1,7,1,'color',mycolor(7,:),'linewidth',linewidth,'DisplayName',"melting");
  plot(V_T7_plot,Gamma_T7_plot,...
            '-.','color',mycolor(7,:),'linewidth',linewidth*2,'DisplayName',"melting");
%          legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
ylim([2.53 2.7]);
xlim([18 25]);
ytickformat('%.2f');
    xlabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
ylabel('\it \gamma','fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
        set(gcf,'position',[500,200,500,340])
        set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../../paper writing/figure/extrapolation/')
        mkdir('../../paper writing/figure/extrapolation/')
    end
    print(gcf,'-dtiff','-r900','../../paper writing/figure/extrapolation/gamma.bmp')       
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 13
    fontsize = 12;    markersize = 7;    linewidth = 1;
    figure(1)
    T1 = 0.0001;
    V_T1 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T1)
        props_T1 = computeThermoPropsTVForResult(T1, V_T1(pp));
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
        props_T2 = computeThermoPropsTVForResult(T2, V_T2(pp));
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
        props_T3 = computeThermoPropsTVForResult(T3, V_T3(pp));
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
        props_T4 = computeThermoPropsTVForResult(T4, V_T4(pp));
        props_T4_check = computeThermoPropsForResult(T4, pmelt(T4));          
        if props_T4(1)>0&&V_T4(pp)<props_T4_check(1)
            p_T4_plot(kk) = props_T4(1);
            V_T4_plot(kk) = V_T4(pp);
            kk = kk + 1;
        end
    end
%  plot(V_T4_plot,p_T4_plot,...
%             ':','color',mycolor(4,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 150 K");  
  plot(V_T4_plot,p_T4_plot,...
            ':','color',mycolor(4,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 150 K");       
     hold on   
     T5 = 300;
    V_T5 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T5)
        props_T5 = computeThermoPropsTVForResult(T5, V_T5(pp));  
        props_T5_check = computeThermoPropsForResult(T5, pmelt(T5));  
        if props_T5(1)>0&&V_T5(pp)<props_T5_check(1)
            p_T5_plot(kk) = props_T5(1);
            
            V_T5_plot(kk) = V_T5(pp);
            kk = kk + 1;
        end
    end

%  dashline(V_T5_plot,p_T5_plot,4,1,2,1,'color',mycolor(5,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 300 K");
  plot(V_T5_plot,p_T5_plot,...
            '-','color',mycolor(5,:),'linewidth',linewidth*2,'DisplayName',"\itT\rm = 300 K"); 
    
 T6 = [0.0001	2	4	6	8	10	12	14	16	18	20	22	24	26	28	30	32	34	36	38	40	42	44	46	48	50	52	54	56	58	60	62	64	66	68	70	72	74	76	78	80	82	83.806];
  for ii = 1:length(T6)

          p6(ii) = psub(T6(ii));
      props_T6 = computeThermoPropsForResult(T6(ii), p6(ii));
      V_T6_plot(ii) = props_T6(1);
  end
%  dashline(V_T6_plot,p6,6,6,6,6,'color',mycolor(6,:),'linewidth',linewidth,'DisplayName',"sublimation"); 
   plot(V_T6_plot,p6,...
            '--','color',mycolor(6,:),'linewidth',linewidth*2,'DisplayName',"sublimation");
 T7 = [83.806	84	86	88	90	92	94	96	98	100	102	104	106	108	110	112	114	116	118	120	122	124	126	128	130	132	134	136	138	140	142	144	146	148	150	152	154	156	158	160	162	164	166	168	170	172	174	176	178	180	182	184	186	188	190	192	194	196	198	200	202	204	206	208	210	212	214	216	218	220	222	224	226	228	230	232	234	236	238	240	242	244	246	248	250	252	254	256	258	260	262	264	266	268	270	272	274	276	278	280	282	284	286	288	290	292	294	296	298	300];
   for ii = 1:length(T7)

          p7(ii) = pmelt(T7(ii));
      props_T7 = computeThermoPropsForResult(T7(ii), p7(ii));
      V_T7_plot(ii) = props_T7(1);
  end
%  dashline(V_T7_plot,p7,7,1,7,1,'color',mycolor(7,:),'linewidth',linewidth,'DisplayName',"melting");
  plot(V_T7_plot,p7,...
            '-.','color',mycolor(7,:),'linewidth',linewidth*2,'DisplayName',"melting");
%          legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
xlim([18 25]);
    xlabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    ylabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth); 
    annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
        set(gcf,'position',[500,200,300,260])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../../paper writing/figure/extrapolation/')
        mkdir('../../paper writing/figure/extrapolation/')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/extrapolation/pV.png')
    
    figure(2)
    T1 = 0.0001;
    V_T1 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T1)
        props_T1 = computeThermoPropsTVForResult(T1, V_T1(pp));
        if props_T1(1)>0
            BetaS_T1_plot(kk) = props_T1(3);
            V_T1_plot(kk) = V_T1(pp);
            kk = kk + 1;
        end
    end
 plot(V_T1_plot,1./BetaS_T1_plot,...
            'color',mycolor(1,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 0 K");
  hold on  
     T2 = 50;
    V_T2 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T2)
        props_T2 = computeThermoPropsTVForResult(T2, V_T2(pp));
        if props_T2(1)>0
            BetaS_T2_plot(kk) = props_T2(3);
            V_T2_plot(kk) = V_T2(pp);
            kk = kk + 1;
        end
    end
 plot(V_T2_plot,1./BetaS_T2_plot,...
            '--','color',mycolor(2,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 50 K"); 
     hold on   
     
     T3 = 83.806;
    V_T3 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T3)
        props_T3 = computeThermoPropsTVForResult(T3, V_T3(pp));
        if props_T3(1)>0
            BetaS_T3_plot(kk) = props_T3(3);
            V_T3_plot(kk) = V_T3(pp);
            kk = kk + 1;
        end
    end
 plot(V_T3_plot,1./BetaS_T3_plot,...
            '-.','color',mycolor(3,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 83.806 K");  
        hold on
        
      T4 = 150;
    V_T4 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T4)
        props_T4 = computeThermoPropsTVForResult(T4, V_T4(pp));
        props_T4_check = computeThermoPropsForResult(T4, pmelt(T4));          
        if props_T4(1)>0&&V_T4(pp)<props_T4_check(1)
            BetaS_T4_plot(kk) = props_T4(3);
            V_T4_plot(kk) = V_T4(pp);
            kk = kk + 1;
        end
    end
 plot(V_T4_plot,1./BetaS_T4_plot,...
            ':','color',mycolor(4,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 150 K");  
     hold on   
     T5 = 300;
    V_T5 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T5)
        props_T5 = computeThermoPropsTVForResult(T5, V_T5(pp));  
        props_T5_check = computeThermoPropsForResult(T5, pmelt(T5));  
        if props_T5(1)>0&&V_T5(pp)<props_T5_check(1)
            BetaS_T5_plot(kk) = props_T5(3);
            
            V_T5_plot(kk) = V_T5(pp);
            kk = kk + 1;
        end
    end

%  dashline(V_T5_plot,BetaS_T5_plot*10^4,4,1,2,1,'color',mycolor(5,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 300 K");
 plot(V_T5_plot,1./BetaS_T5_plot,...
            '-','color',mycolor(5,:),'linewidth',linewidth*2,'DisplayName',"\itT\rm = 300 K"); 
    V_T7_plot = [];
 T6 = [0.0001	2	4	6	8	10	12	14	16	18	20	22	24	26	28	30	32	34	36	38	40	42	44	46	48	50	52	54	56	58	60	62	64	66	68	70	72	74	76	78	80	82	83.806];
  for ii = 1:length(T6)

          p6(ii) = psub(T6(ii));
      props_T6 = computeThermoPropsForResult(T6(ii), p6(ii));
      BetaS_T6_plot(ii) = props_T6(3);
      V_T6_plot(ii) = props_T6(1);
  end
%  dashline(V_T6_plot,BetaS_T6_plot*10^4,6,6,6,6,'color',mycolor(6,:),'linewidth',linewidth,'DisplayName',"sublimation");   
  plot(V_T6_plot,1./BetaS_T6_plot,...
            '--','color',mycolor(6,:),'linewidth',linewidth*2,'DisplayName',"sublimation"); 
 
 
 T7 = [83.806	84	86	88	90	92	94	96	98	100	102	104	106	108	110	112	114	116	118	120	122	124	126	128	130	132	134	136	138	140	142	144	146	148	150	152	154	156	158	160	162	164	166	168	170	172	174	176	178	180	182	184	186	188	190	192	194	196	198	200	202	204	206	208	210	212	214	216	218	220	222	224	226	228	230	232	234	236	238	240	242	244	246	248	250	252	254	256	258	260	262	264	266	268	270	272	274	276	278	280	282	284	286	288	290	292	294	296	298	300];
   for ii = 1:length(T7)

          p7(ii) = pmelt(T7(ii));
      props_T7 = computeThermoPropsForResult(T7(ii), p7(ii));
       BetaS_T7_plot(ii) = props_T7(3);
      V_T7_plot(ii) = props_T7(1);
  end
%  dashline(V_T7_plot,BetaS_T7_plot*10^4,7,1,7,1,'color',mycolor(7,:),'linewidth',linewidth,'DisplayName',"melting");
   plot(V_T7_plot,1./BetaS_T7_plot,...
            '-.','color',mycolor(7,:),'linewidth',linewidth*2,'DisplayName',"melting");
%          legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
xlim([18 25]);
    xlabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
ylabel('\it K_S \rm\rm / MPa','fontsize',fontsize,'linewidth',linewidth);
    annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfb','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
        set(gcf,'position',[500,200,300,260])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    print(gcf,'-dtiff','-r300','../../paper writing/figure/extrapolation/KS.png')
%     print(gcf,'-dtiff','-r300','../solid data fitting/figure/evaluation/KS.png')    

    figure(3)   
     T1 = 0.0001;
    V_T1 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T1)
        props_T1 = computeThermoPropsTVForResult(T1, V_T1(pp));
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
        props_T2 = computeThermoPropsTVForResult(T2, V_T2(pp));
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
        props_T3 = computeThermoPropsTVForResult(T3, V_T3(pp));
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
        props_T4 = computeThermoPropsTVForResult(T4, V_T4(pp));
        props_T4_check = computeThermoPropsForResult(T4, pmelt(T4));          
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
        props_T5 = computeThermoPropsTVForResult(T5, V_T5(pp));  
        props_T5_check = computeThermoPropsForResult(T5, pmelt(T5));  
        if props_T5(1)>0&&V_T5(pp)<props_T5_check(1)
            alpha_T5_plot(kk) = props_T5(4);
            
            V_T5_plot(kk) = V_T5(pp);
            kk = kk + 1;
        end
    end

%  dashline(V_T5_plot,alpha_T5_plot*10^4,4,1,2,1,'color',mycolor(5,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 300 K");
   plot(V_T5_plot,alpha_T5_plot*10^4,...
            '-','color',mycolor(5,:),'linewidth',linewidth*2,'DisplayName',"\itT\rm = 300 K");
    
 T6 = [0.0001	2	4	6	8	10	12	14	16	18	20	22	24	26	28	30	32	34	36	38	40	42	44	46	48	50	52	54	56	58	60	62	64	66	68	70	72	74	76	78	80	82	83.806];
  for ii = 1:length(T6)

          p6(ii) = psub(T6(ii));
      props_T6 = computeThermoPropsForResult(T6(ii), p6(ii));
      alpha_T6_plot(ii) = props_T6(4);
      V_T6_plot(ii) = props_T6(1);
  end
%  dashline(V_T6_plot,alpha_T6_plot*10^4,6,6,6,6,'color',mycolor(6,:),'linewidth',linewidth,'DisplayName',"sublimation");   
  plot(V_T6_plot,alpha_T6_plot*10^4,...
            '--','color',mycolor(6,:),'linewidth',linewidth*2,'DisplayName',"sublimation"); 
 T7 = [83.806	84	86	88	90	92	94	96	98	100	102	104	106	108	110	112	114	116	118	120	122	124	126	128	130	132	134	136	138	140	142	144	146	148	150	152	154	156	158	160	162	164	166	168	170	172	174	176	178	180	182	184	186	188	190	192	194	196	198	200	202	204	206	208	210	212	214	216	218	220	222	224	226	228	230	232	234	236	238	240	242	244	246	248	250	252	254	256	258	260	262	264	266	268	270	272	274	276	278	280	282	284	286	288	290	292	294	296	298	300];
   for ii = 1:length(T7)

          p7(ii) = pmelt(T7(ii));
      props_T7 = computeThermoPropsForResult(T7(ii), p7(ii));
       alpha_T7_plot(ii) = props_T7(4);
      V_T7_plot(ii) = props_T7(1);
  end
%  dashline(V_T7_plot,alpha_T7_plot*10^4,7,1,7,1,'color',mycolor(7,:),'linewidth',linewidth,'DisplayName',"melting");
  plot(V_T7_plot,alpha_T7_plot*10^4,...
            '-.','color',mycolor(7,:),'linewidth',linewidth*2,'DisplayName',"melting");
%          legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
xlim([18 25]);
    xlabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
ylabel(['\it \alpha \rm\cdot10^4\rm / K^',char(hex2dec('2212')),'^1'],'fontsize',fontsize,'linewidth',linewidth);
    annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfc','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
        set(gcf,'position',[500,200,300,260])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
%     if ~isfolder('../solid data fitting/figure/evaluation')
%         mkdir('../solid data fitting/figure/evaluation')
%     end
%     print(gcf,'-dtiff','-r300','../solid data fitting/figure/evaluation/alpha.png') 
    print(gcf,'-dtiff','-r300','../../paper writing/figure/extrapolation/alpha.png')

     figure(4)
    T1 = 0.0001;
    V_T1 = 18:0.01:25;
    kk = 1;
    for pp = 1:length(V_T1)
        props_T1 = computeThermoPropsTVForResult(T1, V_T1(pp));
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
        props_T2 = computeThermoPropsTVForResult(T2, V_T2(pp));
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
        props_T3 = computeThermoPropsTVForResult(T3, V_T3(pp));
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
        props_T4 = computeThermoPropsTVForResult(T4, V_T4(pp));
        props_T4_check = computeThermoPropsForResult(T4, pmelt(T4));          
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
        props_T5 = computeThermoPropsTVForResult(T5, V_T5(pp));  
        props_T5_check = computeThermoPropsForResult(T5, pmelt(T5));  
        if props_T5(1)>0&&V_T5(pp)<props_T5_check(1)
            Beta_T5_plot(kk) = props_T5(4)/props_T5(2);
            
            V_T5_plot(kk) = V_T5(pp);
            kk = kk + 1;
        end
    end

%  dashline(V_T5_plot,Beta_T5_plot,4,1,2,1,'color',mycolor(5,:),'linewidth',linewidth,'DisplayName',"\itT\rm = 300 K");
  plot(V_T5_plot,Beta_T5_plot,...
            '-','color',mycolor(5,:),'linewidth',linewidth*2,'DisplayName',"\itT\rm = 300 K"); 
    
 T6 = [0.0001	2	4	6	8	10	12	14	16	18	20	22	24	26	28	30	32	34	36	38	40	42	44	46	48	50	52	54	56	58	60	62	64	66	68	70	72	74	76	78	80	82	83.806];
  for ii = 1:length(T6)

          p6(ii) = psub(T6(ii));
      props_T6 = computeThermoPropsForResult(T6(ii), p6(ii));
      Beta_T6_plot(ii) = props_T6(4)/props_T6(2);
      V_T6_plot(ii) = props_T6(1);
  end
%  dashline(V_T6_plot,Beta_T6_plot,6,6,6,6,'color',mycolor(6,:),'linewidth',linewidth,'DisplayName',"sublimation");
   plot(V_T6_plot,Beta_T6_plot,...
            '--','color',mycolor(6,:),'linewidth',linewidth*2,'DisplayName',"sublimation");
 T7 = [83.806	84	86	88	90	92	94	96	98	100	102	104	106	108	110	112	114	116	118	120	122	124	126	128	130	132	134	136	138	140	142	144	146	148	150	152	154	156	158	160	162	164	166	168	170	172	174	176	178	180	182	184	186	188	190	192	194	196	198	200	202	204	206	208	210	212	214	216	218	220	222	224	226	228	230	232	234	236	238	240	242	244	246	248	250	252	254	256	258	260	262	264	266	268	270	272	274	276	278	280	282	284	286	288	290	292	294	296	298	300];
   for ii = 1:length(T7)

          p7(ii) = pmelt(T7(ii));
      props_T7 = computeThermoPropsForResult(T7(ii), p7(ii));
       Beta_T7_plot(ii) = props_T7(4)/props_T7(2);
      V_T7_plot(ii) = props_T7(1);
  end
%  dashline(V_T7_plot,Beta_T7_plot,7,1,7,1,'color',mycolor(7,:),'linewidth',linewidth,'DisplayName',"melting");
  plot(V_T7_plot,Beta_T7_plot,...
            '-.','color',mycolor(7,:),'linewidth',linewidth*2,'DisplayName',"melting");
%          legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);

    xlabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
ylabel(['\it \beta \rm / MPa^',char(hex2dec('2212')),'^1'],'fontsize',fontsize,'linewidth',linewidth);
    annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfd','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
        set(gcf,'position',[500,200,300,260])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
%     if ~isfolder('../solid data fitting/figure/evaluation')
%         mkdir('../solid data fitting/figure/evaluation')
%     end
%     print(gcf,'-dtiff','-r300','../solid data fitting/figure/evaluation/thermal pressure coefficient.png')  
   print(gcf,'-dtiff','-r300','../../paper writing/figure/extrapolation/thermal pressure coefficient.png')
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 14
        fontsize = 12;    markersize = 7;    linewidth = 1;

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

figure(1)
%  my_dashline=[0 0 0 0;4 1 1 2;2 0.5 0.5 1;6 6 6 6;3 2 3 2;1 1 1 1;7 1 7 1;13 3 13 3;1 4 1 4];
for ii = 1:a1_Vm_highp
    for jj=1:a2_Vm_highp
        if Nystart_Vm_highp(ii,jj)~=0
            ll1_Vm_highp(ii,jj)=string(strcat(Author_Vm_highp(Nystart_Vm_highp(ii,jj)),Year_Vm_highp(Nystart_Vm_highp(ii,jj)),' - ',{' '},char(string(T_Vm_highp(Nystart_Vm_highp(ii,jj)))),{' '},'K'));            
            plot(p_Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)),Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)),...
            mymarker{jj},'markersize',markersize,'color',mycolor(jj,:),'linewidth',linewidth,'DisplayName',ll1_Vm_highp(ii,jj)); 
        hold on;

            T_plot_Vm_highp = mean(T_Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)));
            pp=1;
            p_plot_Vm_highp=[];
            Vm_p_plot_Vm_highp=[];  
            p_plot_Vm_highp_min=max(min(p_Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)))-2,0);
            p_plot_Vm_highp_max=max(p_Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)))+30;     
            for jjk=p_plot_Vm_highp_min:1:p_plot_Vm_highp_max
                p_plot_Vm_highp(pp,:)=jjk;
                props_Vm_highp(pp,:)=computeThermoPropsForResult(T_plot_Vm_highp, p_plot_Vm_highp(pp,:));
    %             benzeneprops(T_plot,p_plot_Vm_highp(pp,:));
                 Vm_p_plot_Vm_highp(pp,:)=props_Vm_highp(pp,1);
                pp=pp+1;
            end        
            if jj == 1
                plot(p_plot_Vm_highp,Vm_p_plot_Vm_highp,'color',mycolor(jj,:),'linewidth',linewidth);
            elseif jj == 2
                plot(p_plot_Vm_highp,Vm_p_plot_Vm_highp,...
                '--','color',mycolor(jj,:),'linewidth',linewidth);    
            elseif jj == 3
                plot(p_plot_Vm_highp,Vm_p_plot_Vm_highp,...
                '-.','color',mycolor(jj,:),'linewidth',linewidth);   
            elseif jj == 4
                plot(p_plot_Vm_highp,Vm_p_plot_Vm_highp,...
                ':','color',mycolor(jj,:),'linewidth',linewidth);
            elseif jj == 5
                plot(p_plot_Vm_highp,Vm_p_plot_Vm_highp,...
                '-','color',mycolor(jj,:),'linewidth',linewidth*2);      
            elseif jj == 6
                plot(p_plot_Vm_highp,Vm_p_plot_Vm_highp,...
                '--','color',mycolor(jj,:),'linewidth',linewidth*2);  
            elseif jj == 7
                plot(p_plot_Vm_highp,Vm_p_plot_Vm_highp,...
                '-.','color',mycolor(jj,:),'linewidth',linewidth*2);              
            end
        end
    end    
end  
        ytickformat('%.1f');  
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth); 
    ylabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
        set(gcf,'position',[500,200,300,260])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../../paper writing/figure/cell volume/high pressure/')
        mkdir('../../paper writing/figure/cell volume/high pressure/')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/cell volume/high pressure/high pressure_1.png')  
    
    
    figure(2)
 for ii = 1:a1_Vm_highp
    for jj=1:a2_Vm_highp
        Vm_highp_dev = [];
        if Nystart_Vm_highp(ii,jj)~=0
        ll1_Vm_highp(ii,jj)=string(strcat(Author_Vm_highp(Nystart_Vm_highp(ii,jj)),Year_Vm_highp(Nystart_Vm_highp(ii,jj)),' - ',{' '},char(string(T_Vm_highp(Nystart_Vm_highp(ii,jj)))),{' '},'K'));            
        for  kkk = 1:length(p_Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)))
            props_Vm_highp = computeThermoPropsForResult(T_Vm_highp(Nystart_Vm_highp(ii,jj)+kkk-1), p_Vm_highp(Nystart_Vm_highp(ii,jj)+kkk-1));
            Vm_highp_dev(kkk,:) = 100*(Vm_highp(Nystart_Vm_highp(ii,jj)+kkk-1) - props_Vm_highp(1)) ./ Vm_highp(Nystart_Vm_highp(ii,jj)+kkk-1);
        end
        plot(p_Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)),Vm_highp_dev,...
        mymarker{jj},'markersize',markersize,'color',mycolor(jj,:),'linewidth',linewidth,'DisplayName',ll1_Vm_highp(ii,jj)); 
    hold on;
        end
    end    
 end
        c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
        ytickformat('%.1f'); 
        xlim([0 220]);
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth);
        annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfb','FontSize',fontsize,'FontName','Times');

        set(gcf,'position',[500,200,300,260])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
                ylabel(['100\cdot',char(8201),'(\itV\rm_m_,_e_x_p',char(hex2dec('2212')),'\it V\rm_m_,_c_a_l_c) /\it V\rm_m_,_e_x_p'],'fontsize',fontsize-1,'linewidth',linewidth)  

    if ~isfolder('../../paper writing/figure/cell volume/high pressure/')
        mkdir('../../paper writing/figure/cell volume/high pressure/')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/cell volume/high pressure/high pressure_dev_1.png')  
    
    
    
    [Year_Vm_highp2,Author_Vm_highp2,T_Vm_highp2,Vm_highp2,p_Vm_highp2] = textread('../evaluation data/cell volume/Vm_high_pressure_2.txt','%s%s%f%f%f','headerlines',2);
ndata_Vm_highp2=length(T_Vm_highp2);xxcount_Vm_highp2=1; Nxstart_Vm_highp2 = 1;AA_Vm_highp2=Author_Vm_highp2(1);YY_Vm_highp2=Year_Vm_highp2(1);
for ii = 2:ndata_Vm_highp2         
        if ~strcmp(AA_Vm_highp2,Author_Vm_highp2(ii))||~strcmp(YY_Vm_highp2,Year_Vm_highp2(ii))
            Nxend_Vm_highp2(xxcount_Vm_highp2) = ii-1;
            xxcount_Vm_highp2 = xxcount_Vm_highp2 + 1;
            Nxstart_Vm_highp2(xxcount_Vm_highp2) = ii;
            AA_Vm_highp2=Author_Vm_highp2(ii);YY_Vm_highp2=Year_Vm_highp2(ii);
        end
        if ii == ndata_Vm_highp2
            Nxend_Vm_highp2(xxcount_Vm_highp2) = ndata_Vm_highp2;
        end
end


for jj=1:xxcount_Vm_highp2
    ndata_Vm_highp2=length(T_Vm_highp2(Nxstart_Vm_highp2(jj):Nxend_Vm_highp2(jj)));yycount_Vm_highp2=1; Nystart_Vm_highp2(jj,yycount_Vm_highp2) = Nxstart_Vm_highp2(jj);TT_Vm_highp2=T_Vm_highp2(Nxstart_Vm_highp2(jj));
    for ii = 2:ndata_Vm_highp2   

            %if ~strcmp(AA,Author(Nxstart(jj)+ii-1))&&~strcmp(YY,Year(Nxstart(jj)+ii-1))
            if abs(T_Vm_highp2(Nxstart_Vm_highp2(jj)+ii-1) - TT_Vm_highp2)> 0.05
                Nyend_Vm_highp2(jj,yycount_Vm_highp2) = Nxstart_Vm_highp2(jj)+ii-2;
                yycount_Vm_highp2 = yycount_Vm_highp2 + 1;
                Nystart_Vm_highp2(jj,yycount_Vm_highp2) = Nxstart_Vm_highp2(jj)+ii-1;
                TT_Vm_highp2=T_Vm_highp2(Nxstart_Vm_highp2(jj)+ii-1);
            end
            if ii == ndata_Vm_highp2
                Nyend_Vm_highp2(jj,yycount_Vm_highp2) = ndata_Vm_highp2 + Nxstart_Vm_highp2(jj)-1;
            end
    end     
end
[a1_Vm_highp2 a2_Vm_highp2]=size(Nystart_Vm_highp2);

figure(3)
%  my_dashline=[0 0 0 0;4 1 1 2;2 0.5 0.5 1;6 6 6 6;3 2 3 2;1 1 1 1;7 1 7 1;13 3 13 3;1 4 1 4];
for ii = 1:a1_Vm_highp2
    for jj=1:a2_Vm_highp2
        if Nystart_Vm_highp2(ii,jj)~=0
            ll1_Vm_highp2(ii,jj)=string(strcat(Author_Vm_highp2(Nystart_Vm_highp2(ii,jj)),Year_Vm_highp2(Nystart_Vm_highp2(ii,jj)),' - ',{' '},char(string(T_Vm_highp2(Nystart_Vm_highp2(ii,jj)))),{' '},'K'));            
            plot(p_Vm_highp2(Nystart_Vm_highp2(ii,jj):Nyend_Vm_highp2(ii,jj)),Vm_highp2(Nystart_Vm_highp2(ii,jj):Nyend_Vm_highp2(ii,jj)),...
            mymarker{jj},'markersize',markersize,'color',mycolor(jj,:),'linewidth',linewidth,'DisplayName',ll1_Vm_highp2(ii,jj)); 
        hold on;

            T_plot_Vm_highp2 = mean(T_Vm_highp2(Nystart_Vm_highp2(ii,jj):Nyend_Vm_highp2(ii,jj)));
            pp=1;
            p_plot_Vm_highp2=[];
            Vm_p_plot_Vm_highp2=[];  
            p_plot_Vm_highp2_min=max(min(p_Vm_highp2(Nystart_Vm_highp2(ii,jj):Nyend_Vm_highp2(ii,jj)))-2,0);
            p_plot_Vm_highp2_max=max(p_Vm_highp2(Nystart_Vm_highp2(ii,jj):Nyend_Vm_highp2(ii,jj)))+30;     
            for jjk=p_plot_Vm_highp2_min:1:p_plot_Vm_highp2_max
                p_plot_Vm_highp2(pp,:)=jjk;
                props_Vm_highp2(pp,:)=computeThermoPropsForResult(T_plot_Vm_highp2, p_plot_Vm_highp2(pp,:));
    %             benzeneprops(T_plot,p_plot_Vm_highp2(pp,:));
                 Vm_p_plot_Vm_highp2(pp,:)=props_Vm_highp2(pp,1);
                pp=pp+1;
            end        
            if jj == 1
                plot(p_plot_Vm_highp2,Vm_p_plot_Vm_highp2,'color',mycolor(jj,:),'linewidth',linewidth);
            elseif jj == 2
                plot(p_plot_Vm_highp2,Vm_p_plot_Vm_highp2,...
                '--','color',mycolor(jj,:),'linewidth',linewidth);    
            elseif jj == 3
                plot(p_plot_Vm_highp2,Vm_p_plot_Vm_highp2,...
                '-.','color',mycolor(jj,:),'linewidth',linewidth);   
            elseif jj == 4
                plot(p_plot_Vm_highp2,Vm_p_plot_Vm_highp2,...
                ':','color',mycolor(jj,:),'linewidth',linewidth);
            elseif jj == 5
                plot(p_plot_Vm_highp2,Vm_p_plot_Vm_highp2,...
                '-','color',mycolor(jj,:),'linewidth',linewidth*2);      
            elseif jj == 6
                plot(p_plot_Vm_highp2,Vm_p_plot_Vm_highp2,...
                '--','color',mycolor(jj,:),'linewidth',linewidth*2);  
            elseif jj == 7
                plot(p_plot_Vm_highp2,Vm_p_plot_Vm_highp2,...
                '-.','color',mycolor(jj,:),'linewidth',linewidth*2);              
            end
        end
    end    
end  
%         ytickformat('%.1f');  
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth); 
    ylabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfc','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
        set(gcf,'position',[500,200,300,260])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../../paper writing/figure/cell volume/high pressure/')
        mkdir('../../paper writing/figure/cell volume/high pressure/')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/cell volume/high pressure/high pressure_2.png')  
    
    
    figure(4)
 for ii = 1:a1_Vm_highp2
    for jj=1:a2_Vm_highp2
        Vm_highp2_dev = [];
        if Nystart_Vm_highp2(ii,jj)~=0
        ll1_Vm_highp2(ii,jj)=string(strcat(Author_Vm_highp2(Nystart_Vm_highp2(ii,jj)),Year_Vm_highp2(Nystart_Vm_highp2(ii,jj)),' - ',{' '},char(string(T_Vm_highp2(Nystart_Vm_highp2(ii,jj)))),{' '},'K'));            
        for  kkk = 1:length(p_Vm_highp2(Nystart_Vm_highp2(ii,jj):Nyend_Vm_highp2(ii,jj)))
            props_Vm_highp2 = computeThermoPropsForResult(T_Vm_highp2(Nystart_Vm_highp2(ii,jj)+kkk-1), p_Vm_highp2(Nystart_Vm_highp2(ii,jj)+kkk-1));
            Vm_highp2_dev(kkk,:) = 100*(Vm_highp2(Nystart_Vm_highp2(ii,jj)+kkk-1) - props_Vm_highp2(1)) ./ Vm_highp2(Nystart_Vm_highp2(ii,jj)+kkk-1);
        end
        plot(p_Vm_highp2(Nystart_Vm_highp2(ii,jj):Nyend_Vm_highp2(ii,jj)),Vm_highp2_dev,...
        mymarker{jj},'markersize',markersize,'color',mycolor(jj,:),'linewidth',linewidth,'DisplayName',ll1_Vm_highp2(ii,jj)); 
    hold on;
        end
    end    
 end
        c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
        ytickformat('%.1f'); 
        xlim([0 500]);
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth);
%     ylabel(['100\cdot',char(8201),'(\itV\rm_m_,_e_x_p',char(hex2dec('2212')),'\it V\rm_m_,_c_a_l_c) /\it V\rm_m_,_c_a_l_c'],'fontsize',fontsize,'linewidth',linewidth)  
        annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfd','FontSize',fontsize,'FontName','Times');
        set(gcf,'position',[500,200,300,260])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
                    ylabel(['100\cdot',char(8201),'(\itV\rm_m_,_e_x_p',char(hex2dec('2212')),'\it V\rm_m_,_c_a_l_c) /\it V\rm_m_,_e_x_p'],'fontsize',fontsize-1,'linewidth',linewidth)  

    if ~isfolder('../../paper writing/figure/cell volume/high pressure/')
        mkdir('../../paper writing/figure/cell volume/high pressure/')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/cell volume/high pressure/high pressure_dev_2.png')
    [Year_Vm_highp3,Author_Vm_highp3,T_Vm_highp3,Vm_highp3,p_Vm_highp3] = textread('../evaluation data/cell volume/Vm_high_pressure_3.txt','%s%s%f%f%f','headerlines',2);
ndata_Vm_highp3=length(T_Vm_highp3);xxcount_Vm_highp3=1; Nxstart_Vm_highp3 = 1;AA_Vm_highp3=Author_Vm_highp3(1);YY_Vm_highp3=Year_Vm_highp3(1);
for ii = 2:ndata_Vm_highp3         
        if ~strcmp(AA_Vm_highp3,Author_Vm_highp3(ii))||~strcmp(YY_Vm_highp3,Year_Vm_highp3(ii))
            Nxend_Vm_highp3(xxcount_Vm_highp3) = ii-1;
            xxcount_Vm_highp3 = xxcount_Vm_highp3 + 1;
            Nxstart_Vm_highp3(xxcount_Vm_highp3) = ii;
            AA_Vm_highp3=Author_Vm_highp3(ii);YY_Vm_highp3=Year_Vm_highp3(ii);
        end
        if ii == ndata_Vm_highp3
            Nxend_Vm_highp3(xxcount_Vm_highp3) = ndata_Vm_highp3;
        end
end


for jj=1:xxcount_Vm_highp3
    ndata_Vm_highp3=length(T_Vm_highp3(Nxstart_Vm_highp3(jj):Nxend_Vm_highp3(jj)));yycount_Vm_highp3=1; Nystart_Vm_highp3(jj,yycount_Vm_highp3) = Nxstart_Vm_highp3(jj);TT_Vm_highp3=T_Vm_highp3(Nxstart_Vm_highp3(jj));
    for ii = 2:ndata_Vm_highp3   

            %if ~strcmp(AA,Author(Nxstart(jj)+ii-1))&&~strcmp(YY,Year(Nxstart(jj)+ii-1))
            if abs(T_Vm_highp3(Nxstart_Vm_highp3(jj)+ii-1) - TT_Vm_highp3)> 0.05
                Nyend_Vm_highp3(jj,yycount_Vm_highp3) = Nxstart_Vm_highp3(jj)+ii-2;
                yycount_Vm_highp3 = yycount_Vm_highp3 + 1;
                Nystart_Vm_highp3(jj,yycount_Vm_highp3) = Nxstart_Vm_highp3(jj)+ii-1;
                TT_Vm_highp3=T_Vm_highp3(Nxstart_Vm_highp3(jj)+ii-1);
            end
            if ii == ndata_Vm_highp3
                Nyend_Vm_highp3(jj,yycount_Vm_highp3) = ndata_Vm_highp3 + Nxstart_Vm_highp3(jj)-1;
            end
    end     
end
[a1_Vm_highp3 a2_Vm_highp3]=size(Nystart_Vm_highp3);

figure(5)
%  my_dashline=[0 0 0 0;4 1 1 2;2 0.5 0.5 1;6 6 6 6;3 2 3 2;1 1 1 1;7 1 7 1;13 3 13 3;1 4 1 4];
for ii = 1:a1_Vm_highp3
    for jj=1:a2_Vm_highp3
        if Nystart_Vm_highp3(ii,jj)~=0
            ll1_Vm_highp3(ii,jj)=string(strcat(Author_Vm_highp3(Nystart_Vm_highp3(ii,jj)),Year_Vm_highp3(Nystart_Vm_highp3(ii,jj)),' - ',{' '},char(string(T_Vm_highp3(Nystart_Vm_highp3(ii,jj)))),{' '},'K'));            
            plot(p_Vm_highp3(Nystart_Vm_highp3(ii,jj):Nyend_Vm_highp3(ii,jj)),Vm_highp3(Nystart_Vm_highp3(ii,jj):Nyend_Vm_highp3(ii,jj)),...
            mymarker{jj},'markersize',markersize,'color',mycolor(jj,:),'linewidth',linewidth,'DisplayName',ll1_Vm_highp3(ii,jj)); 
        hold on;

            T_plot_Vm_highp3 = mean(T_Vm_highp3(Nystart_Vm_highp3(ii,jj):Nyend_Vm_highp3(ii,jj)));
            pp=1;
            p_plot_Vm_highp3=[];
            Vm_p_plot_Vm_highp3=[];  
            p_plot_Vm_highp3_min=max(min(p_Vm_highp3(Nystart_Vm_highp3(ii,jj):Nyend_Vm_highp3(ii,jj)))-2,0);
            p_plot_Vm_highp3_max=max(p_Vm_highp3(Nystart_Vm_highp3(ii,jj):Nyend_Vm_highp3(ii,jj)))+30;     
            for jjk=p_plot_Vm_highp3_min:1:p_plot_Vm_highp3_max
                p_plot_Vm_highp3(pp,:)=jjk;
                props_Vm_highp3(pp,:)=computeThermoPropsForResult(T_plot_Vm_highp3, p_plot_Vm_highp3(pp,:));
    %             benzeneprops(T_plot,p_plot_Vm_highp3(pp,:));
                 Vm_p_plot_Vm_highp3(pp,:)=props_Vm_highp3(pp,1);
                pp=pp+1;
            end        
            if jj == 1
                plot(p_plot_Vm_highp3,Vm_p_plot_Vm_highp3,'color',mycolor(jj,:),'linewidth',linewidth);
            elseif jj == 2
                plot(p_plot_Vm_highp3,Vm_p_plot_Vm_highp3,...
                '--','color',mycolor(jj,:),'linewidth',linewidth);    
            elseif jj == 3
                plot(p_plot_Vm_highp3,Vm_p_plot_Vm_highp3,...
                '-.','color',mycolor(jj,:),'linewidth',linewidth);   
            elseif jj == 4
                plot(p_plot_Vm_highp3,Vm_p_plot_Vm_highp3,...
                ':','color',mycolor(jj,:),'linewidth',linewidth);
            elseif jj == 5
                plot(p_plot_Vm_highp3,Vm_p_plot_Vm_highp3,...
                '-','color',mycolor(jj,:),'linewidth',linewidth*2);      
            elseif jj == 6
                plot(p_plot_Vm_highp3,Vm_p_plot_Vm_highp3,...
                '--','color',mycolor(jj,:),'linewidth',linewidth*2);  
            elseif jj == 7
                plot(p_plot_Vm_highp3,Vm_p_plot_Vm_highp3,...
                '-.','color',mycolor(jj,:),'linewidth',linewidth*2);              
            end
        end
    end    
end  
        ytickformat('%.1f');  
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth); 
    ylabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
        set(gcf,'position',[500,200,300,260])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../../paper writing/figure/cell volume/high pressure/')
        mkdir('../../paper writing/figure/cell volume/high pressure/')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/cell volume/high pressure/high pressure_3.png')  
    
    
    figure(6)
 for ii = 1:a1_Vm_highp3
    for jj=1:a2_Vm_highp3
        Vm_highp3_dev = [];
        if Nystart_Vm_highp3(ii,jj)~=0
        ll1_Vm_highp3(ii,jj)=string(strcat(Author_Vm_highp3(Nystart_Vm_highp3(ii,jj)),Year_Vm_highp3(Nystart_Vm_highp3(ii,jj)),' - ',{' '},char(string(T_Vm_highp3(Nystart_Vm_highp3(ii,jj)))),{' '},'K'));            
        for  kkk = 1:length(p_Vm_highp3(Nystart_Vm_highp3(ii,jj):Nyend_Vm_highp3(ii,jj)))
            props_Vm_highp3 = computeThermoPropsForResult(T_Vm_highp3(Nystart_Vm_highp3(ii,jj)+kkk-1), p_Vm_highp3(Nystart_Vm_highp3(ii,jj)+kkk-1));
            Vm_highp3_dev(kkk,:) = 100*(Vm_highp3(Nystart_Vm_highp3(ii,jj)+kkk-1) - props_Vm_highp3(1)) ./ Vm_highp3(Nystart_Vm_highp3(ii,jj)+kkk-1);
        end
        plot(p_Vm_highp3(Nystart_Vm_highp3(ii,jj):Nyend_Vm_highp3(ii,jj)),Vm_highp3_dev,...
        mymarker{jj},'markersize',markersize,'color',mycolor(jj,:),'linewidth',linewidth,'DisplayName',ll1_Vm_highp3(ii,jj)); 
    hold on;
        end
    end    
 end
        c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
        ytickformat('%.1f'); 
%         xlim([0 800]);
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth);
%     ylabel(['100\cdot',char(8201),'(\itV\rm_m_,_e_x_p',char(hex2dec('2212')),'\it V\rm_m_,_c_a_l_c) /\it V\rm_m_,_c_a_l_c'],'fontsize',fontsize,'linewidth',linewidth)  
        annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfb','FontSize',fontsize,'FontName','Times');
        set(gcf,'position',[500,200,300,260])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    ylabel(['100\cdot',char(8201),'(\itV\rm_m_,_e_x_p',char(hex2dec('2212')),'\it V\rm_m_,_c_a_l_c) /\it V\rm_m_,_e_x_p'],'fontsize',fontsize-1,'linewidth',linewidth)  

    if ~isfolder('../../paper writing/figure/cell volume/high pressure/')
        mkdir('../../paper writing/figure/cell volume/high pressure/')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/cell volume/high pressure/high pressure_dev_3.png')
[Year_Vm_highp4,Author_Vm_highp4,T_Vm_highp4,Vm_highp4,p_Vm_highp4] = textread('../evaluation data/cell volume/Vm_high_pressure_4.txt','%s%s%f%f%f','headerlines',2);
ndata_Vm_highp4=length(T_Vm_highp4);xxcount_Vm_highp4=1; Nxstart_Vm_highp4 = 1;AA_Vm_highp4=Author_Vm_highp4(1);YY_Vm_highp4=Year_Vm_highp4(1);
for ii = 2:ndata_Vm_highp4         
        if ~strcmp(AA_Vm_highp4,Author_Vm_highp4(ii))||~strcmp(YY_Vm_highp4,Year_Vm_highp4(ii))
            Nxend_Vm_highp4(xxcount_Vm_highp4) = ii-1;
            xxcount_Vm_highp4 = xxcount_Vm_highp4 + 1;
            Nxstart_Vm_highp4(xxcount_Vm_highp4) = ii;
            AA_Vm_highp4=Author_Vm_highp4(ii);YY_Vm_highp4=Year_Vm_highp4(ii);
        end
        if ii == ndata_Vm_highp4
            Nxend_Vm_highp4(xxcount_Vm_highp4) = ndata_Vm_highp4;
        end
end


for jj=1:xxcount_Vm_highp4
    ndata_Vm_highp4=length(T_Vm_highp4(Nxstart_Vm_highp4(jj):Nxend_Vm_highp4(jj)));yycount_Vm_highp4=1; Nystart_Vm_highp4(jj,yycount_Vm_highp4) = Nxstart_Vm_highp4(jj);TT_Vm_highp4=T_Vm_highp4(Nxstart_Vm_highp4(jj));
    for ii = 2:ndata_Vm_highp4   

            %if ~strcmp(AA,Author(Nxstart(jj)+ii-1))&&~strcmp(YY,Year(Nxstart(jj)+ii-1))
            if abs(T_Vm_highp4(Nxstart_Vm_highp4(jj)+ii-1) - TT_Vm_highp4)> 0.05
                Nyend_Vm_highp4(jj,yycount_Vm_highp4) = Nxstart_Vm_highp4(jj)+ii-2;
                yycount_Vm_highp4 = yycount_Vm_highp4 + 1;
                Nystart_Vm_highp4(jj,yycount_Vm_highp4) = Nxstart_Vm_highp4(jj)+ii-1;
                TT_Vm_highp4=T_Vm_highp4(Nxstart_Vm_highp4(jj)+ii-1);
            end
            if ii == ndata_Vm_highp4
                Nyend_Vm_highp4(jj,yycount_Vm_highp4) = ndata_Vm_highp4 + Nxstart_Vm_highp4(jj)-1;
            end
    end     
end
[a1_Vm_highp4 a2_Vm_highp4]=size(Nystart_Vm_highp4);
kkkkk = 1;
figure(7)
%  my_dashline=[0 0 0 0;4 1 1 2;2 0.5 0.5 1;6 6 6 6;3 2 3 2;1 1 1 1;7 1 7 1;13 3 13 3;1 4 1 4];
for ii = 1:a1_Vm_highp4
    for jj=1:a2_Vm_highp4
        if Nystart_Vm_highp4(ii,jj)~=0
            ll1_Vm_highp4(ii,jj)=string(strcat(Author_Vm_highp4(Nystart_Vm_highp4(ii,jj)),Year_Vm_highp4(Nystart_Vm_highp4(ii,jj)),' - ',{' '},char(string(T_Vm_highp4(Nystart_Vm_highp4(ii,jj)))),{' '},'K'));            
            plot(p_Vm_highp4(Nystart_Vm_highp4(ii,jj):Nyend_Vm_highp4(ii,jj)),Vm_highp4(Nystart_Vm_highp4(ii,jj):Nyend_Vm_highp4(ii,jj)),...
            mymarker{kkkkk},'markersize',markersize,'color',mycolor(kkkkk,:),'linewidth',linewidth,'DisplayName',ll1_Vm_highp4(ii,jj)); 
        hold on;

            T_plot_Vm_highp4 = mean(T_Vm_highp4(Nystart_Vm_highp4(ii,jj):Nyend_Vm_highp4(ii,jj)));
            pp=1;
            p_plot_Vm_highp4=[];
            Vm_p_plot_Vm_highp4=[];  
            p_plot_Vm_highp4_min=max(min(p_Vm_highp4(Nystart_Vm_highp4(ii,jj):Nyend_Vm_highp4(ii,jj)))-2,0);
            p_plot_Vm_highp4_max=max(p_Vm_highp4(Nystart_Vm_highp4(ii,jj):Nyend_Vm_highp4(ii,jj)))+30;     
            for jjk=p_plot_Vm_highp4_min:1:p_plot_Vm_highp4_max
                p_plot_Vm_highp4(pp,:)=jjk;
                props_Vm_highp4(pp,:)=computeThermoPropsForResult(T_plot_Vm_highp4, p_plot_Vm_highp4(pp,:));
    %             benzeneprops(T_plot,p_plot_Vm_highp4(pp,:));
                 Vm_p_plot_Vm_highp4(pp,:)=props_Vm_highp4(pp,1);
                pp=pp+1;
            end        
            if jj == 1
                plot(p_plot_Vm_highp4,Vm_p_plot_Vm_highp4,'color',mycolor(kkkkk,:),'linewidth',linewidth);
            elseif jj == 2
                plot(p_plot_Vm_highp4,Vm_p_plot_Vm_highp4,...
                '--','color',mycolor(kkkkk,:),'linewidth',linewidth);    
            elseif jj == 3
                plot(p_plot_Vm_highp4,Vm_p_plot_Vm_highp4,...
                '-.','color',mycolor(kkkkk,:),'linewidth',linewidth);   
            elseif jj == 4
                plot(p_plot_Vm_highp4,Vm_p_plot_Vm_highp4,...
                ':','color',mycolor(kkkkk,:),'linewidth',linewidth);
            elseif jj == 5
                plot(p_plot_Vm_highp4,Vm_p_plot_Vm_highp4,...
                '-','color',mycolor(kkkkk,:),'linewidth',linewidth*2);      
            elseif jj == 6
                plot(p_plot_Vm_highp4,Vm_p_plot_Vm_highp4,...
                '--','color',mycolor(kkkkk,:),'linewidth',linewidth*2);  
            elseif jj == 7
                plot(p_plot_Vm_highp4,Vm_p_plot_Vm_highp4,...
                '-.','color',mycolor(kkkkk,:),'linewidth',linewidth*2);              
            end
            kkkkk = kkkkk + 1;
        end
    end    
end  
%         
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth); 
    ylabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfc','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
        set(gcf,'position',[500,200,300,260])
        ytickformat('%.1f');  
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../../paper writing/figure/cell volume/high pressure/')
        mkdir('../../paper writing/figure/cell volume/high pressure/')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/cell volume/high pressure/high pressure_4.png')  
    
    kkkkk = 1;
    figure(8)
 for ii = 1:a1_Vm_highp4
    for jj=1:a2_Vm_highp4
        Vm_highp4_dev = [];
        if Nystart_Vm_highp4(ii,jj)~=0
        ll1_Vm_highp4(ii,jj)=string(strcat(Author_Vm_highp4(Nystart_Vm_highp4(ii,jj)),Year_Vm_highp4(Nystart_Vm_highp4(ii,jj)),' - ',{' '},char(string(T_Vm_highp4(Nystart_Vm_highp4(ii,jj)))),{' '},'K'));            
        for  kkk = 1:length(p_Vm_highp4(Nystart_Vm_highp4(ii,jj):Nyend_Vm_highp4(ii,jj)))
            props_Vm_highp4 = computeThermoPropsForResult(T_Vm_highp4(Nystart_Vm_highp4(ii,jj)+kkk-1), p_Vm_highp4(Nystart_Vm_highp4(ii,jj)+kkk-1));
            Vm_highp4_dev(kkk,:) = 100*(Vm_highp4(Nystart_Vm_highp4(ii,jj)+kkk-1) - props_Vm_highp4(1)) ./ Vm_highp4(Nystart_Vm_highp4(ii,jj)+kkk-1);
        end
        plot(p_Vm_highp4(Nystart_Vm_highp4(ii,jj):Nyend_Vm_highp4(ii,jj)),Vm_highp4_dev,...
        mymarker{kkkkk},'markersize',markersize,'color',mycolor(kkkkk,:),'linewidth',linewidth,'DisplayName',ll1_Vm_highp4(ii,jj)); 
    hold on;
    kkkkk = kkkkk + 1;
        end
    end    
 end
        c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
        ytickformat('%.1f'); 
        xlim([0 40]);
%         legend show
%         legend('boxoff');
%         legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth);
%     ylabel(['100\cdot',char(8201),'(\itV\rm_m_,_e_x_p',char(hex2dec('2212')),'\it V\rm_m_,_c_a_l_c) /\it V\rm_m_,_c_a_l_c'],'fontsize',fontsize,'linewidth',linewidth)  
        annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfd','FontSize',fontsize,'FontName','Times');
        set(gcf,'position',[500,200,300,260])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    ylabel(['100\cdot',char(8201),'(\itV\rm_m_,_e_x_p',char(hex2dec('2212')),'\it V\rm_m_,_c_a_l_c) /\it V\rm_m_,_e_x_p'],'fontsize',fontsize-1,'linewidth',linewidth)  

    if ~isfolder('../../paper writing/figure/cell volume/high pressure/')
        mkdir('../../paper writing/figure/cell volume/high pressure/')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/cell volume/high pressure/high pressure_dev_4.png')    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 14.5
        fontsize = 12;    markersize = 7;    linewidth = 1;

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

figure(1)
%  my_dashline=[0 0 0 0;4 1 1 2;2 0.5 0.5 1;6 6 6 6;3 2 3 2;1 1 1 1;7 1 7 1;13 3 13 3;1 4 1 4];
for ii = 1:a1_Vm_highp
    for jj=1:a2_Vm_highp
        if Nystart_Vm_highp(ii,jj)~=0
            ll1_Vm_highp(ii,jj)=string(strcat(Author_Vm_highp(Nystart_Vm_highp(ii,jj)),Year_Vm_highp(Nystart_Vm_highp(ii,jj)),' - ',{' '},char(string(T_Vm_highp(Nystart_Vm_highp(ii,jj)))),{' '},'K'));            
            plot(p_Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)),Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)),...
            mymarker{jj},'markersize',markersize,'color',mycolor(jj,:),'linewidth',linewidth,'DisplayName',ll1_Vm_highp(ii,jj)); 
        hold on;

            T_plot_Vm_highp = mean(T_Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)));
            pp=1;
            p_plot_Vm_highp=[];
            Vm_p_plot_Vm_highp=[];  
            p_plot_Vm_highp_min=max(min(p_Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)))-2,0);
            p_plot_Vm_highp_max=max(p_Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)))+30;     
            for jjk=p_plot_Vm_highp_min:1:p_plot_Vm_highp_max
                p_plot_Vm_highp(pp,:)=jjk;
                props_Vm_highp(pp,:)=computeThermoPropsForResult(T_plot_Vm_highp, p_plot_Vm_highp(pp,:));
    %             benzeneprops(T_plot,p_plot_Vm_highp(pp,:));
                 Vm_p_plot_Vm_highp(pp,:)=props_Vm_highp(pp,1);
                pp=pp+1;
            end        
            if jj == 1
                plot(p_plot_Vm_highp,Vm_p_plot_Vm_highp,'color',mycolor(jj,:),'linewidth',linewidth);
            elseif jj == 2
                plot(p_plot_Vm_highp,Vm_p_plot_Vm_highp,...
                '--','color',mycolor(jj,:),'linewidth',linewidth);    
            elseif jj == 3
                plot(p_plot_Vm_highp,Vm_p_plot_Vm_highp,...
                '-.','color',mycolor(jj,:),'linewidth',linewidth);   
            elseif jj == 4
                plot(p_plot_Vm_highp,Vm_p_plot_Vm_highp,...
                ':','color',mycolor(jj,:),'linewidth',linewidth);
            elseif jj == 5
                plot(p_plot_Vm_highp,Vm_p_plot_Vm_highp,...
                '-','color',mycolor(jj,:),'linewidth',linewidth*2);      
            elseif jj == 6
                plot(p_plot_Vm_highp,Vm_p_plot_Vm_highp,...
                '--','color',mycolor(jj,:),'linewidth',linewidth*2);  
            elseif jj == 7
                plot(p_plot_Vm_highp,Vm_p_plot_Vm_highp,...
                '-.','color',mycolor(jj,:),'linewidth',linewidth*2);              
            end
        end
    end    
end  
        ytickformat('%.1f');  
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth); 
    ylabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
      set(gcf,'position',[300,10,800,750])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../../paper writing/figure/cell volume/high pressure/legend/')
        mkdir('../../paper writing/figure/cell volume/high pressure/legend/')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/cell volume/high pressure/legend/high pressure_1.png')  
    
    
    figure(2)
 for ii = 1:a1_Vm_highp
    for jj=1:a2_Vm_highp
        Vm_highp_dev = [];
        if Nystart_Vm_highp(ii,jj)~=0
        ll1_Vm_highp(ii,jj)=string(strcat(Author_Vm_highp(Nystart_Vm_highp(ii,jj)),Year_Vm_highp(Nystart_Vm_highp(ii,jj)),' - ',{' '},char(string(T_Vm_highp(Nystart_Vm_highp(ii,jj)))),{' '},'K'));            
        for  kkk = 1:length(p_Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)))
            props_Vm_highp = computeThermoPropsForResult(T_Vm_highp(Nystart_Vm_highp(ii,jj)+kkk-1), p_Vm_highp(Nystart_Vm_highp(ii,jj)+kkk-1));
            Vm_highp_dev(kkk,:) = 100*(Vm_highp(Nystart_Vm_highp(ii,jj)+kkk-1) - props_Vm_highp(1)) ./ Vm_highp(Nystart_Vm_highp(ii,jj)+kkk-1);
        end
        plot(p_Vm_highp(Nystart_Vm_highp(ii,jj):Nyend_Vm_highp(ii,jj)),Vm_highp_dev,...
        mymarker{jj},'markersize',markersize,'color',mycolor(jj,:),'linewidth',linewidth,'DisplayName',ll1_Vm_highp(ii,jj)); 
    hold on;
        end
    end    
 end
        c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
        ytickformat('%.1f'); 
        xlim([0 220]);
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth);
        annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfb','FontSize',fontsize,'FontName','Times');

        set(gcf,'position',[300,10,800,750])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
                ylabel(['100\cdot',char(8201),'(\itV\rm_m_,_e_x_p',char(hex2dec('2212')),'\it V\rm_m_,_c_a_l_c) /\it V\rm_m_,_c_a_l_c'],'fontsize',fontsize-1,'linewidth',linewidth)  

    if ~isfolder('../../paper writing/figure/cell volume/high pressure/')
        mkdir('../../paper writing/figure/cell volume/high pressure/')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/cell volume/high pressure/legend/high pressure_dev_1.png')  
    
    
    
    [Year_Vm_highp2,Author_Vm_highp2,T_Vm_highp2,Vm_highp2,p_Vm_highp2] = textread('../evaluation data/cell volume/Vm_high_pressure_2.txt','%s%s%f%f%f','headerlines',2);
ndata_Vm_highp2=length(T_Vm_highp2);xxcount_Vm_highp2=1; Nxstart_Vm_highp2 = 1;AA_Vm_highp2=Author_Vm_highp2(1);YY_Vm_highp2=Year_Vm_highp2(1);
for ii = 2:ndata_Vm_highp2         
        if ~strcmp(AA_Vm_highp2,Author_Vm_highp2(ii))||~strcmp(YY_Vm_highp2,Year_Vm_highp2(ii))
            Nxend_Vm_highp2(xxcount_Vm_highp2) = ii-1;
            xxcount_Vm_highp2 = xxcount_Vm_highp2 + 1;
            Nxstart_Vm_highp2(xxcount_Vm_highp2) = ii;
            AA_Vm_highp2=Author_Vm_highp2(ii);YY_Vm_highp2=Year_Vm_highp2(ii);
        end
        if ii == ndata_Vm_highp2
            Nxend_Vm_highp2(xxcount_Vm_highp2) = ndata_Vm_highp2;
        end
end


for jj=1:xxcount_Vm_highp2
    ndata_Vm_highp2=length(T_Vm_highp2(Nxstart_Vm_highp2(jj):Nxend_Vm_highp2(jj)));yycount_Vm_highp2=1; Nystart_Vm_highp2(jj,yycount_Vm_highp2) = Nxstart_Vm_highp2(jj);TT_Vm_highp2=T_Vm_highp2(Nxstart_Vm_highp2(jj));
    for ii = 2:ndata_Vm_highp2   

            %if ~strcmp(AA,Author(Nxstart(jj)+ii-1))&&~strcmp(YY,Year(Nxstart(jj)+ii-1))
            if abs(T_Vm_highp2(Nxstart_Vm_highp2(jj)+ii-1) - TT_Vm_highp2)> 0.05
                Nyend_Vm_highp2(jj,yycount_Vm_highp2) = Nxstart_Vm_highp2(jj)+ii-2;
                yycount_Vm_highp2 = yycount_Vm_highp2 + 1;
                Nystart_Vm_highp2(jj,yycount_Vm_highp2) = Nxstart_Vm_highp2(jj)+ii-1;
                TT_Vm_highp2=T_Vm_highp2(Nxstart_Vm_highp2(jj)+ii-1);
            end
            if ii == ndata_Vm_highp2
                Nyend_Vm_highp2(jj,yycount_Vm_highp2) = ndata_Vm_highp2 + Nxstart_Vm_highp2(jj)-1;
            end
    end     
end
[a1_Vm_highp2 a2_Vm_highp2]=size(Nystart_Vm_highp2);

figure(3)
%  my_dashline=[0 0 0 0;4 1 1 2;2 0.5 0.5 1;6 6 6 6;3 2 3 2;1 1 1 1;7 1 7 1;13 3 13 3;1 4 1 4];
for ii = 1:a1_Vm_highp2
    for jj=1:a2_Vm_highp2
        if Nystart_Vm_highp2(ii,jj)~=0
            ll1_Vm_highp2(ii,jj)=string(strcat(Author_Vm_highp2(Nystart_Vm_highp2(ii,jj)),Year_Vm_highp2(Nystart_Vm_highp2(ii,jj)),' - ',{' '},char(string(T_Vm_highp2(Nystart_Vm_highp2(ii,jj)))),{' '},'K'));            
            plot(p_Vm_highp2(Nystart_Vm_highp2(ii,jj):Nyend_Vm_highp2(ii,jj)),Vm_highp2(Nystart_Vm_highp2(ii,jj):Nyend_Vm_highp2(ii,jj)),...
            mymarker{jj},'markersize',markersize,'color',mycolor(jj,:),'linewidth',linewidth,'DisplayName',ll1_Vm_highp2(ii,jj)); 
        hold on;

            T_plot_Vm_highp2 = mean(T_Vm_highp2(Nystart_Vm_highp2(ii,jj):Nyend_Vm_highp2(ii,jj)));
            pp=1;
            p_plot_Vm_highp2=[];
            Vm_p_plot_Vm_highp2=[];  
            p_plot_Vm_highp2_min=max(min(p_Vm_highp2(Nystart_Vm_highp2(ii,jj):Nyend_Vm_highp2(ii,jj)))-2,0);
            p_plot_Vm_highp2_max=max(p_Vm_highp2(Nystart_Vm_highp2(ii,jj):Nyend_Vm_highp2(ii,jj)))+30;     
            for jjk=p_plot_Vm_highp2_min:1:p_plot_Vm_highp2_max
                p_plot_Vm_highp2(pp,:)=jjk;
                props_Vm_highp2(pp,:)=computeThermoPropsForResult(T_plot_Vm_highp2, p_plot_Vm_highp2(pp,:));
    %             benzeneprops(T_plot,p_plot_Vm_highp2(pp,:));
                 Vm_p_plot_Vm_highp2(pp,:)=props_Vm_highp2(pp,1);
                pp=pp+1;
            end        
            if jj == 1
                plot(p_plot_Vm_highp2,Vm_p_plot_Vm_highp2,'color',mycolor(jj,:),'linewidth',linewidth);
            elseif jj == 2
                plot(p_plot_Vm_highp2,Vm_p_plot_Vm_highp2,...
                '--','color',mycolor(jj,:),'linewidth',linewidth);    
            elseif jj == 3
                plot(p_plot_Vm_highp2,Vm_p_plot_Vm_highp2,...
                '-.','color',mycolor(jj,:),'linewidth',linewidth);   
            elseif jj == 4
                plot(p_plot_Vm_highp2,Vm_p_plot_Vm_highp2,...
                ':','color',mycolor(jj,:),'linewidth',linewidth);
            elseif jj == 5
                plot(p_plot_Vm_highp2,Vm_p_plot_Vm_highp2,...
                '-','color',mycolor(jj,:),'linewidth',linewidth*2);      
            elseif jj == 6
                plot(p_plot_Vm_highp2,Vm_p_plot_Vm_highp2,...
                '--','color',mycolor(jj,:),'linewidth',linewidth*2);  
            elseif jj == 7
                plot(p_plot_Vm_highp2,Vm_p_plot_Vm_highp2,...
                '-.','color',mycolor(jj,:),'linewidth',linewidth*2);              
            end
        end
    end    
end  
%         ytickformat('%.1f');  
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth); 
    ylabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfc','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
        set(gcf,'position',[300,10,800,750])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../../paper writing/figure/cell volume/high pressure/')
        mkdir('../../paper writing/figure/cell volume/high pressure/')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/cell volume/high pressure/legend/high pressure_2.png')  
    
    
    figure(4)
 for ii = 1:a1_Vm_highp2
    for jj=1:a2_Vm_highp2
        Vm_highp2_dev = [];
        if Nystart_Vm_highp2(ii,jj)~=0
        ll1_Vm_highp2(ii,jj)=string(strcat(Author_Vm_highp2(Nystart_Vm_highp2(ii,jj)),Year_Vm_highp2(Nystart_Vm_highp2(ii,jj)),' - ',{' '},char(string(T_Vm_highp2(Nystart_Vm_highp2(ii,jj)))),{' '},'K'));            
        for  kkk = 1:length(p_Vm_highp2(Nystart_Vm_highp2(ii,jj):Nyend_Vm_highp2(ii,jj)))
            props_Vm_highp2 = computeThermoPropsForResult(T_Vm_highp2(Nystart_Vm_highp2(ii,jj)+kkk-1), p_Vm_highp2(Nystart_Vm_highp2(ii,jj)+kkk-1));
            Vm_highp2_dev(kkk,:) = 100*(Vm_highp2(Nystart_Vm_highp2(ii,jj)+kkk-1) - props_Vm_highp2(1)) ./ Vm_highp2(Nystart_Vm_highp2(ii,jj)+kkk-1);
        end
        plot(p_Vm_highp2(Nystart_Vm_highp2(ii,jj):Nyend_Vm_highp2(ii,jj)),Vm_highp2_dev,...
        mymarker{jj},'markersize',markersize,'color',mycolor(jj,:),'linewidth',linewidth,'DisplayName',ll1_Vm_highp2(ii,jj)); 
    hold on;
        end
    end    
 end
        c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
        ytickformat('%.1f'); 
        xlim([0 500]);
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth);
%     ylabel(['100\cdot',char(8201),'(\itV\rm_m_,_e_x_p',char(hex2dec('2212')),'\it V\rm_m_,_c_a_l_c) /\it V\rm_m_,_c_a_l_c'],'fontsize',fontsize,'linewidth',linewidth)  
        annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfd','FontSize',fontsize,'FontName','Times');
        set(gcf,'position',[300,10,800,750])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
                    ylabel(['100\cdot',char(8201),'(\itV\rm_m_,_e_x_p',char(hex2dec('2212')),'\it V\rm_m_,_c_a_l_c) /\it V\rm_m_,_c_a_l_c'],'fontsize',fontsize-1,'linewidth',linewidth)  

    if ~isfolder('../../paper writing/figure/cell volume/high pressure/')
        mkdir('../../paper writing/figure/cell volume/high pressure/')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/cell volume/high pressure/legend/high pressure_dev_2.png')
    [Year_Vm_highp3,Author_Vm_highp3,T_Vm_highp3,Vm_highp3,p_Vm_highp3] = textread('../evaluation data/cell volume/Vm_high_pressure_3.txt','%s%s%f%f%f','headerlines',2);
ndata_Vm_highp3=length(T_Vm_highp3);xxcount_Vm_highp3=1; Nxstart_Vm_highp3 = 1;AA_Vm_highp3=Author_Vm_highp3(1);YY_Vm_highp3=Year_Vm_highp3(1);
for ii = 2:ndata_Vm_highp3         
        if ~strcmp(AA_Vm_highp3,Author_Vm_highp3(ii))||~strcmp(YY_Vm_highp3,Year_Vm_highp3(ii))
            Nxend_Vm_highp3(xxcount_Vm_highp3) = ii-1;
            xxcount_Vm_highp3 = xxcount_Vm_highp3 + 1;
            Nxstart_Vm_highp3(xxcount_Vm_highp3) = ii;
            AA_Vm_highp3=Author_Vm_highp3(ii);YY_Vm_highp3=Year_Vm_highp3(ii);
        end
        if ii == ndata_Vm_highp3
            Nxend_Vm_highp3(xxcount_Vm_highp3) = ndata_Vm_highp3;
        end
end


for jj=1:xxcount_Vm_highp3
    ndata_Vm_highp3=length(T_Vm_highp3(Nxstart_Vm_highp3(jj):Nxend_Vm_highp3(jj)));yycount_Vm_highp3=1; Nystart_Vm_highp3(jj,yycount_Vm_highp3) = Nxstart_Vm_highp3(jj);TT_Vm_highp3=T_Vm_highp3(Nxstart_Vm_highp3(jj));
    for ii = 2:ndata_Vm_highp3   

            %if ~strcmp(AA,Author(Nxstart(jj)+ii-1))&&~strcmp(YY,Year(Nxstart(jj)+ii-1))
            if abs(T_Vm_highp3(Nxstart_Vm_highp3(jj)+ii-1) - TT_Vm_highp3)> 0.05
                Nyend_Vm_highp3(jj,yycount_Vm_highp3) = Nxstart_Vm_highp3(jj)+ii-2;
                yycount_Vm_highp3 = yycount_Vm_highp3 + 1;
                Nystart_Vm_highp3(jj,yycount_Vm_highp3) = Nxstart_Vm_highp3(jj)+ii-1;
                TT_Vm_highp3=T_Vm_highp3(Nxstart_Vm_highp3(jj)+ii-1);
            end
            if ii == ndata_Vm_highp3
                Nyend_Vm_highp3(jj,yycount_Vm_highp3) = ndata_Vm_highp3 + Nxstart_Vm_highp3(jj)-1;
            end
    end     
end
[a1_Vm_highp3 a2_Vm_highp3]=size(Nystart_Vm_highp3);

figure(5)
%  my_dashline=[0 0 0 0;4 1 1 2;2 0.5 0.5 1;6 6 6 6;3 2 3 2;1 1 1 1;7 1 7 1;13 3 13 3;1 4 1 4];
for ii = 1:a1_Vm_highp3
    for jj=1:a2_Vm_highp3
        if Nystart_Vm_highp3(ii,jj)~=0
            ll1_Vm_highp3(ii,jj)=string(strcat(Author_Vm_highp3(Nystart_Vm_highp3(ii,jj)),Year_Vm_highp3(Nystart_Vm_highp3(ii,jj)),' - ',{' '},char(string(T_Vm_highp3(Nystart_Vm_highp3(ii,jj)))),{' '},'K'));            
            plot(p_Vm_highp3(Nystart_Vm_highp3(ii,jj):Nyend_Vm_highp3(ii,jj)),Vm_highp3(Nystart_Vm_highp3(ii,jj):Nyend_Vm_highp3(ii,jj)),...
            mymarker{jj},'markersize',markersize,'color',mycolor(jj,:),'linewidth',linewidth,'DisplayName',ll1_Vm_highp3(ii,jj)); 
        hold on;

            T_plot_Vm_highp3 = mean(T_Vm_highp3(Nystart_Vm_highp3(ii,jj):Nyend_Vm_highp3(ii,jj)));
            pp=1;
            p_plot_Vm_highp3=[];
            Vm_p_plot_Vm_highp3=[];  
            p_plot_Vm_highp3_min=max(min(p_Vm_highp3(Nystart_Vm_highp3(ii,jj):Nyend_Vm_highp3(ii,jj)))-2,0);
            p_plot_Vm_highp3_max=max(p_Vm_highp3(Nystart_Vm_highp3(ii,jj):Nyend_Vm_highp3(ii,jj)))+30;     
            for jjk=p_plot_Vm_highp3_min:1:p_plot_Vm_highp3_max
                p_plot_Vm_highp3(pp,:)=jjk;
                props_Vm_highp3(pp,:)=computeThermoPropsForResult(T_plot_Vm_highp3, p_plot_Vm_highp3(pp,:));
    %             benzeneprops(T_plot,p_plot_Vm_highp3(pp,:));
                 Vm_p_plot_Vm_highp3(pp,:)=props_Vm_highp3(pp,1);
                pp=pp+1;
            end        
            if jj == 1
                plot(p_plot_Vm_highp3,Vm_p_plot_Vm_highp3,'color',mycolor(jj,:),'linewidth',linewidth);
            elseif jj == 2
                plot(p_plot_Vm_highp3,Vm_p_plot_Vm_highp3,...
                '--','color',mycolor(jj,:),'linewidth',linewidth);    
            elseif jj == 3
                plot(p_plot_Vm_highp3,Vm_p_plot_Vm_highp3,...
                '-.','color',mycolor(jj,:),'linewidth',linewidth);   
            elseif jj == 4
                plot(p_plot_Vm_highp3,Vm_p_plot_Vm_highp3,...
                ':','color',mycolor(jj,:),'linewidth',linewidth);
            elseif jj == 5
                plot(p_plot_Vm_highp3,Vm_p_plot_Vm_highp3,...
                '-','color',mycolor(jj,:),'linewidth',linewidth*2);      
            elseif jj == 6
                plot(p_plot_Vm_highp3,Vm_p_plot_Vm_highp3,...
                '--','color',mycolor(jj,:),'linewidth',linewidth*2);  
            elseif jj == 7
                plot(p_plot_Vm_highp3,Vm_p_plot_Vm_highp3,...
                '-.','color',mycolor(jj,:),'linewidth',linewidth*2);              
            end
        end
    end    
end  
        ytickformat('%.1f');  
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth); 
    ylabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfc','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
        set(gcf,'position',[300,10,800,750])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../../paper writing/figure/cell volume/high pressure/')
        mkdir('../../paper writing/figure/cell volume/high pressure/')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/cell volume/high pressure/legend/high pressure_3.png')  
    
    
    figure(6)
 for ii = 1:a1_Vm_highp3
    for jj=1:a2_Vm_highp3
        Vm_highp3_dev = [];
        if Nystart_Vm_highp3(ii,jj)~=0
        ll1_Vm_highp3(ii,jj)=string(strcat(Author_Vm_highp3(Nystart_Vm_highp3(ii,jj)),Year_Vm_highp3(Nystart_Vm_highp3(ii,jj)),' - ',{' '},char(string(T_Vm_highp3(Nystart_Vm_highp3(ii,jj)))),{' '},'K'));            
        for  kkk = 1:length(p_Vm_highp3(Nystart_Vm_highp3(ii,jj):Nyend_Vm_highp3(ii,jj)))
            props_Vm_highp3 = computeThermoPropsForResult(T_Vm_highp3(Nystart_Vm_highp3(ii,jj)+kkk-1), p_Vm_highp3(Nystart_Vm_highp3(ii,jj)+kkk-1));
            Vm_highp3_dev(kkk,:) = 100*(Vm_highp3(Nystart_Vm_highp3(ii,jj)+kkk-1) - props_Vm_highp3(1)) ./ Vm_highp3(Nystart_Vm_highp3(ii,jj)+kkk-1);
        end
        plot(p_Vm_highp3(Nystart_Vm_highp3(ii,jj):Nyend_Vm_highp3(ii,jj)),Vm_highp3_dev,...
        mymarker{jj},'markersize',markersize,'color',mycolor(jj,:),'linewidth',linewidth,'DisplayName',ll1_Vm_highp3(ii,jj)); 
    hold on;
        end
    end    
 end
        c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
        ytickformat('%.1f'); 
%         xlim([0 800]);
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth);
%     ylabel(['100\cdot',char(8201),'(\itV\rm_m_,_e_x_p',char(hex2dec('2212')),'\it V\rm_m_,_c_a_l_c) /\it V\rm_m_,_c_a_l_c'],'fontsize',fontsize,'linewidth',linewidth)  
        annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfd','FontSize',fontsize,'FontName','Times');
        set(gcf,'position',[300,10,800,750])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    ylabel(['100\cdot',char(8201),'(\itV\rm_m_,_e_x_p',char(hex2dec('2212')),'\it V\rm_m_,_c_a_l_c) /\it V\rm_m_,_c_a_l_c'],'fontsize',fontsize-1,'linewidth',linewidth)  

    if ~isfolder('../../paper writing/figure/cell volume/high pressure/')
        mkdir('../../paper writing/figure/cell volume/high pressure/')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/cell volume/high pressure/legend/high pressure_dev_3.png')
    
[Year_Vm_highp4,Author_Vm_highp4,T_Vm_highp4,Vm_highp4,p_Vm_highp4] = textread('../evaluation data/cell volume/Vm_high_pressure_4.txt','%s%s%f%f%f','headerlines',2);
ndata_Vm_highp4=length(T_Vm_highp4);xxcount_Vm_highp4=1; Nxstart_Vm_highp4 = 1;AA_Vm_highp4=Author_Vm_highp4(1);YY_Vm_highp4=Year_Vm_highp4(1);
for ii = 2:ndata_Vm_highp4         
        if ~strcmp(AA_Vm_highp4,Author_Vm_highp4(ii))||~strcmp(YY_Vm_highp4,Year_Vm_highp4(ii))
            Nxend_Vm_highp4(xxcount_Vm_highp4) = ii-1;
            xxcount_Vm_highp4 = xxcount_Vm_highp4 + 1;
            Nxstart_Vm_highp4(xxcount_Vm_highp4) = ii;
            AA_Vm_highp4=Author_Vm_highp4(ii);YY_Vm_highp4=Year_Vm_highp4(ii);
        end
        if ii == ndata_Vm_highp4
            Nxend_Vm_highp4(xxcount_Vm_highp4) = ndata_Vm_highp4;
        end
end


for jj=1:xxcount_Vm_highp4
    ndata_Vm_highp4=length(T_Vm_highp4(Nxstart_Vm_highp4(jj):Nxend_Vm_highp4(jj)));yycount_Vm_highp4=1; Nystart_Vm_highp4(jj,yycount_Vm_highp4) = Nxstart_Vm_highp4(jj);TT_Vm_highp4=T_Vm_highp4(Nxstart_Vm_highp4(jj));
    for ii = 2:ndata_Vm_highp4   

            %if ~strcmp(AA,Author(Nxstart(jj)+ii-1))&&~strcmp(YY,Year(Nxstart(jj)+ii-1))
            if abs(T_Vm_highp4(Nxstart_Vm_highp4(jj)+ii-1) - TT_Vm_highp4)> 0.05
                Nyend_Vm_highp4(jj,yycount_Vm_highp4) = Nxstart_Vm_highp4(jj)+ii-2;
                yycount_Vm_highp4 = yycount_Vm_highp4 + 1;
                Nystart_Vm_highp4(jj,yycount_Vm_highp4) = Nxstart_Vm_highp4(jj)+ii-1;
                TT_Vm_highp4=T_Vm_highp4(Nxstart_Vm_highp4(jj)+ii-1);
            end
            if ii == ndata_Vm_highp4
                Nyend_Vm_highp4(jj,yycount_Vm_highp4) = ndata_Vm_highp4 + Nxstart_Vm_highp4(jj)-1;
            end
    end     
end
[a1_Vm_highp4 a2_Vm_highp4]=size(Nystart_Vm_highp4);
kkkkk = 1;
figure(7)
%  my_dashline=[0 0 0 0;4 1 1 2;2 0.5 0.5 1;6 6 6 6;3 2 3 2;1 1 1 1;7 1 7 1;13 3 13 3;1 4 1 4];
for ii = 1:a1_Vm_highp4
    for jj=1:a2_Vm_highp4
        if Nystart_Vm_highp4(ii,jj)~=0
            ll1_Vm_highp4(ii,jj)=string(strcat(Author_Vm_highp4(Nystart_Vm_highp4(ii,jj)),Year_Vm_highp4(Nystart_Vm_highp4(ii,jj)),' - ',{' '},char(string(T_Vm_highp4(Nystart_Vm_highp4(ii,jj)))),{' '},'K'));            
            plot(p_Vm_highp4(Nystart_Vm_highp4(ii,jj):Nyend_Vm_highp4(ii,jj)),Vm_highp4(Nystart_Vm_highp4(ii,jj):Nyend_Vm_highp4(ii,jj)),...
            mymarker{kkkkk},'markersize',markersize,'color',mycolor(kkkkk,:),'linewidth',linewidth,'DisplayName',ll1_Vm_highp4(ii,jj)); 
        hold on;

            T_plot_Vm_highp4 = mean(T_Vm_highp4(Nystart_Vm_highp4(ii,jj):Nyend_Vm_highp4(ii,jj)));
            pp=1;
            p_plot_Vm_highp4=[];
            Vm_p_plot_Vm_highp4=[];  
            p_plot_Vm_highp4_min=max(min(p_Vm_highp4(Nystart_Vm_highp4(ii,jj):Nyend_Vm_highp4(ii,jj)))-2,0);
            p_plot_Vm_highp4_max=max(p_Vm_highp4(Nystart_Vm_highp4(ii,jj):Nyend_Vm_highp4(ii,jj)))+30;     
            for jjk=p_plot_Vm_highp4_min:1:p_plot_Vm_highp4_max
                p_plot_Vm_highp4(pp,:)=jjk;
                props_Vm_highp4(pp,:)=computeThermoPropsForResult(T_plot_Vm_highp4, p_plot_Vm_highp4(pp,:));
    %             benzeneprops(T_plot,p_plot_Vm_highp4(pp,:));
                 Vm_p_plot_Vm_highp4(pp,:)=props_Vm_highp4(pp,1);
                pp=pp+1;
            end        
            if jj == 1
                plot(p_plot_Vm_highp4,Vm_p_plot_Vm_highp4,'color',mycolor(kkkkk,:),'linewidth',linewidth);
            elseif jj == 2
                plot(p_plot_Vm_highp4,Vm_p_plot_Vm_highp4,...
                '--','color',mycolor(kkkkk,:),'linewidth',linewidth);    
            elseif jj == 3
                plot(p_plot_Vm_highp4,Vm_p_plot_Vm_highp4,...
                '-.','color',mycolor(kkkkk,:),'linewidth',linewidth);   
            elseif jj == 4
                plot(p_plot_Vm_highp4,Vm_p_plot_Vm_highp4,...
                ':','color',mycolor(kkkkk,:),'linewidth',linewidth);
            elseif jj == 5
                plot(p_plot_Vm_highp4,Vm_p_plot_Vm_highp4,...
                '-','color',mycolor(kkkkk,:),'linewidth',linewidth*2);      
            elseif jj == 6
                plot(p_plot_Vm_highp4,Vm_p_plot_Vm_highp4,...
                '--','color',mycolor(kkkkk,:),'linewidth',linewidth*2);  
            elseif jj == 7
                plot(p_plot_Vm_highp4,Vm_p_plot_Vm_highp4,...
                '-.','color',mycolor(kkkkk,:),'linewidth',linewidth*2);              
            end
            kkkkk = kkkkk + 1;
        end
    end    
end  
%         ytickformat('%.1f');  
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth); 
    ylabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfc','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
        set(gcf,'position',[300,10,800,750])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../../paper writing/figure/cell volume/high pressure/')
        mkdir('../../paper writing/figure/cell volume/high pressure/')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/cell volume/high pressure/legend/high pressure_4.png')  
    
    kkkkk = 1;
    figure(8)
 for ii = 1:a1_Vm_highp4
    for jj=1:a2_Vm_highp4
        Vm_highp4_dev = [];
        if Nystart_Vm_highp4(ii,jj)~=0
        ll1_Vm_highp4(ii,jj)=string(strcat(Author_Vm_highp4(Nystart_Vm_highp4(ii,jj)),Year_Vm_highp4(Nystart_Vm_highp4(ii,jj)),' - ',{' '},char(string(T_Vm_highp4(Nystart_Vm_highp4(ii,jj)))),{' '},'K'));            
        for  kkk = 1:length(p_Vm_highp4(Nystart_Vm_highp4(ii,jj):Nyend_Vm_highp4(ii,jj)))
            props_Vm_highp4 = computeThermoPropsForResult(T_Vm_highp4(Nystart_Vm_highp4(ii,jj)+kkk-1), p_Vm_highp4(Nystart_Vm_highp4(ii,jj)+kkk-1));
            Vm_highp4_dev(kkk,:) = 100*(Vm_highp4(Nystart_Vm_highp4(ii,jj)+kkk-1) - props_Vm_highp4(1)) ./ Vm_highp4(Nystart_Vm_highp4(ii,jj)+kkk-1);
        end
        plot(p_Vm_highp4(Nystart_Vm_highp4(ii,jj):Nyend_Vm_highp4(ii,jj)),Vm_highp4_dev,...
        mymarker{kkkkk},'markersize',markersize,'color',mycolor(kkkkk,:),'linewidth',linewidth,'DisplayName',ll1_Vm_highp4(ii,jj)); 
    hold on;
    kkkkk = kkkkk + 1;
        end
    end    
 end
        c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
        ytickformat('%.1f'); 
        xlim([0 40]);
        legend show
        legend('boxoff');
        legend ('Location','bestoutside','NumColumns',1);
    xlabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth);
%     ylabel(['100\cdot',char(8201),'(\itV\rm_m_,_e_x_p',char(hex2dec('2212')),'\it V\rm_m_,_c_a_l_c) /\it V\rm_m_,_c_a_l_c'],'fontsize',fontsize,'linewidth',linewidth)  
        annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfd','FontSize',fontsize,'FontName','Times');
   set(gcf,'position',[300,10,800,750])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    ylabel(['100\cdot',char(8201),'(\itV\rm_m_,_e_x_p',char(hex2dec('2212')),'\it V\rm_m_,_c_a_l_c) /\it V\rm_m_,_c_a_l_c'],'fontsize',fontsize-1,'linewidth',linewidth)  

    if ~isfolder('../../paper writing/figure/cell volume/high pressure/')
        mkdir('../../paper writing/figure/cell volume/high pressure/')
    end
    print(gcf,'-dtiff','-r300','../../paper writing/figure/cell volume/high pressure/legend/high pressure_dev_4.png')    
end
