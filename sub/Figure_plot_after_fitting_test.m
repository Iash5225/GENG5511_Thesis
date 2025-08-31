% if Lplot == 1
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
                TT_Vm_highp=T_Vm_highp(Nxstart_Vm_highp(jj)+ii);
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
xlim([0 85]);
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
    ylabel(['100\cdot',char(8201),'(\it\alpha\rm_,_e_x_p',char(hex2dec('2212')),'\it \alpha\rm_,_c_a_l_c) /\it \alpha\rm_,_c_a_l_c'],'fontsize',fontsize,'linewidth',linewidth)  
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
 V_T7_plot=[];p7 = [];
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

    xlabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
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
    T1 = 0.0001;
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
% end 