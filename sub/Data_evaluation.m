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
RP = py.ctREFPROP.ctREFPROP.REFPROPFunctionLibrary('C:\Refprop\REFPROP');
RP.SETPATHdll('C:\Refprop\REFPROP');
SI = RP.GETENUMdll(int8(0),'MASS SI').iEnum;
iMass = int8(0);% 0: molar fractions; 1: mass fractions
iFlag = int8(1);% 0: don't call SATSPLN; 1: call SATSPLN
RP.FLAGSdll("Reset HMX",int8(1));
RP.FLAGSdll("Peng-Robinson",int8(0));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choice = 1.1;%Vm along sublimation
% choice = 1.2;%Vm along melting
% choice = 1.3;%Vm high pressure
% choice = 2.1;%cp along sublimation
% choice = 3.1;%alpha along sublimation
% choice = 4.1;%Beta T along sublimation
% choice = 4.15;%Beta S along sublimation
choice = 5;%sublimation and fit
% choice = 6;%melting
% choice = 7.1;%enthalpy of sublimation
% choice = 7.2;%enthalpy of fusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 1.1
    [Year,Author,T_Vm,Vm] = textread('../evaluation data/cell volume/Vm_sublimation.txt','%s%s%f%f','headerlines',2);
    AA=Author(1);
NYstart=1;xxcount=1;
ndata=length(Year);
 Author_str=string(Author);
 YY=Year(1);
for i = 1:ndata
    Author_Y(i,:)=strcat(Year(i,:),'-',Author_str(i,:));
end
    for ii=2:ndata
        if ~strcmp(AA,Author(ii))||~strcmp(YY,Year(ii))
            NYend(xxcount)=ii-1;
            xxcount=xxcount+1;
            Legend_Author_Vm(xxcount-1,:)=Author_Y(ii-1,:);
            Author_alpha(xxcount-1,:)=Author_str(ii-1,:);
            NYstart(xxcount)=ii;
            AA=Author(ii);YY=Year(ii);
        end
        if ii==ndata
            NYend(xxcount)=ndata;
            Legend_Author_Vm(xxcount,:)=Author_Y(ii,:); 
            Author_alpha(xxcount,:)=Author_str(ii,:);
        end
    end
for ii = 1:xxcount
    if mymarker{ii}=='s'
            plot(T_Vm(NYstart(ii):NYend(ii)),Vm(NYstart(ii):NYend(ii)),...
                mymarker{ii},'markersize',markersize+2,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_Vm(ii)); %generate data symbol
            hold on;
    else
            plot(T_Vm(NYstart(ii):NYend(ii)),Vm(NYstart(ii):NYend(ii)),...
                mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_Vm(ii)); %generate data symbol
            hold on;
    end
end    
ytickformat('%.1f');
ylim([22.3 25]);
        legend (Legend_Author_Vm,'Location','bestoutside');
        legend('boxoff');
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
    ylabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../evaluation data/figure/cell volume')
        mkdir('../evaluation data/figure/cell volume')
    end
    print(gcf,'-dtiff','-r300','../evaluation data/figure/cell volume/Vm_sublimation.png')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 1.2
    [Year,Author,T_Vm,Vm] = textread('../evaluation data/cell volume/Vm_melting.txt','%s%s%f%f','headerlines',2);
    AA=Author(1);
NYstart=1;xxcount=1;
ndata=length(Year);
 Author_str=string(Author);
 YY=Year(1);
for i = 1:ndata
    Author_Y(i,:)=strcat(Year(i,:),'-',Author_str(i,:));
end
    for ii=2:ndata
        if ~strcmp(AA,Author(ii))||~strcmp(YY,Year(ii))
            NYend(xxcount)=ii-1;
            xxcount=xxcount+1;
            Legend_Author_Vm(xxcount-1,:)=Author_Y(ii-1,:);
            Author_alpha(xxcount-1,:)=Author_str(ii-1,:);
            NYstart(xxcount)=ii;
            AA=Author(ii);YY=Year(ii);
        end
        if ii==ndata
            NYend(xxcount)=ndata;
            Legend_Author_Vm(xxcount,:)=Author_Y(ii,:); 
            Author_alpha(xxcount,:)=Author_str(ii,:);
        end
    end
for ii = 1:xxcount
    if mymarker{ii}=='s'
            plot(T_Vm(NYstart(ii):NYend(ii)),Vm(NYstart(ii):NYend(ii)),...
                mymarker{ii},'markersize',markersize+2,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_Vm(ii)); %generate data symbol
            hold on;
    else
            plot(T_Vm(NYstart(ii):NYend(ii)),Vm(NYstart(ii):NYend(ii)),...
                mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_Vm(ii)); %generate data symbol
            hold on;
    end
end    
% ytickformat('%.1f');
% ylim([22.3 25]);
xlim([50 370]);
ylim([18.7 25]);
        legend (Legend_Author_Vm,'Location','bestoutside');
        legend('boxoff');
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
    ylabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../evaluation data/figure/cell volume')
        mkdir('../evaluation data/figure/cell volume')
    end
    print(gcf,'-dtiff','-r300','../evaluation data/figure/cell volume/Vm_melting.png')    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 1.3
        [Year,Author,T,Vm,p] = textread('../evaluation data/cell volume/Vm_high_pressure.txt','%s%s%f%f%f','headerlines',2);

ndata=length(T);xxcount=1; Nxstart = 1;AA=Author(1);YY=Year(1);
for ii = 2:ndata         
        if ~strcmp(AA,Author(ii))||~strcmp(YY,Year(ii))
            Nxend(xxcount) = ii-1;
            xxcount = xxcount + 1;
            Nxstart(xxcount) = ii;
            AA=Author(ii);YY=Year(ii);
        end
        if ii == ndata
            Nxend(xxcount) = ndata;
        end
end


for jj=1:xxcount
    ndata=length(T(Nxstart(jj):Nxend(jj)));yycount=1; Nystart(jj,yycount) = Nxstart(jj);TT=T(Nxstart(jj));
    for ii = 2:ndata   

            %if ~strcmp(AA,Author(Nxstart(jj)+ii-1))&&~strcmp(YY,Year(Nxstart(jj)+ii-1))
            if abs(T(Nxstart(jj)+ii-1) - TT)> 0.05
                Nyend(jj,yycount) = Nxstart(jj)+ii-2;
                yycount = yycount + 1;
                Nystart(jj,yycount) = Nxstart(jj)+ii-1;
                TT=T(Nxstart(jj)+ii);
            end
            if ii == ndata
                Nyend(jj,yycount) = ndata+Nxstart(jj)-1;
            end
    end     
end
[a1 a2]=size(Nystart);

figure(1)
for ii = 1:a1
    for jj=1:a2
        if Nystart(ii,jj)~=0
        ll1(ii,jj)=string(strcat(Author(Nystart(ii,jj)),Year(Nystart(ii,jj)),' - ',{' '},char(string(T(Nystart(ii,jj)))),{' '},'K'));            
        plot(p(Nystart(ii,jj):Nyend(ii,jj)),Vm(Nystart(ii,jj):Nyend(ii,jj)),...
        mymarker{ii},'markersize',markersize,'color',mycolor(jj,:),'linewidth',linewidth,'DisplayName',ll1(ii,jj)); 
    hold on;
        end
    end
end  
legend show
legend ('Location','southoutside','NumColumns',3);
legend('boxoff');
    xlabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth); 
    ylabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,10,800,750])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../evaluation data/figure/cell volume')
        mkdir('../evaluation data/figure/cell volume')
    end
    print(gcf,'-dtiff','-r300','../evaluation data/figure/cell volume/high pressure.png')  
    
figure(2)
for ii = 1:a1
    for jj=1:a2
        if Nystart(ii,jj)~=0
        ll1(ii,jj)=string(strcat(Author(Nystart(ii,jj)),Year(Nystart(ii,jj)),' - ',{' '},char(string(T(Nystart(ii,jj)))),{' '},'K'));            
        plot(p(Nystart(ii,jj):Nyend(ii,jj)),Vm(Nystart(ii,jj):Nyend(ii,jj)),...
        mymarker{ii},'markersize',markersize,'color',mycolor(jj,:),'linewidth',linewidth,'DisplayName',ll1(ii,jj)); 
    hold on;
        end
    end
end  
legend show
legend ('Location','southoutside','NumColumns',3);
legend('boxoff');
ylim([20 25]);
    xlabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth); 
    ylabel(['\it V\rm_m / (cm^3 mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,10,800,750])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../evaluation data/figure/cell volume')
        mkdir('../evaluation data/figure/cell volume')
    end
    print(gcf,'-dtiff','-r300','../evaluation data/figure/cell volume/high pressure_low_region.png')      
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 2.1
    [Year,Author,T_cp,cp] = textread('../evaluation data/heat capacity/cp_sublimation.txt','%s%s%f%f','headerlines',2);
    AA=Author(1);
NYstart=1;xxcount=1;
ndata=length(Year);
 Author_str=string(Author);
 YY=Year(1);
for i = 1:ndata
    Author_Y(i,:)=strcat(Year(i,:),'-',Author_str(i,:));
end
    for ii=2:ndata
        if ~strcmp(AA,Author(ii))||~strcmp(YY,Year(ii))
            NYend(xxcount)=ii-1;
            xxcount=xxcount+1;
            Legend_Author_cp(xxcount-1,:)=Author_Y(ii-1,:);
            Author_alpha(xxcount-1,:)=Author_str(ii-1,:);
            NYstart(xxcount)=ii;
            AA=Author(ii);YY=Year(ii);
        end
        if ii==ndata
            NYend(xxcount)=ndata;
            Legend_Author_cp(xxcount,:)=Author_Y(ii,:); 
            Author_alpha(xxcount,:)=Author_str(ii,:);
        end
    end
for ii = 1:xxcount
    if mymarker{ii}=='s'
            plot(T_cp(NYstart(ii):NYend(ii)),cp(NYstart(ii):NYend(ii)),...
                mymarker{ii},'markersize',markersize+2,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_cp(ii)); %generate data symbol
            hold on;
    else
            plot(T_cp(NYstart(ii):NYend(ii)),cp(NYstart(ii):NYend(ii)),...
                mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_cp(ii)); %generate data symbol
            hold on;
    end
end    
% ytickformat('%.1f');
% ylim([22.3 25]);
        legend (Legend_Author_cp,'Location','bestoutside');
        legend('boxoff');
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
    ylabel('\it c_p\rm / (J K^-^1 mol^-^1)','fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../evaluation data/figure/heat capacity')
        mkdir('../evaluation data/figure/heat capacity')
    end
    print(gcf,'-dtiff','-r300','../evaluation data/figure/heat capacity/cp_sublimation.png')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 3.1
    [Year,Author,T_alpha,alpha] = textread('../evaluation data/thermal expansion/alpha_sublimation.txt','%s%s%f%f','headerlines',2);
    AA=Author(1);
NYstart=1;xxcount=1;
ndata=length(Year);
 Author_str=string(Author);
 YY=Year(1);
for i = 1:ndata
    Author_Y(i,:)=strcat(Year(i,:),'-',Author_str(i,:));
end
    for ii=2:ndata
        if ~strcmp(AA,Author(ii))||~strcmp(YY,Year(ii))
            NYend(xxcount)=ii-1;
            xxcount=xxcount+1;
            Legend_Author_alpha(xxcount-1,:)=Author_Y(ii-1,:);
            Author_alpha(xxcount-1,:)=Author_str(ii-1,:);
            NYstart(xxcount)=ii;
            AA=Author(ii);YY=Year(ii);
        end
        if ii==ndata
            NYend(xxcount)=ndata;
            Legend_Author_alpha(xxcount,:)=Author_Y(ii,:); 
            Author_alpha(xxcount,:)=Author_str(ii,:);
        end
    end
for ii = 1:xxcount
    if mymarker{ii}=='s'
            plot(T_alpha(NYstart(ii):NYend(ii)),alpha(NYstart(ii):NYend(ii))*10^4,...
                mymarker{ii},'markersize',markersize+2,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_alpha(ii)); %generate data symbol
            hold on;
    else
            plot(T_alpha(NYstart(ii):NYend(ii)),alpha(NYstart(ii):NYend(ii))*10^4,...
                mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_alpha(ii)); %generate data symbol
            hold on;
    end
end    
% ytickformat('%.1f');
% ylim([22.3 25]);
        legend (Legend_Author_alpha,'Location','bestoutside');
        legend('boxoff');
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
ylabel('\it \alpha \rm\cdot10^4\rm / K^-^1','fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../evaluation data/figure/thermal expansion')
        mkdir('../evaluation data/figure/thermal expansion')
    end
    print(gcf,'-dtiff','-r300','../evaluation data/figure/thermal expansion/alpha_sublimation.png')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 4.1%Beta T
    [Year,Author,T_BetaT,BetaT] = textread('../evaluation data/bulk modulus/BetaT_sublimation.txt','%s%s%f%f','headerlines',2);
    AA=Author(1);
NYstart=1;xxcount=1;
ndata=length(Year);
 Author_str=string(Author);
 YY=Year(1);
for i = 1:ndata
    Author_Y(i,:)=strcat(Year(i,:),'-',Author_str(i,:));
end
    for ii=2:ndata
        if ~strcmp(AA,Author(ii))||~strcmp(YY,Year(ii))
            NYend(xxcount)=ii-1;
            xxcount=xxcount+1;
            Legend_Author_BetaT(xxcount-1,:)=Author_Y(ii-1,:);
            Author_BetaT(xxcount-1,:)=Author_str(ii-1,:);
            NYstart(xxcount)=ii;
            AA=Author(ii);YY=Year(ii);
        end
        if ii==ndata
            NYend(xxcount)=ndata;
            Legend_Author_BetaT(xxcount,:)=Author_Y(ii,:); 
            Author_BetaT(xxcount,:)=Author_str(ii,:);
        end
    end
for ii = 1:xxcount
    if mymarker{ii}=='s'
            plot(T_BetaT(NYstart(ii):NYend(ii)),BetaT(NYstart(ii):NYend(ii))*10^4,...
                mymarker{ii},'markersize',markersize+2,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_BetaT(ii)); %generate data symbol
            hold on;
    else
            plot(T_BetaT(NYstart(ii):NYend(ii)),BetaT(NYstart(ii):NYend(ii))*10^4,...
                mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_BetaT(ii)); %generate data symbol
            hold on;
    end
end    
% ytickformat('%.1f');
% ylim([22.3 25]);
        legend (Legend_Author_BetaT,'Location','bestoutside');
        legend('boxoff');
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
ylabel('\it \beta_T \rm\cdot10^4\rm / MPa^-^1','fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../evaluation data/figure/bulk modulus')
        mkdir('../evaluation data/figure/bulk modulus')
    end
    print(gcf,'-dtiff','-r300','../evaluation data/figure/bulk modulus/BetaT_sublimation.png')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 4.15%Beta S
    [Year,Author,T_BetaS,BetaS] = textread('../evaluation data/bulk modulus/BetaS_sublimation.txt','%s%s%f%f','headerlines',2);
    AA=Author(1);
NYstart=1;xxcount=1;
ndata=length(Year);
 Author_str=string(Author);
 YY=Year(1);
for i = 1:ndata
    Author_Y(i,:)=strcat(Year(i,:),'-',Author_str(i,:));
end
    for ii=2:ndata
        if ~strcmp(AA,Author(ii))||~strcmp(YY,Year(ii))
            NYend(xxcount)=ii-1;
            xxcount=xxcount+1;
            Legend_Author_BetaS(xxcount-1,:)=Author_Y(ii-1,:);
            Author_BetaS(xxcount-1,:)=Author_str(ii-1,:);
            NYstart(xxcount)=ii;
            AA=Author(ii);YY=Year(ii);
        end
        if ii==ndata
            NYend(xxcount)=ndata;
            Legend_Author_BetaS(xxcount,:)=Author_Y(ii,:); 
            Author_BetaS(xxcount,:)=Author_str(ii,:);
        end
    end
for ii = 1:xxcount
    if mymarker{ii}=='s'
            plot(T_BetaS(NYstart(ii):NYend(ii)),BetaS(NYstart(ii):NYend(ii))*10^4,...
                mymarker{ii},'markersize',markersize+2,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_BetaS(ii)); %generate data symbol
            hold on;
    else
            plot(T_BetaS(NYstart(ii):NYend(ii)),BetaS(NYstart(ii):NYend(ii))*10^4,...
                mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_BetaS(ii)); %generate data symbol
            hold on;
    end
end    
% ytickformat('%.1f');
% ylim([22.3 25]);
ytickformat('%.1f');
        legend (Legend_Author_BetaS,'Location','bestoutside');
        legend('boxoff');
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
ylabel('\it \beta_S \rm\cdot10^4\rm / MPa^-^1','fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../evaluation data/figure/bulk modulus')
        mkdir('../evaluation data/figure/bulk modulus')
    end
    print(gcf,'-dtiff','-r300','../evaluation data/figure/bulk modulus/BetaS_sublimation.png')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 5%sublimation
%         r1=RP.REFPROPdll("argon","TP","TTRP",int8(2),iMass,iFlag,300,1,{1.0});
%         o1 = double(r1.Output);
        T_triple=83.806 ;   
        r1=RP.REFPROPdll("argon","TP","PTRP",int8(2),iMass,iFlag,300,1,{1.0});
        o1 = double(r1.Output);
        p_triple=6.889100000000001e+04;      %unit Pa        
    [Year,Author,T_sub,p_sub] = textread('../evaluation data/sublimation/sublimation.txt','%s%s%f%f','headerlines',2);
    AA=Author(1);
NYstart=1;xxcount=1;
ndata=length(Year);
 Author_str=string(Author);
 YY=Year(1);
for i = 1:ndata
    Author_Y(i,:)=strcat(Year(i,:),'-',Author_str(i,:));
end
    for ii=2:ndata
        if ~strcmp(AA,Author(ii))||~strcmp(YY,Year(ii))
            NYend(xxcount)=ii-1;
            xxcount=xxcount+1;
            Legend_Author_sub(xxcount-1,:)=Author_Y(ii-1,:);
            Author_sub(xxcount-1,:)=Author_str(ii-1,:);
            NYstart(xxcount)=ii;
            AA=Author(ii);YY=Year(ii);
        end
        if ii==ndata
            NYend(xxcount)=ndata;
            Legend_Author_sub(xxcount,:)=Author_Y(ii,:); 
            Author_sub(xxcount,:)=Author_str(ii,:);
        end
    end
         T = [T_sub',T_triple];
         p = [p_sub',p_triple];
    fun = @(b,T)T_triple./T.*(b(1).*(1-T./T_triple).^1+b(2).*(1-T./T_triple).^1.5+b(3).*(1-T./T_triple).^5);
    options=optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt','TolX', 1e-25, 'TolFun', 1e-25,'MaxFunctionEvaluations',2000000, 'MaxIterations', 200000);
    lb = [-Inf, -5, -Inf];
    ub = []; 
     b0 = [-11.00316274 -0.965384697 -1.471436611];	

    b_fit = lsqcurvefit(fun,b0,T,log(p./p_triple),lb,ub,options);
        p_fit=exp(fun(b_fit,T))*p_triple;
    p_dev = 100 .* (p - p_fit)./p_fit;
for ii = 1:xxcount
    if mymarker{ii}=='s'
            plot(T_sub(NYstart(ii):NYend(ii)),p_sub(NYstart(ii):NYend(ii)),...
                mymarker{ii},'markersize',markersize+2,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_sub(ii)); %generate data symbol
            hold on;
    else
            plot(T_sub(NYstart(ii):NYend(ii)),p_sub(NYstart(ii):NYend(ii)),...
                mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_sub(ii)); %generate data symbol
            hold on;
    end
end    
        plot(T_triple,p_triple,...
            mymarker{19},'markersize',markersize,'color',mycolor(3,:),'linewidth',linewidth,'DisplayName',"Triple point");
        hold on
    jj = 1;        
    for T_plot = 0:0.005:T_triple
           p_plot_fit(jj) =  poly_fit2(b_fit,T_plot,T_triple,p_triple);
        jj = jj + 1;
    end
        T_plot =0:0.005:T_triple;
    plot(T_plot,p_plot_fit,...
           'color',mycolor(1,:),'linewidth',linewidth,'DisplayName',"polynomial fitting");    
% ytickformat('%.1f');
ylim([10^-9 10^5]);
% ytickformat('%.1f');
set(gca, 'YScale', 'log')
        legend (Legend_Author_sub,'Location','bestoutside');
        legend('boxoff');
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
ylabel('\it p \rm/ Pa','fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../evaluation data/figure/sublimation')
        mkdir('../evaluation data/figure/sublimation')
    end
    print(gcf,'-dtiff','-r300','../evaluation data/figure/sublimation/sublimation.png')
    
    figure(2)
for ii = 1:xxcount
    if mymarker{ii}=='s'
            plot(T_sub(NYstart(ii):NYend(ii)),p_sub(NYstart(ii):NYend(ii)),...
                mymarker{ii},'markersize',markersize+2,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_sub(ii)); %generate data symbol
            hold on;
    else
            plot(T_sub(NYstart(ii):NYend(ii)),p_sub(NYstart(ii):NYend(ii)),...
                mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_sub(ii)); %generate data symbol
            hold on;
    end
end    
        plot(T_triple,p_triple,...
            mymarker{19},'markersize',markersize,'color',mycolor(3,:),'linewidth',linewidth,'DisplayName',"Triple point");
        hold on
    jj = 1;        
    for T_plot = 0:0.005:T_triple
           p_plot_fit(jj) =  poly_fit2(b_fit,T_plot,T_triple,p_triple);
        jj = jj + 1;
    end
        T_plot =0:0.005:T_triple;
    plot(T_plot,p_plot_fit,...
           'color',mycolor(1,:),'linewidth',linewidth,'DisplayName',"polynomial fitting");    
% ytickformat('%.1f');
% ylim([22.3 25]);
% ytickformat('%.1f');
set(gca, 'YScale', 'log')
ylim([1000 10^6]);
        legend (Legend_Author_sub,'Location','bestoutside');
        legend('boxoff');
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
ylabel('\it p \rm/ Pa','fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../evaluation data/figure/sublimation')
        mkdir('../evaluation data/figure/sublimation')
    end
    print(gcf,'-dtiff','-r300','../evaluation data/figure/sublimation/sublimation_zoom_in.png') 
    
    figure(3)
 for ii = 1:xxcount
    if mymarker{ii}=='s'
            plot(T_sub(NYstart(ii):NYend(ii)),p_dev(NYstart(ii):NYend(ii)),...
                mymarker{ii},'markersize',markersize+2,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_sub(ii)); %generate data symbol
            hold on;
    else
            plot(T_sub(NYstart(ii):NYend(ii)),p_dev(NYstart(ii):NYend(ii)),...
                mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_sub(ii)); %generate data symbol
            hold on;
    end
 end      
      p_fit_triple=exp(fun(b_fit,T_triple))*p_triple;
      p_dev_triple = 100 * (p_triple-p_fit_triple)/p_fit_triple;
         plot(T_triple,p_dev_triple,...
            mymarker{19},'markersize',markersize,'color',mycolor(3,:),'linewidth',linewidth,'DisplayName',"Triple point"); 
        hold on
            c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
        set( get( get( c(1), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );%not show c(1) legend
        xlim([0 90]);
        xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);  
        ylabel(['10^2\cdot',char(8201),'(\itp\rm_,_e_x_p',char(hex2dec('2212')),'\it p\rm_,_c_a_l_c) /\it p\rm_,_c_a_l_c'],'fontsize',fontsize,'linewidth',linewidth)  
     set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../evaluation data/figure/sublimation')
        mkdir('../evaluation data/figure/sublimation')
    end
    print(gcf,'-dtiff','-r300','../evaluation data/figure/sublimation/dev.png')        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 6%melting
    [Year,Author,T_melt,p_melt] = textread('../evaluation data/melting/melting.txt','%s%s%f%f','headerlines',2);
            T_triple=83.806 ;        p_triple=6.889100000000001e+04/10^6;      %unit MPa  
    AA=Author(1);
NYstart=1;xxcount=1;
ndata=length(Year);
 Author_str=string(Author);
 YY=Year(1);
for i = 1:ndata
    Author_Y(i,:)=strcat(Year(i,:),'-',Author_str(i,:));
end
    for ii=2:ndata
        if ~strcmp(AA,Author(ii))||~strcmp(YY,Year(ii))
            NYend(xxcount)=ii-1;
            xxcount=xxcount+1;
            Legend_Author_melt(xxcount-1,:)=Author_Y(ii-1,:);
            Author_melt(xxcount-1,:)=Author_str(ii-1,:);
            NYstart(xxcount)=ii;
            AA=Author(ii);YY=Year(ii);
        end
        if ii==ndata
            NYend(xxcount)=ndata;
            Legend_Author_melt(xxcount,:)=Author_Y(ii,:); 
            Author_melt(xxcount,:)=Author_str(ii,:);
        end
    end
         T = [T_melt',T_triple];
         p = [p_melt',p_triple];    
     fun = @(b,T)(1+b(1).*(T./T_triple-1).^(b(2))+b(3).*(T./T_triple-1).^b(4));   
%      options=optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt','TolX', 1e-25, 'TolFun', 1e-25,'MaxFunctionEvaluations',2000000, 'MaxIterations', 200000);
%     lb = [];
%     ub = []; 
%      b0 = [3160 1 3040 1];
%     b_fit = lsqcurvefit(fun,b0,T,(p./p_triple),lb,ub,options);
b_fit = [1506.5415	1.73136	4677.1597	0.9849295];

% b_fit = [3251.003 1.473837 2963.07 0.881421];
        p_fit=fun(b_fit,T)*p_triple;
    p_dev = 100 .* (p - p_fit)./p_fit;    
for ii = 1:xxcount
    if mymarker{ii}=='s'
            plot(T_melt(NYstart(ii):NYend(ii)),p_melt(NYstart(ii):NYend(ii)),...
                mymarker{ii},'markersize',markersize+2,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_melt(ii)); %generate data symbol
            hold on;
    else
            plot(T_melt(NYstart(ii):NYend(ii)),p_melt(NYstart(ii):NYend(ii)),...
                mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_melt(ii)); %generate data symbol
            hold on;
    end
end    
% ytickformat('%.1f');
% ylim([22.3 25]);
% ytickformat('%.1f');
% set(gca, 'YScale', 'log')
        plot(T_triple,p_triple,...
            mymarker{19},'markersize',markersize,'color',mycolor(3,:),'linewidth',linewidth,'DisplayName',"Triple point");
        hold on
    jj = 1;        
    for T_plot = T_triple:0.005:3200
           p_plot_fit(jj) =  poly_fit3(b_fit,T_plot,T_triple,p_triple);
        jj = jj + 1;
    end
        T_plot = T_triple:0.005:3200;
    plot(T_plot,p_plot_fit,...
           'color',mycolor(1,:),'linewidth',linewidth,'DisplayName',"polynomial fitting");    
        legend (Legend_Author_melt,'Location','bestoutside');
        legend('boxoff');
        xlim([0 750]);
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
ylabel('\it p \rm/ MPa','fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../evaluation data/figure/melting')
        mkdir('../evaluation data/figure/melting')
    end
    print(gcf,'-dtiff','-r300','../evaluation data/figure/melting/melting.png')
    
%     figure(2)
% for ii = 1:xxcount
%     if mymarker{ii}=='s'
%             plot(T_melt(NYstart(ii):NYend(ii)),p_melt(NYstart(ii):NYend(ii)),...
%                 mymarker{ii},'markersize',markersize+2,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_melt(ii)); %generate data symbol
%             hold on;
%     else
%             plot(T_melt(NYstart(ii):NYend(ii)),p_melt(NYstart(ii):NYend(ii)),...
%                 mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_melt(ii)); %generate data symbol
%             hold on;
%     end
% end    
% % ytickformat('%.1f');
% % ylim([22.3 25]);
% % ytickformat('%.1f');
% % set(gca, 'YScale', 'log')
% ylim([0.05 7000]);
%         legend (Legend_Author_melt,'Location','bestoutside');
%         legend('boxoff');
%     xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
% ylabel('\it p \rm/ MPa','fontsize',fontsize,'linewidth',linewidth);
%     % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
%     % set(gcf,'position',[300,300,500,340])
%     set(gcf,'position',[300,100,650,540])
%     set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
%     if ~isfolder('../evaluation data/figure/melting')
%         mkdir('../evaluation data/figure/melting')
%     end
%     print(gcf,'-dtiff','-r300','../evaluation data/figure/melting/melting_zoom_in.png')   
  
    figure(2)
 for ii = 1:xxcount
    if mymarker{ii}=='s'
            plot(T_melt(NYstart(ii):NYend(ii)),p_dev(NYstart(ii):NYend(ii)),...
                mymarker{ii},'markersize',markersize+2,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_melt(ii)); %generate data symbol
            hold on;
    else
            plot(T_melt(NYstart(ii):NYend(ii)),p_dev(NYstart(ii):NYend(ii)),...
                mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_melt(ii)); %generate data symbol
            hold on;
    end
 end      
      p_fit_triple=fun(b_fit,T_triple)*p_triple;
      p_dev_triple = 100 * (p_triple-p_fit_triple)/p_fit_triple;
         plot(T_triple,p_dev_triple,...
            mymarker{19},'markersize',markersize,'color',mycolor(3,:),'linewidth',linewidth,'DisplayName',"Triple point"); 
        hold on
            c(1)=plot([0,1000],[0,0],'k','linewidth',linewidth-0.2);
        set( get( get( c(1), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );%not show c(1) legend
        legend (Legend_Author_melt,'Location','bestoutside');
        legend('boxoff');
        xlim([0 800]);
        xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth);  
        ylabel(['10^2\cdot',char(8201),'(\itp\rm_,_e_x_p',char(hex2dec('2212')),'\it p\rm_,_c_a_l_c) /\it p\rm_,_c_a_l_c'],'fontsize',fontsize,'linewidth',linewidth)  
     set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../evaluation data/figure/melting')
        mkdir('../evaluation data/figure/melting')
    end
    print(gcf,'-dtiff','-r300','../evaluation data/figure/melting/dev.png') 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 7.1%enthalpy of sublimation

    [Year,Author,T_H_sub,H_sub] = textread('../evaluation data/enthalpy/enthalpy of sublimation.txt','%s%s%f%f','headerlines',2);
    AA=Author(1);
NYstart=1;xxcount=1;
ndata=length(Year);
 Author_str=string(Author);
 YY=Year(1);
for i = 1:ndata
    Author_Y(i,:)=strcat(Year(i,:),'-',Author_str(i,:));
end
    for ii=2:ndata
        if ~strcmp(AA,Author(ii))||~strcmp(YY,Year(ii))
            NYend(xxcount)=ii-1;
            xxcount=xxcount+1;
            Legend_Author_H_sub(xxcount-1,:)=Author_Y(ii-1,:);
            Author_H_sub(xxcount-1,:)=Author_str(ii-1,:);
            NYstart(xxcount)=ii;
            AA=Author(ii);YY=Year(ii);
        end
        if ii==ndata
            NYend(xxcount)=ndata;
            Legend_Author_H_sub(xxcount,:)=Author_Y(ii,:); 
            Author_H_sub(xxcount,:)=Author_str(ii,:);
        end
    end
for ii = 1:xxcount
    if mymarker{ii}=='s'
            plot(T_H_sub(NYstart(ii):NYend(ii)),H_sub(NYstart(ii):NYend(ii)),...
                mymarker{ii},'markersize',markersize+2,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_H_sub(ii)); %generate data symbol
            hold on;
    else
            plot(T_H_sub(NYstart(ii):NYend(ii)),H_sub(NYstart(ii):NYend(ii)),...
                mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_H_sub(ii)); %generate data symbol
            hold on;
    end
end    
% ytickformat('%.1f');
% ylim([1.1 3]);
ytickformat('%.2f');
% set(gca, 'YScale', 'log')
        legend (Legend_Author_H_sub,'Location','bestoutside');
        legend('boxoff');
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
ylabel(['\it H \rm/ kJ mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../evaluation data/figure/enthalpy')
        mkdir('../evaluation data/figure/enthalpy')
    end
    print(gcf,'-dtiff','-r300','../evaluation data/figure/enthalpy/enthalpy of sublimation.png')
     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 7.2%enthalpy of fusion
    [Year,Author,T_H_melt,H_melt] = textread('../evaluation data/enthalpy/enthalpy of fusion.txt','%s%s%f%f','headerlines',2);
    AA=Author(1);
NYstart=1;xxcount=1;
ndata=length(Year);
 Author_str=string(Author);
 YY=Year(1);
for i = 1:ndata
    Author_Y(i,:)=strcat(Year(i,:),'-',Author_str(i,:));
end
    for ii=2:ndata
        if ~strcmp(AA,Author(ii))||~strcmp(YY,Year(ii))
            NYend(xxcount)=ii-1;
            xxcount=xxcount+1;
            Legend_Author_H_melt(xxcount-1,:)=Author_Y(ii-1,:);
            Author_H_melt(xxcount-1,:)=Author_str(ii-1,:);
            NYstart(xxcount)=ii;
            AA=Author(ii);YY=Year(ii);
        end
        if ii==ndata
            NYend(xxcount)=ndata;
            Legend_Author_H_melt(xxcount,:)=Author_Y(ii,:); 
            Author_H_melt(xxcount,:)=Author_str(ii,:);
        end
    end
for ii = 1:xxcount
    if mymarker{ii}=='s'
            plot(T_H_melt(NYstart(ii):NYend(ii)),H_melt(NYstart(ii):NYend(ii)),...
                mymarker{ii},'markersize',markersize+2,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_H_melt(ii)); %generate data symbol
            hold on;
    else
            plot(T_H_melt(NYstart(ii):NYend(ii)),H_melt(NYstart(ii):NYend(ii)),...
                mymarker{ii},'markersize',markersize,'color',mycolor(ii,:),'linewidth',linewidth,'DisplayName',Legend_Author_H_melt(ii)); %generate data symbol
            hold on;
    end
end    
% ytickformat('%.1f');
ylim([1.1 3]);
ytickformat('%.1f');
% set(gca, 'YScale', 'log')
        legend (Legend_Author_H_melt,'Location','bestoutside');
        legend('boxoff');
    xlabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
ylabel(['\it H \rm/ kJ mol^',char(hex2dec('2212')),'^1)'],'fontsize',fontsize,'linewidth',linewidth);
    % annotation('textbox', [0, 1, 0, 0], 'string', '\it\bfa','FontSize',fontsize,'FontName','Times');
    % set(gcf,'position',[300,300,500,340])
    set(gcf,'position',[300,100,650,540])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');    %  print(gcf,'-dtiff','-r600','Volume/Volumelegendoff.png')
    if ~isfolder('../evaluation data/figure/enthalpy')
        mkdir('../evaluation data/figure/enthalpy')
    end
    print(gcf,'-dtiff','-r300','../evaluation data/figure/enthalpy/enthalpy of fusion.png')
     
end