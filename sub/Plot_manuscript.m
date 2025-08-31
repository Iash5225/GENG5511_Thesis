%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;path(path,[pwd,'\..\..\SUB']);
%fontsize = 20; markersize = 11; linewidth = 0.8;
fontsize = 14; markersize = 8; linewidth = 0.9;
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
choice = 1;%phase envelope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 1
    [T_VLE,p_VLE] = textread('../phase envelope data/VLE.txt','%f%f','headerlines',2);
    plot(p_VLE,T_VLE,...
       'color',mycolor(1,:),'linewidth',linewidth); %generate data symbol
   
    set(gca, 'XScale', 'log')

    hold on;
    kk = 1;
    for T_melting = 83.806:1000.806
        p_melting(kk,:) = pmelt(T_melting);
        kk = kk + 1;
    end
    T_melting = 83.806:1000.806;
%     [T_melting,p_melting] = textread('../argon/data/phase envelope/melting.txt','%f%f','headerlines',2);    
    plot(p_melting,T_melting,...
       'color',mycolor(1,:),'linewidth',linewidth); %generate data symbol    
   hold on;

    kk = 1;
    for T_sublimation = 0.806:83.806
        p_sublimation(kk,:) = psub(T_sublimation);
        kk = kk + 1;
    end    
    T_sublimation = 0.806:83.806;
%     [T_sublimation,p_sublimation,~] = textread('../argon/data/phase envelope/sublimation.txt','%f%f%f','headerlines',2);    
    plot(p_sublimation,T_sublimation,...
       'color',mycolor(1,:),'linewidth',linewidth); %generate data symbol       
   hold on;
       
T_crit=150.69;p_crit=4.863;
plot(p_crit,T_crit,mymarker{1},'color',mycolor(1,:),'linewidth',linewidth,'markersize',markersize);
T_triple1=83.806;p_triple1=0.068891;
plot(p_triple1,T_triple1,mymarker{2},'color',mycolor(3,:),'linewidth',linewidth,'markersize',9);
    xlim([10^-10 10^3]);
    ylim([0 250]);
    text(10^-3,38,'Solid','fontsize',14);
    text(10^-7,125,'Vapour','fontsize',14);
    text(2,115,'Liquid','fontsize',14);    
    ylabel('\it T\rm / K','fontsize',fontsize,'linewidth',linewidth); 
    xlabel('\it p\rm / MPa','fontsize',fontsize,'linewidth',linewidth);    
     set(gcf,'position',[300,300,500,340])
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'FontName','Times');
    if ~isfolder('../../figures/phase envelope/')
     mkdir('../../figures/phase envelope/')
    end    
     print(gcf,'-dtiff','-r600','../../figures/phase envelope/argon.png')   
    
end