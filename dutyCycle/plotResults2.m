
clearvars
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork  
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions
ft_defaults
dataDir = '/Users/marleneponcet/Documents/data/dutyCycle/Axx/';


listData = dir([dataDir '*.mat']);
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};


for ss = 1:length(listData)
    clear Axx;
    load([dataDir listData(ss).name]);
    dataSbj(:,ss) = Axx;
end

% keepSbj = [1:size(dataSbj,2)]; 
keepSbj = [1:5 7:11 13:15]; % (there is no S1) reject S7 and S12: just noise?

[avData, proj_Amp] = averageAxxWithStd(dataSbj(:,keepSbj));

col={'b','r','g'};

pickElec = 23; % 23 Oz, 9, B6=38, 16=best SNR


%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Topography
% flicker
figure; hold on;
for cond=1:15 
    subplot(3,5,cond); hold on;
    plotTopo(avData(cond).amp(:,avData(cond).i1f1),cfg.layout);
    colorbar
    if cond < 6
        caxis([0 1.5]);
    else
        caxis([0 2.5]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,'topoFlickF1.jpg')
% motion
position = [11 7 3 9 15 13 8];
figure; hold on;
for cond=16:length(avData) 
    subplot(3,5,position(cond-15)); hold on;
    plotTopo(avData(cond).amp(:,avData(cond).i1f1*2-1),cfg.layout);
    colorbar
    if cond < 6
        caxis([0 1.5]);
    else
        caxis([0 2.5]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,'topoMotionF2.jpg')

%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SSVEP
for cond=1:15 
   avAmp(cond)=avData(cond).amp(pickElec,(avData(cond).i1f1));
%    avProj(cond)=mean(proj_Amp(pickElec,avData(cond).i1f1,cond,:));
   semProj(cond) = std(proj_Amp(pickElec,avData(cond).i1f1,cond,:)) / sqrt(size(dataSbj,2));
end
for cond=16:22 % for motion take the 2f1
   avAmp(cond)=avData(cond).amp(pickElec,(avData(cond).i1f1*2-1));
   semProj(cond) = std(proj_Amp(pickElec,avData(cond).i1f1*2-1,cond,:)) / sqrt(size(dataSbj,2));
end
% % compute noise level per freq
% for cond=1:15 % exclude motion conditions
%    baseline(cond) = (avData(cond).amp(pickElec,avData(cond).i1f1-1)+avData(cond).amp(pickElec,avData(cond).i1f1+1)) / 2;
% end
% compute noise level all cond pooled
for cond = 1:length(avData)
    baseline(cond) = (avData(cond).amp(pickElec,avData(cond).i1f1-1)+avData(cond).amp(pickElec,avData(cond).i1f1+1)) / 2;
end
temp=reshape(baseline(1:20),4,5);
avBaseline = mean(temp);
avBaseline(3) = mean([temp(1:4,3)' baseline(21) baseline(22)]); % add the 50% DC for other freq


figure;hold on;
for freq=1:3
    errorbar(avAmp((freq-1)*5+1:freq*5),semProj((freq-1)*5+1:freq*5),['.-' col{freq}],'MarkerSize',40,'Linewidth',2);
end
plot(avBaseline,'k','Linewidth',1)
errorbar(3,avAmp(18),semProj(18),['^:' col{1}],'MarkerSize',15,'Linewidth',2)
errorbar([2 3 4],avAmp([17 22 19]),semProj([17 22 19]),['^:' col{2}],'MarkerSize',15,'Linewidth',1)
errorbar([1 3 5],avAmp([16 21 20]),semProj([17 22 19]),['^:' col{3}],'MarkerSize',15,'Linewidth',1)
xlim([0.5 5.5]);
xticks([1:1:5]);
xticklabels({'12.5','25','50','75','87.5'})
ylim([0 3])
legend('10','5','2.5','noise level','Location','Best')
xlabel('Duty Cycle')
ylabel('SSVEP amplitude')
set(gca,'FontSize',15)
saveas(gcf,'ampDC.png')
saveas(gcf,'ampDC.eps','epsc')



%%%%%%%%%%%%%%%%%%%%%
%%%% RATINGS
load('fullRatings.mat')
static = mean(tabStatic,3); moving = mean(tabMot,3); 
for fq=1:3
    for dc=1:5
        statSEM(fq,dc) = std(tabStatic(fq,dc,:))/sqrt(6);
        movSEM(fq,dc) = std(tabMot(fq,dc,:))/sqrt(6);
    end
end
figure;hold on;
for fq=1:3
    errorbar(static(fq,:),statSEM(1,:),['.-'  col{fq}],'MarkerSize',40,'LineWidth',2)
end
errorbar(3,moving(7),movSEM(7),['^'  col{1}],'MarkerSize',15,'Linewidth',1)
errorbar([2 3 4],moving([5 8 11]),movSEM([5 8 11]),['^:'  col{2}],'MarkerSize',15,'Linewidth',1)
errorbar([1 3 5],moving([3 9 15]),movSEM([3 9 15]),['^:'  col{3}],'MarkerSize',15,'Linewidth',1)
xlim([0.5 5.5]);
xticks([1:1:5]);
ylim([0 3])
xticklabels({'12.5','25','50','75','87.5'})
legend('10','5','2.5','Location','Best')
xlabel('Duty Cycle')
ylabel('motion ratings')
set(gca,'FontSize',15)
saveas(gcf,'ratingsDC.png')
saveas(gcf,'ratingsDC.eps','epsc')



%%%%%%%%%%%%%%%%%%%%%
%%%%%%% CORREL
figure;hold on; 
y = [static(2,:) static(3,:)];
x = avAmp(6:15);
p = polyfit(x,y,1); 
f = polyval(p,x); 
plot(x,y,'.',x,f,'-') 
y = moving([3 5 7 11 15 9 8]);
x = avAmp(16:22);
p = polyfit(x,y,1); 
f = polyval(p,x); 
plot(x,y,'^',x,f,'-') 
ylim([0 3])
legend('flickering','','moving')
ylabel('motion rating')
xlabel('SSVEP amplitude')
saveas(gcf,'correlSSVEPratings.png')
saveas(gcf,'correlSSVEPratings.eps','epsc')



