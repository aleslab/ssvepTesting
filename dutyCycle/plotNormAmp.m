% normalisation applied on each sbj then average
clearvars

% load individual data
dataDir = '/Volumes/Amrutam/Marlene/JUSTIN/DutyCycle/data/Axx/';
listData = dir([dataDir '*.mat']);
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};

keepSbj = [1:5 7:11 13:15 17:18 20]; 
dcVal = [12.5 25 50 75 87.5];

for ss = 1:length(keepSbj)
    clear Axx;
    load([dataDir listData(keepSbj(ss)).name]);
    dataSbj(:,ss) = Axx;
end

for ss = 1:length(keepSbj)
    for cond = 1:22
        normAmp(:,cond,ss) = dataSbj(cond,ss).normFq;
    end
end
avAmp = mean(normAmp,3);
semAmp = std(normAmp,[],3) / sqrt(size(normAmp,3));

col={'b','r','g'};

chan = 23;
figure;hold on;
for freq=1:3
    errorbar(dcVal,avAmp(chan,(freq-1)*5+1:freq*5),semAmp(chan,(freq-1)*5+1:freq*5),['.-' col{freq}],'MarkerSize',40,'Linewidth',2);
end
% line([dcVal(1) dcVal(5)], [min(baseline) min(baseline)],'Linewidth',2,'Color','k','LineStyle','--')
% line([dcVal(1) dcVal(5)], [mean(baseline) mean(baseline)],'Linewidth',2,'Color','k','LineStyle','-')
% line([dcVal(1) dcVal(5)], [max(baseline) max(baseline)],'Linewidth',2,'Color','k','LineStyle','--')
errorbar(50,avAmp(chan,18),semAmp(chan,18),['^:' col{1}],'MarkerSize',15,'Linewidth',2)
errorbar([25 50 75],avAmp(chan,[17 22 19]),semAmp(chan,[17 22 19]),['^:' col{2}],'MarkerSize',15,'Linewidth',2)
errorbar([12.5 50 87.5],avAmp(chan,[16 21 20]),semAmp(chan,[16 21 20]),['^:' col{3}],'MarkerSize',15,'Linewidth',2)
xlim([0 100]);
xticks([0:12.5:100]);
% ylim([0 3])
legend('10Hz','5Hz','2.5Hz','Location','Best')
xlabel('Duty Cycle')
ylabel('norm SSVEP amplitude')
set(gca,'FontSize',15)
title('Oz')
saveas(gcf,['figures' filesep 'ampDCNormAllFqSbj.jpg'])



% flicker
figure; hold on;
colormap('hot')
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(avAmp(:,cond),cfg.layout);
    colorbar
    if cond>10
        caxis([0 3.5]);
    elseif cond<6
        caxis([0 2]);
    elseif cond>5 && cond<11
        caxis([0 3]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoFlickNormAllFqSbj.jpg'])
% motion
position = [11 7 3 9 15 13 8];
figure; hold on;
colormap('hot')
for cond=16:22
    subplot(3,5,position(cond-15)); hold on;
    plotTopo(avAmp(:,cond),cfg.layout);
    colorbar
    if ismember(cond,[16 21 20])
        caxis([0 3.5]);
    elseif ismember(cond,[17 22 19])
        caxis([0 3]);
    else
        caxis([0 2]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoMotionNormAllFqSbj.jpg'])

