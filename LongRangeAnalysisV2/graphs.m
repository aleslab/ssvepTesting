clearvars;
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
ft_defaults


condition = 10;
channel = 23;

% time-frequency plot
plot(squeeze(gpData(condition).freq),squeeze(gpData(condition).amp(channel,:)))
bar(squeeze(gpData(condition).freq((gpData(condition).freq<30 & gpData(condition).freq>0))),squeeze(gpData(condition).amp(channel,2:60)))

figure;
condition = 9;
bar(squeeze(gpData(condition).freq((gpData(condition).freq<30 & gpData(condition).freq>0))),squeeze(gpData(condition).amp(channel,2:60)),0.8)
hold on;
condition = 10;
bar(squeeze(gpData(condition).freq((gpData(condition).freq<30 & gpData(condition).freq>0))),squeeze(gpData(condition).amp(channel,2:60)),0.3,'FaceColor','r')

ft_topoplotTFR(cfg,squeeze(gpData(condition).amp(2:end,10)))
plotTopo(squeeze(gpData(condition).amp(:,10)),cfg.layout)


ft_topoplotTFR(cfg,squeeze(gpData(condition).amp(2:end,10)))



% average topo for all condition at 5 Hz
freqToPlot1 = find(gpData(1).freq == 2.5);
freqToPlot2 = find(gpData(1).freq == 5);
for cond=2:2:8 % only take the original signal (not the predictions)
    dataCosLR(:,:,cond) = gpData(cond).cos;
    dataSinLR(:,:,cond) = gpData(cond).sin;
    dataCosSR(:,:,cond+8) = gpData(cond+8).cos;
    dataSinSR(:,:,cond+8) = gpData(cond+8).sin;
end
allAmpLR = sqrt((mean(dataCosLR(:,:,:),3)).^2 + (mean(dataSinLR(:,:,:),3)).^2);
allAmpSR = sqrt((mean(dataCosSR(:,:,:),3)).^2 + (mean(dataSinSR(:,:,:),3)).^2);
figure;
subplot(2,2,1); hold on;
plotTopo(allAmpLR(:,freqToPlot1),cfg.layout); colorbar;
caxis([0 0.5])
title('LR 2.5Hz')
subplot(2,2,2); 
plotTopo(allAmpLR(:,freqToPlot2),cfg.layout); colorbar;
caxis([0 0.5])
title('LR 5Hz')
subplot(2,2,3); 
plotTopo(allAmpSR(:,freqToPlot1),cfg.layout); colorbar;
caxis([0 0.5])
title('SR 2.5Hz')
subplot(2,2,4); 
plotTopo(allAmpSR(:,freqToPlot2),cfg.layout); colorbar;
caxis([0 0.5])
title('SR 5Hz')

maxFreq = find(gpData(1).freq == 15);
minFreq = find(gpData(1).freq==0.5);
figure;bar(allAmpLR(23, minFreq: maxFreq))

% choose electrode A23 = Oz
channel = 44;
limFq = 1.5;
limERPmin = -2.5;
limERPmax = 3.5;

limFqDiff = 1;
limERPminDiff = -2;
limERPmaxDiff = 2.5;

% power spectrum & waveform
% LR
figure; hold on;
for condition=1:2:8
    subplot(4,2,condition);hold on;
    plot(squeeze(gpData(condition).time),squeeze(gpData(condition).wave(channel,:)),'b','LineWidth',2)
    plot(squeeze(gpData(condition).time),squeeze(gpData(condition+1).wave(channel,:)),'r','LineWidth',2)
%     plot(squeeze(gpData(condition).time),squeeze((gpData(condition).wave(channel,:))-gpData(condition+1).wave(channel,:)),'g','LineWidth',2)    
    ylim([limERPmin limERPmax])
end
for condition=1:2:8
    subplot(4,2,condition+1);hold on;
    bar(squeeze(gpData(condition).freq((gpData(condition).freq<15.5 & gpData(condition).freq>0))),squeeze(gpData(condition).amp(channel,2:31)),0.8)
    bar(squeeze(gpData(condition).freq((gpData(condition).freq<15.5 & gpData(condition).freq>0))),squeeze(gpData(condition+1).amp(channel,2:31)),0.3,'FaceColor','r')
    ylim([0 limFq])
end

% SR
figure; hold on;
for condition=1:2:8
    subplot(4,2,condition);hold on;
    plot(squeeze(gpData(condition).time),squeeze(gpData(condition+8).wave(channel,:)),'b','LineWidth',2)
    plot(squeeze(gpData(condition).time),squeeze(gpData(condition+1+8).wave(channel,:)),'r','LineWidth',2)
    ylim([limERPmin limERPmax])
end
for condition=1:2:8
    subplot(4,2,condition+1);hold on;
    bar(squeeze(gpData(condition).freq((gpData(condition).freq<15.5 & gpData(condition).freq>0))),squeeze(gpData(condition+8).amp(channel,2:31)),0.8)
    bar(squeeze(gpData(condition).freq((gpData(condition).freq<15.5 & gpData(condition).freq>0))),squeeze(gpData(condition+1+8).amp(channel,2:31)),0.3,'FaceColor','r')
    ylim([0 limFq])
end




% differences
% LR
figure; hold on;
for condition=1:4
    subplot(4,2,condition+(condition-1));hold on;
    plot(squeeze(gpDataDiff(condition).time),squeeze(gpDataDiff(condition).wave(channel,:)),'c','LineWidth',2)
    ylim([limERPminDiff limERPmaxDiff])
end
for condition=1:4
    subplot(4,2,condition+condition);hold on;
    bar(squeeze(gpDataDiff(condition).freq((gpData(condition).freq<15.5 & gpData(condition).freq>0))),squeeze(gpDataDiff(condition).amp(channel,2:31)),0.8,'FaceColor','c')
    ylim([0 limFqDiff])
end

% SR
figure; hold on;
for condition=1:4
    subplot(4,2,condition+(condition-1));hold on;
    plot(squeeze(gpDataDiff(condition).time),squeeze(gpDataDiff(condition+4).wave(channel,:)),'g','LineWidth',2)
    ylim([limERPminDiff limERPmaxDiff])
end
for condition=1:4
    subplot(4,2,condition+condition);hold on;
    bar(squeeze(gpDataDiff(condition).freq((gpData(condition).freq<15.5 & gpData(condition).freq>0))),squeeze(gpDataDiff(condition+4).amp(channel,2:31)),0.8,'FaceColor','g')
    ylim([0 limFqDiff])
end

