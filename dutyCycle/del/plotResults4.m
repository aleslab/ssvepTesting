% sum harmonics per sbj then average then normalise (easier for error bars)
clearvars

addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/commonFunctions
addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/biosemiUpdated
addpath /Volumes/Amrutam/Marlene/Git/fieldtrip-aleslab-fork
addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/svndlCopy
ft_defaults
ft_defaultscfg.layout = 'biosemi128.lay';

% load individual data
dataDir = '/Volumes/Amrutam/Marlene/JUSTIN/DutyCycle/data/Axx/';
listData = dir([dataDir '*.mat']);
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};

keepSbj = [1:5 7:11 13:15 17:18 20]; 
dcVal = [12.5 25 50 75 87.5];

% load the amplitudes for a square function
load('squareFFT.mat')

for ss = 1:length(keepSbj)
    clear Axx;
    load([dataDir listData(keepSbj(ss)).name]);
    
    % sum the harmonics and normalise
    for cond = 1:length(Axx)
        sumSq = sum(sqAllFFT(Axx(cond).harmIdx,cond).^2);
        for ch=1:Axx(cond).nchan
            Axx(cond).powerSumHarm(ch) = sum(Axx(cond).amp(ch,Axx(cond).harmIdx).^2);
            % same as: sum(Axx(cond).cos(ch,Axx(cond).harmIdx).^2 + Axx(cond).sin(ch,Axx(cond).harmIdx).^2)
            Axx(cond).normFq(ch) = Axx(cond).powerSumHarm(ch) / sumSq;
            % for noise level as well
            Axx(cond).powerSumNoise(ch) = sum(Axx(cond).amp(ch,[Axx(cond).harmIdx-1 Axx(cond).harmIdx+1]).^2) /2;
            % = sum(Axx(cond).amp(ch,[Axx(cond).harmIdx-1 Axx(cond).harmIdx+1]).^2 /2)
        end
    end
    dataSbj(:,ss) = Axx;
end

% data dd = [2 4 8 9]; valNorm = 5;
% mean power: mean(dd)
% mean power norm: mean(dd)/valNorm = mean(dd/valNorm)
% std power norm: std(dd)/valNorm = std(dd/valNorm)
% mean amp: sqrt(mean(dd))'RMS' ~= mean(sqrt(dd))
% std amp: std(sqrt(dd))
% mean amp norm: mean(sqrt(dd))/sqrt(valNorm) = mean(sqrt(dd/valNorm))
% std amp: std(sqrt(dd))/sqrt(valNorm) = std(sqrt(dd/valNorm))
% use the same normalisation value for the noise level

[avData, proj_Amp] = averageAxxWithStd(dataSbj);

for cond=1:length(avData)
    % all harmonics
    sumSq = sum(sqAllFFT(avData(cond).harmIdx,cond).^2);
    clear dataTmpPower dataTmpPowerNoise
    for ss=1:size(dataSbj,2)
        dataTmpPower(:,ss) = [dataSbj(cond,ss).powerSumHarm];
        dataTmpPowerNoise(:,ss) = [dataSbj(cond,ss).powerSumNoise];
    end
    avData(cond).powerSumHarm = mean(dataTmpPower,2);
    avData(cond).powerNorm = mean(dataTmpPower,2) / sumSq;
    avData(cond).powerNormStd = std(dataTmpPower,[],2) / sumSq;
    avData(cond).ampSumHarm = mean(sqrt(dataTmpPower),2);
    avData(cond).ampSumHarmStd = std(sqrt(dataTmpPower),[],2);
    avData(cond).ampNorm = mean(sqrt(dataTmpPower/ sumSq),2);
    avData(cond).ampNormStd = std(sqrt(dataTmpPower/ sumSq),[],2);
    % noise level
    avData(cond).noise = mean(sqrt(dataTmpPowerNoise),2);
    avData(cond).noiseStd = std(sqrt(dataTmpPowerNoise),[],2); 
    avData(cond).noiseNorm = mean(sqrt(dataTmpPowerNoise/ sumSq),2);
    avData(cond).noiseNormStd = std(sqrt(dataTmpPowerNoise/ sumSq),[],2);    
end

   
save('16sbjDC','avData','proj_Amp','cfg')

col={'b','r','g'};
chan = 23; % Oz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plot only 1f1 at Oz
for cond=1:length(avData)
    avAmp(cond)=avData(cond).amp(chan,(avData(cond).harmIdx(1)));
    semProj(cond) = std(proj_Amp(chan,avData(cond).harmIdx(1),cond,:)) / sqrt(length(keepSbj));
    % compute noise level
    baseline(cond) = (avData(cond).amp(chan,avData(cond).harmIdx(1)-1)+avData(cond).amp(chan,avData(cond).harmIdx(1)+1)) / 2;
end

figure;hold on;
for freq=1:3
    errorbar(dcVal,avAmp((freq-1)*5+1:freq*5),semProj((freq-1)*5+1:freq*5),['.-' col{freq}],'MarkerSize',30,'Linewidth',2);
end
line([dcVal(1) dcVal(5)], [min(baseline) min(baseline)],'Linewidth',2,'Color','k','LineStyle','--')
line([dcVal(1) dcVal(5)], [mean(baseline) mean(baseline)],'Linewidth',2,'Color','k','LineStyle','-')
line([dcVal(1) dcVal(5)], [max(baseline) max(baseline)],'Linewidth',2,'Color','k','LineStyle','--')
errorbar(50,avAmp(18),semProj(18),['^:' col{1}],'MarkerSize',20,'Linewidth',2)
errorbar([25 50 75],avAmp([17 22 19]),semProj([17 22 19]),['^:' col{2}],'MarkerSize',20,'Linewidth',2)
errorbar([12.5 50 87.5],avAmp([16 21 20]),semProj([16 21 20]),['^:' col{3}],'MarkerSize',20,'Linewidth',2)
xlim([0 100]);
xticks([0:12.5:100]);
ylim([0 3])
legend('10Hz','5Hz','2.5Hz','noise level','Location','Best')
xlabel('Duty Cycle')
ylabel('1f1 SSVEP amplitude')
set(gca,'FontSize',15)
title('Oz')
saveas(gcf,['figures' filesep 'amp1f1.jpg'])





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plot all harmonics at Oz

% powerNorm = [avData.powerNorm];
% powerNormSem = [avData.powerNormStd]/sqrt(length(keepSbj));
ampNorm = [avData.ampNorm];
ampNormSem = [avData.ampNormStd]/sqrt(length(keepSbj));
% mean noise across frequencies and static/motion -> get one value per DC
% noiseNorm = [avData.noiseNorm];
% noiseSem = [avData.ampNoiseStd];
for dd=1:5
    if dd==3
        noiseNorm(:,dd) = mean([avData([3 8 13 18 22 21]).noiseNorm],2);
        noiseNormSem(:,dd) = sqrt(sum([avData([3 8 13 18 22 21]).noiseNormStd].^2,2)/7);
        noise(:,dd) = mean([avData([3 8 13 18 22 21]).noise],2);
        noiseSem(:,dd) = sqrt(sum([avData([3 8 13 18 22 21]).noiseStd].^2,2)/7);
    else
        noiseNorm(:,dd) = mean([avData([0+dd 5+dd 10+dd 15+dd]).noiseNorm],2);
        noiseNormSem(:,dd) = sqrt(sum([avData([0+dd 5+dd 10+dd 15+dd]).noiseNormStd].^2,2)/4);
        noise(:,dd) = mean([avData([0+dd 5+dd 10+dd 15+dd]).noise],2);
        noiseSem(:,dd) = sqrt(sum([avData([0+dd 5+dd 10+dd 15+dd]).noiseStd].^2,2)/4);
    end
end

%%%%%%%%%% normalised amplitudes
figure;hold on;
for freq=1:3
    errorbar(dcVal,ampNorm(chan,(freq-1)*5+1:freq*5),ampNormSem(chan,(freq-1)*5+1:freq*5),['.-' col{freq}],'MarkerSize',30,'Linewidth',2);
end
errorbar(dcVal,noiseNorm(chan,:),noiseNormSem(chan,:),['--' 'k'],'Linewidth',2);
errorbar(50,ampNorm(chan,18),ampNormSem(chan,18),['^:' col{1}],'MarkerSize',20,'Linewidth',2)
errorbar([25 50 75],ampNorm(chan,[17 22 19]),ampNormSem(chan,[17 22 19]),['^:' col{2}],'MarkerSize',20,'Linewidth',2)
errorbar([12.5 50 87.5],ampNorm(chan,[16 21 20]),ampNormSem(chan,[16 21 20]),['^:' col{3}],'MarkerSize',20,'Linewidth',2)
xlim([0 100]);
xticks([0:12.5:100]);
ylim([0 4])
legend('10Hz','5Hz','2.5Hz','noise level','Location','Best')
xlabel('Duty Cycle')
ylabel('norm SSVEP amplitude')
set(gca,'FontSize',15)
title('Oz')
saveas(gcf,['figures' filesep 'ampDCNormAllFqSbj.jpg'])

%%%%%%%%%% NOT norm amplitudes
ampNorm = [avData.ampSumHarm];
ampNormSem = [avData.ampSumHarmStd]/sqrt(length(keepSbj));
figure;hold on;
for freq=1:3
    errorbar(dcVal,ampNorm(chan,(freq-1)*5+1:freq*5),ampNormSem(chan,(freq-1)*5+1:freq*5),['.-' col{freq}],'MarkerSize',30,'Linewidth',2);
end
errorbar(dcVal,noise(chan,:),noiseSem(chan,:),['--' 'k'],'Linewidth',2);
errorbar(50,ampNorm(chan,18),ampNormSem(chan,18),['^:' col{1}],'MarkerSize',20,'Linewidth',2)
errorbar([25 50 75],ampNorm(chan,[17 22 19]),ampNormSem(chan,[17 22 19]),['^:' col{2}],'MarkerSize',20,'Linewidth',2)
errorbar([12.5 50 87.5],ampNorm(chan,[16 21 20]),ampNormSem(chan,[16 21 20]),['^:' col{3}],'MarkerSize',20,'Linewidth',2)
xlim([0 100]);
xticks([0:12.5:100]);
ylim([0 4])
legend('10Hz','5Hz','2.5Hz','noise level','Location','Best')
xlabel('Duty Cycle')
ylabel('sumHarm SSVEP amplitude')
set(gca,'FontSize',15)
title('Oz')
saveas(gcf,['figures' filesep 'ampDCSumHarm.jpg'])


% % instead of noise level per DC just use min max and mean?
% baseline = noiseNorm(chan,:);
% chan = 23;
% figure;hold on;
% for freq=1:3
%     errorbar(dcVal,ampNorm(chan,(freq-1)*5+1:freq*5),ampNormSem(chan,(freq-1)*5+1:freq*5),['.-' col{freq}],'MarkerSize',30,'Linewidth',2);
% end
% line([dcVal(1) dcVal(5)], [mean(baseline) mean(baseline)],'Linewidth',2,'Color','k','LineStyle','-')
% patch([dcVal(1) dcVal(5) dcVal(5) dcVal(1)],[min(baseline) min(baseline) max(baseline) max(baseline)],'k','facealpha',0.1,'EdgeColor','none')
% errorbar(50,ampNorm(chan,18),ampNormSem(chan,18),['^:' col{1}],'MarkerSize',20,'Linewidth',2)
% errorbar([25 50 75],ampNorm(chan,[17 22 19]),ampNormSem(chan,[17 22 19]),['^:' col{2}],'MarkerSize',20,'Linewidth',2)
% errorbar([12.5 50 87.5],ampNorm(chan,[16 21 20]),ampNormSem(chan,[16 21 20]),['^:' col{3}],'MarkerSize',20,'Linewidth',2)
% xlim([0 100]);
% xticks([0:12.5:100]);
% ylim([0 4])
% legend('10Hz','5Hz','2.5Hz','noise level','Location','Best')
% xlabel('Duty Cycle')
% ylabel('norm SSVEP amplitude')
% set(gca,'FontSize',15)
% title('Oz')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plot all harmonics at Oz but using the average data for summing the
%%%%% harmonics
for cond=1:length(avData)
    avAmp(:,cond) = sqrt(sum(avData(cond).amp(:,avData(cond).harmIdx).^2,2));
    sumSq = sum(sqAllFFT(avData(cond).harmIdx,cond).^2);
    avAmpNorm(:,cond) = avAmp(:,cond) / sqrt(sumSq);
end
figure;hold on;
for freq=1:3
    plot(dcVal,avAmp(chan,(freq-1)*5+1:freq*5),['.-' col{freq}],'MarkerSize',30,'Linewidth',2);
end
plot(50,avAmp(chan,18),['^:' col{1}],'MarkerSize',20,'Linewidth',2)
plot([25 50 75],avAmp(chan,[17 22 19]),['^:' col{2}],'MarkerSize',20,'Linewidth',2)
plot([12.5 50 87.5],avAmp(chan,[16 21 20]),['^:' col{3}],'MarkerSize',20,'Linewidth',2)
xlim([0 100]);
xticks([0:12.5:100]);
% ylim([0 4])
legend('10Hz','5Hz','2.5Hz','Location','Best')
xlabel('Duty Cycle')
ylabel('sumHarm SSVEP amplitude')
set(gca,'FontSize',15)
title('Oz')
saveas(gcf,['figures' filesep 'ampDCSumHarmAv.jpg'])

figure;hold on;
for freq=1:3
    plot(dcVal,avAmpNorm(chan,(freq-1)*5+1:freq*5),['.-' col{freq}],'MarkerSize',30,'Linewidth',2);
end
plot(50,avAmpNorm(chan,18),['^:' col{1}],'MarkerSize',20,'Linewidth',2)
plot([25 50 75],avAmpNorm(chan,[17 22 19]),['^:' col{2}],'MarkerSize',20,'Linewidth',2)
plot([12.5 50 87.5],avAmpNorm(chan,[16 21 20]),['^:' col{3}],'MarkerSize',20,'Linewidth',2)
xlim([0 100]);
xticks([0:12.5:100]);
% ylim([0 4])
legend('10Hz','5Hz','2.5Hz','Location','Best')
xlabel('Duty Cycle')
ylabel('norm SSVEP amplitude')
set(gca,'FontSize',15)
title('Oz')
saveas(gcf,['figures' filesep 'ampDCNormAv.jpg'])






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% topographies
ampNorm = [avData.ampNorm];
ampNormSem = [avData.ampNormStd]/sqrt(length(keepSbj));
ampNorm = [avData.ampSumHarm];
ampNormSem = [avData.ampSumHarmStd]/sqrt(length(keepSbj));

% flicker
figure; hold on;
colormap('hot')
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(ampNorm(:,cond),cfg.layout);
    colorbar
    if cond < 6
        caxis([0 2]);
    elseif cond > 5 && cond < 11
        caxis([0 3.5]);
    elseif cond > 10
        caxis([0 4]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
% saveas(gcf,['figures' filesep 'topoFlickNormAll.jpg'])
saveas(gcf,['figures' filesep 'topoFlickSumHarm.jpg'])

% motion
position = [11 7 3 9 15 13 8];
figure; hold on;
colormap('hot')
for cond=16:length(avData)
    %     max(avData(cond).amp(:,avData(cond).i1f1*2-1))
    subplot(3,5,position(cond-15)); hold on;
    plotTopo(ampNorm(:,cond),cfg.layout);
    colorbar
    if position(cond-15) < 6
        caxis([0 2]);
    elseif position(cond-15) > 5 && position(cond-15) < 11
        caxis([0 3.5]);
    elseif position(cond-15) > 10
        caxis([0 4]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
% saveas(gcf,['figures' filesep 'topoMotionNormAll.jpg'])
saveas(gcf,['figures' filesep 'topoMotionSumHarm.jpg'])



% flicker
figure; hold on;
colormap('hot')
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(avAmp(:,cond),cfg.layout);
    colorbar
    if cond < 6
        caxis([0 2]);
    elseif cond > 5 && cond < 11
        caxis([0 3.5]);
    elseif cond > 10
        caxis([0 4]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
% saveas(gcf,['figures' filesep 'topoFlickNormAll.jpg'])
saveas(gcf,['figures' filesep 'topoFlickSumHarmAv.jpg'])

% motion
position = [11 7 3 9 15 13 8];
figure; hold on;
colormap('hot')
for cond=16:length(avData)
    %     max(avData(cond).amp(:,avData(cond).i1f1*2-1))
    subplot(3,5,position(cond-15)); hold on;
    plotTopo(avAmp(:,cond),cfg.layout);
    colorbar
    if position(cond-15) < 6
        caxis([0 2]);
    elseif position(cond-15) > 5 && position(cond-15) < 11
        caxis([0 3.5]);
    elseif position(cond-15) > 10
        caxis([0 4]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
% saveas(gcf,['figures' filesep 'topoMotionNormAll.jpg'])
saveas(gcf,['figures' filesep 'topoMotionSumHarmAv.jpg'])
