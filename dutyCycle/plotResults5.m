% sum of power computed on sbj averaged harmonics (using averaged cos and
% sin), which corresponds to the average data that will be run through the
% areas model
% Different from sum of power per sbj then average
clearvars

addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/commonFunctions
addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/biosemiUpdated
addpath /Volumes/Amrutam/Marlene/Git/fieldtrip-aleslab-fork
addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/svndlCopy
ft_defaults

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

[avData, proj_Amp] = averageAxxWithStd(dataSbj);

% normalisation term
% load the amplitudes for a square function
load('squareFFT.mat')
for cond=1:length(avData)
    sumSq(cond) = sum(sqAllFFT(avData(cond).harmIdx,cond).^2);
end
% compute sum of power for signal and noise
for cond=1:length(dataSbj)
    avData(cond).sumPower = sum(avData(cond).amp(:,avData(cond).harmIdx).^2,2);
    avData(cond).sumPowerNoise = sum(avData(cond).amp(:,[avData(cond).harmIdx-1 avData(cond).harmIdx+1]).^2,2) /2;
end
% bootstrap for getting confidence intervals
nbBoot = 1000;
tic
for cond=1:length(dataSbj)
    power = zeros(nbBoot,avData(1).nchan);powerNoise = zeros(nbBoot,avData(1).nchan);
    for bb=1:nbBoot
        if mod(bb,100)==0
            fprintf('condition%d boostrap%d \n' ,cond,bb)
        end
        % pick randomly 16 sbj with replacement
        pickSS = randi(length(keepSbj),1,length(keepSbj));
        tmpSin = zeros(16,128,577); tmpCos = zeros(16,128,577);
        for ss=1:length(keepSbj)
            tmpSin(ss,:,:) = dataSbj(cond,pickSS(ss)).sin;
            tmpCos(ss,:,:) = dataSbj(cond,pickSS(ss)).cos;
        end
        meanAmp = sqrt(squeeze(mean(tmpSin)).^2 + squeeze(mean(tmpCos)).^2);
        power(bb,:) = sum(meanAmp(:,avData(cond).harmIdx).^2,2);
        powerNoise(bb,:) = sum(meanAmp(:,[avData(cond).harmIdx-1 avData(cond).harmIdx+1]).^2,2) /2;
    end
    avData(cond).ciPower = prctile(power,[2.5 97.5])'; 
    avData(cond).ciPowerNoise = prctile(powerNoise,[2.5 97.5])';
    avData(cond).stdPower = std(power)'; 
    avData(cond).stdPowerNoise = std(powerNoise)';
    avData(cond).stdAmp = std(sqrt(power))';
    avData(cond).stdAmpNoise = std(sqrt(powerNoise))';
    avData(cond).ciAmp = prctile(sqrt(power),[2.5 97.5])'; 
    avData(cond).ciAmpNoise = prctile(sqrt(powerNoise),[2.5 97.5])';
    % note that std(power)/sumSq = std(power/sumSq) so normalisation can be
    % done later on the stdPower
end
toc

save('16sbjDC','avData','proj_Amp','cfg')



load('16sbjDC.mat')

col={'b','r','g'};
chan = 23; % Oz



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plot only 1f1 at Oz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean noise across frequencies and static/motion -> get one value per DC
% for dd=1:5
%     if dd==3
%         noise(:,dd) = mean([avData([3 8 13 18 22 21]).sumPowerNoise],2);
%         noiseNorm(:,dd) = mean([avData([3 8 13 18 22 21]).sumPowerNoise]./sqrt(sumSq([3 8 13 18 22 21])),2);
%         noiseStd(:,dd) = sqrt(sum([avData([3 8 13 18 22 21]).stdPowerNoise].^2,2)/6);
%         noiseStdAmp(:,dd) = sqrt(sum([avData([3 8 13 18 22 21]).stdAmpNoise].^2,2)/6);
%         noiseStdAmpNorm(:,dd) = sqrt(sum([avData([3 8 13 18 22 21]).stdAmpNoise]./sqrt(sumSq([3 8 13 18 22 21])).^2,2)/6);
%         tmpCI2 = [avData([3 8 13 18 22 21]).ciPowerNoise]./sqrt(sumSq([3 3 8 8 13 13 18 18 22 22 21 21])).^2;
%         noiseCIpowerLow(:,dd) = sqrt(sum(tmpCI2(:,1:2:end),2)/6);
%         noiseCIpowerHigh(:,dd) = sqrt(sum(tmpCI2(:,2:2:end),2)/6);
%     else
%         noise(:,dd) = mean([avData([0+dd 5+dd 10+dd 15+dd]).sumPowerNoise],2);
%         noiseNorm(:,dd) = mean([avData([0+dd 5+dd 10+dd 15+dd]).sumPowerNoise]./sqrt(sumSq([0+dd 5+dd 10+dd 15+dd])),2);
%         noiseStd(:,dd) = sqrt(sum([avData([0+dd 5+dd 10+dd 15+dd]).stdPowerNoise].^2,2)/4);
%         noiseStdAmp(:,dd) = sqrt(sum([avData([0+dd 5+dd 10+dd 15+dd]).stdAmpNoise].^2,2)/4);
%         noiseStdAmpNorm(:,dd) = sqrt(sum([avData([0+dd 5+dd 10+dd 15+dd]).stdAmpNoise]./sqrt(sumSq([0+dd 5+dd 10+dd 15+dd])).^2,2)/4);
%         tmpCI = [avData([0+dd 5+dd 10+dd 15+dd]).ciPowerNoise]./sqrt(sumSq([0+dd 0+dd 5+dd 5+dd 10+dd 10+dd 15+dd 15+dd])).^2;
%         noiseCIpowerLow(:,dd) = sqrt(sum(tmpCI(:,1:2:end),2)/4);
%         noiseCIpowerHigh(:,dd) = sqrt(sum(tmpCI(:,2:2:end),2)/4);
%     end
% end

% power
clear toPlot toPlotSTD toPlotCI
toPlot = [avData.sumPower]; toPlotSTD = [avData.stdPower];
toPlotCI = [avData.ciPower];  
toPlotNoise = [avData.sumPowerNoise];
toPlotCINoise = [avData.ciPowerNoise];  
toPlotSTDNoise = [avData.stdPowerNoise];
% toPlotCI = [avData.stdPower] * 1.96;
figure;hold on
for freq=1:3
    errorbar(dcVal,toPlot(chan,(freq-1)*5+1:freq*5),...
        toPlot(chan,(freq-1)*5+1:freq*5) - toPlotCI(chan,(freq-1)*10+1:2:freq*10),...
        toPlot(chan,(freq-1)*5+1:freq*5) - toPlotCI(chan,(freq-1)*10+2:2:freq*10),...
        ['.-' col{freq}],'MarkerSize',30,'Linewidth',2);
end
for freq=1:3
    errorbar(dcVal,toPlotNoise(chan,(freq-1)*5+1:freq*5),...
        toPlotNoise(chan,(freq-1)*5+1:freq*5) - toPlotCINoise(chan,(freq-1)*10+1:2:freq*10),...
        toPlotNoise(chan,(freq-1)*5+1:freq*5) - toPlotCINoise(chan,(freq-1)*10+2:2:freq*10),...
        ['.-' col{freq}],'MarkerSize',15,'Linewidth',1);
end
errorbar(50,toPlot(chan,18),toPlot(chan,18)-toPlotCI(chan,35),toPlot(chan,18)-toPlotCI(chan,36),['^:' col{1}],'MarkerSize',20,'Linewidth',2)
errorbar([25 50 75],toPlot(chan,[17 22 19]),toPlot(chan,[17 22 19])-toPlotCI(chan,[17*2-1 22*2-1 19*2-1]),toPlot(chan,[17 22 19])-toPlotCI(chan,[17*2 22*2 19*2]),['^:' col{2}],'MarkerSize',20,'Linewidth',2)
errorbar([12.5 50 87.5],toPlot(chan,[16 21 20]),toPlot(chan,[16 21 20])-toPlotCI(chan,[16*2-1 21*2-1 20*2-1]),toPlot(chan,[16 21 20])-toPlotCI(chan,[16*2 21*2 20*2]),['^:' col{3}],'MarkerSize',20,'Linewidth',2)
xlim([0 100]);
xticks([0:12.5:100]);
% ylim([0 4])
legend('10Hz','5Hz','2.5Hz','noise level','Location','Best')
xlabel('Duty Cycle')
ylabel('SSVEP power harmonics')
set(gca,'FontSize',15)
title('Oz bootstrapCI')
saveas(gcf,['figures' filesep 'powerDCallFqOz2.jpg'])

figure;hold on
for freq=1:3
    errorbar(dcVal,toPlot(chan,(freq-1)*5+1:freq*5),...
        toPlotSTD(chan,(freq-1)*5+1:freq*5),...
        ['.-' col{freq}],'MarkerSize',30,'Linewidth',2);
end
for freq=1:3
    errorbar(dcVal,toPlotNoise(chan,(freq-1)*5+1:freq*5),...
        toPlotSTDNoise(chan,(freq-1)*5+1:freq*5),...
        ['.-' col{freq}],'MarkerSize',15,'Linewidth',1);
end
errorbar(50,toPlot(chan,18),toPlotSTD(chan,18),['^:' col{1}],'MarkerSize',20,'Linewidth',2)
errorbar([25 50 75],toPlot(chan,[17 22 19]),toPlotSTD(chan,[17 22 19]),['^:' col{2}],'MarkerSize',20,'Linewidth',2)
errorbar([12.5 50 87.5],toPlot(chan,[16 21 20]),toPlotSTD(chan,[16 21 20]),['^:' col{3}],'MarkerSize',20,'Linewidth',2)
xlim([0 100]);
xticks([0:12.5:100]);
% ylim([0 4])
legend('10Hz','5Hz','2.5Hz','noise level','Location','Best')
xlabel('Duty Cycle')
ylabel('SSVEP power harmonics')
set(gca,'FontSize',15)
title('Oz bootstrapStd')
saveas(gcf,['figures' filesep 'powerDCallFqOz.jpg'])



% normalised amplitudes
clear toPlot toPlotSTD toPlotCI toPlotNoise toPlotSTDNoise toPlotCINoise
toPlot = sqrt([avData.sumPower] ./ sumSq); 
toPlotSTD = [avData.stdAmp] ./ sqrt(sumSq); 
toPlotNoise = sqrt([avData.sumPowerNoise] ./ sumSq); 
toPlotSTDNoise = [avData.stdAmpNoise] ./ sqrt(sumSq); 

count = 0;
for aa=1:2:size([avData.ciAmp],2)
    count = count+1;
    sumSqforCI(aa) = sumSq(count);
    sumSqforCI(aa+1) = sumSq(count);
end
toPlotCI = [avData.ciAmp] ./ sqrt(sumSqforCI);  
toPlotCINoise = [avData.ciAmpNoise] ./ sqrt(sumSqforCI); 

% CI
figure;hold on
for freq=1:3
    errorbar(dcVal,toPlot(chan,(freq-1)*5+1:freq*5),...
        toPlot(chan,(freq-1)*5+1:freq*5) - toPlotCI(chan,(freq-1)*10+1:2:freq*10),...
        toPlot(chan,(freq-1)*5+1:freq*5) - toPlotCI(chan,(freq-1)*10+2:2:freq*10),...
        ['.-' col{freq}],'MarkerSize',30,'Linewidth',2);
end
for freq=1:3
    errorbar(dcVal,toPlotNoise(chan,(freq-1)*5+1:freq*5),...
        toPlotNoise(chan,(freq-1)*5+1:freq*5) - toPlotCINoise(chan,(freq-1)*10+1:2:freq*10),...
        toPlotNoise(chan,(freq-1)*5+1:freq*5) - toPlotCINoise(chan,(freq-1)*10+2:2:freq*10),...
        ['.-' col{freq}],'MarkerSize',15,'Linewidth',1);
end
errorbar(50,toPlot(chan,18),toPlot(chan,18)-toPlotCI(chan,35),toPlot(chan,18)-toPlotCI(chan,36),['^:' col{1}],'MarkerSize',20,'Linewidth',2)
errorbar([25 50 75],toPlot(chan,[17 22 19]),toPlot(chan,[17 22 19])-toPlotCI(chan,[17*2-1 22*2-1 19*2-1]),toPlot(chan,[17 22 19])-toPlotCI(chan,[17*2 22*2 19*2]),['^:' col{2}],'MarkerSize',20,'Linewidth',2)
errorbar([12.5 50 87.5],toPlot(chan,[16 21 20]),toPlot(chan,[16 21 20])-toPlotCI(chan,[16*2-1 21*2-1 20*2-1]),toPlot(chan,[16 21 20])-toPlotCI(chan,[16*2 21*2 20*2]),['^:' col{3}],'MarkerSize',20,'Linewidth',2)
errorbar(50,toPlotNoise(chan,18),toPlotNoise(chan,18)-toPlotCINoise(chan,35),toPlotNoise(chan,18)-toPlotCINoise(chan,36),['^:' col{1}],'MarkerSize',10,'Linewidth',1)
errorbar([25 50 75],toPlotNoise(chan,[17 22 19]),toPlotNoise(chan,[17 22 19])-toPlotCINoise(chan,[17*2-1 22*2-1 19*2-1]),toPlotNoise(chan,[17 22 19])-toPlotCINoise(chan,[17*2 22*2 19*2]),['^:' col{2}],'MarkerSize',10,'Linewidth',1)
errorbar([12.5 50 87.5],toPlotNoise(chan,[16 21 20]),toPlotNoise(chan,[16 21 20])-toPlotCINoise(chan,[16*2-1 21*2-1 20*2-1]),toPlotNoise(chan,[16 21 20])-toPlotCINoise(chan,[16*2 21*2 20*2]),['^:' col{3}],'MarkerSize',10,'Linewidth',1)
xlim([0 100]);
xticks([0:12.5:100]);
% ylim([0 4])
legend('10Hz','5Hz','2.5Hz','noise level','Location','Best')
xlabel('Duty Cycle')
ylabel('norm SSVEP amplitude harmonics')
set(gca,'FontSize',15)
title('Oz bootstrapCI')
saveas(gcf,['figures' filesep 'ampDCallFqOzCI.jpg'])

% STD = SEM
figure;hold on
for freq=1:3
    errorbar(dcVal,toPlot(chan,(freq-1)*5+1:freq*5),...
        toPlotSTD(chan,(freq-1)*5+1:freq*5),...
        ['.-' col{freq}],'MarkerSize',30,'Linewidth',2);
end
errorbar(50,toPlot(chan,18),toPlotSTD(chan,18),['^:' col{1}],'MarkerSize',20,'Linewidth',2)
errorbar([25 50 75],toPlot(chan,[17 22 19]),toPlotSTD(chan,[17 22 19]),['^:' col{2}],'MarkerSize',20,'Linewidth',2)
errorbar([12.5 50 87.5],toPlot(chan,[16 21 20]),toPlotSTD(chan,[16 21 20]),['^:' col{3}],'MarkerSize',20,'Linewidth',2)
for freq=1:3
    errorbar(dcVal,toPlotNoise(chan,(freq-1)*5+1:freq*5),...
        toPlotSTDNoise(chan,(freq-1)*5+1:freq*5),...
        ['.-' col{freq}],'MarkerSize',15,'Linewidth',1);
end
errorbar(50,toPlotNoise(chan,18),toPlotSTDNoise(chan,18),['^:' col{1}],'MarkerSize',10,'Linewidth',1)
errorbar([25 50 75],toPlotNoise(chan,[17 22 19]),toPlotSTDNoise(chan,[17 22 19]),['^:' col{2}],'MarkerSize',10,'Linewidth',1)
errorbar([12.5 50 87.5],toPlotNoise(chan,[16 21 20]),toPlotSTDNoise(chan,[16 21 20]),['^:' col{3}],'MarkerSize',10,'Linewidth',1)
xlim([0 100]);
xticks([0:12.5:100]);
% ylim([0 3.5])
legend('10Hz','5Hz','2.5Hz','noise level','Location','Best')
xlabel('Duty Cycle')
ylabel('norm SSVEP amplitude harmonics')
set(gca,'FontSize',15)
title('Oz bootstrapStd')
saveas(gcf,['figures' filesep 'ampDCallFqOz.jpg'])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plot topographies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% power
% flicker
figure; hold on;
colormap('hot')
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo([avData(cond).sumPower],cfg.layout);
    colorbar
    if cond < 6
        caxis([0 1.5]);
    elseif cond > 5 && cond < 11
        caxis([0 6.5]);
    elseif cond > 10
        caxis([0 7.5]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoFlickPower.jpg'])
% motion
position = [11 7 3 9 15 13 8];
figure; hold on;
colormap('hot')
for cond=16:length(avData)
    subplot(3,5,position(cond-15)); hold on;
    plotTopo(avData(cond).sumPower,cfg.layout);
    colorbar
    if position(cond-15) < 6
        caxis([0 1.5]);
    elseif position(cond-15) > 5 && position(cond-15) < 11
        caxis([0 6.5]);
    elseif position(cond-15) > 10
        caxis([0 7.5]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoMotionPower.jpg'])


%%%% norm amp
% flicker
figure; hold on;
colormap('hot')
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(sqrt([avData(cond).sumPower]),cfg.layout);
    colorbar
    if cond < 6
        caxis([0 1]);
    elseif cond > 5 && cond < 11
        caxis([0 3]);
    elseif cond > 10
        caxis([0 3]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoFlickAmp.jpg'])
% motion
position = [11 7 3 9 15 13 8];
figure; hold on;
colormap('hot')
for cond=16:length(avData)
    subplot(3,5,position(cond-15)); hold on;
    plotTopo(sqrt([avData(cond).sumPower]),cfg.layout);
    colorbar
    if position(cond-15) < 6
        caxis([0 1]);
    elseif position(cond-15) > 5 && position(cond-15) < 11
        caxis([0 3]);
    elseif position(cond-15) > 10
        caxis([0 3]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoMotionAmp.jpg'])

%%%% norm amp
% flicker
figure; hold on;
colormap('hot')
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(sqrt([avData(cond).sumPower] ./ sumSq(cond)),cfg.layout);
    colorbar
    if cond < 6
        caxis([0 1]);
    elseif cond > 5 && cond < 11
        caxis([0 3]);
    elseif cond > 10
        caxis([0 3]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoFlickAmpNorm.jpg'])
% motion
position = [11 7 3 9 15 13 8];
figure; hold on;
colormap('hot')
for cond=16:length(avData)
    subplot(3,5,position(cond-15)); hold on;
    plotTopo(sqrt([avData(cond).sumPower] ./ sumSq(cond)),cfg.layout);
    colorbar
    if position(cond-15) < 6
        caxis([0 1]);
    elseif position(cond-15) > 5 && position(cond-15) < 11
        caxis([0 3]);
    elseif position(cond-15) > 10
        caxis([0 3]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoMotionAmpNorm.jpg'])











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% RATINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('fullRatings9.mat')
static = mean(tabStatic,3); moving = mean(tabMot,3);
for fq=1:3
    for dc=1:5
        statSEM(fq,dc) = std(tabStatic(fq,dc,:))/sqrt(6);
        movSEM(fq,dc) = std(tabMot(fq,dc,:))/sqrt(6);
    end
end

figure;hold on;
for fq=1:3
    errorbar(dcVal,static(fq,:),statSEM(1,:),['.-'  col{fq}],'MarkerSize',40,'LineWidth',2)
end
errorbar(dcVal(3),moving(7),movSEM(7),['^'  col{1}],'MarkerSize',15,'Linewidth',2)
errorbar([25 50 75],moving([5 8 11]),movSEM([5 8 11]),['^:'  col{2}],'MarkerSize',15,'Linewidth',2)
errorbar([12.5 50 87.5],moving([3 9 15]),movSEM([3 9 15]),['^:'  col{3}],'MarkerSize',15,'Linewidth',2)
xlim([0 100]);
xticks([0:12.5:100]);
ylim([0 3])
% xticklabels({'12.5','25','50','75','87.5'})
legend('10','5','2.5','Location','Best')
xlabel('Duty Cycle')
ylabel('motion ratings')
set(gca,'FontSize',15)
saveas(gcf,['figures' filesep 'ratingsDC.png'])
% saveas(gcf,'ratingsDC.eps','epsc')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% correl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% 1f1
figure;hold on;
scatter(avAmp(6:10),static(2,:),80,'filled','MarkerEdgeColor','none');
scatter(avAmp(11:15),static(3,:),80,'filled','MarkerEdgeColor','none');
R = corrcoef(avAmp(6:15),[static(2,:) static(3,:)]);
RsqS = R(1,2).^2;
R = corrcoef(avAmp([16:17 19:22]),moving([3 5 11 15 9 8]));
RsqM = R(1,2).^2;
lsline
title(['Oz S=' num2str(RsqS,2) ' M=' num2str(RsqM,2) ])
ylabel('motion rating')
xlabel('SSVEP amplitude 1f1 Oz')


%%% all harmonics
power = [avData.sumPower];
amp = sqrt([avData.sumPower]);
ampNorm = sqrt([avData.sumPower] ./ sumSq);

figure;hold on;
scatter(ampNorm(chan,1:5),static(1,:),80,'filled','MarkerEdgeColor','none');
scatter(ampNorm(chan,6:10),static(2,:),80,'filled','MarkerEdgeColor','none');
scatter(ampNorm(chan,11:15),static(3,:),80,'filled','MarkerEdgeColor','none');
ylabel('motion rating')
xlabel('SSVEP amplitude')
saveas(gcf,['figures' filesep 'correlAll.png'])

y = [static(2,:) static(3,:)];
yE = [statSEM(2,:) statSEM(3,:)];
x = ampNorm(chan,6:15);
stdAmp = [avData.stdAmp] ./ sqrt(sumSq);
% ciAmp = [avData.ciAmp] ./ sqrt(sumSqforCI);
xE = stdAmp(chan,6:15);
y2 = moving([3 5 11 15 9 8]); % do not include 10Hz
yE2 = movSEM([3 5 11 15 9 8]);
x2 = ampNorm(chan,[16:17 19:22]);% do not include 10Hz
xE2 = stdAmp(chan,[16:17 19:22]);

figure;hold on;
scatter(x, y,80,'filled','MarkerEdgeColor','none');
scatter(x2, y2,80,'filled','MarkerEdgeColor','none');
eb(1) = errorbar(x,y,xE, 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(x,y,yE, 'vertical', 'LineStyle', 'none');
set(eb, 'color', 'b', 'LineWidth', 1)
eb2(1) = errorbar(x2,y2,xE2, 'horizontal', 'LineStyle', 'none');
eb2(2) = errorbar(x2,y2,yE2, 'vertical', 'LineStyle', 'none');
set(eb2, 'color', 'r', 'LineWidth', 1)
lsline
ylim([0 3])
legend('flickering','moving')
ylabel('motion rating')
xlabel('norm SSVEP amplitude')
R = corrcoef(ampNorm(chan,6:15),[static(2,:) static(3,:)]);
RsqS = R(1,2).^2;
R = corrcoef(ampNorm(chan,[16:17 19:22]),moving([3 5 11 15 9 8]));
RsqM = R(1,2).^2;
title(['Oz S=' num2str(RsqS,2) ' M=' num2str(RsqM,2) ]) %only 2 digits
saveas(gcf,['figures' filesep 'correlSSVEPratingsAllFqOz.png'])


figure;hold on;
scatter(log(x), log(y),80,'filled','MarkerEdgeColor','none');
scatter(log(x2), log(y2),80,'filled','MarkerEdgeColor','none');
lsline
legend('flickering','moving')
ylabel('log(motion rating)')
xlabel('log(norm SSVEP amplitude)')
