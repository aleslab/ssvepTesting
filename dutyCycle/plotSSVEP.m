%%% Plot results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plot only 1f1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load dataDCavRef.mat
chan = 23; % Oz
dcVal = [12.5 25 50 75 87.5];
col={'b','r','g'};

%%%%%%%% Oz
% get values to plot
avAmp = zeros(1,length(avData));semProj = avAmp; baseline = avAmp; semBase = avAmp; 
for cond=1:length(avData)
    avAmp(cond)=avData(cond).amp(chan,(avData(cond).harmIdx(1)));
    semProj(cond) = std(proj_Amp(chan,avData(cond).harmIdx(1),cond,:)) / sqrt(size(proj_Amp,4));
    % compute noise level
    baseline(cond) = (avData(cond).amp(chan,avData(cond).harmIdx(1)-1)+avData(cond).amp(chan,avData(cond).harmIdx(1)+1)) / 2;
    % get an average of the noise amplitude between F-1 F+1 then std & SEM (/sqrtN)
    semBase(cond) = std((proj_Amp(chan,avData(cond).harmIdx(1)-1,cond,:) + proj_Amp(chan,avData(cond).harmIdx(1)+1,cond,:) ) / 2) / sqrt(size(proj_Amp,4));
end
% pool noise?
for cond=1:length(avData)
    allNoise(cond,:) = proj_Amp(chan,avData(cond).harmIdx(1)-1,cond,:) + proj_Amp(chan,avData(cond).harmIdx(1)+1,cond,:) ;
end
stdNoise = std(mean(allNoise))/sqrt(size(allNoise,2));
meanNoise = mean(mean(allNoise));

figure;hold on;
for freq=1:3
    errorbar(dcVal,avAmp((freq-1)*5+1:freq*5),semProj((freq-1)*5+1:freq*5),['.-' col{freq}],'MarkerSize',30,'Linewidth',2);
end
line([dcVal(1) dcVal(5)], [meanNoise-stdNoise meanNoise-stdNoise],'Linewidth',2,'Color','k','LineStyle','--')
line([dcVal(1) dcVal(5)], [meanNoise meanNoise],'Linewidth',2,'Color','k','LineStyle','-')
line([dcVal(1) dcVal(5)], [meanNoise+stdNoise meanNoise+stdNoise],'Linewidth',2,'Color','k','LineStyle','--')
errorbar(50,avAmp(18),semProj(18),['^:' col{1}],'MarkerSize',20,'Linewidth',2)
errorbar([25 50 75],avAmp([17 22 19]),semProj([17 22 19]),['^:' col{2}],'MarkerSize',20,'Linewidth',2)
errorbar([12.5 50 87.5],avAmp([16 21 20]),semProj([16 21 20]),['^:' col{3}],'MarkerSize',20,'Linewidth',2)
xlim([0 100]);
xticks([0:12.5:100]);
ylim([0 2])
legend('10Hz','5Hz','2.5Hz','noise level','Location','Best')
xlabel('Duty Cycle')
ylabel('1f1 SSVEP amplitude')
set(gca,'FontSize',15)
title('Oz')
saveas(gcf,['figures' filesep 'amp1f1.png'])


%%%%%%%% Topographies
% flicker
figure; hold on;
colormap('hot')
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(avData(cond).amp(:,(avData(cond).harmIdx(1))),cfg.layout);
    colorbar
    if cond < 6
        caxis([0 1]);
    else
        caxis([0 2]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoFlickAvgRef.png'])
% motion
position = [11 7 3 9 15 13 8];
figure; hold on;
colormap('hot')
for cond=16:length(avData)
    subplot(3,5,position(cond-15)); hold on;
    plotTopo(avData(cond).amp(:,(avData(cond).harmIdx(1))),cfg.layout);
    colorbar
    if position(cond-15) < 6
        caxis([0 1]);
    elseif position(cond-15) > 5 && position(cond-15) < 11
        caxis([0 2]);
    elseif position(cond-15) > 10
        caxis([0 2]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoMotionAvgRef.png'])


%%%%%%%% Amplitude per Area

addpath /Users/marleneponcet/Documents/Git/eegSourceTemplateMatching
load('biosemi128.mat');
% pool cond sin and cos reg parameter SAME for all conditions
poolData = [];
for cond = 1:length(avData)
    poolData = [poolData avData(cond).cos(:,avData(cond).harmIdx(1)) avData(cond).sin(:,avData(cond).harmIdx(1))];
end
% match with templates
[beta, betaUnscaled] = fitEEGTemplates(poolData,biosemi128);
% compute amplitude (power=cos^2+sin^2)
betaAmp = sqrt(betaUnscaled(:,1:2:end).^2 + betaUnscaled(:,2:2:end).^2); % odd betas = cos, even betas = sin

col={'b','r','g'};
dcVal = [12.5 25 50 75 87.5];
figure;set(gcf, 'Position', [0 0 1600 700])
for roi = 1:18
    subplot(3,6,roi); hold on;
    plot(dcVal,betaAmp(roi,1:5),['.-'  col{1}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,betaAmp(roi,6:10),['.-'  col{2}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,betaAmp(roi,11:15),['.-'  col{3}],'MarkerSize',40,'LineWidth',2)
    plot(50,betaAmp(roi,18),['^:' col{1}],'MarkerSize',15,'Linewidth',2)
    plot([25 50 75],betaAmp(roi,[17 22 19]),['^:' col{2}],'MarkerSize',15,'Linewidth',2)
    plot([12.5 50 87.5],betaAmp(roi,[16 21 20]),['^:' col{3}],'MarkerSize',15,'Linewidth',2)    
    xticks([0:12.5:100]);
    title(biosemi128.listROIs(roi))
    ylim([0 max(betaAmp(:))]);
    ylabel('amplitude (a.u.)')
end
xlabel('Duty Cycle')
legend('10','5','2.5','10mov','5mov','2.5mov','Location','Best')
saveas(gcf,'figures/beta1f1','png')

% pie chart?
figure;
for cond=1:15
subplot(3,5,cond); 
pie(betaAmp(:,cond)/sum(betaAmp(:,cond)),biosemi128.listROIs)
end
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,'figures/beta1f1pie','png')
% % Create legend
% lgd = legend(biosemi128.listROIs);
% lgd.Layout.Tile = 'east';

% motion
position = [11 7 3 9 15 13 8];
figure;
for cond=16:22
subplot(3,5,position(cond-15)); 
pie(betaAmp(:,cond)/sum(betaAmp(:,cond)))
end
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,'figures/beta1f1motionpie','png')

% check that I retrieve the right topographies
retrieve = biosemi128.weights * betaUnscaled;
retrieveAmp = sqrt(retrieve(:,1:2:end).^2 + retrieve(:,2:2:end).^2); % odd betas = cos, even betas = sin
figure; hold on;
colormap('hot')
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(retrieveAmp(:,cond),cfg.layout)
    colorbar
    if cond < 6
        caxis([0 1]);
    else
        caxis([0 2]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoFlickAvgRefretrieved.png'])
% motion
position = [11 7 3 9 15 13 8];
figure; hold on;
colormap('hot')
for cond=16:length(avData)
    subplot(3,5,position(cond-15)); hold on;
    plotTopo(retrieveAmp(:,cond),cfg.layout);
    colorbar
    if position(cond-15) < 6
        caxis([0 1]);
    else 
        caxis([0 2]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoMotionAvgRefretrieved.png'])





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% All harmonics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need to normalise amplitudes
% load the amplitudes for a square function
load('squareFFTnew.mat')
% compute sum of power
for cond=1:length(avData)
    sumSq(cond) = sum(sqAnalytic(avData(cond).harmIdx,cond).^2);
end
% normalise amplitudes
toPlot = sqrt([avData.sumPower] ./ sumSq); 
toPlotNoise = sqrt([avData.sumPowerNoise] ./ sumSq); 
% toPlotSTD = [avData.stdAmp] ./ sqrt(sumSq); 
% toPlotSTDNoise = [avData.stdAmpNoise] ./ sqrt(sumSq); 
% need copied sumSq for the CI
sumSqforCI = repmat(sumSq,2,1); sumSqforCI = reshape(sumSqforCI,[1 44]);
toPlotCI = [avData.ciAmp] ./ sqrt(sumSqforCI);  
toPlotCINoise = [avData.ciAmpNoise] ./ sqrt(sumSqforCI); 

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
ylabel('norm SSVEP amplitude (sum of harmonics)')
set(gca,'FontSize',15)
title('Oz bootstrapCI')
saveas(gcf,['figures' filesep 'allFqNormAmpOz.png'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Correl RATINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('ratings10.mat')
for fq=1:3
    for dc=1:5
        statSEM(fq,dc) = std(tabStatic(fq,dc,:))/sqrt(length(tabStatic));
        movSEM(fq,dc) = std(tabMot(fq,dc,:))/sqrt(length(tabMot));
    end
end

ampNorm = toPlot;

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


figure;
subplot(1,3,1);hold on;
scatter(x, y,80,'filled','MarkerEdgeColor','none');
scatter(x2, y2,80,'filled','MarkerEdgeColor','none');
lsline
legend('flickering','moving')
ylabel('(motion rating)')
xlabel('(norm SSVEP amplitude)')
subplot(1,3,2);hold on;
scatter(log(x), y,80,'filled','MarkerEdgeColor','none');
scatter(log(x2), y2,80,'filled','MarkerEdgeColor','none');
lsline
legend('flickering','moving')
ylabel('log(motion rating)')
xlabel('norm SSVEP amplitude')
subplot(1,3,3);hold on;
scatter(log(x), y,80,'filled','MarkerEdgeColor','none');
scatter(log(x2), y2,80,'filled','MarkerEdgeColor','none');
lsline
legend('flickering','moving')
ylabel('log(motion rating)')
xlabel('log(norm SSVEP amplitude)')
saveas(gcf,['figures' filesep 'correlTestLog.png'])


%%%% norm amp
% flicker
figure; hold on;
colormap('hot')
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(toPlot(:,cond),cfg.layout);
    colorbar
    if cond < 6
        caxis([0 1]);
    elseif cond > 5 && cond < 11
        caxis([0 2]);
    elseif cond > 10
        caxis([0 2.5]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoFlickSumAmpNormAvgRef.png'])
% motion
position = [11 7 3 9 15 13 8];
figure; hold on;
colormap('hot')
for cond=16:length(avData)
    subplot(3,5,position(cond-15)); hold on;
    plotTopo(toPlot(:,cond),cfg.layout);
    colorbar
    if position(cond-15) < 6
        caxis([0 1]);
    elseif position(cond-15) > 5 && position(cond-15) < 11
        caxis([0 2]);
    elseif position(cond-15) > 10
        caxis([0 2.5]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoMotionSumAmpNormAvgRef.png'])


%%%%%%%% Amplitude per Area
addpath /Users/marleneponcet/Documents/Git/eegSourceTemplateMatching
load('biosemi128.mat');
% do on normalised sum of amplitudes NOOOOOOO! Need cos and sin!!
% Amplitudes (or power) are not negative!!
[beta, betaUnscaled,lambda] = fitEEGTemplates(toPlot,biosemi128,1);
col={'b','r','g'};
dcVal = [12.5 25 50 75 87.5];
figure;set(gcf, 'Position', [0 0 1600 700])
for roi = 1:18
    subplot(3,6,roi); hold on;
    plot(dcVal,betaUnscaled(roi,1:5),['.-'  col{1}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,betaUnscaled(roi,6:10),['.-'  col{2}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,betaUnscaled(roi,11:15),['.-'  col{3}],'MarkerSize',40,'LineWidth',2)
    plot(50,betaUnscaled(roi,18),['^:' col{1}],'MarkerSize',15,'Linewidth',2)
    plot([25 50 75],betaUnscaled(roi,[17 22 19]),['^:' col{2}],'MarkerSize',15,'Linewidth',2)
    plot([12.5 50 87.5],betaUnscaled(roi,[16 21 20]),['^:' col{3}],'MarkerSize',15,'Linewidth',2)    
%     plot(dcVal,ampNormNoise(roi,1:5),['.-'  col{1}],'MarkerSize',20,'LineWidth',1)
%     plot(dcVal,ampNormNoise(roi,6:10),['.-'  col{2}],'MarkerSize',20,'LineWidth',1)
%     plot(dcVal,ampNormNoise(roi,11:15),['.-'  col{3}],'MarkerSize',20,'LineWidth',1)
%     plot(50,ampNormNoise(roi,18),['^:' col{1}],'MarkerSize',10,'Linewidth',1)
%     plot([25 50 75],ampNormNoise(roi,[17 22 19]),['^:' col{2}],'MarkerSize',10,'Linewidth',1)
%     plot([12.5 50 87.5],ampNormNoise(roi,[16 21 20]),['^:' col{3}],'MarkerSize',10,'Linewidth',1)    
    xticks([0:12.5:100]);
    title(biosemi128.listROIs(roi))
%     ylim([min(betaUnscaled(:)) max(betaUnscaled(:))]);
    ylabel('amplitude (a.u.)')
end
xlabel('Duty Cycle')
legend('10','5','2.5','10mov','5mov','2.5mov','Location','Best')

% check that I retrieve the right topographies
retrieveAmp = biosemi128.weights * betaUnscaled;
figure; hold on;
colormap('hot')
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(retrieveAmp(:,cond),cfg.layout)
    colorbar
%     if cond < 6
%         caxis([0 1]);
%     elseif cond > 5 && cond < 11
%         caxis([0 2]);
%     elseif cond > 10
%         caxis([0 2.5]);
%     end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoFlickSumFqretrieved.png'])



% pool cond sin and cos reg parameter SAME for all conditions
poolData = [];poolNoise=[];
for cond = 1:length(avData)
    poolData = [poolData avData(cond).cos(:,avData(cond).harmIdx) avData(cond).sin(:,avData(cond).harmIdx)];
    poolNoise = [poolNoise (avData(cond).cos(:,avData(cond).harmIdx-1)+ avData(cond).cos(:,avData(cond).harmIdx+1)/2) (avData(cond).sin(:,avData(cond).harmIdx-1)+ avData(cond).sin(:,avData(cond).harmIdx+1)/2)];
end
% match with templates
[beta, betaUnscaled,lambda] = fitEEGTemplates(poolData,biosemi128,1);
% [betaNoise, betaUnscaledNoise,lambda] = fitEEGTemplates(poolNoise,biosemi128,0,1000);
% [beta, betaUnscaled,lambda] = fitEEGTemplates(poolData,biosemi128,0,1000);

% cos and sin are stacked across conditions
% need to get their indexes for computing power per condition
nbHarm = arrayfun(@(x) numel([avData(x).harmIdx]),1:22);
% get the range: cumulative sum of elements
rangeHarm = [0 cumsum(nbHarm*2)]; % *2 to account for sin and cos

% compute amplitude (power=cos^2+sin^2)
for cond = 1:22
    power(:,cond) = sum(betaUnscaled(:,rangeHarm(cond)+1:2:rangeHarm(cond+1)).^2 ...
        + betaUnscaled(:,rangeHarm(cond)+2:2:rangeHarm(cond+1)).^2,2);% odd betas = cos, even betas = sin    
    % normalised amplitudes
    ampNorm(:,cond) = sqrt( power(:,cond) / sumSq(cond));
%     % same for noise
%     powerNoise(:,cond) = sum(betaUnscaledNoise(:,rangeHarm(cond)+1:2:rangeHarm(cond+1)).^2 ...
%         + betaUnscaledNoise(:,rangeHarm(cond)+2:2:rangeHarm(cond+1)).^2,2);% odd betas = cos, even betas = sin    
%     % normalised amplitudes
%     ampNormNoise(:,cond) = sqrt( powerNoise(:,cond) / sumSq(cond));
end

col={'b','r','g'};
dcVal = [12.5 25 50 75 87.5];
figure;set(gcf, 'Position', [0 0 1600 700])
for roi = 1:18
    subplot(3,6,roi); hold on;
    plot(dcVal,ampNorm(roi,1:5),['.-'  col{1}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,ampNorm(roi,6:10),['.-'  col{2}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,ampNorm(roi,11:15),['.-'  col{3}],'MarkerSize',40,'LineWidth',2)
    plot(50,ampNorm(roi,18),['^:' col{1}],'MarkerSize',15,'Linewidth',2)
    plot([25 50 75],ampNorm(roi,[17 22 19]),['^:' col{2}],'MarkerSize',15,'Linewidth',2)
    plot([12.5 50 87.5],ampNorm(roi,[16 21 20]),['^:' col{3}],'MarkerSize',15,'Linewidth',2)    
%     plot(dcVal,ampNormNoise(roi,1:5),['.-'  col{1}],'MarkerSize',20,'LineWidth',1)
%     plot(dcVal,ampNormNoise(roi,6:10),['.-'  col{2}],'MarkerSize',20,'LineWidth',1)
%     plot(dcVal,ampNormNoise(roi,11:15),['.-'  col{3}],'MarkerSize',20,'LineWidth',1)
%     plot(50,ampNormNoise(roi,18),['^:' col{1}],'MarkerSize',10,'Linewidth',1)
%     plot([25 50 75],ampNormNoise(roi,[17 22 19]),['^:' col{2}],'MarkerSize',10,'Linewidth',1)
%     plot([12.5 50 87.5],ampNormNoise(roi,[16 21 20]),['^:' col{3}],'MarkerSize',10,'Linewidth',1)    
    xticks([0:12.5:100]);
    title(biosemi128.listROIs(roi))
    ylim([0 max(ampNorm(:))]);
    ylabel('amplitude (a.u.)')
end
xlabel('Duty Cycle')
legend('10','5','2.5','10mov','5mov','2.5mov','Location','Best')
saveas(gcf,'figures/betaSumFqNorm','png')


% pie chart?
figure;
for cond=1:15
subplot(3,5,cond); 
pie(ampNorm(:,cond)/sum(ampNorm(:,cond)),biosemi128.listROIs) % normalise so the total is 1
end
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,'figures/betaSumFqpie','png')
% % Create legend
% lgd = legend(biosemi128.listROIs);
% lgd.Layout.Tile = 'east';

% motion
position = [11 7 3 9 15 13 8];
figure;
for cond=16:22
subplot(3,5,position(cond-15)); 
pie(ampNorm(:,cond)/sum(ampNorm(:,cond)),biosemi128.listROIs)
end
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,'figures/betaSumFqmotionpie','png')

% check that I retrieve the right topographies
retrieve = biosemi128.weights * betaUnscaled;
for cond=1:22
retrieveAmp(:,cond) = sqrt(sum(retrieve(:,rangeHarm(cond)+1:2:rangeHarm(cond+1)).^2 ...
        + retrieve(:,rangeHarm(cond)+2:2:rangeHarm(cond+1)).^2,2));% odd betas = cos, even betas = sin    
end
figure; hold on;
colormap('hot')
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(retrieveAmp(:,cond),cfg.layout)
    colorbar
    if cond < 6
        caxis([0 1]);
    elseif cond > 5 && cond < 11
        caxis([0 2]);
    elseif cond > 10
        caxis([0 2.5]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoFlickSumFqretrieved.png'])
% motion
position = [11 7 3 9 15 13 8];
figure; hold on;
colormap('hot')
for cond=16:length(avData)
    subplot(3,5,position(cond-15)); hold on;
    plotTopo(retrieveAmp(:,cond),cfg.layout);
    colorbar
    if cond < 6
        caxis([0 1]);
    elseif cond > 5 && cond < 11
        caxis([0 2]);
    elseif cond > 10
        caxis([0 2.5]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoMotionSumFqretrieved.png'])



% pool cond sin and cos reg parameter DIFFERENT for all conditions
for cond = 1:length(avData)
    poolData = [];clear beta betaUnscaled
    poolData = [avData(cond).cos(:,avData(cond).harmIdx) avData(cond).sin(:,avData(cond).harmIdx)];
    % match with templates
    [beta, betaUnscaled,lambda(cond)] = fitEEGTemplates(poolData,biosemi128,1);
%     saveas(gcf,['figures' filesep 'lambdaCond' num2str(cond) '.png'])
    % compute amplitude (power=cos^2+sin^2)
    betaPower(:,cond) = sum(betaUnscaled(:,1:2:end).^2 + betaUnscaled(:,2:2:end).^2,2); % odd betas = cos, even betas = sin
    % compute normalised amplitude 
    ampNorm2(:,cond) = sqrt( betaPower(:,cond) / sumSq(cond));
    % check that I retrieve the right topographies
    retrieve2 = biosemi128.weights * betaUnscaled;
    retrieveAmp2(:,cond) = sqrt(sum(retrieve2(:,1:2:end).^2 + retrieve2(:,2:2:end).^2,2)); % odd betas = cos, even betas = sin
end

col={'b','r','g'};
dcVal = [12.5 25 50 75 87.5];
figure;set(gcf, 'Position', [0 0 1600 700])
for roi = 1:18
    subplot(3,6,roi); hold on;
    plot(dcVal,ampNorm2(roi,1:5),['.-'  col{1}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,ampNorm2(roi,6:10),['.-'  col{2}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,ampNorm2(roi,11:15),['.-'  col{3}],'MarkerSize',40,'LineWidth',2)
    plot(50,ampNorm2(roi,18),['^:' col{1}],'MarkerSize',15,'Linewidth',2)
    plot([25 50 75],ampNorm2(roi,[17 22 19]),['^:' col{2}],'MarkerSize',15,'Linewidth',2)
    plot([12.5 50 87.5],ampNorm2(roi,[16 21 20]),['^:' col{3}],'MarkerSize',15,'Linewidth',2)    
    xticks([0:12.5:100]);
    title(biosemi128.listROIs(roi))
%     ylim([0 max(ampNorm2(:))]);
    ylabel('amplitude (a.u.)')
end
xlabel('Duty Cycle')
legend('10','5','2.5','10mov','5mov','2.5mov','Location','Best')

figure; hold on;
colormap('hot')
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(retrieveAmp2(:,cond),cfg.layout)
    colorbar
    if cond < 6
        caxis([0 1]);
    elseif cond > 5 && cond < 11
        caxis([0 2]);
    elseif cond > 10
        caxis([0 2.5]);
    end
end
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoFlickSumFqretrieved2.png'])





% try to do fitEEGtemplates per fq
poolCond = [1:5; 6:10 ; 11:15];
ampNorm3 = [];
for pooling=1:3
    poolData = [];clear beta betaUnscaled power
    for cond=poolCond(pooling,:)
        poolData = [poolData avData(cond).cos(:,avData(cond).harmIdx) avData(cond).sin(:,avData(cond).harmIdx)];
    end
    [beta, betaUnscaled,lambda(pooling)] = fitEEGTemplates(poolData,biosemi128);
    % get the range of harmonics: cumulative sum of elements
    rangeHarm = [0 cumsum(repmat(length(avData(cond).harmIdx),1,5)*2)]; % *2 to account for sin and cos
    % compute amplitude (power=cos^2+sin^2)
    for cond=1:5
        power(:,cond) = sum(betaUnscaled(:,rangeHarm(cond)+1:2:rangeHarm(cond+1)).^2 ...
        + betaUnscaled(:,rangeHarm(cond)+2:2:rangeHarm(cond+1)).^2,2);% odd betas = cos, even betas = sin
    end
    % compute normalised amplitude 
    ampNorm3 = [ampNorm3 sqrt( power ./ sumSq(poolCond(pooling,:)))];
%     % check that I retrieve the right topographies
%     retrieve3 = biosemi128.weights * betaUnscaled;
%     retrieveAmp3(:,pooling) = sqrt(sum(retrieve3(:,1:2:end).^2 + retrieve3(:,2:2:end).^2,2)); % odd betas = cos, even betas = sin
end

col={'b','r','g'};
dcVal = [12.5 25 50 75 87.5];
figure;set(gcf, 'Position', [0 0 1600 700])
for roi = 1:18
    subplot(3,6,roi); hold on;
    plot(dcVal,ampNorm3(roi,1:5),['.-'  col{1}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,ampNorm3(roi,6:10),['.-'  col{2}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,ampNorm3(roi,11:15),['.-'  col{3}],'MarkerSize',40,'LineWidth',2)
%     plot(50,ampNorm3(roi,18),['^:' col{1}],'MarkerSize',15,'Linewidth',2)
%     plot([25 50 75],ampNorm3(roi,[17 22 19]),['^:' col{2}],'MarkerSize',15,'Linewidth',2)
%     plot([12.5 50 87.5],ampNorm3(roi,[16 21 20]),['^:' col{3}],'MarkerSize',15,'Linewidth',2)    
    xticks([0:12.5:100]);
    title(biosemi128.listROIs(roi))
    ylim([0 max(ampNorm3(:))]);
    ylabel('amplitude (a.u.)')
end
xlabel('Duty Cycle')
legend('10','5','2.5','10mov','5mov','2.5mov','Location','Best')





%%% Try to do fitEEGtemplates on the normalised data 
normPower = [avData.sumPower] ./ sumSq; 
cos = real(test);
sin = -imag(test);

% transform to complex numbers to be able to use ifft
% power=cos^2+sin^2. 
% cos=real=amplitude
% sin =-imag 
% amp = abs(dftData);
% cos = real(dftData);
% sin = -imag(dftData);
% phase is angle(complex number: dftData)
abs(complex(avData(2).cos(2,2), avData(2).sin(2,2)))
avData(2).amp(2,2)
pow = avData(2).cos(2,2)^2 + avData(2).sin(2,2)^2;
cosN = sqrt(pow - avData(2).sin(2,2)^2 );

% test on one condition
testPower = (avData(5).cos.^2 + avData(5).sin.^2) ;
cosN = sqrt(testPower - avData(5).sin.^2 ); % can be - or + !! 
wAmp = sqrt(testPower);
wPhase = complex(avData(5).cos,avData(5).sin);
complexBack = cos(wPhase).*wAmp + 1i.*sin(wPhase).*wAmp;
fftBack = fft(complexBack);

nT = 577;
fourierBasis = dftmtx(nT);
test=avData(5).amp*avData(5).freq';

invFourier  = conj(fourierBasis)./nT;
waveForm = complexBack'.*fourierBasis;


wPhase = (-1*wList)+repmat([0*pi/2 0*pi/2]',nF/2,1);
wAmp = cos(wPhase).*wAmp + 1i.*sin(wPhase).*wAmp;
waveForm = wAmp*fourierBasis;
waveForm = real(waveForm); 
    
wTest = complex(avData(5).cos,avData(5).sin);
wTest(:,2:end) = wTest(:,2:end)/2;
wTest2 = complex(avData(5).cos,-avData(5).sin);
wTest3 = fliplr(wTest2(:,2:end)/2);
wTest4 = [wTest wTest3];
soWhat = ifft(wTest4,[],2);
check = real(soWhat);


%%%%%%%% try fitEEGtemplates on waveforms
whichFq = [1:5; 6:10; 11:15];
for freq=1:size(whichFq,1)
    clear beta betaUnscaled rangeHarm
    [beta, betaUnscaled,lambda(freq)] = fitEEGTemplates([avData(whichFq(freq,:)).wave],biosemi128,0,158);
    rangeHarm = [0 cumsum(repmat(size(avData(whichFq(freq,1)).wave,2),1,5))];
    figure;
    for roi=1:18
        subplot(3,6,roi);hold on
        for cond=1:5
            plot(betaUnscaled(roi,rangeHarm(cond)+1:rangeHarm(cond+1)) )
        end
        title(biosemi128.listROIs(roi))
        ylim([min(betaUnscaled(:)) max(betaUnscaled(:))]);
    end
    legend({'12.5' '25' '50' '75' '87.5'})
end
[beta, betaUnscaled,lambda2] = fitEEGTemplates([avData(:).wave],biosemi128,1);
