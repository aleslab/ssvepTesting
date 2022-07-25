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
% pool cond sin and cos reg parameter SAME for all conditions
poolData = [];
for cond = 1:length(avData)
    poolData = [poolData avData(cond).cos(:,avData(cond).harmIdx) avData(cond).sin(:,avData(cond).harmIdx)];
end
% match with templates
[beta, betaUnscaled] = fitEEGTemplates(poolData,biosemi128,1);
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
    xticks([0:12.5:100]);
    title(biosemi128.listROIs(roi))
%     ylim([0 max(ampNorm(:))]);
    ylabel('amplitude (a.u.)')
end
xlabel('Duty Cycle')
legend('10','5','2.5','10mov','5mov','2.5mov','Location','Best')
saveas(gcf,'figures/betaSumFq','png')


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
pie(ampNorm(:,cond)/sum(ampNorm(:,cond)))
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
    saveas(gcf,['figures' filesep 'lambdaCond' num2str(cond) '.png'])
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
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoFlickSumFqretrieved2.png'])
