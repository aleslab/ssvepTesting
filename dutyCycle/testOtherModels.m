% regression with the area loc model fitting the harmonics separately
addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/biosemiUpdated
addpath /Volumes/Amrutam/Marlene/Git/fieldtrip-aleslab-fork
ft_defaults


clearvars

modelName = 'model-501.mat';

%%% 502: p = columns 
%%% 503: p = columns 
%%% 504: p < columns

% load models
load(['/Volumes/Amrutam/Marlene/JUSTIN/topoMapForBiosemi/' modelName])
% load data
load('16sbjDC.mat')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do the regression 
% this is without normalisation (useful for plotting the predicted maps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for cond = 1:length(avData)
    for hh = 1 : length(avData(cond).harmIdx)
        fitData(cond).betaCosRaw(:,hh)  = regress (avData(cond).cos(:,avData(cond).harmIdx(hh)), minModel.cos);
        fitData(cond).betaSinRaw(:,hh)  = regress (avData(cond).sin(:,avData(cond).harmIdx(hh)), minModel.cos);
    end
end


%%% create topographies from betas and model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('squareFFT.mat') % load the amplitudes for a square function
for cond = 1: length(avData)
    for hh = 1 : length(avData(cond).harmIdx)
        pred(cond).cos(hh,:) = minModel.cos * fitData(cond).betaCosRaw(:,hh);
        pred(cond).sin(hh,:) = minModel.cos * fitData(cond).betaSinRaw(:,hh);
    end
    % for summing harmonics do sum of power, NOT amplitudes!
    % same as: 
    % test = sqrt(pred(cond).cos(:,:).^2 + pred(cond).sin(:,:).^2);
    % test2 = sqrt(sum(test.^2));
    pred(cond).power = sum(pred(cond).cos.^2 + pred(cond).sin.^2);
    pred(cond).sumAmp = sqrt(sum(pred(cond).cos.^2 + pred(cond).sin.^2));
    % normalise the amplitudes
    sumSq = sum(sqAllFFT(avData(cond).harmIdx,cond).^2);
    pred(cond).normAmp = pred(cond).sumAmp / sqrt(sumSq);     
end

figure('position', [200, 0, 1500, 1000])
colormap('hot')
for cond = 1:15
    subplot(3,5,cond)
    plotTopo(pred(cond).normAmp,cfg.layout)
    if cond < 6
        caxis([0 1]);
    elseif cond > 5 && cond < 11
        caxis([0 3]);
    elseif cond > 10
        caxis([0 3]);
    end
    colorbar
end
saveas(gcf,['figures' filesep modelName(1:end-4) filesep 'predTopoFlickAmpNorm'],'png')

figure('position', [200, 0, 1500, 1000])
colormap('hot')
for cond = 1:15
    subplot(3,5,cond)
    plotTopo(pred(cond).power,cfg.layout)
    if cond < 6
        caxis([0 1.5]);
    elseif cond > 5 && cond < 11
        caxis([0 6.5]);
    elseif cond > 10
        caxis([0 7.5]);
    end
    colorbar
end
saveas(gcf,['figures' filesep modelName(1:end-4) filesep 'predTopoFlickPower'],'png')

figure('position', [200, 0, 1500, 1000])
colormap('hot')
position = [11 7 3 9 15 13 8];
for cond = 16:22
    subplot(3,5,position(cond-15))
    plotTopo(pred(cond).normAmp,cfg.layout)
    if position(cond-15) < 6
        caxis([0 1]);
    elseif position(cond-15) > 5 && position(cond-15) < 11
        caxis([0 3]);
    elseif position(cond-15) > 10
        caxis([0 3]);
    end
    colorbar
end
saveas(gcf,['figures' filesep modelName(1:end-4) 'predTopoMotionAmpNorm'],'png')

figure('position', [200, 0, 1500, 1000])
colormap('hot')
position = [11 7 3 9 15 13 8];
for cond = 16:22
    subplot(3,5,position(cond-15))
    plotTopo(pred(cond).power,cfg.layout)
    if position(cond-15) < 6
        caxis([0 1.5]);
    elseif position(cond-15) > 5 && position(cond-15) < 11
        caxis([0 6.5]);
    elseif position(cond-15) > 10
        caxis([0 7.5]);
    end
    colorbar
end
saveas(gcf,['figures' filesep modelName(1:end-4) 'predTopoMotionPower'],'png')










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do the regression with a normalised model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do a normalisation for each ROI -> unit norming
% so that the betas represent microVolts (instead of microVolts/area size
% as it is now)
% unit norming is: all electrodes are squared and summed. These values are
% then divided so that the total of the electrodes for each ROI (power) is
% equal to 1
model = minModel.cos;
modelNorm = sqrt(sum(minModel.cos.^2,1));
model = bsxfun(@rdivide,model,modelNorm);

for cond = 1:length(avData)
    for hh = 1 : length(avData(cond).harmIdx)
        fitData(cond).betaCos(:,hh)  = regress (avData(cond).cos(:,avData(cond).harmIdx(hh)), model);
        fitData(cond).betaSin(:,hh)  = regress (avData(cond).sin(:,avData(cond).harmIdx(hh)), model);
    end
end

load('squareFFT.mat') % load the amplitudes for a square function

for cond = 1:length(avData)
    % sum power harmonics then amplitude
    fitData(cond).sumAmp =  sqrt(sum(fitData(cond).betaCos.^2 + fitData(cond).betaSin.^2,2));
    % normalise the amplitudes
    sumSq = sum(sqAllFFT(avData(cond).harmIdx,cond));
    fitData(cond).normAmp = fitData(cond).sumAmp / sqrt(sumSq);    
end

% easier format for plotting
allNormAmp = [fitData.normAmp];

% plot sources * 22 conditions
col={'b','r','g'};
dcVal = [12.5 25 50 75 87.5];
figure;set(gcf, 'Position', [0 0 1200 800])
for area = 1:size(allNormAmp,1)
    subplot(2,round(size(allNormAmp,1)/2),area); hold on;
    plot(dcVal,allNormAmp(area,1:5),['.-'  col{1}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,allNormAmp(area,6:10),['.-'  col{2}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,allNormAmp(area,11:15),['.-'  col{3}],'MarkerSize',40,'LineWidth',2)
    plot(50,allNormAmp(area,18),['^:' col{1}],'MarkerSize',15,'Linewidth',2)
    plot([25 50 75],allNormAmp(area,[17 22 19]),['^:' col{2}],'MarkerSize',15,'Linewidth',2)
    plot([12.5 50 87.5],allNormAmp(area,[16 21 20]),['^:' col{3}],'MarkerSize',15,'Linewidth',2)
    title(minModel.locNames(area))
    xticks([0:12.5:100]);
%     ylim([0 12])
    legend('10','5','2.5','Location','Best')
    xlabel('Duty Cycle')
    ylabel('betaNormAllFq')
end
saveas(gcf,['figures' filesep modelName(1:end-4) filesep 'betaNormAllFqPerArea'],'png')



%%%% correlation with behavioural percept
load('fullRatings9.mat')
static = mean(tabStatic,3); moving = mean(tabMot,3);

flickRating = static(2:3,:)';
movRating = moving([3 5 11 15 9 8]);

figure;set(gcf, 'Position', [0 0 1000 600])
for area = 1:size(allNormAmp,1)
    subplot(2,round(size(allNormAmp,1)/2),area); hold on;
    scatter(allNormAmp(area,6:15),flickRating(:), 80,'filled','MarkerEdgeColor','none');
    R = corrcoef(allNormAmp(area,6:15),flickRating(:));
    RsqS = R(1,2).^2;
    scatter(allNormAmp(area,[16:17 19:22]),movRating, 80,'filled','MarkerEdgeColor','none');
    R = corrcoef(allNormAmp(area,[16:17 19:22]),movRating);
    RsqM = R(1,2).^2;
%     R = corrcoef([allNormAmp(area,6:15) allNormAmp(area,[16:17 19:22])],[flickRating(:)' movRating]);
%     RsqF = R(1,2).^2;
    ylim([0 3]);
    xlabel('betaNormECC')
    ylabel('motion rating')
    title([minModel.locNames{area} ' S=' num2str(RsqS,2) ' M=' num2str(RsqM,2)]) %only 2 digits
    lsline
end
% Rsquare reported for static and moving 
saveas(gcf,['figures' filesep modelName(1:end-4) filesep 'correlWithBetaNorm'],'png')



