%%% load the sbj predictions
load('sbjprediction.mat')
% do the rectification for the time 
for aa=1:size(sbj,1)
    for bb=1:size(sbj,2)
        sbj(aa,bb).data.wave(:,:) = sbj(aa,bb).data.wave(:,[42:204 1:41]);
    end
end

% do the average
avPredictions = averageSbj(sbj(:,:));

%%%%
avPredictions(1).condLabel = 'LR_originalmotion';
avPredictions(2).condLabel = 'LR_linearPred';
avPredictions(3).condLabel = 'LR_spatialPred';
avPredictions(4).condLabel = 'LR_tempPred';
avPredictions(5).condLabel = 'LR_spat&tempPred';
avPredictions(6).condLabel = 'SR_originalmotion';
avPredictions(7).condLabel = 'SR_linearPred';
avPredictions(8).condLabel = 'SR_spatialPred';
avPredictions(9).condLabel = 'SR_tempPred';
avPredictions(10).condLabel = 'SR_spat&tempPred';
avPredictions(11).condLabel = 'LR_STnl';
avPredictions(12).condLabel = 'SR_STnl';

% create a new field "filteredWave" after low-pass filter
% warning about filter is NORMAL
for condIdx=1:length(avPredictions)
    filtIdx = determineFilterIndices( 'low49', avPredictions(condIdx).freq, avPredictions(condIdx).i1f1 );
    
    %Create a logical matrix selecting frequency components.
    filtMat = false(size(avPredictions(condIdx).amp));
    filtMat(:,filtIdx) = true;
    
    %Combine the filter and sig vaules with a logical AND.
    % condData(condIdx).activeFreq = (filtMat.*sigFreqs)>=1; %Store the filtered coefficients for the spec plot
    avPredictions(condIdx).activeFreq = (filtMat)>=1; %Store the filtered coefficients for the spec plot
    cfg.activeFreq =  avPredictions(condIdx).activeFreq;
    
    [ avPredictions(condIdx).filteredWave ] = filterSteadyState( cfg, avPredictions(condIdx) );
end


%%% Compute baseline for the diff between each prediction and motion
for ss = 1:size(sbj,1)
    for cond=1:size(sbj,2)
        dataWave(:,:,cond,ss) = sbj(ss,cond).data.wave;
    end
end
condCompare = 2:5; fixCond = 1;
condCompare = 7:10; fixCond = 6;
[baseDiff_LR, baseCI_LR]= compute_baseline(dataWave,fixCond,condCompare); %data,motion condition, conditions to compare to

%%% plot the topo and ERP differences between motion and 1/L+R 2/L+R+S
%%% 3/L+R+T 4/L+R+TS
elec = 23;

figure;hold on;
for cond=1:4
    subplot(2,2,cond); hold on;
    plot(avPredictions(fixCond).time,avPredictions(fixCond).filteredWave(elec,:),'LineWidth',2);
    plot(avPredictions(fixCond).time,avPredictions(condCompare(cond)).filteredWave(elec,:),'LineWidth',2)
    plot(avPredictions(fixCond).time,avPredictions(fixCond).filteredWave(elec,:) - avPredictions(condCompare(cond)).filteredWave(elec,:),'LineWidth',4)
%     plot(avPredictions(fixCond).time,baseDiff_LR(elec,:,cond),'k')
%     fill([avPredictions(fixCond).time avPredictions(fixCond).time(end:-1:1) avPredictions(fixCond).time(1)],...
%         [baseDiff_LR(elec,:,cond)-baseCI_LR(elec,:,cond) ...
%         baseDiff_LR(elec,end:-1:1,cond)+baseCI_LR(elec,end:-1:1,cond) ...
%         baseDiff_LR(elec,1,cond)-baseCI_LR(elec,1,cond)], 'k','EdgeColor','none','facealpha',0.1);
    ylim([-2.5 3])
    xlabel('time (ms)')
    ylabel('Oz amplitude')
    legend(avPredictions(fixCond).condLabel,avPredictions(condCompare(cond)).condLabel,'difference')
end



figure;plotTopo(mean(avPredictions(1).filteredWave(:,:)-avPredictions(2).filteredWave(:,:),2),cfg.layout);colorbar;
plot(avPredictions(1).time,avPredictions(1).filteredWave(elec,[42:204 1:41]))
compCond = 1;
figure; hold on;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Topographies to determine which electrodes to include
condition = [1 6]; condName = {'long range' 'short range'};
harm = avPredictions(1).freq(avPredictions(1).i1f1).*[1:6];
for hh=1:length(harm)
    harmIndex(hh) = find(avPredictions(1).freq==harm(hh));
end

% for the short and long range conditions = noisy
for cc=1:length(condition)
    figure;
    for hh=1:length(harmIndex)
        subplot(2,length(harmIndex)/2,hh); hold on;
        plotTopo(avPredictions(condition(cc)).amp(:,harmIndex(hh)),cfg.layout);
        colorbar;
%         caxis([0 0.8]);
        title([num2str(harm(hh)) 'Hz'])
    end
    set(gcf,'PaperPositionMode','auto')
    set(gcf, 'Position', [0 0 2000 2500])
    saveas(gcf,['figures/topo' condName{cc} '.png'])
end

% for left stim at short and long
condition = [2 8];
for cc=1:length(condition)
    figure;
    for hh=1:length(harmIndex)
        subplot(2,length(harmIndex)/2,hh); hold on;
        plotTopo(avData(condition(cc)).amp(:,harmIndex(hh)),cfg.layout);
        colorbar;
%         caxis([0 0.8]);
        title([num2str(harm(hh)) 'Hz'])
    end
    set(gcf,'PaperPositionMode','auto')
    set(gcf, 'Position', [0 0 2000 2500])
end