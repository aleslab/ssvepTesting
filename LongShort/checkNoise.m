% check noise calculation

clearvars;
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};

dataPath = {'/Users/marleneponcet/Documents/data/LongRangeV2/', '/Users/marleneponcet/Documents/data/LRshortDC/V2/'};
pickElec = [23 9 38];


%% plot noise waveform
figure;hold on;
chan=1;
for ee=1:2 % which experiment
    clear avPredictions
    load([dataPath{ee} 'sbjprediction.mat'])

    %%%%%% compute average
    avPredictions = averageSbj(sbj);
    condRange = {'LR','SR'};
    for ss=1:2
%         for chan=1:length(pickElec)
            subplot(2,2,ee+2*(ss-1)); hold on;
            plot(avPredictions(1+5*(ss-1)).time,avPredictions(1+5*(ss-1)).filteredWave(pickElec(chan),:),'LineWidth',2);
            plot(avPredictions(1+5*(ss-1)).time,avPredictions(2+5*(ss-1)).filteredWave(pickElec(chan),:) - avPredictions(1+5*(ss-1)).filteredWave(pickElec(chan),:),'LineWidth',2);
            plot(avPredictions(1+5*(ss-1)).time,avPredictions(1+5*(ss-1)).noiseWave(pickElec(chan),:),'LineWidth',2);
            line([0 400],[0 0],'Color','k','LineStyle','--')
            legend('AM','diff','noise','Location','Best')
            title([condRange{ss} num2str(pickElec(chan))])
%         end
    end 
end
saveas(gcf,['figures' filesep 'noise'],'png')

%% noise for each condition
    figure; hold on;
for ee=1:2 % which experiment
    clear avPredictions
    load([dataPath{ee} 'sbjprediction.mat'])
    avPredictions = averageSbj(sbj);
    subplot(2,2,ee+1*(ee-1)); hold on;
    for ss=1:5
        plot(avPredictions(ss).time,avPredictions(ss).noiseWave(pickElec(chan),:),'LineWidth',2);
    end
    subplot(2,2,ee+1+1*(ee-1)); hold on;
    for ss=6:10
        plot(avPredictions(ss).time,avPredictions(ss).noiseWave(pickElec(chan),:),'LineWidth',2);
    end
end
saveas(gcf,['figures' filesep 'noiseEachCond'],'png')
